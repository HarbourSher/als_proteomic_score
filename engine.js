// Tree path dependent SHAP recursion, following Lundberg et al. and SHAP's BSD implementation.
// Only tree aggregate node statistics are used; no participant background rows are shipped.
export function forestProbability(model, x) {
  let sum=0;
  for(const t of model.trees) {let n=0;while(t.feature[n]>=0)n=x[t.feature[n]]<=t.threshold[n]?t.left[n]:t.right[n];sum+=t.value[n];}
  return sum/model.trees.length;
}
function extend(p,d,z,o,f) {
  p.f[d]=f;p.z[d]=z;p.o[d]=o;p.w[d]=d===0?1:0;
  for(let i=d-1;i>=0;i--){p.w[i+1]+=o*p.w[i]*(i+1)/(d+1);p.w[i]=z*p.w[i]*(d-i)/(d+1);}
}
function unwind(p,d,k) {
  const o=p.o[k],z=p.z[k];let next=p.w[d];
  for(let i=d-1;i>=0;i--){if(o!==0){const old=p.w[i];p.w[i]=next*(d+1)/((i+1)*o);next=old-p.w[i]*z*(d-i)/(d+1);}else p.w[i]=p.w[i]*(d+1)/(z*(d-i));}
  for(let i=k;i<d;i++){p.f[i]=p.f[i+1];p.z[i]=p.z[i+1];p.o[i]=p.o[i+1];}
}
function pathSum(p,d,k) {
  const o=p.o[k],z=p.z[k];let next=p.w[d],sum=0;
  for(let i=d-1;i>=0;i--){if(o!==0){const v=next*(d+1)/((i+1)*o);sum+=v;next=p.w[i]-v*z*(d-i)/(d+1);}else sum+=(p.w[i]/z)/((d-i)/(d+1));}
  return sum;
}
function visit(t,x,phi,n,d,path,z,o,f) {
  const p={f:path.f.slice(),z:path.z.slice(),o:path.o.slice(),w:path.w.slice()};extend(p,d,z,o,f);
  const split=t.feature[n];
  if(split<0){for(let i=1;i<=d;i++)phi[p.f[i]]+=pathSum(p,d,i)*(p.o[i]-p.z[i])*t.value[n];return;}
  const hot=x[split]<=t.threshold[n]?t.left[n]:t.right[n],cold=hot===t.left[n]?t.right[n]:t.left[n];
  const hz=t.cover[hot]/t.cover[n],cz=t.cover[cold]/t.cover[n];let iz=1,io=1;
  let k=0;while(k<=d&&p.f[k]!==split)k++;
  if(k<=d){iz=p.z[k];io=p.o[k];unwind(p,d,k);d--;}
  if(hz*iz!==0||io!==0)visit(t,x,phi,hot,d+1,p,hz*iz,io,split);
  if(cz*iz!==0)visit(t,x,phi,cold,d+1,p,cz*iz,0,split);
}
function treeExpected(t,n=0) {return t.feature[n]<0?t.value[n]:(t.cover[t.left[n]]*treeExpected(t,t.left[n])+t.cover[t.right[n]]*treeExpected(t,t.right[n]))/t.cover[n];}
export function forestBaseline(model) {return model.trees.reduce((s,t)=>s+treeExpected(t),0)/model.trees.length;}
export function treeShap(model,x) {
  const phi=Array(model.proteins.length).fill(0);
  for(const t of model.trees)visit(t,x,phi,0,0,{f:[],z:[],o:[],w:[]},1,1,-1);
  return phi.map(v=>v/model.trees.length);
}
export function prepareValues(model,values) {
  if(!values||typeof values!=='object'||Array.isArray(values))throw new Error('invalid_values');
  const x=[],observed=[],raw=[];
  for(let i=0;i<model.proteins.length;i++) {
    const p=model.proteins[i],v=values[p.gene];
    const missing=v===undefined||v===null||(typeof v==='string'&&v.trim()==='');
    if(!missing&&!(typeof v==='number'||typeof v==='string'))throw new Error(`invalid_value:${p.gene}`);
    const n=missing?null:Number(v);
    if(!missing&&(!Number.isFinite(n)||(model.platform==='somascan'&&n<0)||(model.platform==='serum'&&n<=0)))throw new Error(`invalid_value:${p.gene}`);
    observed.push(!missing);raw.push(n);
    if(model.platform==='somascan')x.push(((missing?model.preprocessing.median[i]:Math.log2(n+1))-model.preprocessing.center[i])/model.preprocessing.scale[i]);
    else x.push(missing?0:((model.platform==='serum'?Math.log2(n):n)-p.center)/p.scale);
  }
  if(!observed.some(Boolean))throw new Error('no_applicable_proteins');
  return {x,observed,raw};
}
export function predictLocal(model,values,referenceId,{explain=true}={}) {
  if(referenceId!==model.referenceId)throw new Error('incompatible_reference');
  const {x,observed,raw}=prepareValues(model,values),n=observed.filter(Boolean).length;
  let probability,riskScore,classification,contributions,modelOutput;
  if(model.platform==='somascan') {
    probability=forestProbability(model,x);modelOutput=probability;
    const min=model.score.probabilityMin,max=model.score.probabilityMax;
    riskScore=Math.min(100,Math.max(0,100*(probability-min)/(max-min)));
    classification=probability>=model.threshold?'ALS':'non-ALS';contributions=explain?treeShap(model,x):Array(x.length).fill(null);
  } else {
    contributions=x.map((z,i)=>z*model.proteins[i].direction*model.proteins[i].weight);
    modelOutput=contributions.reduce((s,v)=>s+v,0);
    probability=1/(1+Math.exp(-(model.intercept+model.slope*modelOutput)));riskScore=100*probability;
    classification=model.threshold==null?null:(modelOutput>=model.threshold?'ALS':'non-ALS');
  }
  const importance=model.platform==='somascan'?model.importance:model.proteins.map(p=>Math.abs(p.weight));
  if(!importance)throw new Error('model_not_finalized');
  const sum=importance.reduce((s,v)=>s+v,0),coverage=importance.reduce((s,v,i)=>s+(observed[i]?v:0),0)/sum;
  const partial=n<model.proteins.length;
  return {status:'ok',platform:model.platform,version:model.version,referenceId,classification,probability,modelOutput,riskScore,threshold:model.threshold,observedCount:n,totalCount:model.proteins.length,importanceCoverage:coverage,partial,warnings:partial?['partial_input_exploratory','combination_not_individually_validated']:[],contributionType:model.platform==='somascan'?'TreeSHAP (probability scale)':'signed standardized signature contribution',proteins:model.proteins.map((p,i)=>({gene:p.gene,raw:raw[i],observed:observed[i],imputed:!observed[i],inputUsed:model.platform==='somascan'?(observed[i]?raw[i]:2**model.preprocessing.median[i]-1):(observed[i]?raw[i]:(model.platform==='serum'?2**p.center:p.center)),importance:importance[i]/sum,contribution:contributions[i]}))};
}
export function parseCSV(text) {
  const rows=[];let row=[],cell='',quoted=false;
  text=text.replace(/^\uFEFF/,'');
  for(let i=0;i<text.length;i++){const c=text[i];if(quoted){if(c==='"'&&text[i+1]==='"'){cell+='"';i++;}else if(c==='"')quoted=false;else cell+=c;}else if(c==='"'){if(cell!=='')throw new Error('malformed_csv');quoted=true;}else if(c===','){row.push(cell);cell='';}else if(c==='\n'||c==='\r'){if(c==='\r'&&text[i+1]==='\n')i++;row.push(cell);if(row.some(v=>v!==''))rows.push(row);row=[];cell='';}else cell+=c;}
  if(quoted)throw new Error('malformed_csv');row.push(cell);if(row.some(v=>v!==''))rows.push(row);return rows;
}
export function stageTable(model,table) {
  if(!Array.isArray(table)||table.length<2)throw new Error('empty_file');
  if(table.length>10001)throw new Error('too_many_rows');
  const aliases=new Map();for(const p of model.proteins)for(const s of [p.gene,p.seqId,p.olinkId,p.uniprot].filter(Boolean))aliases.set(s.toLowerCase(),p.gene);
  const headers=table[0].map(v=>String(v??'').trim()),mapped=headers.map(h=>aliases.get(h.toLowerCase())??null);
  const cols=mapped.filter(Boolean);if(!cols.length)throw new Error('no_applicable_columns');if(new Set(cols).size!==cols.length)throw new Error('duplicate_protein_columns');
  const ids=headers.map(h=>h.toLowerCase()).map((h,i)=>h==='anonymous_id'?i:-1).filter(i=>i>=0);if(ids.length>1)throw new Error('duplicate_id_columns');
  const refCols=headers.map((h,i)=>h.toLowerCase()==='reference_id'?i:-1).filter(i=>i>=0),platformCols=headers.map((h,i)=>h.toLowerCase()==='platform'?i:-1).filter(i=>i>=0);
  if(refCols.length>1||platformCols.length>1)throw new Error('duplicate_metadata_columns');
  const known=new Set(),rows=[];
  for(let i=1;i<table.length;i++){const cells=table[i];if(!cells.some(v=>v!==null&&v!==undefined&&String(v).trim()!==''))continue;
    const id=ids.length?String(cells[ids[0]]??'').trim()||`S${String(i).padStart(4,'0')}`:`S${String(i).padStart(4,'0')}`;
    const values={};mapped.forEach((g,j)=>{if(g)values[g]=cells[j]??null;});let error=null;
    if(known.has(id))error='duplicate_anonymous_id';known.add(id);
    if(refCols.length&&String(cells[refCols[0]]??'').trim()!==model.referenceId)error='incompatible_reference';
    if(platformCols.length&&String(cells[platformCols[0]]??'').trim()!==model.platform)error='incompatible_platform';
    if(cells.length>headers.length&&cells.slice(headers.length).some(v=>String(v??'').trim()!==''))error='extra_cells';
    try{prepareValues(model,values);}catch(e){error??=e.message;}
    rows.push({row:i+1,id,values,error});
  }
  return {headers,mapping:headers.map((h,i)=>({header:h,gene:mapped[i]})),ignored:headers.filter((_,i)=>!mapped[i]&&i!==ids[0]&&i!==refCols[0]&&i!==platformCols[0]),rows,validCount:rows.filter(r=>!r.error).length};
}
export function encodeCSV(table) {return '\uFEFF'+table.map(row=>row.map(v=>{const s=String(v??'');const safe=/^[=+@\t\r]/.test(s)||(/^-(?!\d|\.)/.test(s))?'\''+s:s;return '"'+safe.replaceAll('"','""')+'"';}).join(',')).join('\r\n');}
