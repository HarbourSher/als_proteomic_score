import {predictLocal} from './engine.js';
let model=null;
self.onmessage=async ({data})=>{
  if(data.type==='init'){model=data.model;self.postMessage({type:'ready'});return;}
  if(data.type==='single'){try{self.postMessage({type:'single',requestId:data.requestId,result:predictLocal(model,data.values,data.referenceId)});}catch(e){self.postMessage({type:'error',requestId:data.requestId,error:e.message});}return;}
  if(data.type==='batch'){
    const results=[];
    for(let i=0;i<data.rows.length;i++) {const row=data.rows[i];try{results.push({id:row.id,row:row.row,...predictLocal(model,row.values,data.referenceId,{explain:false})});}catch(e){results.push({id:row.id,row:row.row,status:'error',error:e.message});}if((i+1)%10===0)self.postMessage({type:'progress',jobId:data.jobId,done:i+1,total:data.rows.length});}
    self.postMessage({type:'batch',jobId:data.jobId,results});
  }
};
