export const localValidation = {
  "schemaVersion": 1,
  "analysisDate": "2026-10-06",
  "evidenceType": "five-protein cross-platform signature evaluation",
  "proteins": [
    "CKM",
    "ANTXR2",
    "COMP",
    "CHAD",
    "SFRP4"
  ],
  "cohort": {
    "restrictedN": 154,
    "completeN": 153,
    "alsN": 90,
    "controlN": 63,
    "excludedOtherDiagnoses": 4,
    "missingSFRP4": 1
  },
  "controlGroups": [
    {
      "label": "AD",
      "n": 18
    },
    {
      "label": "HD",
      "n": 20
    },
    {
      "label": "Peripheral neuropathy",
      "n": 25
    }
  ],
  "performance": [
    {
      "estimate": 0.783950617283951,
      "ci_low": 0.703693576388889,
      "ci_high": 0.86093963569224,
      "metric": "roc_auc"
    },
    {
      "estimate": 0.797248774171951,
      "ci_low": 0.722426552658976,
      "ci_high": 0.883693292081399,
      "metric": "average_precision"
    }
  ],
  "subgroups": [
    {
      "comparison": "ALS vs Peripheral neuropathy",
      "n": 115,
      "auc": 0.649333333333333
    },
    {
      "comparison": "ALS vs AD",
      "n": 108,
      "auc": 0.885802469135803
    },
    {
      "comparison": "ALS vs HD",
      "n": 110,
      "auc": 0.860555555555555
    }
  ],
  "missingSensitivity": {
    "n": 154,
    "auc": 0.786,
    "ap": 0.797,
    "method": "single missing SFRP4 standardized value set to zero"
  },
  "directionConcordance": {
    "concordant": 5,
    "tested": 5,
    "spearman": 0.9
  },
  "bootstrap": {
    "reps": 2000,
    "seed": 20261006,
    "reestimateControlReference": true
  },
  "limitations": [
    "Cohort-specific disease-control reference",
    "Not full 20-protein random-forest validation",
    "No validated clinical threshold or absolute-risk calibration"
  ],
  "aggregateSources": [
    "signature_performance.csv",
    "subgroup_performance.csv",
    "missingness_sensitivity.csv",
    "local_serum_external_validation_report.md"
  ]
};
