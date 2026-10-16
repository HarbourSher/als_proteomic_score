# ALS Proteomic Score

A browser-based research scoring tool for ALS protein panels. Supports a 20-protein SomaScan model, a 13-protein Olink signature and a 5-protein SAHZU serum mass-spectrometry signature, with English and Chinese interfaces.

[Open the calculator](https://harboursher.github.io/als_proteomic_score/) · [Alternative website](https://als-proteomic-score.harboursher.chatgpt.site)

## Use and interpretation

- Enter a single sample or import a CSV, XLSX or XLS file. Missing proteins use the frozen reference values; partial panels are exploratory.
- Review protein contributions, export results as CSV, or print a report using the browser's Print / Save PDF function.
- Use only measurements matching the selected assay and reference scale. New laboratory batches require calibration.
- Scores are for research, not independent ALS diagnosis or exclusion, and are not absolute disease probabilities. Scores from different platforms are not directly comparable.

SomaScan uses a separately trained protein-only probability forest. Olink and serum MS use fixed signatures; their cross-platform evaluations do not validate the full 20-protein forest. Olink and SAHZU reference mean/SD were estimated from controls in their respective evaluation cohorts. Serum scoring has no established diagnostic classification threshold. The website's expandable methods section provides model definitions, performance summaries and limitations.

## Run locally

This repository contains the static website and frozen model files. No installation or model training is needed. With Python 3 installed, run from the repository directory:

```sh
python -m http.server 4173 --bind 127.0.0.1
```

Open `http://127.0.0.1:4173/` in a browser. Serve the files over HTTP rather than opening `index.html` directly.

## Data handling and reproducibility

Protein inputs, file parsing, scoring and exports are processed in the browser. The application has no input-upload API, account requirement, tracking scripts or persistent input storage. Batch input supports up to 10,000 rows and 25 MB; use anonymous sample identifiers.

The `models/` directory contains frozen parameters and aggregate verification summaries, not participant-level background records. Original research data and restricted validation fixtures are not distributed in this static repository. Serving the site reproduces the deployed calculations; it does not retrain the models or independently reproduce the clinical evaluation.

## Third-party software

SheetJS is bundled for spreadsheet parsing. TreeSHAP follows the SHAP algorithm implementation. Their licenses are included in `vendor/`. The website does not require runtime CDN access.
