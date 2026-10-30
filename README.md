# Body weight modifies the efficacy of lifespan-extending interventions in *ad libitum*-fed mice

Code and data supporting the manuscript analyses and figures.

## Reproduction

Tested with Python 3.12.14 and 3.13.14, R 4.5.1 and `survival` 3.8-3.

```bash
python -m pip install -r requirements.txt
python run_all.py
```

The runner refits the analyses and verifies the released results, supplement, figures and checksums. Use `--python` or `--rscript` to specify an executable.

## Contents

- `data/inputs`: source measurements, publication mappings and external-study assessment
- `data/outputs`: manuscript results
- `scripts/analysis`: analysis, figure, table and verification scripts
- `figures`: Figures 1–5 and Extended Data Fig. 1 in PNG, TIFF and PDF
- `supplementary/Supplementary_Data_1-10.xlsx`: the ten supplementary datasets
- `checksums.sha256`: publication-file checksums

## Data sources

ITP measurements are from the [Mouse Phenome Database](https://phenome.jax.org/centers/ITP) (ITP2004–ITP2021; survival project ITP1) and [GeneNetwork HET3-ITP](https://genenetwork.org/show_trait?dataset=HET3-ITPPublish&trait_id=50601). External datasets, strain surveys and publication sources are identified in the supplement and `DATA_LICENSE.md`.

Copyright (c) 2026 the authors. Code is licensed under AGPL-3.0-or-later. Source data retain their providers' terms; see `DATA_LICENSE.md`.
