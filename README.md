# Body weight modifies the efficacy of lifespan-extending interventions in *ad libitum*-fed mice

## Running

Python 3.13, R 4.5.1 and the R package `survival` 3.8-3 were used.

```bash
python -m pip install -r requirements.txt
python run_all.py
```

`run_all.py` runs the scripts in `scripts` in order, writes the result tables to `data/outputs` and
the figures to `figures`, checks the values reported in the manuscript and compares all files with
`checksums.sha256`. A full run takes about five minutes.

## Contents

- `data/inputs`: mouse-level and external input data
- `data/outputs`: result tables
- `scripts`: analysis steps 01–11 and shared functions (`common.py`)
- `figures`: Figures 1–5 and Extended Data Fig. 1 (PNG, TIFF and PDF), with source data in `figures/source_data`
- `supplementary`: Supplementary Data 1–10. Sheets 2, 5, 9 and 10 are analysis or figure inputs.

## Data sources

ITP body weights and survival are from the
[Mouse Phenome Database](https://phenome.jax.org/centers/ITP) (projects ITP1 and ITP2004–ITP2021)
and [GeneNetwork](https://genenetwork.org/show_trait?dataset=HET3-ITPPublish&trait_id=50601)
(HET3-ITP), harmonized in `data/inputs/ITP_bodyweight_all_cohorts.xlsx` (Supplementary Data 2).
External sources are the published Diversity Outbred dietary-restriction release and the source data
of Garratt et al. 2022 (*eLife* 11, e84060). The cross-study corpus is from
[DrugAge](https://genomics.senescence.info/drugs/) (Build 5, accessed 3 September 2026) and subsequent ITP reports. Strain data are
from Mouse Phenome Database measures 32603, 32680 and 10331 and from Yuan et al. 2020 (*Aging Cell* 19, e13263).

## Licence

Code is licensed under AGPL-3.0-or-later (`LICENSE`); data terms are in `DATA_LICENSE.md`.
