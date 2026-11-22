# Lifespan extension is weight-dependent in *ad libitum*-fed mice

## Running

Install Python, R and the R package `survival` before running the commands below.
Tested on 64-bit Windows 11 Pro with Python 3.13.14, R 4.5.1 and `survival` 3.8-3.
No special hardware is required. Installing the Python dependencies took a few
minutes with cached packages. A full run took about 15 minutes.

```bash
python -m pip install -r requirements.txt
python run_all.py
```

If `Rscript` is outside `PATH`, pass its full path with `--rscript`. The pipeline
regenerates the results and figures, checks manuscript values and compares files
against `checksums.sha256`.

## Contents

- `data/inputs` and `data/outputs` hold the input data and result tables.
- `figures` holds Figures 1–5 and their source workbooks.
- `supplementary` holds Supplementary Data 1–10.

## Data sources

ITP records come from the [Mouse Phenome Database](https://phenome.jax.org/centers/ITP)
and [GeneNetwork](https://genenetwork.org/show_trait?dataset=HET3-ITPPublish&trait_id=50601).
The 42-day interaction analysis standardizes weights against controls alive at
treatment initiation within each experiment, site and sex.

External inputs come from the Diversity Outbred dietary-restriction study, the
olfactory-exposure study, Mouse Phenome Database strain surveys and
[DrugAge Build 5](https://genomics.senescence.info/drugs/), supplemented by later ITP reports.

## Licence

Code is licensed under AGPL-3.0-or-later. Data terms are in `DATA_LICENSE.md`.
