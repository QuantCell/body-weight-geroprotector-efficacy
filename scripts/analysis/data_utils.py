"""Apply the adult ITP weight range."""
from pathlib import Path

import numpy as np
import pandas as pd

MIN_WEIGHT_G = 10.0
MAX_WEIGHT_G = 100.0
WEIGHT_COLUMNS = ("bw_6", "bw_12", "bw_18", "bw_24")


def clean_adult_weights(data, columns=WEIGHT_COLUMNS):
    data = data.copy()
    for column in columns:
        if column in data:
            values = pd.to_numeric(data[column], errors="coerce")
            data[column] = values.where(values.between(MIN_WEIGHT_G, MAX_WEIGHT_G))
    return data


def write_weight_qc():
    root = Path(__file__).absolute().parents[2]
    data = pd.read_excel(root / "data/inputs/ITP_bodyweight_all_cohorts.xlsx")
    rows = []
    for column in WEIGHT_COLUMNS:
        values = pd.to_numeric(data[column], errors="coerce")
        # List out-of-range measurements and retained upper-tail weights.
        flagged = values.notna() & ((values < MIN_WEIGHT_G) | (values > 80))
        for index in data.index[flagged]:
            value = values.loc[index]
            excluded = not (MIN_WEIGHT_G <= value <= MAX_WEIGHT_G)
            rows.append({
                **data.loc[index, ["cohort", "site", "id", "sex", "group"]].to_dict(),
                "weight_age_months": int(column.split("_")[1]),
                "source_weight_g": value,
                "analysis_weight_g": np.nan if excluded else value,
                "excluded": excluded,
                "rule": "outside 10-100 g" if excluded else "within 10-100 g",
            })
    result = pd.DataFrame(rows)
    result.to_csv(root / "data/derived/adult_weight_qc.csv", index=False)
    print(f"Adult-weight QC: {result.excluded.sum()} excluded measurements; "
          f"{(~result.excluded).sum()} reviewed upper-tail measurements retained.")




def initiation_months(values):
    return pd.to_numeric(values.astype(str).str.extract(r"(\d+(?:\.\d+)?)")[0],
                         errors="coerce")
