"""Shared paths, the adult-weight plausibility range, and reading and writing the result tables."""
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).absolute().parents[1]
DERIVED = ROOT / "data" / "derived"
OUTPUTS = ROOT / "data" / "outputs"
MIN_WEIGHT_G = 10.0
MAX_WEIGHT_G = 100.0


def within_range(values):
    """Return the values that fall inside the plausibility range, others as missing."""
    numeric = pd.to_numeric(values, errors="coerce")
    return numeric.where(numeric.between(MIN_WEIGHT_G, MAX_WEIGHT_G))


def write_results(tables: dict, combined: dict | None = None) -> None:
    """Write each table to data/outputs/<name>.csv."""
    combined = combined or {}
    grouped = {name for names in combined.values() for name in names}
    for name, table in tables.items():
        if name not in grouped:
            table.to_csv(OUTPUTS / f"{name}.csv", index=False, lineterminator="\n")
    for file, names in combined.items():
        stacked = pd.concat([tables[name].assign(model_set=name)[["model_set", *tables[name].columns]]
                             for name in names], ignore_index=True, sort=False)
        stacked.to_csv(OUTPUTS / f"{file}.csv", index=False, lineterminator="\n")


def _restore(column: pd.Series) -> pd.Series:
    """Undo the type changes that sharing a file imposes."""
    if column.dtype == object:
        if set(column.map(str)) <= {"True", "False", "TRUE", "FALSE"}:
            return column.map(lambda value: str(value).lower() == "true")
        try:
            column = pd.to_numeric(column)
        except (ValueError, TypeError):
            return column
    if column.dtype == float and column.notna().all() and (column == column.round()).all():
        return column.astype("int64")
    return column


def read_result(name: str) -> pd.DataFrame:
    """One result table."""
    if (OUTPUTS / f"{name}.csv").exists():
        return pd.read_csv(OUTPUTS / f"{name}.csv", low_memory=False).apply(_restore)
    for path in sorted(OUTPUTS.glob("*.csv")):
        if "model_set" not in pd.read_csv(path, nrows=0).columns:
            continue
        stacked = pd.read_csv(path, low_memory=False)
        rows = stacked[stacked.model_set.eq(name)]
        if len(rows):
            rows = rows.drop(columns="model_set").dropna(axis=1, how="all").reset_index(drop=True)
            return rows.apply(_restore)
    raise KeyError(f"No result table named {name}")
