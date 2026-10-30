"""Prepare analysis datasets from the deposited source files."""
from __future__ import annotations

from pathlib import Path
import warnings

from lifelines import CoxPHFitter
from scipy.stats import ttest_ind
import numpy as np
import pandas as pd

from data_utils import clean_adult_weights, initiation_months, write_weight_qc


def prepare_experiment_contrasts():
    ROOT = Path(__file__).absolute().parents[2]

    def derive(data):
        data = data.copy()
        data["age_initiation"] = initiation_months(data["age_initiation"])
        rows = []
        for (cohort, sex), block in data.groupby(["cohort", "sex"]):
            controls = block[block.group.eq("Control")]
            if controls.empty:
                continue
            cw = controls.bw_18.dropna()
            cd = controls.loc[controls.dead.eq(1), "lifespan_days"].dropna()
            for group, experiment in block.groupby("group"):
                if group == "Control":
                    continue
                init = experiment.age_initiation.dropna()
                tw = experiment.bw_18.dropna()
                td = experiment.loc[experiment.dead.eq(1), "lifespan_days"].dropna()
                rows.append({
                    "cohort": cohort, "sex": sex, "group": group,
                    "age_init": init.iloc[0] if len(init) else np.nan,
                    "n_tw": len(tw), "n_cw": len(cw), "n_td": len(td), "n_cd": len(cd),
                    "wt_pct": 100 * (tw.mean() - cw.mean()) / cw.mean()
                    if len(tw) and len(cw) else np.nan,
                    "life_pct": 100 * (td.median() - cd.median()) / cd.median()
                    if len(td) and len(cd) else np.nan,
                })
        result = pd.DataFrame(rows)
        return result.loc[
            result.age_init.lt(18) & result[["n_tw", "n_cw", "n_td", "n_cd"]].ge(8).all(axis=1)
        ].dropna(subset=["wt_pct", "life_pct"]).reset_index(drop=True)

    mice = clean_adult_weights(pd.read_excel(ROOT / "data/inputs/ITP_bodyweight_all_cohorts.xlsx"))

    result = derive(mice)

    result.to_csv(ROOT / "data/derived/experiment_weight_changes.csv", index=False)

    print(f"Supplementary experiment contrasts: {len(result)} rows.")


def prepare_experiment_models():
    ROOT = Path(__file__).absolute().parents[2]

    SOURCE = ROOT / "data" / "inputs" / "ITP_bodyweight_all_cohorts.xlsx"

    OUTPUT = ROOT / "data" / "derived"

    DAYS_PER_MONTH = 30.4

    WEIGHTS = {6: "bw_6", 12: "bw_12", 18: "bw_18", 24: "bw_24"}

    def standardize(values: pd.Series) -> pd.Series:
        return (values - values.mean()) / values.std(ddof=1)

    def load_mouse_data() -> pd.DataFrame:
        mice = clean_adult_weights(pd.read_excel(SOURCE))
        mice["init"] = initiation_months(mice["age_initiation"])
        return mice[mice["dead"].notna() & mice["lifespan_days"].notna()].copy()

    def derive_weight_effects(mice: pd.DataFrame) -> pd.DataFrame:
        rows: list[dict] = []
        treated = mice[mice["group"] != "Control"]
        for (cohort, group), experiment in treated.groupby(["cohort", "group"]):
            initiation = experiment["init"].iloc[0]
            controls = mice[(mice["cohort"] == cohort) & (mice["group"] == "Control")]
            for age, column in WEIGHTS.items():
                tx_all = experiment[column].dropna()
                ctl_all = controls[column].dropna()
                if len(tx_all) < 40 or len(ctl_all) < 40:
                    continue
                differences, sample_sizes = [], []
                for sex in ["m", "f"]:
                    tx = experiment.loc[experiment["sex"] == sex, column].dropna()
                    ctl = controls.loc[controls["sex"] == sex, column].dropna()
                    if len(tx) < 20 or len(ctl) < 20:
                        continue
                    differences.append(tx.mean() - ctl.mean())
                    sample_sizes.append(len(tx))
                if not differences:
                    continue
                difference = float(np.average(differences, weights=sample_sizes))
                rows.append({
                    "cohort": cohort,
                    "group": group,
                    "init": initiation,
                    "weight_age": age,
                    "post_dosing": bool(not np.isnan(initiation) and age >= initiation),
                    "n_tx": len(tx_all),
                    "diff_g": difference,
                    "pct": 100 * difference / ctl_all.mean(),
                    "p": ttest_ind(tx_all, ctl_all, equal_var=False).pvalue,
                })
        return pd.DataFrame(rows)

    def fit_experiment_interaction(
        mice: pd.DataFrame, cohort: str, group: str, age: int
    ) -> list[dict]:
        column = WEIGHTS[age]
        initiation = mice.loc[(mice["cohort"] == cohort) & (mice["group"] == group), "init"].iloc[0]
        landmark = max(age, initiation if np.isfinite(initiation) else age) * DAYS_PER_MONTH
        experiment = mice[(mice["cohort"] == cohort) & mice["group"].isin([group, "Control"])].copy()
        experiment = experiment[experiment[column].notna() & (experiment["lifespan_days"] > landmark)]
        experiment["Tx"] = (experiment["group"] == group).astype(int)
        rows: list[dict] = []
        for sex in ["m", "f"]:
            frame = experiment[experiment["sex"] == sex].copy()
            if frame["Tx"].nunique() < 2 or len(frame) < 60 or frame["dead"].sum() < 40:
                continue
            frame["wz"] = frame.groupby(["site", "Tx"])[column].transform(standardize)
            frame = frame.dropna(subset=["wz"])
            frame["entry"] = landmark
            frame["Tx_x_wz"] = frame["Tx"] * frame["wz"]
            try:
                model = CoxPHFitter().fit(
                    frame[["lifespan_days", "dead", "entry", "Tx", "wz", "Tx_x_wz", "site"]],
                    "lifespan_days", "dead", entry_col="entry", strata=["site"],
                )
            except Exception:
                continue
            result = model.summary.loc["Tx_x_wz"]
            rows.append({
                "cohort": cohort,
                "group": group,
                "sex": sex,
                "weight_age": age,
                "within_experiment": True,
                "n": len(frame),
                "int_hr": result["exp(coef)"],
                "int_lo": result["exp(coef) lower 95%"],
                "int_hi": result["exp(coef) upper 95%"],
                "int_p": result["p"],
                "coef": result["coef"],
                "se": result["se(coef)"],
            })
        return rows

    def derive_withinexperiment_interactions(mice: pd.DataFrame) -> pd.DataFrame:
        rows: list[dict] = []
        experiments = mice[mice["group"] != "Control"][["cohort", "group"]].drop_duplicates()
        for experiment in experiments.itertuples(index=False):
            rows.extend(fit_experiment_interaction(mice, experiment.cohort, experiment.group, 6))
        return pd.DataFrame(rows)

    def derive_comparator_cells(mice: pd.DataFrame) -> pd.DataFrame:
        rows: list[dict] = []
        treated = mice[mice["group"] != "Control"]
        landmark = 6 * DAYS_PER_MONTH
        for (cohort, group), experiment in treated.groupby(["cohort", "group"]):
            for site in experiment["site"].dropna().unique():
                for sex in ["m", "f"]:
                    tx = experiment[(experiment["site"] == site) & (experiment["sex"] == sex)]
                    ctl = mice[
                        (mice["cohort"] == cohort)
                        & (mice["group"] == "Control")
                        & (mice["site"] == site)
                        & (mice["sex"] == sex)
                    ]
                    if len(tx) < 25 or len(ctl) < 25:
                        continue
                    frame = pd.concat([tx.assign(Tx=1), ctl.assign(Tx=0)])
                    frame = frame[frame["lifespan_days"] > landmark].copy()
                    if frame["dead"].sum() < 30:
                        continue
                    frame["entry"] = landmark
                    try:
                        model = CoxPHFitter().fit(
                            frame[["lifespan_days", "dead", "entry", "Tx"]],
                            "lifespan_days", "dead", entry_col="entry",
                        )
                    except Exception:
                        continue
                    result = model.summary.loc["Tx"]
                    rows.append({
                        "cohort": cohort,
                        "group": group,
                        "site": site,
                        "sex": sex,
                        "ctl_bw6": ctl["bw_6"].mean(),
                        "ctl_median": ctl["lifespan_days"].median(),
                        "n_ctl": len(ctl),
                        "n_tx": len(tx),
                        "log_hr": result["coef"],
                        "se": result["se(coef)"],
                        "hr": result["exp(coef)"],
                    })
        return pd.DataFrame(rows).dropna(subset=["ctl_bw6", "ctl_median"])

    warnings.filterwarnings("ignore")
    OUTPUT.mkdir(parents=True, exist_ok=True)
    mice = load_mouse_data()
    tables = {
        "experiment_weight_effects.csv": derive_weight_effects(mice),
        "experiment_interactions.csv": derive_withinexperiment_interactions(mice),
        "comparator_effect_cells.csv": derive_comparator_cells(mice),
    }
    for name, table in tables.items():
        table.to_csv(OUTPUT / name, index=False)
        print(f"{name}: {len(table):,} rows")


def prepare_primary_mice():
    PACKAGE = Path(__file__).absolute().parents[2]

    SUPPLEMENT = PACKAGE / "supplementary" / "Supplementary_Data_1-10.xlsx"

    OUTPUTS = PACKAGE / "data" / "derived"

    DAYS_PER_MONTH = 30.4

    mice = clean_adult_weights(pd.read_excel(SUPPLEMENT, sheet_name="SD2_mouse_level"))
    meta = pd.read_csv(OUTPUTS / "experiment_interactions.csv")
    mice = mice[mice.dead.notna() & mice.lifespan_days.notna()].copy()
    mice["init"] = initiation_months(mice["age_initiation"])

    parts = []
    for row in meta[["cohort", "group", "sex", "weight_age"]].drop_duplicates().itertuples(index=False):
        age = int(row.weight_age)
        weight_col = f"bw_{age}"
        initiation = mice.loc[(mice.cohort == row.cohort) & (mice.group == row.group), "init"].iloc[0]
        entry = max(age, initiation if pd.notna(initiation) else age) * DAYS_PER_MONTH
        z = mice[
            (mice.cohort == row.cohort)
            & (mice.sex == row.sex)
            & mice.group.isin(["Control", row.group])
            & mice[weight_col].notna()
            & (mice.lifespan_days > entry)
        ].copy()
        z["Tx"] = (z.group == row.group).astype(int)
        if z.Tx.nunique() < 2 or len(z) < 60 or z.dead.sum() < 40:
            continue
        z["wz"] = z.groupby(["site", "Tx"])[weight_col].transform(
            lambda x: (x - x.mean()) / x.std(ddof=1)
        )
        z = z.dropna(subset=["wz"])
        z["entry"] = entry
        z["experiment_sex"] = f"{row.cohort}|{row.group}|{row.sex}"
        z["stratum"] = z.site.astype(str) + "|" + z.experiment_sex
        z["mouse"] = z.cohort.astype(str) + "|" + z.id.astype(str)
        z["interaction"] = z.Tx * z.wz
        parts.append(z[["sex", "experiment_sex", "stratum", "mouse", "Tx", "wz", "entry", "lifespan_days", "dead", "interaction"]])

    out = pd.concat(parts, ignore_index=True)
    path = OUTPUTS / "primary_model_mice.csv"
    out.to_csv(path, index=False)
    summary = out.groupby("sex").agg(experiments=("experiment_sex", "nunique"), rows=("experiment_sex", "size"), unique_mice=("mouse", "nunique"))
    expected = {"m": (92, 17943), "f": (88, 15837)}
    for sex, (experiments, unique_mice) in expected.items():
        got = summary.loc[sex]
        if int(got.experiments) != experiments or int(got.unique_mice) != unique_mice:
            raise RuntimeError(f"Primary reconciliation failed for {sex}: {got.to_dict()}")
    print(summary)
    print(path)


def prepare_predosing_mice():
    PACKAGE = Path(__file__).absolute().parents[2]

    SUPPLEMENT = PACKAGE / "supplementary" / "Supplementary_Data_1-10.xlsx"

    INPUTS = PACKAGE / "data" / "inputs"

    OUTPUTS = PACKAGE / "data" / "derived"

    DAYS_PER_MONTH = 30.4375

    def quartile_from_rank(series: pd.Series) -> pd.Series:
        rank = series.rank(method="first", pct=True)
        return np.ceil(rank * 4).clip(1, 4).astype(int)

    mice = clean_adult_weights(pd.read_excel(SUPPLEMENT, sheet_name="SD2_mouse_level"))
    weight_effects = pd.read_csv(OUTPUTS / "experiment_weight_effects.csv")
    experiment_meta = pd.read_csv(OUTPUTS / "experiment_interactions.csv")

    post = weight_effects[
        weight_effects["post_dosing"].astype(str).str.lower().eq("true")
    ].copy()
    post["weighted_diff"] = post["diff_g"] * post["n_tx"]
    classification = (
        post.groupby(["cohort", "group"], as_index=False)
        .agg(weighted_diff=("weighted_diff", "sum"), n_tx=("n_tx", "sum"))
    )
    classification["mean_diff_g"] = (
        classification["weighted_diff"] / classification["n_tx"]
    )
    classification["weight_lowering"] = classification["mean_diff_g"] < -1

    init = post[["cohort", "group", "init"]].drop_duplicates()
    eligible = experiment_meta.merge(init, on=["cohort", "group"], how="left")
    eligible = eligible.merge(
        classification[["cohort", "group", "mean_diff_g", "weight_lowering"]],
        on=["cohort", "group"],
        how="left",
    )
    eligible = eligible[
        (eligible["sex"] == "m")
        & eligible["weight_lowering"].eq(True)
        & (eligible["weight_age"] < eligible["init"])
    ].drop_duplicates(["cohort", "group"])

    parts = []
    for row in eligible.itertuples(index=False):
        weight_col = f"bw_{int(row.weight_age)}"
        entry = float(row.weight_age) * DAYS_PER_MONTH
        z = mice[
            (mice["cohort"] == row.cohort)
            & (mice["sex"] == "m")
            & mice["group"].isin(["Control", row.group])
            & mice[weight_col].notna()
            & (mice["lifespan_days"] > entry)
        ].copy()
        z["experiment"] = f"{row.cohort}|{row.group}"
        z["Tx"] = (z["group"] == row.group).astype(int)
        z["weight_age"] = int(row.weight_age)
        z["entry"] = entry
        z["body_weight"] = pd.to_numeric(z[weight_col], errors="coerce")
        z["mean_treatment_weight_effect_g"] = row.mean_diff_g
        parts.append(z)

    data = pd.concat(parts, ignore_index=True)
    data["mouse"] = data["cohort"].astype(str) + "|" + data["id"].astype(str)
    data["stratum"] = data["experiment"] + "|" + data["site"].astype(str)
    data["matched_cell"] = data["stratum"] + "|Tx" + data["Tx"].astype(str)
    data["matched_quartile"] = data.groupby("matched_cell", group_keys=False)[
        "body_weight"
    ].apply(quartile_from_rank)
    data["weight_z_matched"] = data.groupby("matched_cell")["body_weight"].transform(
        lambda x: (x - x.mean()) / x.std(ddof=1)
    )

    # Equal experiment-site-treatment weights within each quartile.
    n_cell_q = data.groupby(["matched_cell", "matched_quartile"])["mouse"].transform("size")
    data["standardization_weight"] = 1.0 / n_cell_q

    # Reused controls are identifiable and handled by mouse-clustered variance.
    data["comparison_reuse_count"] = data.groupby("mouse")["experiment"].transform("nunique")

    cols = [
        "cohort", "site", "id", "mouse", "experiment", "stratum", "group", "Tx",
        "weight_age", "entry", "body_weight", "weight_z_matched",
        "matched_quartile", "standardization_weight", "comparison_reuse_count",
        "lifespan_days", "dead", "mean_treatment_weight_effect_g",
    ]
    OUTPUTS.mkdir(parents=True, exist_ok=True)
    out = OUTPUTS / "predosing_comparison_mice.csv"
    data[cols].to_csv(out, index=False)
    eligible[["cohort", "group", "weight_age", "init", "mean_diff_g"]].to_csv(
        OUTPUTS / "predosing_eligibility.csv", index=False
    )

    summary = {
        "experiments": int(data["experiment"].nunique()),
        "comparison_rows": int(len(data)),
        "unique_mice": int(data["mouse"].nunique()),
        "treated_mice": int(data.loc[data.Tx.eq(1), "mouse"].nunique()),
        "control_mice": int(data.loc[data.Tx.eq(0), "mouse"].nunique()),
    }
    expected = {"experiments": 30, "comparison_rows": 13044, "unique_mice": 8137}
    for key, value in expected.items():
        if summary[key] != value:
            raise RuntimeError(f"Reconciliation failed for {key}: {summary[key]} != {value}")
    print(summary)
    print(out)


def prepare_nonlowering_mice():
    PACKAGE = Path(__file__).absolute().parents[2]

    SUPPLEMENT = PACKAGE / "supplementary" / "Supplementary_Data_1-10.xlsx"

    INPUTS = PACKAGE / "data" / "inputs"

    OUTPUTS = PACKAGE / "data" / "derived"

    mice = clean_adult_weights(pd.read_excel(SUPPLEMENT, sheet_name="SD2_mouse_level"))
    weight_effects = pd.read_csv(OUTPUTS / "experiment_weight_effects.csv")
    experiment_meta = pd.read_csv(OUTPUTS / "experiment_interactions.csv")

    post = weight_effects[
        weight_effects["post_dosing"].astype(str).str.lower().eq("true")
    ].copy()
    post["weighted_diff"] = post["diff_g"] * post["n_tx"]
    classification = post.groupby(["cohort", "group"], as_index=False).agg(
        weighted_diff=("weighted_diff", "sum"), n_tx=("n_tx", "sum")
    )
    classification["mean_diff_g"] = (
        classification["weighted_diff"] / classification["n_tx"]
    )
    initiation = post[["cohort", "group", "init"]].drop_duplicates()
    eligible = experiment_meta.merge(
        initiation, on=["cohort", "group"], how="left"
    ).merge(
        classification[["cohort", "group", "mean_diff_g"]],
        on=["cohort", "group"], how="left",
    )
    eligible = eligible[
        (eligible["sex"] == "m")
        & (eligible["mean_diff_g"] >= -1)
        & (eligible["weight_age"] < eligible["init"])
    ].drop_duplicates(["cohort", "group"])

    parts = []
    for row in eligible.itertuples(index=False):
        weight_col = f"bw_{int(row.weight_age)}"
        entry = float(row.weight_age) * 30.4375
        subset = mice[
            (mice["cohort"] == row.cohort)
            & (mice["sex"] == "m")
            & mice["group"].isin(["Control", row.group])
            & mice[weight_col].notna()
            & (mice["lifespan_days"] > entry)
        ].copy()
        subset["Tx"] = (subset["group"] == row.group).astype(int)
        subset["entry"] = entry
        subset["body_weight"] = pd.to_numeric(subset[weight_col], errors="coerce")
        parts.append(subset)

    data = pd.concat(parts, ignore_index=True)
    data["mouse"] = data["cohort"].astype(str) + "|" + data["id"].astype(str)
    data = data.drop_duplicates("mouse")
    if len(eligible) != 30 or len(data) != 8252:
        raise RuntimeError(
            f"Expected 30 experiments/8,252 mice; found {len(eligible)}/{len(data)}"
        )

    OUTPUTS.mkdir(parents=True, exist_ok=True)
    output = OUTPUTS / "nonlowering_comparison_mice.csv"
    columns = [
        "cohort", "site", "id", "mouse", "group", "Tx", "entry",
        "body_weight", "lifespan_days", "dead",
    ]
    data[columns].to_csv(output, index=False)
    print({"experiments": len(eligible), "unique_mice": len(data), "output": str(output)})


def prepare_landmark_mice():
    PACKAGE = Path(__file__).absolute().parents[2]

    SUPPLEMENT = PACKAGE / "supplementary" / "Supplementary_Data_1-10.xlsx"

    INPUTS = PACKAGE / "data" / "inputs"

    OUTPUTS = PACKAGE / "data" / "derived"

    DAYS = 30.4375

    mice = clean_adult_weights(pd.read_excel(SUPPLEMENT, sheet_name="SD2_mouse_level"))
    specs = pd.read_csv(INPUTS / "published_lifespan_extenders.csv")
    rows = []
    for spec in specs.itertuples(index=False):
        for landmark in (6, 12, 18, 24):
            if landmark <= spec.start_age:
                continue
            col = f"bw_{landmark}"
            z = mice[
                (mice.cohort == spec.cohort)
                & mice.group.isin(["Control", spec.group])
                & mice[col].notna()
                & (mice.lifespan_days > landmark * DAYS)
            ].copy()
            z["Tx"] = (z.group == spec.group).astype(int)
            z["bw"] = pd.to_numeric(z[col], errors="coerce")
            z = z[z.bw.notna()]
            z["weight_z"] = z.groupby(["site", "sex", "Tx"])["bw"].transform(
                lambda x: (x - x.mean()) / x.std(ddof=1)
            )
            z = z.dropna(subset=["weight_z"])
            z["entry"] = landmark * DAYS
            z["experiment"] = f"{spec.cohort}|{spec.group}"
            z["experiment_stage"] = z.experiment + f"|L{landmark}"
            z["landmark"] = landmark
            z["positive_m"] = bool(spec.positive_m)
            z["positive_f"] = bool(spec.positive_f)
            z["mouse"] = z.cohort.astype(str) + "|" + z.id.astype(str)
            rows.append(z)
    out = pd.concat(rows, ignore_index=True)
    cols = ["cohort", "site", "sex", "id", "mouse", "group", "Tx", "dead", "lifespan_days", "entry", "bw", "weight_z", "experiment", "experiment_stage", "landmark", "positive_m", "positive_f"]
    path = OUTPUTS / "landmark_model_mice.csv"
    out[cols].to_csv(path, index=False)
    print(f"{out.experiment.nunique()} experiments; {out.mouse.nunique()} unique mice; {len(out)} rows")
    print(path)


def prepare_adult_weight_change():
    ROOT = Path(__file__).absolute().parents[2]

    def prepare(data):
        data = clean_adult_weights(data)
        for column in ("lifespan_days", "dead"):
            data[column] = pd.to_numeric(data[column], errors="coerce")
        data["sex"] = data["sex"].astype(str).str[0].str.lower()
        data["is_control"] = data["group"].astype(str).str.casefold().eq("control")
        data["init_months"] = initiation_months(data["age_initiation"])
        data["cell"] = (
            data["cohort"].astype(str) + "|" + data["site"].astype(str)
            + "|" + data["group"].astype(str)
            + "|" + data["dose"].fillna("").astype(str)
            + "|" + data["sex"].astype(str)
        )
        data["pct_6_12"] = 100.0 * (data["bw_12"] - data["bw_6"]) / data["bw_6"]
        columns = ["cohort", "site", "sex", "group", "dose", "cell", "is_control",
                   "init_months", "bw_6", "bw_12", "pct_6_12", "lifespan_days", "dead"]
        return data.loc[
            data["lifespan_days"].gt(365) & data["bw_6"].notna()
            & data["pct_6_12"].notna() & data["dead"].notna(), columns
        ].copy()

    source = pd.read_excel(ROOT / "data/inputs/ITP_bodyweight_all_cohorts.xlsx")

    result = prepare(source)

    result.to_csv(ROOT / "data/derived/adult_weight_change_mice.csv", index=False)

    print(f"Weight-change dataset: {len(result)} mice in {result.cohort.nunique()} cohorts.")


def prepare_restriction_mice():
    ROOT = Path(__file__).absolute().parents[2]

    INPUT = ROOT / "data" / "inputs" / "drido_bodyweights_survival.csv"

    OUTPUT = ROOT / "data" / "derived" / "drido_180_landmark_pooled_cr.csv"

    LANDMARK = 180

    WINDOW = 45

    DIETS = {"AL", "20", "40"}

    def to_long(frame: pd.DataFrame) -> pd.DataFrame:
        parts: list[pd.DataFrame] = []
        for order, weight_column in enumerate(frame.columns):
            if not weight_column.endswith("_BWTest"):
                continue
            age_column = weight_column[:-7] + "_AgeInDays"
            if age_column not in frame.columns:
                continue
            part = frame[["MouseID", "Diet", "SurvDays", "Died", age_column, weight_column]].copy()
            part = part.rename(columns={age_column: "age_days", weight_column: "weight_g"})
            part["measurement_order"] = order
            parts.append(part.dropna(subset=["age_days", "weight_g"]))
        return pd.concat(parts, ignore_index=True)

    measurements = to_long(pd.read_csv(INPUT))
    measurements = measurements[
        measurements["Diet"].astype(str).isin(DIETS)
        & ((measurements["age_days"] - LANDMARK).abs() <= WINDOW)
    ].copy()
    measurements["gap"] = (measurements["age_days"] - LANDMARK).abs()
    landmark = (
        measurements.sort_values(["gap", "measurement_order", "MouseID"])
        .groupby("MouseID", as_index=False)
        .first()
    )
    landmark = landmark[landmark["SurvDays"] > LANDMARK].copy()
    landmark["group"] = np.where(landmark["Diet"].astype(str).eq("AL"), "AL", "CR")

    # Quarters are formed within the ad libitum and pooled caloric-restriction groups.
    landmark["matched_quarter"] = (
        landmark.groupby("group")["weight_g"]
        .transform(lambda values: pd.qcut(values, 4, labels=[1, 2, 3, 4]))
        .astype(int)
    )
    landmark["weight_z"] = landmark.groupby("group")["weight_g"].transform(
        lambda values: (values - values.mean()) / values.std(ddof=1)
    )
    landmark["event"] = landmark["Died"].astype(int)
    landmark["entry"] = LANDMARK
    landmark["treatment"] = (landmark["group"] == "CR").astype(int)

    prepared = landmark.rename(
        columns={"MouseID": "mouse_id", "Diet": "diet", "age_days": "weight_age_days"}
    )[[
        "mouse_id",
        "diet",
        "group",
        "SurvDays",
        "event",
        "entry",
        "weight_age_days",
        "weight_g",
        "weight_z",
        "matched_quarter",
        "treatment",
    ]]
    prepared.to_csv(OUTPUT, index=False, lineterminator="\n")


def prepare_olfactory_mice():
    root = Path(__file__).absolute().parents[2]

    inputs, out = root / "data" / "inputs", root / "data" / "derived"

    survival = pd.read_excel(inputs / "elife-84060-fig1-data1.xlsx", sheet_name="RAM")

    weights = pd.read_excel(inputs / "elife-84060-fig3-data1.xlsx", sheet_name="Covariates")

    assert survival.UMN.is_unique and weights.UMI.is_unique

    assert set(survival.UMN) == set(weights.UMI)

    mice = survival.merge(weights, left_on="UMN", right_on="UMI", how="inner", suffixes=("", "_weights"), validate="one_to_one")

    assert len(mice) == 360 and mice.UMN.is_unique, "every mouse should link exactly once"

    for column in ("Sex", "Rx", "Gene"):
        assert (mice[column] == mice[f"{column}_weights"]).all(), f"{column} disagrees between source files"

    mice = mice.drop(columns=["UMI", "Sex_weights", "Rx_weights", "Gene_weights"])

    mice.to_csv(out / "olfactory_cohort_mice.csv", index=False)

    print(f"olfactory_cohort_mice.csv: {len(mice)} mice; six-month weights {int(mice.W6m.notna().sum())}")


def main():
    root = Path(__file__).absolute().parents[2]
    (root / "data/derived").mkdir(parents=True, exist_ok=True)
    write_weight_qc()
    columns = ["cohort", "site", "sex", "group", "dead", "lifespan_days", "bw_6", "bw_12", "bw_18", "bw_24"]
    pd.read_excel(root / "data/inputs/ITP_bodyweight_all_cohorts.xlsx")[columns].to_csv(root / "data/derived/control_weight_gradient_input.csv", index=False)
    prepare_experiment_contrasts()
    prepare_experiment_models()
    prepare_primary_mice()
    prepare_predosing_mice()
    prepare_nonlowering_mice()
    prepare_landmark_mice()
    prepare_adult_weight_change()
    prepare_restriction_mice()
    prepare_olfactory_mice()

if __name__ == "__main__":
    main()
