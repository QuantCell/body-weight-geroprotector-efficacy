"""Prepare the analysis-ready tables in data/derived from the inputs and the supplement."""
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter
from scipy.stats import ttest_ind

from common import DERIVED, MAX_WEIGHT_G, MIN_WEIGHT_G, OUTPUTS, within_range


def list_adult_weight_exclusions():
    """List the adult weights excluded by the plausibility range."""
    ROOT = Path(__file__).absolute().parents[1]
    SOURCE = ROOT / "data" / "inputs" / "ITP_bodyweight_all_cohorts.xlsx"
    OUTPUT = OUTPUTS / "adult_weight_qc.csv"
    REVIEW_ABOVE_G = 80.0
    WEIGHT_COLUMNS = ("bw_6", "bw_12", "bw_18", "bw_24")

    def main():
        mice = pd.read_excel(SOURCE)
        records = []
        for column in WEIGHT_COLUMNS:
            values = pd.to_numeric(mice[column], errors="coerce")
            notable = values.notna() & ((values < MIN_WEIGHT_G) | (values > REVIEW_ABOVE_G))
            for index in mice.index[notable]:
                value = float(values.loc[index])
                excluded = not MIN_WEIGHT_G <= value <= MAX_WEIGHT_G
                records.append({
                    **mice.loc[index, ["cohort", "site", "id", "sex", "group"]].to_dict(),
                    "weight_age_months": int(column.split("_")[1]),
                    "source_weight_g": value,
                    "analysis_weight_g": np.nan if excluded else value,
                    "excluded": excluded,
                })
        table = pd.DataFrame(records).sort_values(["weight_age_months", "source_weight_g"])
        table.to_csv(OUTPUT, index=False, lineterminator="\n")

        excluded = int(table.excluded.sum())
        print(f"Adult weights outside {MIN_WEIGHT_G:.0f}-{MAX_WEIGHT_G:.0f} g: {excluded} excluded "
              f"of {len(mice) * len(WEIGHT_COLUMNS)} measurements; "
              f"{len(table) - excluded} retained above {REVIEW_ABOVE_G:.0f} g.")
        for row in table[table.excluded].itertuples(index=False):
            print(f"  {row.cohort} {row.site} {row.sex} {row.group} "
                  f"{row.weight_age_months} months: {row.source_weight_g:g} g")
    main()


def prepare_early_weight_change():
    """Build the six-to-twelve-month weight-change table used for the developmental models."""
    ROOT = Path(__file__).absolute().parents[1]
    SOURCE = ROOT / "data" / "inputs" / "ITP_bodyweight_all_cohorts.xlsx"
    OUTPUT = DERIVED / "early_weight_change.csv"
    COLUMNS = ["cohort", "site", "sex", "group", "dose", "cell", "is_control", "init_months",
               "bw_6", "bw_12", "pct_6_12", "lifespan_days", "dead"]

    def prepare(mice):
        mice = mice.copy()
        for column in ("bw_6", "bw_12"):
            mice[column] = within_range(mice[column])
        for column in ("lifespan_days", "dead"):
            mice[column] = pd.to_numeric(mice[column], errors="coerce")
        mice["sex"] = mice["sex"].astype(str).str[0].str.lower()
        mice["is_control"] = mice["group"].astype(str).str.casefold().eq("control")
        mice["init_months"] = pd.to_numeric(
            mice["age_initiation"].astype(str).str.extract(r"(\d+(?:\.\d+)?)")[0], errors="coerce")
        mice["cell"] = (mice["cohort"].astype(str) + "|" + mice["site"].astype(str) + "|"
                        + mice["group"].astype(str) + "|" + mice["dose"].fillna("").astype(str)
                        + "|" + mice["sex"].astype(str))
        mice["pct_6_12"] = 100.0 * (mice["bw_12"] - mice["bw_6"]) / mice["bw_6"]
        eligible = (mice["lifespan_days"].gt(365) & mice["bw_6"].notna()
                    & mice["pct_6_12"].notna() & mice["dead"].notna())
        return mice.loc[eligible, COLUMNS].copy()
    table = prepare(pd.read_excel(SOURCE))
    table.to_csv(OUTPUT, index=False, lineterminator="\n")
    print(f"Weight-change table: {len(table)} mice across {table.cohort.nunique()} cohorts "
              f"and {table.cell.nunique()} cells.")


def prepare_experiment_contrasts():
    """Derive the cohort by sex by experiment contrasts used in the unselected weight-change analysis."""
    ROOT = Path(__file__).absolute().parents[1]
    SOURCE = ROOT / "data" / "inputs" / "ITP_bodyweight_all_cohorts.xlsx"
    WEIGHT = "bw_18"
    MAX_INITIATION_MONTHS = 18
    MIN_ANIMALS = 8
    EXPECTED_EXPERIMENTS = 163

    def load_mice():
        """Mouse-level records with implausible eighteen-month weights set aside."""
        mice = pd.read_excel(SOURCE)
        mice[WEIGHT] = within_range(mice[WEIGHT])
        mice["init"] = pd.to_numeric(
            mice.age_initiation.astype(str).str.extract(r"(\d+(?:\.\d+)?)")[0], errors="coerce")
        return mice

    def experiment_contrasts(mice):
        rows = []
        for (cohort, sex, group), experiment in mice[mice.group.ne("Control")].groupby(
                ["cohort", "sex", "group"]):
            initiation = experiment.init.iloc[0]
            if not np.isfinite(initiation) or initiation >= MAX_INITIATION_MONTHS:
                continue
            controls = mice[mice.cohort.eq(cohort) & mice.sex.eq(sex) & mice.group.eq("Control")]
            treated_weight = experiment[WEIGHT].dropna()
            control_weight = controls[WEIGHT].dropna()
            treated_death = experiment.loc[experiment.dead.eq(1), "lifespan_days"]
            control_death = controls.loc[controls.dead.eq(1), "lifespan_days"]
            if min(len(treated_weight), len(control_weight), len(treated_death),
                   len(control_death)) < MIN_ANIMALS:
                continue
            rows.append({
                "cohort": cohort, "sex": sex, "group": group, "age_init": initiation,
                "n_tw": len(treated_weight), "n_cw": len(control_weight),
                "n_td": len(treated_death), "n_cd": len(control_death),
                "wt_pct": 100 * (treated_weight.mean() - control_weight.mean())
                / control_weight.mean(),
                "life_pct": 100 * (treated_death.median() - control_death.median())
                / control_death.median(),
            })
        return pd.DataFrame(rows).sort_values(["cohort", "sex", "group"], ignore_index=True)

    def main():
        contrasts = experiment_contrasts(load_mice())
        if len(contrasts) != EXPECTED_EXPERIMENTS:
            raise RuntimeError(f"{len(contrasts)} experiment contrasts, expected {EXPECTED_EXPERIMENTS}")
        contrasts.to_csv(DERIVED / "experiment_contrasts.csv", index=False, lineterminator="\n")
        print(f"experiment_contrasts_prepared.csv: {len(contrasts)} experiments in "
              f"{contrasts.cohort.nunique()} cohorts")
    main()


def experiment_level_tables():
    """Derive the experiment-level tables from the harmonized mouse-level workbook."""
    PACKAGE = Path(__file__).absolute().parents[1]
    SOURCE = PACKAGE / "data" / "inputs" / "ITP_bodyweight_all_cohorts.xlsx"
    DAYS_PER_MONTH = 30.4
    WEIGHTS = {6: "bw_6", 12: "bw_12", 18: "bw_18", 24: "bw_24"}

    def standardize(values):
        return (values - values.mean()) / values.std(ddof=1)

    def load_mice():
        mice = pd.read_excel(SOURCE)
        for column in WEIGHTS.values():
            mice[column] = pd.to_numeric(mice[column], errors="coerce")
            mice.loc[~mice[column].between(10, 100), column] = np.nan
        mice["init"] = pd.to_numeric(
            mice.age_initiation.astype(str).str.extract(r"(\d+(?:\.\d+)?)")[0], errors="coerce")
        return mice[mice.dead.notna() & mice.lifespan_days.notna()].copy()

    def weight_effects(mice):
        rows = []
        for (cohort, group), experiment in mice[mice.group != "Control"].groupby(["cohort", "group"]):
            initiation = experiment.init.iloc[0]
            controls = mice[(mice.cohort == cohort) & (mice.group == "Control")]
            for age, column in WEIGHTS.items():
                treated_all, control_all = experiment[column].dropna(), controls[column].dropna()
                if len(treated_all) < 40 or len(control_all) < 40:
                    continue
                differences, sizes = [], []
                for sex in ("m", "f"):
                    treated = experiment.loc[experiment.sex == sex, column].dropna()
                    control = controls.loc[controls.sex == sex, column].dropna()
                    if len(treated) < 20 or len(control) < 20:
                        continue
                    differences.append(treated.mean() - control.mean())
                    sizes.append(len(treated))
                if not differences:
                    continue
                difference = float(np.average(differences, weights=sizes))
                rows.append({"cohort": cohort, "group": group, "init": initiation, "weight_age": age,
                             "post_dosing": bool(not np.isnan(initiation) and age >= initiation),
                             "n_tx": len(treated_all), "diff_g": difference,
                             "pct": 100 * difference / control_all.mean(),
                             "p": ttest_ind(treated_all, control_all, equal_var=False).pvalue})
        return pd.DataFrame(rows)

    def experiment_interaction(mice, cohort, group, age, entry_months):
        column = WEIGHTS[age]
        landmark = entry_months * DAYS_PER_MONTH
        experiment = mice[(mice.cohort == cohort) & mice.group.isin([group, "Control"])].copy()
        experiment = experiment[experiment[column].notna() & (experiment.lifespan_days > landmark)]
        experiment["Tx"] = (experiment.group == group).astype(int)
        rows = []
        for sex in ("m", "f"):
            frame = experiment[experiment.sex == sex].copy()
            if frame.Tx.nunique() < 2 or frame.dead.sum() < 40:
                continue
            frame["wz"] = frame.groupby(["site", "Tx"])[column].transform(standardize)
            frame = frame.dropna(subset=["wz"])
            frame["entry"] = landmark
            frame["Tx_x_wz"] = frame.Tx * frame.wz
            columns = ["lifespan_days", "dead", "entry", "Tx", "wz", "Tx_x_wz", "site"]
            try:
                model = CoxPHFitter().fit(frame[columns], "lifespan_days", "dead",
                                          entry_col="entry", strata=["site"])
            except Exception:
                continue
            result = model.summary.loc["Tx_x_wz"]
            rows.append({"cohort": cohort, "group": group, "sex": sex, "weight_age": age,
                         "entry_months": entry_months, "within_experiment": True, "n": len(frame), "int_hr": result["exp(coef)"],
                         "int_lo": result["exp(coef) lower 95%"], "int_hi": result["exp(coef) upper 95%"],
                         "int_p": result["p"], "coef": result["coef"], "se": result["se(coef)"]})
        return rows

    def within_experiment_interactions(mice):
        rows = []
        experiments = mice[mice.group != "Control"][["cohort", "group"]].drop_duplicates()
        for experiment in experiments.itertuples(index=False):
            initiation = mice.loc[(mice.cohort == experiment.cohort) & (mice.group == experiment.group), "init"].iloc[0]
            entry_months = 6.0 if np.isnan(initiation) else max(6.0, float(initiation))
            rows.extend(experiment_interaction(mice, experiment.cohort, experiment.group, 6, entry_months))
        return pd.DataFrame(rows)

    def comparator_cells(mice):
        rows = []
        landmark = 6 * DAYS_PER_MONTH
        for (cohort, group), experiment in mice[mice.group != "Control"].groupby(["cohort", "group"]):
            for site in experiment.site.dropna().unique():
                for sex in ("m", "f"):
                    treated = experiment[(experiment.site == site) & (experiment.sex == sex)]
                    control = mice[(mice.cohort == cohort) & (mice.group == "Control")
                                   & (mice.site == site) & (mice.sex == sex)]
                    if len(treated) < 25 or len(control) < 25:
                        continue
                    frame = pd.concat([treated.assign(Tx=1), control.assign(Tx=0)])
                    frame = frame[frame.lifespan_days > landmark].copy()
                    if frame.dead.sum() < 30:
                        continue
                    frame["entry"] = landmark
                    try:
                        model = CoxPHFitter().fit(frame[["lifespan_days", "dead", "entry", "Tx"]],
                                                  "lifespan_days", "dead", entry_col="entry")
                    except Exception:
                        continue
                    result = model.summary.loc["Tx"]
                    rows.append({"cohort": cohort, "group": group, "site": site, "sex": sex,
                                 "ctl_bw6": control.bw_6.mean(),
                                 "ctl_median": control.lifespan_days.median(),
                                 "n_ctl": len(control), "n_tx": len(treated),
                                 "log_hr": result["coef"], "se": result["se(coef)"],
                                 "hr": result["exp(coef)"]})
        return pd.DataFrame(rows).dropna(subset=["ctl_bw6", "ctl_median"])

    def main():
        warnings.filterwarnings("ignore")
        mice = load_mice()
        tables = {"experiment_weight_effects.csv": weight_effects(mice),
                  "experiment_interactions.csv": within_experiment_interactions(mice),
                  "comparator_cells.csv": comparator_cells(mice)}
        expected = {"experiment_weight_effects.csv": 368, "experiment_interactions.csv": 180,
                    "comparator_cells.csv": 534}
        for name, table in tables.items():
            if len(table) != expected[name]:
                raise RuntimeError(f"{name}: {len(table)} rows, expected {expected[name]}")
            table.to_csv(DERIVED / name, index=False)
            print(f"{name}: {len(table):,} rows")
    main()


def prepare_within_experiment_models():
    """Prepare the primary all-experiment treatment-by-weight dataset from the mouse-level file."""
    PACKAGE = Path(__file__).absolute().parents[1]
    SUPPLEMENT = PACKAGE / "supplementary" / "Supplementary_Data.xlsx"
    DAYS_PER_MONTH = 30.4

    def main():
        mice = pd.read_excel(SUPPLEMENT, sheet_name="SD2_mouse_level")
        meta = pd.read_csv(DERIVED / "experiment_interactions.csv")
        for age in (6, 12, 18, 24):
            col = f"bw_{age}"
            mice[col] = pd.to_numeric(mice[col], errors="coerce")
            mice.loc[~mice[col].between(10, 100), col] = np.nan
        mice = mice[mice.dead.notna() & mice.lifespan_days.notna()].copy()

        parts = []
        for row in meta[["cohort", "group", "sex", "weight_age", "entry_months"]].drop_duplicates().itertuples(index=False):
            age = int(row.weight_age)
            weight_col = f"bw_{age}"
            entry = row.entry_months * DAYS_PER_MONTH
            z = mice[
                (mice.cohort == row.cohort)
                & (mice.sex == row.sex)
                & mice.group.isin(["Control", row.group])
                & mice[weight_col].notna()
                & (mice.lifespan_days > entry)
            ].copy()
            z["Tx"] = (z.group == row.group).astype(int)
            if z.Tx.nunique() < 2 or z.dead.sum() < 40:
                continue
            z["wz"] = z.groupby(["site", "Tx"])[weight_col].transform(
                lambda x: (x - x.mean()) / x.std(ddof=1)
            )
            z = z.dropna(subset=["wz"])
            z["entry"] = entry
            z["experimentsex"] = f"{row.cohort}|{row.group}|{row.sex}"
            z["stratum"] = z.site.astype(str) + "|" + z.experimentsex
            z["mouse"] = z.cohort.astype(str) + "|" + z.id.astype(str)
            z["interaction"] = z.Tx * z.wz
            parts.append(z[["sex", "experimentsex", "stratum", "mouse", "Tx", "wz", "entry", "lifespan_days", "dead", "interaction"]])

        out = pd.concat(parts, ignore_index=True)
        path = DERIVED / "within_experiment_long.csv"
        out.to_csv(path, index=False)
        summary = out.groupby("sex").agg(experiments=("experimentsex", "nunique"), rows=("experimentsex", "size"), unique_mice=("mouse", "nunique"))
        expected = {"m": (92, 17943), "f": (88, 15837)}
        for sex, (experiments, unique_mice) in expected.items():
            got = summary.loc[sex]
            if int(got.experiments) != experiments or int(got.unique_mice) != unique_mice:
                raise RuntimeError(f"Primary reconciliation failed for {sex}: {got.to_dict()}")
        print(summary)
        print(path)
    main()


def prepare_matched_quartiles():
    """Build the pre-dosing weight-lowering comparison dataset."""
    PACKAGE = Path(__file__).absolute().parents[1]
    SUPPLEMENT = PACKAGE / "supplementary" / "Supplementary_Data.xlsx"
    DAYS_PER_MONTH = 30.4375

    def quartile_from_rank(series: pd.Series) -> pd.Series:
        rank = series.rank(method="first", pct=True)
        return np.ceil(rank * 4).clip(1, 4).astype(int)

    def main():
        mice = pd.read_excel(SUPPLEMENT, sheet_name="SD2_mouse_level")
        weight_effects = pd.read_csv(DERIVED / "experiment_weight_effects.csv")
        experiment_meta = pd.read_csv(DERIVED / "experiment_interactions.csv")

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
                & pd.to_numeric(mice[weight_col], errors="coerce").between(10, 100)
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

        n_cell_q = data.groupby(["matched_cell", "matched_quartile"])["mouse"].transform("size")
        data["standardization_weight"] = 1.0 / n_cell_q

        data["comparison_reuse_count"] = data.groupby("mouse")["experiment"].transform("nunique")

        cols = [
            "cohort", "site", "id", "mouse", "experiment", "stratum", "group", "Tx",
            "weight_age", "entry", "body_weight", "weight_z_matched",
            "matched_quartile", "standardization_weight", "comparison_reuse_count",
            "lifespan_days", "dead", "mean_treatment_weight_effect_g",
        ]
        DERIVED.mkdir(parents=True, exist_ok=True)
        out = DERIVED / "quartile_comparisons.csv"
        data[cols].to_csv(out, index=False)
        eligible[["cohort", "group", "weight_age", "init", "mean_diff_g"]].to_csv(
            DERIVED / "quartile_eligibility.csv", index=False
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
    main()


def prepare_nonlowering_quartiles():
    """Build the unique-mouse 29-experiment pre-dosing non-lowering sensitivity set."""
    PACKAGE = Path(__file__).absolute().parents[1]
    SUPPLEMENT = PACKAGE / "supplementary" / "Supplementary_Data.xlsx"

    def main():
        mice = pd.read_excel(SUPPLEMENT, sheet_name="SD2_mouse_level")
        weight_effects = pd.read_csv(DERIVED / "experiment_weight_effects.csv")
        experiment_meta = pd.read_csv(DERIVED / "experiment_interactions.csv")
        post = weight_effects[
            weight_effects["post_dosing"].astype(str).str.lower().eq("true")
        ].copy()
        post["weighted_diff"] = post["diff_g"] * post["n_tx"]
        classification = post.groupby(["cohort", "group"], as_index=False).agg(
            weighted_diff=("weighted_diff", "sum"), n_tx=("n_tx", "sum")
        )
        classification["mean_diff_g"] = classification["weighted_diff"] / classification["n_tx"]
        init = post[["cohort", "group", "init"]].drop_duplicates()
        eligible = experiment_meta.merge(init, on=["cohort", "group"], how="left").merge(
            classification[["cohort", "group", "mean_diff_g"]],
            on=["cohort", "group"], how="left"
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
            z = mice[
                (mice["cohort"] == row.cohort)
                & (mice["sex"] == "m")
                & mice["group"].isin(["Control", row.group])
                & pd.to_numeric(mice[weight_col], errors="coerce").between(10, 100)
                & (mice["lifespan_days"] > entry)
            ].copy()
            z["Tx"] = (z["group"] == row.group).astype(int)
            z["entry"] = entry
            z["body_weight"] = pd.to_numeric(z[weight_col], errors="coerce")
            parts.append(z)

        data = pd.concat(parts, ignore_index=True)
        data["mouse"] = data["cohort"].astype(str) + "|" + data["id"].astype(str)
        data = data.drop_duplicates("mouse")
        if len(eligible) != 30 or len(data) != 8252:
            raise RuntimeError(f"Expected 30 experiments/8,252 mice; found {len(eligible)}/{len(data)}")
        DERIVED.mkdir(parents=True, exist_ok=True)
        out = DERIVED / "quartile_nonlowering.csv"
        cols = ["cohort", "site", "id", "mouse", "group", "Tx", "entry",
                "body_weight", "lifespan_days", "dead"]
        data[cols].to_csv(out, index=False)
        print({"experiments": len(eligible), "unique_mice": len(data), "output": str(out)})
    main()


def prepare_landmark_models():
    """Prepare 6/12/18/24-month landmarks for published lifespan extenders."""
    PACKAGE = Path(__file__).absolute().parents[1]
    SUPPLEMENT = PACKAGE / "supplementary" / "Supplementary_Data.xlsx"
    INPUTS = PACKAGE / "data" / "inputs"
    DAYS = 30.4375

    def main():
        mice = pd.read_excel(SUPPLEMENT, sheet_name="SD2_mouse_level")
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
                z = z[z.bw.between(10, 100)]
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
        path = DERIVED / "landmark_long.csv"
        out[cols].to_csv(path, index=False)
        print(f"{out.experiment.nunique()} experiments; {out.mouse.nunique()} unique mice; {len(out)} rows")
        print(path)
    main()


def prepare_drido_landmark():
    """Reshape the Diversity Outbred release and select measurements near the 180-day landmark."""
    ROOT = Path(__file__).absolute().parents[1]
    INPUTS = ROOT / "data" / "inputs"
    LANDMARK = 180
    WINDOW = 45
    DIETS = ["AL", "20", "40"]

    source = pd.read_csv(INPUTS / "drido_bodyweights_survival.csv")
    measurements = []
    for column in source.columns:
        if not column.endswith("_BWTest"):
            continue
        age_column = column[:-len("_BWTest")] + "_AgeInDays"
        if age_column not in source.columns:
            continue
        visit = source[["MouseID", "Diet", "SurvDays", "Died", age_column, column]].rename(
            columns={age_column: "age_days", column: "weight_g"})
        measurements.append(visit.dropna(subset=["age_days", "weight_g"]))

    long = pd.concat(measurements, ignore_index=True)
    long = long[long.Diet.astype(str).isin(DIETS)]
    long = long[(long.age_days - LANDMARK).abs() <= WINDOW]
    long = long.rename(columns={"MouseID": "mouse_id", "Diet": "diet",
                                "SurvDays": "survival_days", "Died": "died"})
    long["died"] = long.died.astype(int)
    long["age_days"] = long.age_days.astype(int)
    long = (long[["mouse_id", "diet", "survival_days", "died", "age_days", "weight_g"]]
            .drop_duplicates()
            .sort_values(["mouse_id", "age_days"])
            .reset_index(drop=True))
    long.to_csv(DERIVED / "drido_landmark.csv", index=False, lineterminator="\n")

    print(f"Diversity Outbred release: {source.MouseID.nunique()} mice, "
          f"{sum(c.endswith('_BWTest') for c in source.columns)} weighing visits.")
    print(f"  within {WINDOW} days of the {LANDMARK}-day landmark and on an eligible diet: "
          f"{len(long)} measurements on {long.mouse_id.nunique()} mice "
          f"({', '.join(sorted(long.diet.astype(str).unique()))}).")


def prepare_olfactory_cohort():
    """Link the external olfactory-exposure lifespan cohort (Garratt et al., eLife 2022) by mouse identifier."""
    root = Path(__file__).absolute().parents[1]
    inputs = root / "data" / "inputs"

    survival = pd.read_excel(inputs / "elife-84060-fig1-data1.xlsx", sheet_name="RAM")
    weights = pd.read_excel(inputs / "elife-84060-fig3-data1.xlsx", sheet_name="Covariates")
    mice = survival.merge(weights, left_on="UMN", right_on="UMI", how="inner", suffixes=("", "_weights"))
    assert len(mice) == 360 and mice.UMN.is_unique, "every mouse should link exactly once"
    for column in ("Sex", "Rx", "Gene"):
        assert (mice[column] == mice[f"{column}_weights"]).all(), f"{column} disagrees between source files"
    mice = mice.drop(columns=["UMI", "Sex_weights", "Rx_weights", "Gene_weights"])
    mice.to_csv(DERIVED / "olfactory_mice.csv", index=False)
    print(f"olfactory_cohort_mice.csv: {len(mice)} mice; six-month weights {int(mice.W6m.notna().sum())}")


def main():
    DERIVED.mkdir(parents=True, exist_ok=True)
    list_adult_weight_exclusions()
    prepare_early_weight_change()
    prepare_experiment_contrasts()
    experiment_level_tables()
    prepare_within_experiment_models()
    prepare_matched_quartiles()
    prepare_nonlowering_quartiles()
    prepare_landmark_models()
    prepare_drido_landmark()
    prepare_olfactory_cohort()


if __name__ == "__main__":
    main()
