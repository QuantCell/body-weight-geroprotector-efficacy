"""Survival summaries."""
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter, KaplanMeierFitter

from common import DERIVED, write_results

TABLES = {}


def comparator_state():
    """Descriptive comparator-state tertiles for Table 3."""
    ROOT = Path(__file__).absolute().parents[1]
    LIFE = ROOT / "data/inputs/ITP_bodyweight_all_cohorts.xlsx"
    CELLS = DERIVED / "comparator_cells.csv"
    Z95 = 1.959963984540054
    LANDMARK_DAYS = 6 * 30.4

    def km_median(frame):
        """Six-month-left-truncated Kaplan-Meier median (lifelines convention)."""
        z = frame.loc[frame["lifespan_days"].gt(LANDMARK_DAYS),
                      ["lifespan_days", "dead"]].dropna()
        if z.empty:
            return np.nan
        km = KaplanMeierFitter()
        km.fit(z["lifespan_days"], z["dead"],
               entry=np.full(len(z), LANDMARK_DAYS))
        return float(km.median_survival_time_)

    life = pd.read_excel(LIFE)
    for column in ("lifespan_days", "dead"):
        life[column] = pd.to_numeric(life[column], errors="coerce")
    life = life[life.lifespan_days.notna() & life.dead.notna()].copy()
    cells = pd.read_csv(CELLS)
    if len(cells) != 534:
        raise AssertionError(f"Expected 534 comparator cells, found {len(cells)}")

    medians = []
    for row in cells.itertuples(index=False):
        treated = life[(life.cohort == row.cohort) & (life.group == row.group) &
                       (life.site == row.site) & (life.sex == row.sex)]
        control = life[(life.cohort == row.cohort) & (life.group == "Control") &
                       (life.site == row.site) & (life.sex == row.sex)]
        cmed, tmed = km_median(control), km_median(treated)
        medians.append((cmed, tmed, 100 * (tmed - cmed) / cmed))

    cells[["control_km_median", "treated_km_median", "median_gain_pct"]] = pd.DataFrame(
        medians, index=cells.index
    )
    cells["control_cluster"] = cells.cohort.astype(str) + "|" + cells.site.astype(str)
    cells["sex_specific_control"] = cells.control_cluster + "|" + cells.sex.astype(str)
    cells.to_csv(DERIVED / "comparator_cell_medians.csv", index=False)

    rows = []
    for sex in ("m", "f"):
        subset = cells[cells.sex == sex].copy()
        for split in (3, 4):
            subset["bin"] = pd.qcut(
                subset.ctl_bw6, split, labels=[f"Q{i}" for i in range(1, split + 1)]
            )
            for label, group in subset.groupby("bin", observed=True):
                weights = 1 / group.se.pow(2)
                beta = np.average(group.log_hr, weights=weights)
                se = np.sqrt(1 / weights.sum())
                rows.append({
                    "sex": sex,
                    "split": split,
                    "bin": str(label),
                    "cells": len(group),
                    "mean_control_weight_g": group.ctl_bw6.mean(),
                    "pooled_hr": np.exp(beta),
                    "lower_95_ci": np.exp(beta - Z95 * se),
                    "upper_95_ci": np.exp(beta + Z95 * se),
                    "hazard_reduction_pct": 100 * (1 - np.exp(beta)),
                    "mean_cell_km_median_gain_pct": group.median_gain_pct.mean(),
                })

    summary = pd.DataFrame(rows)
    TABLES["comparator_strata"] = summary

    print(f"{len(cells)} cells; {cells.control_cluster.nunique()} cohort-site control clusters; "
          f"{cells.sex_specific_control.nunique()} sex-specific control means")
    print(summary.to_string(index=False))


def stratum_consistency():
    """Stratum-level consistency of the six-month weight-lifespan relation."""
    PACKAGE = Path(__file__).absolute().parents[1]
    INPUT = PACKAGE / "supplementary" / "Supplementary_Data.xlsx"

    def main() -> None:
        data = pd.read_excel(INPUT, sheet_name="SD2_mouse_level")
        data["bw_6"] = pd.to_numeric(data["bw_6"], errors="coerce")
        data.loc[~data["bw_6"].between(10, 100), "bw_6"] = np.nan
        data = data[
            data["dead"].eq(1)
            & data["bw_6"].notna()
            & data["lifespan_days"].notna()
        ].copy()

        rows = []
        for (cohort, sex, group), stratum in data.groupby(
            ["cohort", "sex", "group"], sort=True
        ):
            site_estimates = []
            for _, cell in stratum.groupby("site", sort=True):
                if len(cell) < 12:
                    continue
                rho = cell[["bw_6", "lifespan_days"]].corr(method="spearman").iloc[0, 1]
                if not np.isfinite(rho):
                    continue
                site_estimates.append(
                    (np.arctanh(np.clip(rho, -0.999999, 0.999999)), len(cell) - 3, len(cell))
                )
            if not site_estimates:
                continue
            z = np.asarray([value[0] for value in site_estimates])
            weights = np.asarray([value[1] for value in site_estimates], dtype=float)
            rows.append(
                {
                    "cohort": cohort,
                    "sex": sex,
                    "group": group,
                    "nsite": len(site_estimates),
                    "n": int(sum(value[2] for value in site_estimates)),
                    "rho": float(np.tanh(np.average(z, weights=weights))),
                }
            )

        result = pd.DataFrame(rows)
        result.to_csv(DERIVED / "bw6_strata.csv", index=False)

        expected = {"m": (109, 108), "f": (105, 87)}
        for sex, (n_expected, negative_expected) in expected.items():
            subset = result[result["sex"].eq(sex)]
            observed = (len(subset), int(subset["rho"].lt(0).sum()))
            if observed != (n_expected, negative_expected):
                raise RuntimeError(
                    f"Reconciliation failed for sex={sex}: {observed} != "
                    f"{(n_expected, negative_expected)}"
                )
            print(
                f"{sex}: {observed[1]}/{observed[0]} negative; "
                f"median rho={subset['rho'].median():.3f}"
            )
    main()


def external_do_cohort():
    """The external Diversity Outbred analysis reported in Extended Data Table 4."""
    LANDMARK = 180
    WINDOW = 45

    measurements = pd.read_csv(DERIVED / "drido_landmark.csv")
    measurements = measurements[(measurements.age_days - LANDMARK).abs() <= WINDOW].copy()
    measurements["gap"] = (measurements.age_days - LANDMARK).abs()
    at_landmark = (measurements.sort_values(["gap", "mouse_id"])
                   .groupby("mouse_id", as_index=False).first())
    at_landmark = at_landmark[at_landmark.survival_days > LANDMARK].copy()
    at_landmark["experiment"] = np.where(at_landmark.diet.astype(str).eq("AL"), "AL", "CR")
    at_landmark["restriction_40"] = at_landmark.diet.astype(str).eq("40").astype(int)
    at_landmark["quarter"] = (at_landmark.groupby("experiment").weight_g
                              .transform(lambda s: pd.qcut(s, 4, labels=[1, 2, 3, 4])).astype(int))
    at_landmark["entry"] = LANDMARK

    def km_median(frame):
        fitter = KaplanMeierFitter().fit(frame.survival_days, frame.died,
                                         entry=np.full(len(frame), LANDMARK))
        return float(fitter.median_survival_time_)

    quarters = []
    for quarter in (1, 2, 3, 4):
        block = at_landmark[at_landmark.quarter.eq(quarter)]
        al = block[block.experiment.eq("AL")]
        cr = block[block.experiment.eq("CR")]
        al_median, cr_median = km_median(al), km_median(cr)
        quarters.append({"quarter": quarter,
                         "mean_weight_g": block.weight_g.mean(),
                         "al_n": len(al), "al_km_median_days": al_median,
                         "cr_n": len(cr), "cr_km_median_days": cr_median,
                         "gain_days": cr_median - al_median,
                         "gain_pct": 100 * (cr_median - al_median) / al_median})
    TABLES["restriction_quartiles"] = pd.DataFrame(quarters)

    def weight_gradient(frame, adjust_for_level):
        frame = frame.copy()
        frame["weight_z"] = ((frame.weight_g - frame.weight_g.mean())
                             / frame.weight_g.std(ddof=1))
        columns = ["survival_days", "died", "entry", "weight_z"]
        if adjust_for_level:
            columns.append("restriction_40")
        fit = CoxPHFitter().fit(frame[columns], "survival_days", "died", entry_col="entry")
        return fit.summary.loc["weight_z"], len(frame)

    rows = []
    for label, frame, adjust in [("ad_libitum_weight_gradient", at_landmark[at_landmark.experiment.eq("AL")],
                                  False),
                                 ("pooled_restriction_weight_gradient",
                                  at_landmark[at_landmark.experiment.eq("CR")], True)]:
        summary, n = weight_gradient(frame, adjust)
        rows.append({"analysis": label, "n": n, "estimate": summary["exp(coef)"],
                     "lo": summary["exp(coef) lower 95%"], "hi": summary["exp(coef) upper 95%"],
                     "p": summary["p"]})

    trend = at_landmark.copy()
    trend["restricted"] = trend.experiment.eq("CR").astype(int)
    trend["restricted_by_quarter"] = trend.restricted * trend.quarter
    fit = CoxPHFitter().fit(
        trend[["survival_days", "died", "entry", "restricted", "quarter",
               "restricted_by_quarter", "restriction_40"]],
        "survival_days", "died", entry_col="entry")
    summary = fit.summary.loc["restricted_by_quarter"]
    rows.append({"analysis": "treatment_by_quarter_trend", "n": len(trend),
                 "estimate": summary["exp(coef)"], "lo": summary["exp(coef) lower 95%"],
                 "hi": summary["exp(coef) upper 95%"], "p": summary["p"]})

    ad_libitum = at_landmark[at_landmark.experiment.eq("AL")]
    fitter = KaplanMeierFitter().fit(ad_libitum.survival_days, ad_libitum.died)
    rows.append({"analysis": "ad_libitum_day500_mortality_pct", "n": len(ad_libitum),
                 "estimate": 100 * (1 - float(fitter.predict(500))),
                 "lo": np.nan, "hi": np.nan, "p": np.nan})
    all_events = KaplanMeierFitter().fit(ad_libitum.survival_days, np.ones(len(ad_libitum)))
    rows.append({"analysis": "ad_libitum_day500_mortality_all_events_pct", "n": len(ad_libitum),
                 "estimate": 100 * (1 - float(all_events.predict(500))),
                 "lo": np.nan, "hi": np.nan, "p": np.nan})

    TABLES["restriction_models"] = pd.DataFrame(rows)

    print(f"Diversity Outbred cohort: {len(at_landmark)} females at the {LANDMARK}-day landmark "
          f"({len(ad_libitum)} ad libitum).")
    for row in rows:
        if np.isnan(row["lo"]):
            print(f"  {row['analysis']}: {row['estimate']:.1f}% (n = {row['n']})")
        else:
            print(f"  {row['analysis']}: {row['estimate']:.2f} "
                  f"({row['lo']:.2f}-{row['hi']:.2f}); P = {row['p']:.3f}")


def main():
    comparator_state()
    stratum_consistency()
    external_do_cohort()
    write_results(TABLES)


if __name__ == "__main__":
    main()
