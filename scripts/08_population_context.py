"""Population context."""
import csv
import math
from pathlib import Path

import numpy as np
import pandas as pd
from lifelines import KaplanMeierFitter
from scipy import stats

from common import DERIVED, read_result, within_range, write_results

TABLES = {}


def sex_survival_context():
    """Extended Data Table 6."""
    ROOT = Path(__file__).absolute().parents[1]
    INPUTS = ROOT / "data" / "inputs"
    SUPPLEMENTARY = ROOT / "supplementary" / "Supplementary_Data.xlsx"
    CHENG_COHORTS = ["C2004", "C2005", "C2006", "C2007", "C2009", "C2010"]

    mice = pd.read_excel(SUPPLEMENTARY, sheet_name="SD2_mouse_level")
    controls = mice[mice.group.eq("Control") & mice.lifespan_days.notna()].copy()
    controls["event"] = controls.dead.astype(int)

    def km(frame, event):
        return KaplanMeierFitter().fit(frame.lifespan_days, event)

    def summarise(frame, label):
        rows = []
        for sex, name in (("m", "male"), ("f", "female")):
            block = frame[frame.sex.eq(sex)]
            censored = km(block, block.event)
            all_events = km(block, pd.Series(1, index=block.index))
            rows.append({"cohort_set": label, "sex": name, "n": len(block),
                         "cohorts": block.cohort.nunique(),
                         "km_median_days": float(censored.median_survival_time_),
                         "day500_mortality_pct": 100 * (1 - float(censored.predict(500))),
                         "day500_mortality_all_events_pct":
                             100 * (1 - float(all_events.predict(500)))})
        female, male = rows[1], rows[0]
        rows.append({"cohort_set": label, "sex": "female_minus_male", "n": None, "cohorts": None,
                     "km_median_days": female["km_median_days"] - male["km_median_days"],
                     "day500_mortality_pct":
                         female["day500_mortality_pct"] - male["day500_mortality_pct"],
                     "day500_mortality_all_events_pct": None})
        return rows

    summary = pd.DataFrame(summarise(controls, "all_17_cohorts")
                           + summarise(controls[controls.cohort.isin(CHENG_COHORTS)],
                                       "published_six_cohorts"))
    TABLES["sex_survival_summary"] = summary

    by_cohort = []
    for cohort, block in controls.groupby("cohort"):
        medians = {}
        for sex in ("m", "f"):
            sub = block[block.sex.eq(sex)]
            if len(sub):
                medians[sex] = float(km(sub, sub.event).median_survival_time_)
        if len(medians) == 2:
            by_cohort.append({"cohort": cohort, "n": len(block),
                              "male_km_median_days": medians["m"],
                              "female_km_median_days": medians["f"],
                              "female_advantage_days": medians["f"] - medians["m"]})
    cohort_table = pd.DataFrame(by_cohort).sort_values("cohort")
    TABLES["sex_survival_by_cohort"] = cohort_table

    external = pd.read_csv(INPUTS / "external_sex_survival_context.csv")
    external["female_minus_male_days"] = (external.female_survival_days
                                          - external.male_survival_days)
    TABLES["external_survival_context"] = external

    indexed = summary.set_index(["cohort_set", "sex"])
    print("ITP untreated controls:")
    for label in ("all_17_cohorts", "published_six_cohorts"):
        male, female = indexed.loc[(label, "male")], indexed.loc[(label, "female")]
        print(f"  {label}: male {male.km_median_days:.0f} d (n = {male.n:.0f}), "
              f"female {female.km_median_days:.0f} d (n = {female.n:.0f}), "
              f"advantage {indexed.loc[(label, 'female_minus_male'), 'km_median_days']:.0f} d")
    male = indexed.loc[("all_17_cohorts", "male")]
    female = indexed.loc[("all_17_cohorts", "female")]
    print(f"  day-500 mortality: male {male.day500_mortality_pct:.1f}%, "
          f"female {female.day500_mortality_pct:.1f}% "
          f"({female.day500_mortality_all_events_pct:.1f}% counting every animal as an event)")
    print(f"  female advantage by cohort: {cohort_table.female_advantage_days.min():.0f} to "
          f"{cohort_table.female_advantage_days.max():.0f} d across {len(cohort_table)} cohorts")
    for row in external.itertuples(index=False):
        print(f"  {row.population} ({row.study}): male {row.male_survival_days} d, "
              f"female {row.female_survival_days} d, "
              f"difference {row.female_minus_male_days:+.0f} d [{row.metric}]")


def control_early_mortality():
    """Day-500 mortality by six-month weight quartile in all male controls, and when the quartile lifespan gap arises."""
    ROOT = Path(__file__).absolute().parents[1]
    SIX_MONTHS, LANDMARK = 182.625, 500

    mice = pd.read_excel(ROOT / "supplementary" / "Supplementary_Data.xlsx", sheet_name="SD2_mouse_level")
    mice = mice[mice.dead.notna() & mice.lifespan_days.notna()].copy()
    mice["bw_6"] = within_range(mice.bw_6)
    males = mice[mice.group.eq("Control") & mice.sex.eq("m") & (mice.lifespan_days > SIX_MONTHS)
                 & mice.bw_6.notna()].copy()
    males["quartile"] = males.groupby(["cohort", "site"]).bw_6.transform(
        lambda s: np.ceil(s.rank(method="first", pct=True) * 4).clip(1, 4)).astype(int)
    days = np.arange(183, int(mice.lifespan_days.max()) + 2)
    rows, curves = [{"quantity": "male_controls", "value": len(males)}], {}
    for quartile, block in males.groupby("quartile"):
        fit = KaplanMeierFitter().fit(block.lifespan_days, block.dead, entry=np.full(len(block), SIX_MONTHS))
        curves[quartile] = fit.survival_function_at_times(days).values
        rows += [{"quantity": f"q{quartile}_mice", "value": len(block)},
                 {"quantity": f"q{quartile}_dead_by_day_500_pct", "value": 100 * (1 - float(fit.predict(LANDMARK)))}]

    light, heavy = curves[1], curves[4]
    k = int(np.flatnonzero(days == LANDMARK)[0])

    def remaining(curve):
        return curve[k:].sum() / curve[k]

    difference = light.sum() - heavy.sum()
    # exchange survival before day 500 and conditional survival after it, in both orders, and average
    early = ((light[:k].sum() + light[k] * remaining(heavy) - heavy.sum())
             + (light.sum() - heavy[:k].sum() - heavy[k] * remaining(light))) / 2
    rows += [{"quantity": "mean_lifespan_difference_days", "value": difference},
             {"quantity": "difference_from_mortality_before_day_500_days", "value": early},
             {"quantity": "difference_from_mortality_before_day_500_pct", "value": 100 * early / difference}]
    table = pd.DataFrame(rows)
    TABLES["control_early_mortality"] = table
    values = table.set_index("quantity").value
    print(f"Male controls ({len(males)}): dead by day 500 in the lightest quartile {values['q1_dead_by_day_500_pct']:.1f}%, "
          f"heaviest {values['q4_dead_by_day_500_pct']:.1f}%; {early:.0f} of the {difference:.0f}-day mean lifespan "
          f"difference arose from mortality before day 500")


def inbred_strain_comparison():
    """Compare the ITP female survival advantage with the sex difference across inbred strains."""
    ROOT = Path(__file__).absolute().parents[1]
    INPUTS = ROOT / "data" / "inputs"

    strains = pd.read_csv(INPUTS / "inbred_strain_median_lifespan.csv")
    itp = TABLES["sex_survival_summary"]
    itp = itp[itp.cohort_set.eq("all_17_cohorts")].set_index("sex").km_median_days
    difference = strains.female_median_days - strains.male_median_days
    summary = pd.DataFrame([{
        "strains": len(strains),
        "strain_median_female_minus_male_days": difference.median(),
        "strains_female_longer": int((difference > 0).sum()),
        "itp_female_median_days": itp["female"],
        "itp_female_minus_male_days": itp["female"] - itp["male"],
        "strains_above_itp_female_median": int((strains.female_median_days > itp["female"]).sum()),
    }])
    TABLES["inbred_strain_comparison"] = summary
    row = summary.iloc[0]
    print(f"{int(row.strains)} inbred strains: median female-male difference {row.strain_median_female_minus_male_days:.0f} d "
          f"(ITP {row.itp_female_minus_male_days:.0f} d); {int(row.strains_above_itp_female_median)} strains above the "
          f"ITP female median of {row.itp_female_median_days:.0f} d")


def control_lifespan_ceiling():
    """The control-survival comparison reported in Extended Data Table 7."""
    ROOT = Path(__file__).absolute().parents[1]
    SUPPLEMENTARY = ROOT / "supplementary" / "Supplementary_Data.xlsx"

    experiments = pd.read_excel(SUPPLEMENTARY, sheet_name="SD10_control_lifespan")
    experiments = experiments[experiments.sex.isin(["Female", "Male"])].copy()
    experiments["is_itp"] = experiments.itp_experiment.astype(str).str.upper().eq("ITP")

    def control_groups(frame):
        return frame.drop_duplicates(["pmid", "lifespan_control_days"])

    def study_medians(frame):
        return frame.groupby("pmid", as_index=False).lifespan_control_days.median()

    SPECIFICATIONS = [
        ("unique_control_groups", control_groups),
        ("study_level_median", study_medians),
        ("median_endpoints_only",
         lambda frame: control_groups(frame[frame.lifespan_measure.eq("median")])),
        ("published_2000_or_later", lambda frame: control_groups(frame[frame.year.ge(2000)])),
    ]

    rows = []
    for label, select in SPECIFICATIONS:
        for sex in ("Female", "Male"):
            itp = select(experiments[experiments.is_itp & experiments.sex.eq(sex)]).lifespan_control_days
            field = select(experiments[~experiments.is_itp & experiments.sex.eq(sex)]).lifespan_control_days
            itp_median = itp.median()
            rows.append({"specification": label, "sex": sex,
                         "itp_n": len(itp), "itp_median_days": itp_median,
                         "non_itp_n": len(field), "non_itp_median_days": field.median(),
                         "non_itp_reaching_itp_median": int((field >= itp_median).sum()),
                         "p_mannwhitney": stats.mannwhitneyu(itp, field,
                                                             alternative="two-sided").pvalue})

    table = pd.DataFrame(rows)
    TABLES["control_lifespan_comparison"] = table

    field_female = control_groups(experiments[~experiments.is_itp & experiments.sex.eq("Female")])
    correlation = stats.spearmanr(field_female.year, field_female.lifespan_control_days,
                                  nan_policy="omit")
    era = pd.DataFrame([{"quantity": "non_itp_female_year_lifespan_spearman_rho",
                         "value": correlation.statistic},
                        {"quantity": "non_itp_female_year_lifespan_p",
                         "value": correlation.pvalue},
                        {"quantity": "non_itp_female_median_before_2000",
                         "value": field_female[field_female.year < 2000]
                         .lifespan_control_days.median()},
                        {"quantity": "non_itp_female_median_2000_or_later",
                         "value": field_female[field_female.year >= 2000]
                         .lifespan_control_days.median()}])
    TABLES["control_lifespan_era"] = era

    for row in rows:
        print(f"  {row['specification']:<24}{row['sex']:<8}"
              f"ITP {row['itp_median_days']:.1f} (n = {row['itp_n']})  "
              f"non-ITP {row['non_itp_median_days']:.1f} (n = {row['non_itp_n']})  "
              f"reaching {row['non_itp_reaching_itp_median']}  P = {row['p_mannwhitney']:.3g}")
    print(f"  publication year against non-ITP female control lifespan: "
          f"rho = {correlation.statistic:.3f}; P = {correlation.pvalue:.3f}")


def founder_strain_weights():
    """The founder-strain comparison reported in Extended Data Table 8."""
    ROOT = Path(__file__).absolute().parents[1]
    INPUTS = ROOT / "data" / "inputs"

    def read_rows(path):
        with path.open(newline="", encoding="utf-8-sig") as handle:
            return list(csv.DictReader(handle))

    founders = read_rows(INPUTS / "founder_strain_mpd_reed1.csv")
    founder_mean = sum(float(row["reported_weight_g"]) for row in founders) / len(founders)

    mice = pd.read_excel(INPUTS / "ITP_bodyweight_all_cohorts.xlsx")
    mice["bw_6"] = within_range(mice["bw_6"])
    controls = mice[mice.group.eq("Control") & mice.sex.eq("m") & mice.bw_6.notna()
                    & (mice.lifespan_days > 6 * 30.4)]
    itp_male_controls = controls.bw_6.tolist()
    itp_mean = sum(itp_male_controls) / len(itp_male_controls)

    sensitivities = read_result("matched_quartile_survival").astype(str).to_dict("records")
    lowest = [
        row for row in sensitivities
        if row["quartile"] == "1"
    ]
    if len(lowest) != 1:
        raise RuntimeError(f"Expected one lowest-quarter row; found {len(lowest)}")
    lowest_mean = float(lowest[0]["mean_control_weight"])

    output = [
        {"quantity": "four_founder_mean_weight_g", "value": founder_mean},
        {"quantity": "itp_male_control_6m_n", "value": len(itp_male_controls)},
        {"quantity": "itp_male_control_6m_mean_weight_g", "value": itp_mean},
        {"quantity": "itp_minus_founder_mean_g", "value": itp_mean - founder_mean},
        {"quantity": "itp_percent_above_founder_mean", "value": 100 * (itp_mean / founder_mean - 1)},
        {"quantity": "lowest_itp_quarter_mean_weight_g", "value": lowest_mean},
        {"quantity": "lowest_itp_quarter_percent_above_founder_mean",
         "value": 100 * (lowest_mean / founder_mean - 1)},
        {"quantity": "heaviest_founder_reported_mean_weight_g",
         "value": max(float(row["reported_weight_g"]) for row in founders)},
    ]

    TABLES["founder_weight_comparison"] = pd.DataFrame(output)


def cross_study_summaries():
    """Cross-study and body-composition summaries from the supplement."""
    PACKAGE = Path(__file__).absolute().parents[1]
    SOURCE = PACKAGE / "supplementary" / "Supplementary_Data.xlsx"

    def clustered_ols(X, y, cluster):
        X = np.asarray(X, dtype=float)
        y = np.asarray(y, dtype=float)
        beta = np.linalg.solve(X.T @ X, X.T @ y)
        resid = y - X @ beta
        groups = np.unique(cluster)
        meat = np.zeros((X.shape[1], X.shape[1]))
        for group in groups:
            score = X[cluster == group].T @ resid[cluster == group]
            meat += np.outer(score, score)
        bread = np.linalg.inv(X.T @ X)
        n, k, g = len(y), X.shape[1], len(groups)
        correction = (g / (g - 1)) * ((n - 1) / (n - k))
        return beta, correction * bread @ meat @ bread

    Z = stats.norm.ppf(0.975)

    def clustered_summary(estimate, se, studies):
        """Interval and P value on a t reference with studies - 1 degrees of freedom."""
        df = studies - 1
        critical = stats.t.ppf(0.975, df)
        return (estimate - critical * se, estimate + critical * se,
                2 * stats.t.sf(abs(estimate / se), df))

    def fisher_ci(r, n):
        z, se = math.atanh(r), 1 / math.sqrt(n - 3)
        return math.tanh(z - Z * se), math.tanh(z + Z * se)

    rows = []

    d = pd.read_excel(SOURCE, sheet_name="SD9_cross_study_experiments")
    d["source_group"] = d.source_group.replace({"non-ITP": "non-ITP pharmacological"})
    d = d[(d.weight_change_percent < 0) & (d.lifespan_change_percent > 0)]
    d = d[d.sex.isin(["Male", "Female"])].copy()
    d["x"] = -d.weight_change_percent.astype(float)
    d["y"] = d.lifespan_change_percent.astype(float)
    d["study"] = d.pubmed_id.astype(str)

    for (source, sex), z in d.groupby(["source_group", "sex"]):
        beta, vcov = clustered_ols(z[["x"]], z.y, z.study.to_numpy())
        se = math.sqrt(vcov[0, 0])
        rows.append([f"slope_{source}_{sex}", len(z), z.study.nunique(), beta[0],
                     *clustered_summary(beta[0], se, z.study.nunique())])
        rows.append([f"median_ratio_{source}_{sex}", len(z), z.study.nunique(), (z.y / z.x).median(), np.nan, np.nan, np.nan])

    for source in ["ITP", "non-ITP pharmacological"]:
        z = d[d.source_group == source].copy()
        female = (z.sex == "Female").astype(float).to_numpy()
        X = np.column_stack([z.x.to_numpy(), z.x.to_numpy() * female])
        beta, vcov = clustered_ols(X, z.y, z.study.to_numpy())
        se = math.sqrt(vcov[1, 1])
        rows.append([f"female_minus_male_{source}", len(z), z.study.nunique(), beta[1],
                     *clustered_summary(beta[1], se, z.study.nunique())])

    z = d[d.source_group.isin(["ITP", "non-ITP pharmacological"])].copy()
    female = (z.sex == "Female").astype(float).to_numpy()
    non_itp = (z.source_group == "non-ITP pharmacological").astype(float).to_numpy()
    x = z.x.to_numpy()
    X = np.column_stack([x, x * female, x * non_itp, x * female * non_itp])
    beta, vcov = clustered_ols(X, z.y, z.study.to_numpy())
    se = math.sqrt(vcov[3, 3])
    rows.append(["pharmacological_sex_by_source", len(z), z.study.nunique(), beta[3],
                 *clustered_summary(beta[3], se, z.study.nunique())])

    z = d[(d.x >= 5) & d.source_group.isin(["ITP", "non-ITP pharmacological"])].copy()
    female = (z.sex == "Female").astype(float).to_numpy()
    non_itp = (z.source_group == "non-ITP pharmacological").astype(float).to_numpy()
    x = z.x.to_numpy()
    X = np.column_stack([x, x * female, x * non_itp, x * female * non_itp])
    beta, vcov = clustered_ols(X, z.y, z.study.to_numpy())
    se = math.sqrt(vcov[3, 3])
    rows.append(["pharmacological_sex_by_source_weight_loss_ge_5pct", len(z), z.study.nunique(),
                 beta[3], *clustered_summary(beta[3], se, z.study.nunique())])

    out = pd.DataFrame(rows, columns=["analysis", "experiments", "studies", "estimate", "ci_low", "ci_high", "p"])
    TABLES["cross_study_weight_response"] = out
    print(out.to_string(index=False))


def weight_change_response():
    """The within-cohort weight-change analysis reported in Extended Data Table 2."""
    MIN_EXPERIMENTS_PER_COHORT = 4

    experiments = pd.read_csv(DERIVED / "experiment_contrasts.csv")
    experiments = experiments.dropna(subset=["wt_pct", "life_pct"]).copy()

    def cluster_inference(estimate, se, groups):
        """Two-sided P value and 95% interval from a t reference with groups - 1 degrees of freedom."""
        df = max(groups - 1, 1)
        critical = stats.t.ppf(0.975, df)
        return (2 * stats.t.sf(abs(estimate / se), df),
                estimate - critical * se, estimate + critical * se)

    def clustered_fit(response, design, cluster):
        """Least squares with CR1 cluster-robust variance."""
        design = np.column_stack([np.ones(len(response)), design])
        bread = np.linalg.pinv(design.T @ design)
        beta = bread @ design.T @ response
        resid = response - design @ beta
        meat = np.zeros((design.shape[1], design.shape[1]))
        for key in pd.unique(cluster):
            rows = np.asarray(cluster == key)
            block = design[rows].T @ resid[rows]
            meat += np.outer(block, block)
        groups = len(pd.unique(cluster))
        n_obs, n_par = design.shape
        var = bread @ meat @ bread * (groups / max(groups - 1, 1)) * ((n_obs - 1) / (n_obs - n_par))
        se = np.sqrt(np.diag(var))
        return beta[1], se[1], groups

    def quartile_table():
        """Lifespan change by within-cohort quartile of weight change, for each sex."""
        rows, trends = [], []
        for sex in ("m", "f"):
            block = experiments[experiments.sex.eq(sex)].copy()
            counts = block.groupby("cohort").size()
            block = block[block.cohort.isin(counts[counts >= MIN_EXPERIMENTS_PER_COHORT].index)].copy()
            block["quartile"] = (block.groupby("cohort").wt_pct
                                 .transform(lambda s: pd.qcut(s, 4, labels=[1, 2, 3, 4],
                                                              duplicates="drop"))
                                 .astype(float))
            block = block.dropna(subset=["quartile"])
            for quartile, cell in block.groupby("quartile"):
                rows.append({"sex": sex, "quartile": int(quartile), "experiments": len(cell),
                             "median_weight_change_pct": cell.wt_pct.median(),
                             "mean_lifespan_change_pct": cell.life_pct.mean(),
                             "median_lifespan_change_pct": cell.life_pct.median(),
                             "share_with_longer_lifespan": (cell.life_pct > 0).mean()})
            slope, se, groups = clustered_fit(block.life_pct.to_numpy(float),
                                              block.quartile.to_numpy(float), block.cohort)
            p, lo, hi = cluster_inference(slope, se, groups)
            trends.append({"sex": sex, "experiments": len(block), "cohorts": block.cohort.nunique(),
                           "trend_per_quartile_pp": slope, "lo": lo, "hi": hi, "p": p})
        return pd.DataFrame(rows), pd.DataFrame(trends)

    def continuous_slopes():
        """Pooled and within-cohort slopes of lifespan change on weight change."""
        rows = []
        for sex in ("m", "f"):
            block = experiments[experiments.sex.eq(sex)].copy()
            pooled = clustered_fit(block.life_pct.to_numpy(float),
                                   block.wt_pct.to_numpy(float), block.cohort)
            centred = block.assign(
                wt=block.wt_pct - block.groupby("cohort").wt_pct.transform("mean"),
                life=block.life_pct - block.groupby("cohort").life_pct.transform("mean"))
            within = clustered_fit(centred.life.to_numpy(float),
                                   centred.wt.to_numpy(float), centred.cohort)
            for label, (slope, se, groups) in (("pooled", pooled), ("within_cohort", within)):
                p, lo, hi = cluster_inference(slope, se, groups)
                rows.append({"sex": sex, "specification": label, "experiments": len(block),
                             "cohorts": block.cohort.nunique(),
                             "slope_pct_lifespan_per_pct_weight": slope, "se": se,
                             "lo": lo, "hi": hi, "p": p})
        return pd.DataFrame(rows)

    def against_comparator():
        """Treatment weight change and control weight entered together on the cell hazard ratio."""
        cells = pd.read_csv(DERIVED / "comparator_cells.csv")
        joined = cells.merge(experiments[["cohort", "sex", "group", "wt_pct"]],
                             on=["cohort", "sex", "group"], how="inner")
        joined = joined.dropna(subset=["log_hr", "ctl_bw6", "wt_pct", "se"])
        joined["cluster"] = (joined.cohort.astype(str) + "|" + joined.site.astype(str)
                             + "|" + joined.sex)
        rows = []
        for sex in ("m", "f"):
            block = joined[joined.sex.eq(sex)]
            response = block.log_hr.to_numpy(float)
            weights = 1.0 / block.se.to_numpy(float) ** 2
            for terms, label in ((["wt_pct"], "weight_change_alone"),
                                 (["ctl_bw6"], "control_weight_alone"),
                                 (["wt_pct", "ctl_bw6"], "both_together")):
                design = np.column_stack([np.ones(len(block))]
                                         + [block[t].to_numpy(float) for t in terms])
                scaled = np.sqrt(weights)[:, None]
                bread = np.linalg.pinv((design * scaled).T @ (design * scaled))
                beta = bread @ (design * scaled).T @ (response * np.sqrt(weights))
                resid = response - design @ beta
                meat = np.zeros((design.shape[1], design.shape[1]))
                for key in pd.unique(block.cluster):
                    sel = (block.cluster == key).to_numpy()
                    piece = design[sel].T @ (weights[sel] * resid[sel])
                    meat += np.outer(piece, piece)
                groups = block.cluster.nunique()
                n_obs, n_par = design.shape
                var = bread @ meat @ bread * (groups / max(groups - 1, 1)) * ((n_obs - 1) / (n_obs - n_par))
                se = np.sqrt(np.diag(var))
                for position, term in enumerate(terms, start=1):
                    p, lo, hi = cluster_inference(beta[position], se[position], groups)
                    rows.append({"sex": sex, "model": label, "term": term,
                                 "cells": len(block), "control_groups": groups,
                                 "estimate": beta[position], "se": se[position],
                                 "lo": lo, "hi": hi, "p": p})
        correlations = [{"sex": sex, "model": "predictor_correlation", "term": "spearman_rho",
                         "cells": int((joined.sex == sex).sum()),
                         "control_groups": joined[joined.sex.eq(sex)].cluster.nunique(),
                         "estimate": stats.spearmanr(joined[joined.sex.eq(sex)].wt_pct,
                                                     joined[joined.sex.eq(sex)].ctl_bw6).statistic,
                         "se": np.nan, "lo": np.nan, "hi": np.nan,
                         "p": stats.spearmanr(joined[joined.sex.eq(sex)].wt_pct,
                                              joined[joined.sex.eq(sex)].ctl_bw6).pvalue}
                        for sex in ("m", "f")]
        return pd.DataFrame(rows + correlations)

    quartiles, trends = quartile_table()
    TABLES["weight_response_quartiles"] = quartiles
    TABLES["weight_change_trend"] = trends
    TABLES["weight_change_slopes"] = continuous_slopes()
    TABLES["weight_change_vs_comparator"] = against_comparator()

    for row in trends.itertuples(index=False):
        name = "males" if row.sex == "m" else "females"
        print(f"{name}: {row.experiments} experiments in {row.cohorts} cohorts; trend "
              f"{row.trend_per_quartile_pp:+.2f} percentage points per quartile toward less weight "
              f"loss (95% CI {row.lo:+.2f} to {row.hi:+.2f}; P = {row.p:.4f})")
    for row in quartiles.itertuples(index=False):
        name = "male" if row.sex == "m" else "female"
        print(f"  {name:<7} Q{row.quartile}: {row.experiments:>2} experiments, median weight "
              f"{row.median_weight_change_pct:+5.1f}%, mean lifespan "
              f"{row.mean_lifespan_change_pct:+5.1f}%")


def main():
    sex_survival_context()
    control_early_mortality()
    inbred_strain_comparison()
    control_lifespan_ceiling()
    founder_strain_weights()
    cross_study_summaries()
    weight_change_response()
    write_results(TABLES, combined={"weight_response_models": [
        "weight_change_trend", "weight_change_slopes", "weight_change_vs_comparator"]})


if __name__ == "__main__":
    main()
