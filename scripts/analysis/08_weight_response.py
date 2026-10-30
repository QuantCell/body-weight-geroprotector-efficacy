"""Analyze intervention-level weight change within ITP and across studies."""
from __future__ import annotations

from pathlib import Path
import math

from scipy import stats
import numpy as np
import pandas as pd


def analyze_itp_weight_response():
    ROOT = Path(__file__).absolute().parents[2]

    OUTPUTS = ROOT / "data" / "derived"

    def clustered_fit(frame, response, predictors, groups, weights=None):
        """OLS/WLS with an intercept, clustered sandwich variance and small-sample correction."""
        x = np.column_stack([np.ones(len(frame)), frame[predictors].to_numpy(float)])
        y = frame[response].to_numpy(float)
        w = np.ones(len(y)) if weights is None else np.asarray(weights, dtype=float)
        cluster = np.asarray(groups)
        n, k = x.shape
        g = len(pd.unique(cluster))
        if g < 2 or n <= k or np.linalg.matrix_rank(x) < k:
            raise ValueError("Insufficient clusters, residual degrees of freedom, or design rank")
        if not np.isfinite(x).all() or not np.isfinite(y).all() or not (w > 0).all():
            raise ValueError("Non-finite model data or non-positive weights")
        bread = np.linalg.pinv(x.T @ (w[:, None] * x))
        beta = bread @ x.T @ (w * y)
        residual = y - x @ beta
        meat = np.zeros((k, k))
        for group in pd.unique(cluster):
            selected = cluster == group
            score = x[selected].T @ (w[selected] * residual[selected])
            meat += np.outer(score, score)
        correction = g / (g - 1) * (n - 1) / (n - k)
        variance = correction * bread @ meat @ bread
        se = np.sqrt(np.diag(variance))
        critical = stats.t.ppf(0.975, g - 1)
        return [dict(term=term, estimate=beta[i], se=se[i],
                     lo=beta[i] - critical * se[i], hi=beta[i] + critical * se[i],
                     p=2 * stats.t.sf(abs(beta[i] / se[i]), g - 1),
                     clusters=g, df=g - 1, n=n)
                for i, term in enumerate(["intercept", *predictors]) if i]

    def analyze(experiments, cells):
        if experiments.duplicated(["cohort", "sex", "group"]).any():
            raise ValueError("Duplicate experiment contrast key")
        quartiles, trends, slopes, joint = [], [], [], []
        for sex in ("m", "f"):
            block = experiments[experiments.sex.eq(sex)].copy()
            for label in ("pooled", "within_cohort"):
                z = block.copy()
                if label == "within_cohort":
                    for term in ("wt_pct", "life_pct"):
                        z[term] -= z.groupby("cohort")[term].transform("mean")
                row = clustered_fit(z, "life_pct", ["wt_pct"], z.cohort)[0]
                slopes.append(dict(sex=sex, specification=label, **row))
            eligible = block[block.groupby("cohort").cohort.transform("size") >= 4].copy()
            eligible["quartile"] = eligible.groupby("cohort").wt_pct.transform(
                lambda s: pd.qcut(s, 4, labels=[1, 2, 3, 4])).astype(int)
            row = clustered_fit(eligible, "life_pct", ["quartile"], eligible.cohort)[0]
            trends.append(dict(sex=sex, **row))
            for q, z in eligible.groupby("quartile"):
                quartiles.append(dict(sex=sex, quartile=q, experiments=len(z),
                                      median_weight_change_pct=z.wt_pct.median(),
                                      mean_lifespan_change_pct=z.life_pct.mean(),
                                      median_lifespan_change_pct=z.life_pct.median(),
                                      share_with_longer_lifespan=z.life_pct.gt(0).mean()))
        joined = cells.merge(experiments[["cohort", "sex", "group", "wt_pct"]],
                             on=["cohort", "sex", "group"], validate="many_to_one")
        joined = joined.dropna(subset=["log_hr", "ctl_bw6", "wt_pct", "se"])
        for sex in ("m", "f"):
            z = joined[joined.sex.eq(sex)].copy()
            groups = z.cohort.astype(str) + "|" + z.site.astype(str)
            for label, terms in [("weight_change_alone", ["wt_pct"]),
                                 ("control_weight_alone", ["ctl_bw6"]),
                                 ("both_together", ["wt_pct", "ctl_bw6"])]:
                for row in clustered_fit(z, "log_hr", terms, groups, 1 / z.se ** 2):
                    joint.append(dict(sex=sex, model=label, **row))
            rho = stats.spearmanr(z.wt_pct, z.ctl_bw6)
            joint.append(dict(sex=sex, model="predictor_correlation", term="spearman_rho",
                              estimate=rho.statistic, n=len(z), clusters=groups.nunique(),
                              se=np.nan, lo=np.nan, hi=np.nan, p=np.nan, df=np.nan))
        return {"weight_change_quartiles.csv": pd.DataFrame(quartiles),
                "weight_change_trend.csv": pd.DataFrame(trends),
                "weight_change_slopes.csv": pd.DataFrame(slopes),
                "weight_change_vs_comparator.csv": pd.DataFrame(joint)}

    experiments = pd.read_csv(OUTPUTS / "experiment_weight_changes.csv")
    cells = pd.read_csv(OUTPUTS / "comparator_effect_cells.csv")
    for filename, result in analyze(experiments, cells).items():
        result.to_csv(OUTPUTS / filename, index=False, lineterminator="\n")
        print(f"{filename}: {len(result)} rows")


def validate_semaglutide_extraction():
    """Check the curated comparison against the published lifespan-cohort source data."""
    root = Path(__file__).absolute().parents[2]
    inputs = root / "data/inputs"
    survival = pd.read_excel(inputs / "41586_2026_10940_MOESM6_ESM.xlsx", header=None)
    control = pd.to_numeric(survival.iloc[1, 1:], errors="coerce").dropna()
    treated = pd.to_numeric(survival.iloc[2, 1:], errors="coerce").dropna()
    assert (len(control), len(treated)) == (39, 40)
    assert (control.median(), treated.median()) == (742, 834)
    weights = pd.read_excel(inputs / "41586_2026_10940_MOESM11_ESM.xlsx", header=None)
    # Use the observation nearest three months after treatment started at day 600.
    days = pd.to_numeric(weights.iloc[2:33, 0], errors="raise")
    row = weights.loc[(days - 690).abs().idxmin()]
    assert row.iloc[0] == 691
    control_weight = pd.to_numeric(row.iloc[1:40], errors="coerce").dropna()
    treated_weight = pd.to_numeric(row.iloc[40:80], errors="coerce").dropna()
    assert (len(control_weight), len(treated_weight)) == (30, 35)
    weight_change = 100 * (treated_weight.mean() / control_weight.mean() - 1)
    lifespan_change = 100 * (treated.median() / control.median() - 1)
    supplement = root / "supplementary/Supplementary_Data_1-10.xlsx"
    corpus = pd.read_excel(supplement, sheet_name="SD9_cross_study_experiments")
    record = corpus[pd.to_numeric(corpus.pubmed_id, errors="coerce").eq(42686906)]
    assert len(record) == 1
    assert np.isclose(record.iloc[0].weight_change_percent, weight_change, atol=1e-12)
    assert np.isclose(record.iloc[0].lifespan_change_percent, lifespan_change, atol=1e-12)
    controls = pd.read_excel(supplement, sheet_name="SD10_control_lifespan")
    record = controls[controls.pmid.eq(42686906)]
    assert len(record) == 1 and record.iloc[0].lifespan_control_days == control.median()


def analyze_cross_study_weight_response():
    """Analyze recorded DrugAge entries and source-supported curated additions."""


    PACKAGE = Path(__file__).absolute().parents[2]

    SOURCE = PACKAGE / "supplementary" / "Supplementary_Data_1-10.xlsx"

    OUTPUT = PACKAGE / "data" / "derived" / "supplementary_analyses.csv"

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

    def cluster_p(z, groups):
        return 2 * stats.t.sf(abs(float(z)), groups - 1)

    def critical(groups):
        return stats.t.ppf(0.975, groups - 1)

    rows = []

    d = pd.read_excel(SOURCE, sheet_name="SD9_cross_study_experiments")

    d = d[(d.weight_change_percent < 0) & (d.lifespan_change_percent > 0)]

    d = d[d.sex.isin(["Male", "Female"])].copy()

    d["x"] = -d.weight_change_percent.astype(float)

    d["y"] = d.lifespan_change_percent.astype(float)

    d["study"] = d.pubmed_id.astype(str)

    for (source, sex), z in d.groupby(["source_group", "sex"]):
        beta, vcov = clustered_ols(z[["x"]], z.y, z.study.to_numpy())
        se = math.sqrt(vcov[0, 0])
        Z = critical(z.study.nunique())
        rows.append([f"slope_{source}_{sex}", len(z), z.study.nunique(), beta[0], beta[0] - Z * se, beta[0] + Z * se, cluster_p(beta[0] / se, z.study.nunique())])
        rows.append([f"median_ratio_{source}_{sex}", len(z), z.study.nunique(),
                     (z.y / z.x).median(), np.nan, np.nan, np.nan])

    for source in ["ITP", "non-ITP pharmacological", "dietary restriction"]:
        z = d[d.source_group == source].copy()
        female = (z.sex == "Female").astype(float).to_numpy()
        X = np.column_stack([z.x.to_numpy(), z.x.to_numpy() * female])
        beta, vcov = clustered_ols(X, z.y, z.study.to_numpy())
        se = math.sqrt(vcov[1, 1])
        Z = critical(z.study.nunique())
        rows.append([f"female_minus_male_{source}", len(z), z.study.nunique(), beta[1], beta[1] - Z * se, beta[1] + Z * se, cluster_p(beta[1] / se, z.study.nunique())])

    z = d[d.source_group.isin(["ITP", "non-ITP pharmacological"])].copy()

    female = (z.sex == "Female").astype(float).to_numpy()

    non_itp = (z.source_group == "non-ITP pharmacological").astype(float).to_numpy()

    x = z.x.to_numpy()

    X = np.column_stack([x, x * female, x * non_itp, x * female * non_itp])

    beta, vcov = clustered_ols(X, z.y, z.study.to_numpy())

    se = math.sqrt(vcov[3, 3])

    Z = critical(z.study.nunique())

    rows.append(["pharmacological_sex_by_source", len(z), z.study.nunique(), beta[3], beta[3] - Z * se, beta[3] + Z * se, cluster_p(beta[3] / se, z.study.nunique())])

    z = d[(d.x >= 5) & d.source_group.isin(["ITP", "non-ITP pharmacological"])].copy()

    female = (z.sex == "Female").astype(float).to_numpy()

    non_itp = (z.source_group == "non-ITP pharmacological").astype(float).to_numpy()

    x = z.x.to_numpy()

    X = np.column_stack([x, x * female, x * non_itp, x * female * non_itp])

    beta, vcov = clustered_ols(X, z.y, z.study.to_numpy())

    se = math.sqrt(vcov[3, 3])

    Z = critical(z.study.nunique())

    rows.append(["pharmacological_sex_by_source_weight_loss_ge_5pct", len(z), z.study.nunique(), beta[3], beta[3] - Z * se, beta[3] + Z * se, cluster_p(beta[3] / se, z.study.nunique())])

    out = pd.DataFrame(rows, columns=["analysis", "experiments", "studies", "estimate", "ci_low", "ci_high", "p_t"])

    out.to_csv(OUTPUT, index=False)

    print(out.to_string(index=False))


def main():
    analyze_itp_weight_response()
    validate_semaglutide_extraction()
    analyze_cross_study_weight_response()

if __name__ == "__main__":
    main()
