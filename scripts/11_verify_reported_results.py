"""Check regenerated outputs against the values reported in the manuscript."""
from pathlib import Path
import math
import pandas as pd
from scipy.stats import norm as normal, t as t_dist

from common import DERIVED, read_result

root = Path(__file__).absolute().parents[1]
out = root / "data" / "outputs"

primary = read_result("primary_interaction").set_index("sex")
assert int(primary.loc["m", "experiments"]) == 92 and int(primary.loc["m", "unique_mice"]) == 17943
assert int(primary.loc["f", "experiments"]) == 88 and int(primary.loc["f", "unique_mice"]) == 15837
assert math.isclose(primary.loc["m", "HR"], 0.9448510519, rel_tol=0, abs_tol=1e-8)

split = read_result("weight_class").set_index("analysis")
assert int(split.loc["male_weight_lowering", "experiments"]) == 42
assert math.isclose(split.loc["male_weight_lowering", "HR"], 0.896835, rel_tol=0, abs_tol=1e-6)
assert int(split.loc["male_predosing_weight_lowering", "experiments"]) == 30
eligible = pd.read_csv(DERIVED / "quartile_eligibility.csv")
assert ((eligible.cohort == "C2011") & (eligible.group == "Prot")).sum() == 1, "Protandim dropped"
for name, (experiments, mice, hr, lo, hi, p) in {
        "male_nonlowering": (50, 11878, 0.98, 0.94, 1.03, 0.47),
        "female_weight_lowering": (40, 9289, 1.03, 0.98, 1.07, 0.22),
        "female_nonlowering": (48, 10504, 1.03, 0.99, 1.07, 0.16),
        "male_predosing_all": (60, 12689, 0.95, 0.91, 0.99, 0.018),
        "male_predosing_weight_lowering": (30, 7458, 0.89, 0.85, 0.94, 1.9e-05),
        "female_predosing_all": (58, 11703, 1.02, 0.99, 1.07, 0.22)}.items():
    row = split.loc[name]
    assert (int(row.experiments), int(row.unique_mice)) == (experiments, mice), f"{name}: n"
    assert (round(row.HR, 2), round(row.lo, 2), round(row.hi, 2)) == (hr, lo, hi), name
    assert float(f"{row.p:.2g}") == p, f"{name}: P {row.p:.4g}"
assert round(split.loc["male_weight_lowering", "hi"], 2) == 0.94
for name, p in (("male_weight_lowering", 4.2e-06), ("male_class_contrast", 0.00022),
                ("male_predosing_class_contrast", 0.00016)):
    assert float(f"{split.loc[name, 'p']:.2g}") == p, f"{name}: P"
nonlowering = split.loc["male_predosing_nonlowering"]
assert int(nonlowering.experiments) == 30 and round(nonlowering.HR, 2) == 1.01
assert round(nonlowering.lo, 2) == 0.96 and round(nonlowering.hi, 2) == 1.07
assert round(nonlowering.p, 2) == 0.64

q = read_result("matched_quartile_survival")
assert q["gain_pct"].round(1).tolist() == [6.2, 8.0, 12.3, 20.2]
trend = read_result("quartile_trend").iloc[0]
assert int(trend.experiments) == 30 and int(trend.unique_mice) == 8137
assert math.isclose(trend.HR_per_quartile, 0.9087365722, rel_tol=0, abs_tol=1e-8)
assert float(f"{trend.p:.2g}") == 4.7e-05
assert q.standardized_control_median.tolist() == [890, 840, 766, 673]
assert q.standardized_treated_median.tolist() == [945, 907, 860, 809]
day500 = read_result("survival_beyond_day500").set_index(["follow_up", "quartile"])
assert [round(day500.loc[("from six months", k), "control_dead_by_day_500_pct"], 1) for k in ("1", "4")] == [6.3, 27.6]
for follow_up, days, pct in (("from six months", 85, 10.6), ("beyond day 500", 67, 7.9)):
    row = day500.loc[(follow_up, "all")]
    assert (row.gain_days, round(row.gain_pct, 1)) == (days, pct), follow_up
early = read_result("control_early_mortality").set_index("quantity").value
assert int(early["male_controls"]) == 4983
assert [round(early[f"q{k}_dead_by_day_500_pct"], 1) for k in (1, 4)] == [6.1, 27.4]
assert round(early["mean_lifespan_difference_days"]) == 195
assert round(early["difference_from_mortality_before_day_500_days"]) == 101
trend500 = read_result("quartile_trend_beyond_day500").iloc[0]
assert int(trend500.experiments) == 30 and int(trend500.unique_mice) == 6884
assert float(f"{trend500.p:.2g}") == 0.00036

rmst = read_result("restricted_mean_survival").set_index(["grouping", "quartile"])
assert [round(rmst.loc[("separate ranks", k), "gain_days"], 1) for k in (1, 4)] == [37.3, 88.9]
rmst_contrast = read_result("restricted_mean_contrasts").set_index("grouping")
separate = rmst_contrast.loc["separate ranks"]
assert int(separate.replicates) == 2000
assert (round(separate.heaviest_minus_lightest_days, 1), round(separate.lo, 1), round(separate.hi, 1)) == (51.6, 26.7, 76.9)
cutoffs = read_result("quartile_trend_control_cutoffs").iloc[0]
assert int(cutoffs.experiments) == 30 and int(cutoffs.experiment_site_cells) == 88
assert (round(cutoffs.HR_per_quartile, 2), round(cutoffs.lo, 2), round(cutoffs.hi, 2)) == (0.92, 0.87, 0.96)
assert float(f"{cutoffs.p:.2g}") == 0.00018

strata = read_result("quartile_sensitivities")
for name, (light, heavy) in {"weight-lowering standardized thirds": (7.3, 20.4),
                             "weight-lowering standardized quarters": (6.2, 20.2),
                             "weight-lowering standardized fifths": (4.0, 22.2)}.items():
    gains = strata[strata.analysis == name].sort_values("bin").gain_pct
    assert (round(gains.iloc[0], 1), round(gains.iloc[-1], 1)) == (light, heavy), name
above = strata[strata.analysis.str.startswith("weight-lowering") & (strata.control_median > 900)]
assert above[["analysis", "bin", "control_median", "treated_median"]].values.tolist() == [
    ["weight-lowering standardized fifths", 1, 907, 943]]
standardized = strata[strata.analysis == "weight-lowering standardized quarters"].sort_values("bin")
assert standardized.gain_pct.round(6).tolist() == q.gain_pct.round(6).tolist()

comp = read_result("compound_survival")
assert len(comp) == 15 and int(comp.iloc[-1].experiments) == 26
assert int(comp.iloc[-1].unique_mice) == 7238
assert round(comp.iloc[-1].lowest_gain_pct, 1) == 5.0
assert round(comp.iloc[-1].highest_gain_pct, 1) == 20.2
assert "Protandim" not in set(comp.compound)
assert int((comp.iloc[:-1].highest_gain_pct > comp.iloc[:-1].lowest_gain_pct).sum()) == 13
named = comp.set_index("compound")
for compound, lowest, highest in [("Acarbose", 7.9, 24.8),
                                  ("17alpha-estradiol", 5.9, 22.9),
                                  ("Canagliflozin", 6.1, 23.2),
                                  ("Rapamycin", 7.3, 21.0),
                                  ("Glycine", 8.3, 1.5)]:
    row = named.loc[compound]
    assert round(row.lowest_gain_pct, 1) == lowest, f"{compound}: {row.lowest_gain_pct:.3f}"
    assert round(row.highest_gain_pct, 1) == highest, f"{compound}: {row.highest_gain_pct:.3f}"
pooled = comp.iloc[-1]
assert round(pooled.overall_gain_pct, 1) == 9.7
assert round(pooled.overall_gain_pct / pooled.lowest_gain_pct, 1) == 1.9
common = read_result("compound_common_interaction").iloc[0]
assert int(common.compounds) == 14 and int(common.experiments) == 26
assert math.isclose(common.estimate, 0.8964465, rel_tol=0, abs_tol=1e-7)
assert float(f"{common.p:.2g}") == 0.00011

marks = read_result("landmark_models").set_index(["landmark", "sex"])
for (landmark, sex), (hr, lo, hi, p) in {(18, "m"): (1.01, 0.95, 1.06, 0.81),
                                         (24, "m"): (1.11, 1.03, 1.18, 0.0031),
                                         (12, "m"): (0.9, 0.85, 0.96, 0.00055)}.items():
    row = marks.loc[(landmark, sex)]
    assert (round(row.HR, 2), round(row.lo, 2), round(row.hi, 2)) == (hr, lo, hi), landmark
    assert float(f"{row.p:.2g}") == p, f"landmark {landmark}: P {row.p:.3g}"
assert round(marks.loc[(6, "m"), "HR"], 2) == 0.8 and round(marks.loc[(6, "m"), "p"], 3) == 0.003
assert marks.xs("f", level="sex").p.min() > 0.05
lag = read_result("landmark_lag90").set_index(["landmark", "sex"])
assert round(lag.loc[(12, "m"), "HR"], 2) == 0.91 and round(lag.loc[(12, "m"), "p"], 3) == 0.002
loo = read_result("landmark_leave_one_cohort_out")
assert round(loo.HR.min(), 2) == 0.88 and round(loo.HR.max(), 2) == 0.92 and loo.p.max() < 0.05

cells = pd.read_csv(DERIVED / "comparator_cell_medians.csv")
assert len(cells) == 534
assert len(cells.query("cohort == 'C2021' and site == 'TJL'")) == 12
comparator = read_result("comparator_strata")
tertiles = comparator[comparator.split.eq(3)].set_index(["sex", "bin"])
assert round(tertiles.loc[("m", "Q1"), "mean_cell_km_median_gain_pct"], 1) == 4.0
assert round(tertiles.loc[("m", "Q3"), "mean_cell_km_median_gain_pct"], 1) == 8.3
assert round(tertiles.loc[("f", "Q1"), "mean_cell_km_median_gain_pct"], 1) == 1.9
assert round(tertiles.loc[("f", "Q3"), "mean_cell_km_median_gain_pct"], 1) == 2.8
for (sex, bin_), reduction in {("m", "Q1"): 10.9, ("m", "Q3"): 21.5,
                               ("f", "Q1"): 5.3, ("f", "Q3"): 12.0}.items():
    got = tertiles.loc[(sex, bin_), "hazard_reduction_pct"]
    assert round(got, 1) == reduction, f"{sex} {bin_} hazard reduction {got:.2f}"
gradient = read_result("comparator_gradients").set_index("sex")
assert (gradient.analysis == "within_cohort_clustered").all() and (gradient.control_clusters == 50).all()
for sex, (hr, lo, hi, p) in {"m": (0.96, 0.93, 1.0, 0.034), "f": (1.02, 0.98, 1.06, 0.27)}.items():
    row = gradient.loc[sex]
    assert (round(row.hr_per_sd, 2), round(row.ci_low, 2), round(row.ci_high, 2)) == (hr, lo, hi), sex
    assert round(row.p, 3 if sex == "m" else 2) == p, f"{sex}: gradient P {row.p:.4f}"
assert round(gradient.loc["m", "ci_high"], 3) == 0.997

# Extended Data Table 1, attained-weight rows: n, HR, both confidence limits and P.
age = read_result("control_weight_gradients").set_index(["weight_age", "sex"])
EXPECTED = {
    (6, "m"): (4983, 1.39, 1.35, 1.43, 7.1e-100),
    (6, "f"): (4553, 1.08, 1.04, 1.11, 1.5e-06),
    (12, "m"): (4572, 1.36, 1.32, 1.41, 2.1e-76),
    (12, "f"): (4537, 1.1, 1.07, 1.13, 7.9e-10),
    (18, "m"): (3900, 1.13, 1.09, 1.17, 1.2e-11),
    (18, "f"): (4450, 1.05, 1.02, 1.09, 0.00088),
    (24, "m"): (2670, 0.87, 0.84, 0.91, 1.3e-09),
    (24, "f"): (3508, 0.97, 0.94, 1.0, 0.085),
}
for key, (n, hr, lo, hi, p) in EXPECTED.items():
    row = age.loc[key]
    assert int(row.n) == n, f"{key}: n {int(row.n)} != {n}"
    assert round(row.HR_per_SD, 2) == hr, f"{key}: HR {row.HR_per_SD:.4f} != {hr}"
    assert round(row.lo, 2) == lo and round(row.hi, 2) == hi, f"{key}: CI"
    assert float(f"{row.p:.1e}") == p, f"{key}: P {row.p:.3g} != {p:.3g}"

joint = read_result("developmental_estimates").set_index("label")
assert int(joint.loc["6-month adjusted", "n"]) == 539 and round(joint.loc["6-month adjusted", "HR"], 2) == 1.46
# Results: 1.4-month weight on its own in male controls
early_alone = read_result("developmental_estimates")
early_alone = early_alone[(early_alone.panel == "age") & (early_alone.label.astype(str) == "1.4") & (early_alone.sex == "m")].iloc[0]
assert int(early_alone.n) == 562
assert (round(early_alone.HR, 2), round(early_alone.lo, 2), round(early_alone.hi, 2)) == (1.15, 1.05, 1.25)
assert float(f"{early_alone.p:.1e}") == 2.1e-03
maternal = read_result("maternal_adjustment")
assert (maternal.HR_weight - maternal.HR_adjusted_maternal).abs().max() < 0.01
puberty = read_result("puberty_moderation").set_index("sex")
assert int(puberty.loc["m", "experiments"]) == 4 and int(puberty.loc["m", "unique_mice"]) == 730
assert round(puberty.loc["m", "HR"], 2) == 1.00 and round(puberty.loc["m", "p"], 2) == 0.97
assert int(puberty.loc["f", "experiments"]) == 2 and int(puberty.loc["f", "unique_mice"]) == 553
assert round(puberty.loc["f", "HR"], 2) == 1.03 and round(puberty.loc["f", "p"], 2) == 0.77

trajectory = read_result("adult_weight_change")
early = trajectory[trajectory.timing.eq("primary") & trajectory.adjusted_for_w6].set_index("sex")
for sex, n_mice, n_cells, hr, lo, hi, p in [
    ("m", 14118, 257, 1.1, 1.07, 1.12, 3.1e-16),
    ("f", 13167, 251, 1.04, 1.03, 1.06, 5.5e-07),
]:
    row = early.loc[sex]
    assert int(row.n_mice) == n_mice and int(row.cells) == n_cells
    assert round(row.hr, 2) == hr and round(row.ci_low, 2) == lo and round(row.ci_high, 2) == hi
    assert float(f"{row.p:.1e}") == p, f"{sex}: P {row.p:.3g}"

treated = trajectory[trajectory.timing.eq("treated_primary")].iloc[0]
assert int(treated.n_mice) == 9705 and int(treated.cells) == 208
assert round(treated.hr, 2) == 1.1
assert round(treated.ci_low, 2) == 1.07 and round(treated.ci_high, 2) == 1.13
assert float(f"{treated.p:.1e}") == 5.4e-11, f"early-treated P {treated.p:.3g}"

trajectory_input = pd.read_csv(DERIVED / "early_weight_change.csv")
males = trajectory_input[trajectory_input.sex.eq("m")]
controls = males[males.is_control.astype(bool)]
dosed = males[~males.is_control.astype(bool) & males.init_months.notna()]
counts = (len(controls), int((dosed.init_months < 12).sum()),
          int((dosed.init_months >= 12).sum()))
assert counts == (4413, 9705, 2697), f"male trajectory counts {counts}"
assert counts[0] + counts[1] == 14118, "primary male set is not controls plus early-treated"
assert counts[0] + counts[2] == 7110, "pre-treatment set is not controls plus late-treated"

fat = read_result("cross_study_weight_response").set_index("analysis")

for source, sex, experiments, studies, slope, lo, hi, ratio in [
    ("ITP", "Female", 25, 12, 0.43, 0.09, 0.77, 0.40),
    ("ITP", "Male", 30, 13, 1.00, 0.71, 1.29, 1.24),
    ("non-ITP pharmacological", "Female", 21, 15, 0.86, 0.56, 1.16, 1.00),
    ("non-ITP pharmacological", "Male", 26, 21, 0.70, 0.33, 1.07, 1.27),
]:
    row = fat.loc[f"slope_{source}_{sex}"]
    assert int(row.experiments) == experiments and int(row.studies) == studies, f"{source} {sex}: n"
    assert round(row.estimate, 2) == slope, f"{source} {sex}: slope {row.estimate:.4f}"
    assert round(row.ci_low, 2) == lo and round(row.ci_high, 2) == hi, f"{source} {sex}: CI"
    assert round(fat.loc[f"median_ratio_{source}_{sex}", "estimate"], 2) == ratio, \
        f"{source} {sex}: median ratio"
for label, difference, lo, hi, p in [("ITP", -0.57, -1.02, -0.12, 0.017),
                                     ("non-ITP pharmacological", 0.16, -0.28, 0.59, 0.46)]:
    row = fat.loc[f"female_minus_male_{label}"]
    assert round(row.estimate, 2) == difference, f"{label}: {row.estimate:.4f}"
    assert round(row.ci_low, 2) == lo and round(row.ci_high, 2) == hi, \
        f"{label}: CI {row.ci_low:.4f} to {row.ci_high:.4f}"
    assert float(f"{row.p:.2g}") == p, f"{label}: P {row.p:.4g}"
assert round(fat.loc["pharmacological_sex_by_source", "p"], 3) == 0.018
assert round(fat.loc["pharmacological_sex_by_source_weight_loss_ge_5pct", "p"], 3) == 0.014

strain = pd.read_csv(out / "strain_context_summary.csv").set_index("quantity")
assert int(strain.loc["age_matched_strains", "value"]) == 27
assert round(float(strain.loc["survey_mean_weight_g", "value"]), 1) == 30.3
assert int(strain.loc["male_weight_strains", "value"]) == 27
assert int(strain.loc["male_weight_at_or_above_lightest_quarter", "value"]) == 7
assert int(strain.loc["male_body_fat_strains", "value"]) == 40
assert int(strain.loc["female_body_fat_strains", "value"]) == 40
assert int(strain.loc["male_body_fat_at_or_above_itp", "value"]) == 2
assert int(strain.loc["female_body_fat_at_or_above_itp", "value"]) == 1

founder = read_result("founder_weight_comparison").set_index("quantity")
for quantity, expected in [("four_founder_mean_weight_g", 29.4), ("itp_male_control_6m_n", 4983),
                           ("itp_male_control_6m_mean_weight_g", 40.6),
                           ("itp_percent_above_founder_mean", 38.0),
                           ("lowest_itp_quarter_mean_weight_g", 34.4),
                           ("heaviest_founder_reported_mean_weight_g", 32.0)]:
    got = float(founder.loc[quantity, "value"])
    assert round(got, 1) == expected, f"{quantity}: {got:.3f} != {expected}"

# Extended Data Table 4, the external Diversity Outbred cohort.
do_quarters = read_result("restriction_quartiles").set_index("quarter")
assert do_quarters.al_n.tolist() == [47, 47, 47, 47]
assert do_quarters.cr_n.tolist() == [93, 93, 92, 93]
assert do_quarters.mean_weight_g.round(1).tolist() == [24.8, 29.8, 34.6, 42.1]
assert do_quarters.al_km_median_days.tolist() == [877, 825, 779, 714]
assert do_quarters.cr_km_median_days.tolist() == [973, 974, 931, 946]
# Discussion: lightest-to-heaviest fall in median survival under each feeding regimen
assert (int(do_quarters.al_km_median_days[1] - do_quarters.al_km_median_days[4]),
        int(do_quarters.cr_km_median_days[1] - do_quarters.cr_km_median_days[4])) == (163, 27)
assert do_quarters.gain_pct.round(1).tolist() == [10.9, 18.1, 19.5, 32.5]
do = read_result("restriction_models").set_index("analysis")
for analysis, n, estimate, lo, hi in [
    ("ad_libitum_weight_gradient", 188, 1.21, 1.04, 1.41),
    ("pooled_restriction_weight_gradient", 371, 1.05, 0.95, 1.17),
    ("treatment_by_quarter_trend", 559, 0.90, 0.76, 1.05),
]:
    row = do.loc[analysis]
    assert int(row.n) == n, f"{analysis}: n {int(row.n)} != {n}"
    assert (round(row.estimate, 2) == estimate
            and round(row.lo, 2) == lo
            and round(row.hi, 2) == hi), (
        f"{analysis}: {row.estimate:.3f} ({row.lo:.3f}-{row.hi:.3f})")
assert round(do.loc["ad_libitum_weight_gradient", "p"], 3) == 0.015
assert round(do.loc["treatment_by_quarter_trend", "p"], 2) == 0.18
assert int(do.loc["ad_libitum_day500_mortality_pct", "n"]) == 188
assert round(do.loc["ad_libitum_day500_mortality_pct", "estimate"], 1) == 11.7
assert round(do.loc["ad_libitum_day500_mortality_all_events_pct", "estimate"], 1) == 12.2

# Extended Data Table 6, sex-specific control survival within the ITP.
context = read_result("sex_survival_summary").set_index(["cohort_set", "sex"])
full = "all_17_cohorts"
assert int(context.loc[(full, "male"), "n"]) == 5352
assert int(context.loc[(full, "female"), "n"]) == 4833
assert int(context.loc[(full, "male"), "cohorts"]) == 17
assert context.loc[(full, "male"), "km_median_days"] == 800
assert context.loc[(full, "female"), "km_median_days"] == 886
assert context.loc[(full, "female_minus_male"), "km_median_days"] == 86
assert round(context.loc[(full, "male"), "day500_mortality_pct"], 1) == 16.0
assert round(context.loc[(full, "female"), "day500_mortality_pct"], 1) == 2.4
assert round(context.loc[(full, "female"), "day500_mortality_all_events_pct"], 1) == 3.0
six = "published_six_cohorts"
assert int(context.loc[(six, "male"), "cohorts"]) == 6
assert int(context.loc[(six, "male"), "n"]) == 1993
assert int(context.loc[(six, "female"), "n"]) == 1697
assert context.loc[(six, "male"), "km_median_days"] == 803
assert context.loc[(six, "female"), "km_median_days"] == 887
assert context.loc[(six, "female_minus_male"), "km_median_days"] == 84

external = read_result("external_survival_context").set_index("population")
for population, male, female, male_n, female_n, difference in [
    ("Independent UM-HET3", 831, 836, 14, 16, 5),
    ("Independent 5LG6 cohort", 804, 865, 40, 45, 61),
]:
    row = external.loc[population]
    assert int(row.male_survival_days) == male and int(row.female_survival_days) == female
    assert int(row.male_n) == male_n and int(row.female_n) == female_n
    assert int(row.female_minus_male_days) == difference, population
cohorts = read_result("sex_survival_by_cohort")
assert len(cohorts) == 17 and (cohorts.female_advantage_days > 0).all()
assert cohorts.female_advantage_days.min() == 52
assert cohorts.female_advantage_days.max() == 135

# Extended Data Table 7, ITP control survival against the published field.
ceiling = read_result("control_lifespan_comparison").set_index(["specification", "sex"])
for specification, sex, itp_n, itp_median, field_n, field_median, reaching, p in [
    ("unique_control_groups", "Female", 18, 891.0, 17, 815, 0, 2.5e-06),
    ("unique_control_groups", "Male", 18, 800.5, 20, 814, 10, 5.1e-01),
    ("study_level_median", "Female", 13, 893.0, 15, 724, 0, 1.5e-05),
    ("study_level_median", "Male", 14, 795.5, 18, 794, 9, 8.6e-01),
    ("median_endpoints_only", "Female", 18, 891.0, 9, 847, 0, 2.9e-04),
    ("median_endpoints_only", "Male", 18, 800.5, 9, 912, 5, 2.9e-01),
    ("published_2000_or_later", "Female", 18, 891.0, 12, 842, 0, 3.6e-05),
    ("published_2000_or_later", "Male", 18, 800.5, 14, 835, 8, 1.5e-01),
]:
    row = ceiling.loc[(specification, sex)]
    assert int(row.itp_n) == itp_n and int(row.non_itp_n) == field_n, f"{specification} {sex}: n"
    assert row.itp_median_days == itp_median, f"{specification} {sex}: ITP median"
    assert round(row.non_itp_median_days) == field_median, f"{specification} {sex}: field"
    assert int(row.non_itp_reaching_itp_median) == reaching, f"{specification} {sex}: reaching"
    assert float(f"{row.p_mannwhitney:.1e}") == p, f"{specification} {sex}: P"
# Results: site consistency of the primary estimate and the inbred-strain sex comparison
sites = read_result("site_consistency")
female_sites = sites[(sites.sex == "f") & (sites.site != "heterogeneity")]
assert sorted(female_sites.HR.round(2)) == [1.01, 1.02, 1.05]
assert round(float(sites[(sites.sex == "f") & (sites.site == "heterogeneity")].p.iloc[0]), 2) == 0.62
male_sites = sites[(sites.sex == "m") & (sites.site != "heterogeneity")]
assert sorted(male_sites.HR.round(2)) == [0.94, 0.95, 0.95]
assert round(float(sites[(sites.sex == "m") & (sites.site == "heterogeneity")].p.iloc[0]), 2) == 0.99
strain_sex = read_result("inbred_strain_comparison").iloc[0]
assert int(strain_sex.strains) == 31 and strain_sex.strain_median_female_minus_male_days == 9
assert int(strain_sex.strains_above_itp_female_median) == 2 and strain_sex.itp_female_median_days == 886
assert strain_sex.itp_female_minus_male_days == 86
era = read_result("control_lifespan_era").set_index("quantity")
assert round(era.loc["non_itp_female_year_lifespan_spearman_rho", "value"], 2) == 0.71
assert round(era.loc["non_itp_female_year_lifespan_p", "value"], 3) == 0.001

# Extended Data Table 2, the unselected within-cohort weight-change analysis.
contrasts = pd.read_csv(DERIVED / "experiment_contrasts.csv")
assert len(contrasts) == 163, f"{len(contrasts)} experiment contrasts, expected 163"
assert contrasts.cohort.nunique() == 16
assert contrasts.age_init.max() < 18, "an experiment initiated at or after eighteen months was kept"
enalapril = contrasts.set_index(["cohort", "sex", "group"]).loc[("C2005", "m", "Enal")]
assert round(enalapril.wt_pct, 2) == 0.66, (
    f"C2005 male enalapril weight change {enalapril.wt_pct:.3f}; 0.17 means the implausible "
    "five-gram eighteen-month weight was retained")
quartiles = read_result("weight_response_quartiles").set_index(["sex", "quartile"])
for sex, quartile, experiments, weight, lifespan in [
    ("m", 1, 24, -4.8, 10.9), ("m", 2, 16, -0.5, 4.1),
    ("m", 3, 15, 0.4, 4.4), ("m", 4, 19, 1.5, 2.7),
    ("f", 1, 22, -11.3, 4.6), ("f", 2, 15, -4.0, 1.4),
    ("f", 3, 14, -0.5, 1.5), ("f", 4, 17, 2.1, 0.6),
]:
    row = quartiles.loc[(sex, quartile)]
    assert int(row.experiments) == experiments, f"{sex} Q{quartile}: experiments {int(row.experiments)} != {experiments}"
    assert round(row.median_weight_change_pct, 1) == weight, f"{sex} Q{quartile}: weight"
    assert round(row.mean_lifespan_change_pct, 1) == lifespan, f"{sex} Q{quartile}: lifespan"
trend = read_result("weight_change_trend").set_index("sex")
assert int(trend.loc["m", "experiments"]) == 74 and int(trend.loc["m", "cohorts"]) == 13
assert int(trend.loc["f", "experiments"]) == 68 and int(trend.loc["f", "cohorts"]) == 12
assert round(trend.loc["m", "trend_per_quartile_pp"], 1) == -2.6
assert round(trend.loc["m", "lo"], 1) == -4.6 and round(trend.loc["m", "hi"], 1) == -0.5
assert round(trend.loc["m", "p"], 3) == 0.018
assert round(trend.loc["f", "trend_per_quartile_pp"], 1) == -1.3
assert round(trend.loc["f", "lo"], 1) == -2.8 and round(trend.loc["f", "hi"], 1) == 0.3
assert round(trend.loc["f", "p"], 3) == 0.098

slopes = read_result("weight_change_slopes").set_index(["sex", "specification"])
assert round(slopes.loc[("m", "pooled"), "slope_pct_lifespan_per_pct_weight"], 2) == -0.71
assert round(slopes.loc[("m", "within_cohort"), "slope_pct_lifespan_per_pct_weight"], 2) == -0.65
assert round(slopes.loc[("m", "within_cohort"), "p"], 3) == 0.008
assert round(slopes.loc[("m", "within_cohort"), "lo"], 2) == -1.1
assert round(slopes.loc[("m", "within_cohort"), "hi"], 2) == -0.2
assert slopes.loc[("f", "pooled"), "p"] > 0.05 and slopes.loc[("f", "within_cohort"), "p"] > 0.05

joint = read_result("weight_change_vs_comparator")
joint = joint.set_index(["sex", "model", "term"])
alone = joint.loc[("m", "weight_change_alone", "wt_pct"), "estimate"]
both = joint.loc[("m", "both_together", "wt_pct"), "estimate"]
assert round(alone, 3) == 0.022 and round(both, 3) == 0.021
assert joint.loc[("m", "both_together", "wt_pct"), "p"] < 1e-4
assert round(joint.loc[("m", "predictor_correlation", "spearman_rho"), "estimate"], 2) == -0.12
assert int(joint.loc[("m", "both_together", "wt_pct"), "cells"]) == 246
female = joint.loc[("f", "both_together", "wt_pct")]
assert int(female.cells) == 237 and 246 + int(female.cells) == 483
assert round(female.estimate, 3) == 0.011
assert round(female.lo, 3) == 0.004 and round(female.hi, 3) == 0.019
assert round(female.p, 3) == 0.004, f"female adjusted P {female.p:.4f}"

row = slopes.loc[("m", "within_cohort")]
df = int(row.cohorts) - 1
ratio = abs(row.slope_pct_lifespan_per_pct_weight / row.se)
assert abs(row.p - 2 * t_dist.sf(ratio, df)) < 1e-9, "slope P is not on a t reference"
assert abs(row.p - 2 * normal.sf(ratio)) > 1e-3, "slope P looks like a normal approximation"
assert abs(row.hi - row.slope_pct_lifespan_per_pct_weight
           - t_dist.ppf(0.975, df) * row.se) < 1e-9, "slope interval is not on a t reference"

workbook = root / "supplementary" / "Supplementary_Data.xlsx"
for sheet, output, key in (("SD6_experiment_contrasts", "experiment_contrasts.csv", ["cohort", "sex", "group"]),
                           ("SD4_stratum_correlations", "bw6_strata.csv",
                            ["cohort", "sex", "group"])):
    shipped = pd.read_excel(workbook, sheet_name=sheet).set_index(key).sort_index()
    regenerated = pd.read_csv(DERIVED / output).set_index(key).sort_index()
    assert shipped.index.equals(regenerated.index), f"{sheet}: rows differ from {output}"
    for column in shipped.columns:
        gap = (shipped[column] - regenerated[column]).abs().max()
        assert gap < 1e-9, f"{sheet}.{column} differs from {output} by {gap}"
quartile_sheet = pd.read_excel(workbook, sheet_name="SD3_matched_quartiles", header=None)
matched = read_result("matched_quartile_survival")
matched_trend = read_result("quartile_trend").iloc[0]
gap = abs(quartile_sheet.iloc[4:8, 7].astype(float).to_numpy() - matched.gain_pct.to_numpy()).max()
assert gap < 1e-9, "matched-quartile gains (Supplementary Data 3) are stale"
assert int(quartile_sheet.iloc[10, 1]) == int(matched_trend.experiments)
assert int(quartile_sheet.iloc[10, 3]) == int(matched_trend.unique_mice)

# Extended Data Table 5 and Results: the external olfactory-exposure cohort.
olfactory = read_result("olfactory_weight_models")
olfactory = olfactory.set_index(["analysis", "sex"])
for key, (n, estimate, lo, hi, p) in {
        ("control_6m", "F"): (59, 1.52, 1.14, 2.02, 0.0042),
        ("control_6m", "M"): (60, 1.18, 0.91, 1.52, 0.21),
        ("females_6w_adjusted_for_6m", "F"): (170, 1.01, 0.82, 1.24, 0.92),
        ("females_6m_adjusted_for_6w", "F"): (170, 1.23, 1.01, 1.50, 0.041)}.items():
    row = olfactory.loc[key]
    assert int(row.n) == n, f"{key}: n"
    assert (round(row.HR_per_SD, 2), round(row.lo, 2), round(row.hi, 2)) == (estimate, lo, hi), key
    assert float(f"{row.p:.2g}") == p, f"{key}: P {row.p:.4g}"
assert round(olfactory.loc[("control_dead_by_day_500", "F"), "dead_by_day_500_pct"], 1) == 32.1
assert round(olfactory.loc[("control_dead_by_day_500", "M"), "dead_by_day_500_pct"], 1) == 13.3

# Supplementary Data 1: every treated group, with its publication; 8: olfactory-cohort interactions
experiments = pd.read_csv(out / "experiment_table.csv")
primary_rows = experiments[experiments.in_analysis == "yes"]
counts = primary_rows.groupby("sex").agg(n=("code", "size"), lowering=("weight_lowering", lambda s: (s == "yes").sum()),
                                         predosing=("weight_before_dosing", lambda s: (s == "yes").sum()),
                                         matched=("matched_quartile_analysis", lambda s: (s == "yes").sum()))
assert counts.loc["male"].tolist() == [92, 42, 60, 30] and counts.loc["female", ["n", "lowering"]].tolist() == [88, 40]
assert len(experiments) == 182 and experiments.itp_publication.notna().all() and experiments.pmid.nunique() == 17
olf_int = read_result("olfactory_interactions")
assert len(olf_int) == 12 and (olf_int.interaction_p > 0.05).all()
six = olf_int[(olf_int.weight_age == "6 months") & (olf_int.sex == "Female") & (olf_int.exposure == "Adult-female odour")].iloc[0]
assert (round(six.interaction_hr, 2), round(six.interaction_lo, 2), round(six.interaction_hi, 2)) == (0.71, 0.49, 1.04)
book = root / "supplementary" / "Supplementary_Data.xlsx"
# Supplementary Data 8 shows P to two significant figures
shown = olf_int.assign(interaction_p=olf_int.interaction_p.map(lambda p: float(f"{p:.2g}")))
# SD1 uses short citations; the experiment table retains the full references.
citation = r"^(.+? et al\.) .*\. ([^.,]+?) \d+, [^(]+ (\(\d{4}\))\.$"
sd1 = experiments.assign(itp_publication=experiments.itp_publication.str.replace(citation, r"\1 \2 \3", regex=True))
for sheet, frame in (("SD1_ITP_experiments", sd1), ("SD8_olfactory_interactions", shown)):
    s = pd.read_excel(book, sheet_name=sheet)
    assert s.shape == frame.shape, f"{sheet} is stale"
    frame = frame[list(s.columns)]
    assert (s.round(3).astype(str).replace({"nan": "", "<NA>": ""}).values ==
            frame.round(3).astype(str).replace({"nan": "", "<NA>": ""}).values).all(), f"{sheet} is stale"

print("Verified the primary, class, quartile, compound, comparator, control-gradient, "
      "developmental, strain-context, external-cohort, control-survival and "
      "weight-change results.")
