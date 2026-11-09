"""One row per ITP experiment (Supplementary Data 1)."""
from pathlib import Path

import pandas as pd

from common import DERIVED

root = Path(__file__).absolute().parents[1]
inputs, out = root / "data" / "inputs", root / "data" / "outputs"

interactions = pd.read_csv(DERIVED / "experiment_interactions.csv")
long = pd.read_csv(DERIVED / "within_experiment_long.csv",
                   usecols=["sex", "experimentsex", "stratum", "mouse", "Tx", "dead"])
weights = pd.read_csv(DERIVED / "experiment_weight_effects.csv")
matched = pd.read_csv(DERIVED / "quartile_eligibility.csv")
references = pd.read_csv(inputs / "itp_experiment_references.csv")
mice = pd.read_excel(inputs / "ITP_bodyweight_all_cohorts.xlsx",
                     usecols=["cohort", "site", "sex", "group", "dose", "age_initiation", "bw_6"])
mice["sex"] = mice.sex.str.lower().str[0]

treated = mice[mice.group != "Control"]
groups = treated.groupby(["cohort", "group", "sex"], as_index=False).agg(
    adult_weights=("bw_6", lambda s: int(s.notna().sum())))
t = groups.merge(interactions, on=["cohort", "group", "sex"], how="left", validate="one_to_one")
t["in_primary"] = t.int_hr.notna()
t["note"] = ""
t.loc[~t.in_primary & t.adult_weights.eq(0), "note"] = "No recorded six-month weight"
assert (t.in_primary | t.note.ne("")).all(), "a treated group is outside the primary analysis for an unrecorded reason"

long["site"] = long.stratum.str.split("|").str[0]
counts = long.groupby("experimentsex").apply(lambda g: pd.Series({
    "sites": ", ".join(sorted(g.site.unique())),
    "treated_mice": int((g.Tx == 1).sum()), "control_mice": int((g.Tx == 0).sum()),
    "treated_deaths": int(g.dead[g.Tx == 1].sum()), "control_deaths": int(g.dead[g.Tx == 0].sum())}),
    include_groups=False).reset_index()
t["experimentsex"] = t.cohort + "|" + t.group + "|" + t.sex
t = t.merge(counts, on="experimentsex", how="left", validate="one_to_one")
site_list = treated.groupby(["cohort", "group", "sex"]).site.agg(lambda s: ", ".join(sorted(s.unique()))).rename("all_sites")
t = t.merge(site_list.reset_index(), on=["cohort", "group", "sex"], how="left")
t["sites"] = t.sites.fillna(t.all_sites)

post = weights[weights.post_dosing.astype(str).str.lower().eq("true")].copy()
post["weighted_diff"] = post.diff_g * post.n_tx
lowering = post.groupby(["cohort", "group"], as_index=False).agg(weighted_diff=("weighted_diff", "sum"),
                                                                n_tx=("n_tx", "sum"))
lowering["weight_lowering"] = lowering.weighted_diff / lowering.n_tx < -1
init = weights[["cohort", "group", "init"]].drop_duplicates()
dose = treated.groupby(["cohort", "group"], as_index=False).agg(
    dose_ppm=("dose", lambda s: "; ".join(sorted(set(map(str, s.dropna()))))),
    starting_age_months=("age_initiation", lambda s: "; ".join(sorted(set(map(str, s.dropna()))))))
t = (t.merge(lowering[["cohort", "group", "weight_lowering"]], on=["cohort", "group"], how="left")
      .merge(init, on=["cohort", "group"], how="left", validate="many_to_one")
      .merge(dose, on=["cohort", "group"], how="left", validate="many_to_one")
      .merge(references, on=["cohort", "group"], how="left", validate="many_to_one"))
primary = t.in_primary
assert (t.treated_mice[primary] + t.control_mice[primary] == t.n[primary]).all(), "model counts do not add up"
assert t.reference.notna().all() and t.dose_ppm.notna().all(), "missing reference or dose"
yes_no = lambda flag: flag.map({True: "yes", False: "no"})
matched_keys = matched.set_index(["cohort", "group"]).index


def blank(series):
    return series.where(primary)


table = pd.DataFrame({
    "sex": t.sex.map({"m": "male", "f": "female"}), "cohort": t.cohort, "intervention": t.intervention,
    "code": t.group, "dose_ppm": t.dose_ppm, "starting_age_months": t.starting_age_months, "sites": t.sites,
    "in_analysis": yes_no(primary), "note": t.note,
    "treated_mice": blank(t.treated_mice).astype("Int64"), "control_mice": blank(t.control_mice).astype("Int64"),
    "treated_deaths": blank(t.treated_deaths).astype("Int64"), "control_deaths": blank(t.control_deaths).astype("Int64"),
    "follow_up_start_months": blank(t.entry_months).round().astype("Int64"),
    "weight_lowering": yes_no(t.weight_lowering.eq(True)).where(primary, ""),
    "weight_before_dosing": yes_no(t.init.notna() & (t.weight_age < t.init)).where(primary, ""),
    "matched_quartile_analysis": yes_no((t.sex == "m") & t.set_index(["cohort", "group"]).index.isin(matched_keys)).where(primary, ""),
    "interaction_hr": t.int_hr.round(3), "interaction_lo": t.int_lo.round(3), "interaction_hi": t.int_hi.round(3),
    "interaction_p": t.int_p.map(lambda p: float(f"{p:.2g}") if pd.notna(p) else p),
    "itp_publication": t.reference, "pmid": t.pmid, "mpd_project": t.mpd_project})
table = table.sort_values(["in_analysis", "sex", "cohort", "intervention", "code"],
                          ascending=[False, False, True, True, True])
table.to_csv(out / "experiment_table.csv", index=False, lineterminator="\n")
inn = table[table.in_analysis == "yes"].sex.value_counts().to_dict()
print(f"experiment_table.csv: {len(table)} treated groups; primary analysis {inn.get('male', 0)} male and "
      f"{inn.get('female', 0)} female experiments; {table.pmid.nunique()} publications")
