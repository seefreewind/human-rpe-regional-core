#!/usr/bin/env python3
from pathlib import Path
import re
import numpy as np
import pandas as pd
import statsmodels.api as sm
from statsmodels.stats.multitest import multipletests

PROJECT = Path(__file__).resolve().parents[2]
META = PROJECT / "data/metadata/CELLXGENE_OBS_METADATA_RPE.parquet"
OUT = PROJECT / "results/pivot_phase2"
OUT.mkdir(parents=True, exist_ok=True)

def parse_age(x):
    m = re.search(r"(\d+(?:\.\d+)?)", str(x))
    return float(m.group(1)) if m else np.nan

obs = pd.read_parquet(META)
obs = obs[obs["reported_diseases"].astype(str).str.lower().eq("none")].copy()
obs["region"] = obs["Region"].map({"Macular": "MACULA", "Peripheral": "PERIPHERY"}).fillna(obs["Region"].astype(str).str.upper())
obs["age"] = obs["donor_age"].map(parse_age)
obs = obs[obs["region"].isin(["MACULA", "PERIPHERY"]) & obs["author_cell_type"].astype(str).str.startswith("RPE_")].copy()

cell_counts = obs.groupby(["donor_id", "region"]).size().rename("n_RPE_cells")
keep = cell_counts[cell_counts >= 30].index
obs = obs.set_index(["donor_id", "region"]).loc[keep].reset_index()
sample_meta = obs.groupby(["donor_id", "region"], as_index=False).agg(
    age=("age", "first"),
    sex=("sex", "first"),
    ancestry_ethnicity=("self_reported_ethnicity", "first"),
    n_RPE_cells=("author_cell_type", "size"),
)
counts = pd.crosstab([obs["donor_id"], obs["region"]], obs["author_cell_type"])
counts = counts.reindex(index=pd.MultiIndex.from_frame(sample_meta[["donor_id", "region"]]), fill_value=0)
labels = sorted(set(obs["author_cell_type"]))
for label in labels:
    if label not in counts.columns:
        counts[label] = 0
counts = counts[labels]
props = counts.div(counts.sum(axis=1), axis=0)
props.columns = [f"prop_{x}" for x in props.columns]
donor_region = sample_meta.set_index(["donor_id", "region"]).join(counts).join(props).reset_index()
donor_region["age_centered"] = donor_region["age"] - donor_region["age"].median()
donor_region.to_csv(OUT / "RPE_SUBTYPE_COMPOSITION_DONOR_REGION.tsv", sep="\t", index=False)

rows = []
for label in labels:
    prop_col = f"prop_{label}"
    d = donor_region[[prop_col, "region", "age_centered"]].copy()
    d["region_macula"] = (d["region"] == "MACULA").astype(int)
    X = sm.add_constant(d[["region_macula", "age_centered"]])
    X["region_macula_age_centered"] = d["region_macula"] * d["age_centered"]
    fit = sm.OLS(d[prop_col].astype(float), X).fit()
    rows.append({
        "subtype": label,
        "n_donor_regions": len(d),
        "macula_mean_proportion": float(d.loc[d.region_macula == 1, prop_col].mean()),
        "periphery_mean_proportion": float(d.loc[d.region_macula == 0, prop_col].mean()),
        "beta_region_at_median_age": float(fit.params.get("region_macula", np.nan)),
        "p_region": float(fit.pvalues.get("region_macula", np.nan)),
        "beta_age_periphery": float(fit.params.get("age_centered", np.nan)),
        "p_age_periphery": float(fit.pvalues.get("age_centered", np.nan)),
        "beta_age_region_interaction": float(fit.params.get("region_macula_age_centered", np.nan)),
        "p_age_region_interaction": float(fit.pvalues.get("region_macula_age_centered", np.nan)),
        "r_squared": float(fit.rsquared),
    })
summary = pd.DataFrame(rows)
for pcol, qcol in [("p_region", "fdr_region"), ("p_age_periphery", "fdr_age_periphery"), ("p_age_region_interaction", "fdr_age_region_interaction")]:
    valid = summary[pcol].notna()
    summary.loc[valid, qcol] = multipletests(summary.loc[valid, pcol], method="fdr_bh")[1]
summary["descriptive_region_shift_fdr05"] = summary["fdr_region"] < 0.05
summary["descriptive_age_shift_fdr05"] = summary["fdr_age_periphery"] < 0.05
summary["descriptive_age_region_shift_fdr05"] = summary["fdr_age_region_interaction"] < 0.05
summary.to_csv(OUT / "RPE_SUBTYPE_COMPOSITION_AUDIT.tsv", sep="\t", index=False)
print("donor-region rows:", len(donor_region))
print(summary.to_string(index=False))

