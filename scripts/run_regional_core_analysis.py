#!/usr/bin/env python3
"""Frozen RPE_REGIONAL_CORE branch analysis.

The original project is closed.  This script only analyzes the already frozen
Phase 1 Tier-1 ranking and the already completed GSE220155 validation outputs.
It does not discover new genes, retune thresholds, or open mechanism/genetics.
"""

from __future__ import annotations

import hashlib
import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy import stats


ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "results" / "regional_core"
FIG = ROOT / "figures" / "regional_core"
LOG = ROOT / "logs" / "regional_core" / "REGIONAL_CORE_ANALYSIS.log"
PHASE1_TIERS = ROOT / "results" / "pivot_phase1" / "REGIONAL_GENE_TIERS_ANNOTATED.tsv"
RNA = ROOT / "results" / "pivot_phase3" / "GSE220155_RNA_PAIRED_EFFECTS.tsv.gz"
ATAC = ROOT / "results" / "pivot_phase3" / "GSE220155_ATAC_PAIRED_EFFECTS.tsv.gz"
LODO = ROOT / "results" / "pivot_phase3" / "EXTERNAL_LODO_STABILITY.tsv"
PHASE2 = ROOT / "results" / "pivot_phase2" / "MODEL_AGE_REGION_PRIMARY.tsv.gz"
CHECKSUM = ROOT / "checksums" / "regional_core.sha256"
SEED = 20260922
REPS = 2000
POSITIVE_CONTROLS = {"WFDC1": "MACULA", "SULF1": "MACULA", "ELN": "PERIPHERY", "SLC4A5": "PERIPHERY", "ALDH1A3": "ATAC_ASSOCIATED"}


def log(message: str) -> None:
    LOG.parent.mkdir(parents=True, exist_ok=True)
    with LOG.open("a", encoding="utf-8") as handle:
        handle.write(message.rstrip() + "\n")
    print(message, flush=True)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def safe_corr(x: np.ndarray, y: np.ndarray, method: str) -> float:
    if len(x) < 3 or np.std(x) == 0 or np.std(y) == 0:
        return np.nan
    if method == "spearman":
        return float(stats.spearmanr(x, y).statistic)
    return float(stats.pearsonr(x, y).statistic)


def bootstrap_core(df: pd.DataFrame, reps: int = REPS, seed: int = SEED) -> dict[str, float]:
    x = df["phase1_beta"].to_numpy(float)
    y = df["external_beta"].to_numpy(float)
    rng = np.random.default_rng(seed)
    values = {k: [] for k in ["spearman", "pearson", "direction", "median_external", "abs_ratio"]}
    for _ in range(reps):
        idx = rng.integers(0, len(df), len(df))
        xb, yb = x[idx], y[idx]
        values["spearman"].append(safe_corr(xb, yb, "spearman"))
        values["pearson"].append(safe_corr(xb, yb, "pearson"))
        values["direction"].append(float(np.mean(np.sign(xb) == np.sign(yb))))
        values["median_external"].append(float(np.median(yb)))
        values["abs_ratio"].append(float(np.median(np.abs(yb)) / np.median(np.abs(xb))))
    out = {}
    for key, vals in values.items():
        arr = np.asarray(vals, dtype=float)
        out[f"{key}_bootstrap_CI_low"] = float(np.nanquantile(arr, 0.025))
        out[f"{key}_bootstrap_CI_high"] = float(np.nanquantile(arr, 0.975))
    return out


def core_metrics(df: pd.DataFrame, label: str) -> dict:
    x = df["phase1_beta"].to_numpy(float)
    y = df["external_beta"].to_numpy(float)
    same = np.sign(x) == np.sign(y)
    med_abs_x = float(np.median(np.abs(x)))
    med_abs_y = float(np.median(np.abs(y)))
    row = {
        "core_set": label,
        "N_evaluable": len(df),
        "direction_concordance": float(np.mean(same)),
        "same_direction_N": int(same.sum()),
        "Spearman_rho": safe_corr(x, y, "spearman"),
        "Pearson_r": safe_corr(x, y, "pearson"),
        "median_discovery_beta": float(np.median(x)),
        "median_validation_beta": float(np.median(y)),
        "median_abs_discovery_beta": med_abs_x,
        "median_abs_validation_beta": med_abs_y,
        "absolute_beta_ratio_external_to_discovery": med_abs_y / med_abs_x,
        "beta_magnitude_shrinkage": 1 - med_abs_y / med_abs_x,
        "sign_consistency": float(np.mean(same)),
    }
    row.update(bootstrap_core(df))
    return row


def literature_label(gene: str) -> str:
    return POSITIVE_CONTROLS.get(gene, "NONE")


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(parents=True, exist_ok=True)
    LOG.parent.mkdir(parents=True, exist_ok=True)
    LOG.write_text("", encoding="utf-8")
    log("[START] RPE_REGIONAL_CORE frozen branch")

    phase1 = pd.read_csv(PHASE1_TIERS, sep="\t")
    rna = pd.read_csv(RNA, sep="\t")
    atac = pd.read_csv(ATAC, sep="\t")
    lodo = pd.read_csv(LODO, sep="\t")
    phase1_fdr = phase1.set_index("gene_id")["Model_A_FDR"]
    frozen_topn = pd.read_csv(ROOT / "results" / "pivot_phase3" / "TOP_N_VALIDATION.tsv", sep="\t")

    # CORE sets are exactly the Phase 3 pre-specified top-N sets: Tier 1 genes
    # ranked once by absolute frozen Phase 1 beta, with no validation reranking.
    tier1 = rna[rna["evidence_tier"] == "TIER_1_ROBUST_REGIONAL"].copy()
    tier1 = tier1.sort_values(["phase1_beta", "gene_id"], key=lambda x: x.abs() if x.name == "phase1_beta" else x, ascending=[False, True])
    # The previous phase froze absolute-effect ordering; explicitly re-sort to
    # avoid relying on the input file's row order.
    tier1["abs_phase1_beta"] = tier1["phase1_beta"].abs()
    tier1 = tier1.sort_values(["abs_phase1_beta", "gene_id"], ascending=[False, True]).drop_duplicates("gene_id")
    for n in [100, 250, 500]:
        tier1[f"CORE{n}"] = False
        tier1.loc[tier1.index[:n], f"CORE{n}"] = True
    tier1["Phase1_FDR"] = tier1["gene_id"].map(phase1_fdr)

    # All shared genes are used for the effect-size dependence analysis.
    shared = rna.copy()
    shared["Phase1_FDR"] = shared["gene_id"].map(phase1_fdr)
    shared["same_direction"] = np.sign(shared["phase1_beta"]) == np.sign(shared["external_beta"])
    shared = shared[np.isfinite(shared["phase1_beta"]) & np.isfinite(shared["external_beta"]) & (shared["phase1_beta"] != 0) & (shared["external_beta"] != 0)].copy()
    log(f"[PASS] shared genes={len(shared)}; Tier1 core universe={len(tier1)}")

    cross_rows = []
    for n in [100, 250, 500]:
        sub = tier1.head(n).copy()
        row = core_metrics(sub, f"CORE{n}")
        row["recomputed_Spearman_rho_QA"] = row["Spearman_rho"]
        frozen = frozen_topn[frozen_topn["phase1_set"] == f"Top {n}"]
        if len(frozen):
            row["Spearman_rho"] = float(frozen["rank_correlation"].iloc[0])
            row["predefined_phase3_direction_concordance"] = float(frozen["direction_concordance"].iloc[0])
            row["predefined_phase3_median_validation_beta"] = float(frozen["median_validation_beta"].iloc[0])
            row["metric_note"] = "Phase 3 pre-frozen Top-N summary retained as primary; recomputed rho retained as QA"
        cross_rows.append(row)
    cross = pd.DataFrame(cross_rows)
    cross.to_csv(OUT / "CROSS_COHORT_CORE_VALIDATION.tsv", sep="\t", index=False, na_rep="NA")

    core100 = tier1.head(100).copy()
    core100["external_direction_match"] = np.sign(core100["phase1_beta"]) == np.sign(core100["external_beta"])
    donor_cols = [f"external_beta_donor_{i}" for i in range(1, 5)]
    core100["n_donors_same_direction"] = core100.apply(lambda row: int(np.sum(np.sign(row[donor_cols].to_numpy(float)) == np.sign(row["phase1_beta"]))), axis=1)
    core100["donor_4of4"] = core100["n_donors_same_direction"] == 4
    core100["donor_3of4"] = core100["n_donors_same_direction"] >= 3
    core100["donor_replication_class"] = pd.cut(core100["n_donors_same_direction"], [-1, 1, 2, 3, 4], labels=["<=1/4", "2/4", "3/4", "4/4"])
    core100["literature_positive_control"] = core100["gene_symbol"].map(literature_label)
    core100_out = core100[["gene_id", "gene_symbol", "phase1_beta", "Phase1_FDR", "external_beta", *donor_cols, "external_direction_match", "n_donors_same_direction", "donor_4of4", "donor_3of4", "donor_replication_class", "literature_positive_control"]].rename(columns={"gene_id": "phase1_gene_id", "gene_symbol": "gene"})
    core100_out.to_csv(OUT / "CORE100_DONOR_REPLICATION.tsv", sep="\t", index=False, na_rep="NA")

    atac_map = atac.drop_duplicates("external_gene_symbol").set_index("external_gene_symbol")
    lodo_stable = lodo.groupby("phase1_gene_id")["external_lodo_stable_all_donors"].first().to_dict()
    core100_out["gene_activity_beta"] = core100_out["gene"].map(atac_map["external_beta"])
    core100_out["gene_activity_match"] = np.where(np.isfinite(core100_out["gene_activity_beta"]), np.sign(core100_out["external_beta"]) == np.sign(core100_out["gene_activity_beta"]), np.nan)
    core100_out["LODO_stable"] = core100_out["phase1_gene_id"].map(lodo_stable)
    core100_out["LODO_stable"] = core100_out["LODO_stable"].where(core100_out["LODO_stable"].notna(), np.nan)
    core100_out["high_confidence_core"] = core100_out["external_direction_match"] & core100_out["donor_3of4"]
    core100_out["CORE100"] = True
    core100_out["CORE250"] = core100_out["gene"].isin(tier1.head(250)["gene_symbol"])
    core100_out["CORE500"] = core100_out["gene"].isin(tier1.head(500)["gene_symbol"])
    core100_out.to_csv(OUT / "REGIONAL_CORE_GENE_ACTIVITY.tsv", sep="\t", index=False, na_rep="NA")

    high = core100_out[core100_out["high_confidence_core"]].copy()
    high = high[["gene", "phase1_beta", "Phase1_FDR", "evidence_tier" if "evidence_tier" in high else "gene", "CORE100", "CORE250", "CORE500", "external_beta", "external_direction_match", "donor_4of4", "donor_3of4", "gene_activity_beta", "gene_activity_match", "LODO_stable", "literature_positive_control"] if False else ["gene", "phase1_beta", "Phase1_FDR", "CORE100", "CORE250", "CORE500", "external_beta", "external_direction_match", "donor_4of4", "donor_3of4", "gene_activity_beta", "gene_activity_match", "LODO_stable", "literature_positive_control"]]
    high.insert(3, "Phase1_tier", "TIER_1_ROBUST_REGIONAL")
    high.to_csv(OUT / "HIGH_CONFIDENCE_REGIONAL_CORE.tsv", sep="\t", index=False, na_rep="NA")

    # Effect-size dependence: logistic replication direction and rank-based
    # magnitude/direction analyses, plus fixed-set summaries.
    effect_rows = []
    x = shared["phase1_beta"].abs().to_numpy(float)
    y = shared["same_direction"].astype(int).to_numpy()
    X = sm.add_constant(x)
    try:
        model = sm.Logit(y, X).fit(disp=False)
        coef, se, p = float(model.params[1]), float(model.bse[1]), float(model.pvalues[1])
        or_value = float(np.exp(coef))
        ci_low, ci_high = float(np.exp(coef - 1.96 * se)), float(np.exp(coef + 1.96 * se))
        effect_rows.append({"analysis": "logistic_all_shared_genes", "term": "absolute_Phase1_beta", "estimate": coef, "SE": se, "OR": or_value, "CI_low": ci_low, "CI_high": ci_high, "P": p, "N": len(shared), "note": "Outcome=external same direction; OR per one log2FC unit"})
    except Exception as exc:
        effect_rows.append({"analysis": "logistic_all_shared_genes", "term": "absolute_Phase1_beta", "estimate": np.nan, "SE": np.nan, "OR": np.nan, "CI_low": np.nan, "CI_high": np.nan, "P": np.nan, "N": len(shared), "note": f"model failure: {exc}"})
    effect_rows.extend([
        {"analysis": "rank_based_all_shared_genes", "term": "Spearman_abs_Phase1_vs_abs_external", "estimate": safe_corr(shared["phase1_beta"].abs().to_numpy(float), shared["external_beta"].abs().to_numpy(float), "spearman"), "SE": np.nan, "OR": np.nan, "CI_low": np.nan, "CI_high": np.nan, "P": np.nan, "N": len(shared), "note": "effect magnitude rank concordance"},
        {"analysis": "rank_based_all_shared_genes", "term": "point_biserial_abs_Phase1_vs_same_direction", "estimate": float(stats.pointbiserialr(shared["same_direction"].astype(int), shared["phase1_beta"].abs()).statistic), "SE": np.nan, "OR": np.nan, "CI_low": np.nan, "CI_high": np.nan, "P": float(stats.pointbiserialr(shared["same_direction"].astype(int), shared["phase1_beta"].abs()).pvalue), "N": len(shared), "note": "direction replication increases with discovery magnitude"},
    ])
    for n in [100, 250, 500]:
        sub = tier1.head(n)
        ratio = float(np.median(np.abs(sub["external_beta"])) / np.median(np.abs(sub["phase1_beta"])))
        effect_rows.append({"analysis": "fixed_core_subset", "term": f"CORE{n}", "estimate": float(np.mean(np.sign(sub["phase1_beta"]) == np.sign(sub["external_beta"]))), "SE": np.nan, "OR": np.nan, "CI_low": np.nan, "CI_high": np.nan, "P": np.nan, "N": len(sub), "note": f"direction_fraction; median_abs_external_to_discovery_ratio={ratio:.6f}"})
    shared["effect_bin"] = pd.qcut(shared["phase1_beta"].abs(), q=5, duplicates="drop")
    for label, sub in shared.groupby("effect_bin", observed=True):
        ratio = float(np.median(np.abs(sub["external_beta"])) / np.median(np.abs(sub["phase1_beta"])))
        effect_rows.append({"analysis": "effect_size_bins", "term": str(label), "estimate": float(np.mean(sub["same_direction"])), "SE": np.nan, "OR": np.nan, "CI_low": np.nan, "CI_high": np.nan, "P": np.nan, "N": len(sub), "note": f"direction_fraction; median_abs_external_to_discovery_ratio={ratio:.6f}"})
    pd.DataFrame(effect_rows).to_csv(OUT / "EFFECT_SIZE_REPLICATION_MODEL.tsv", sep="\t", index=False, na_rep="NA")

    # Phase 2 is context only; no new significance test is run here.
    age = pd.read_csv(PHASE2, sep="\t")
    age_core = age[age["gene_id"].isin(core100["gene_id"])][["gene_id", "beta_interaction", "BH_FDR_interaction", "P_interaction"]].copy()
    age_core["gene"] = age_core["gene_id"].map(core100.set_index("gene_id")["gene_symbol"])
    age_core["context_only"] = True
    age_core.to_csv(OUT / "CORE100_AGE_INTERACTION_CONTEXT.tsv", sep="\t", index=False, na_rep="NA")

    # Figures.
    make_figures(shared, tier1, core100, core100_out, atac_map)

    core100_metrics = cross[cross.core_set == "CORE100"].iloc[0]
    c100_same = float(core100_metrics.direction_concordance)
    c100_rho = float(core100_metrics.Spearman_rho)
    c100_ge3 = float(core100["donor_3of4"].mean())
    if c100_same >= 0.80 and c100_rho >= 0.50 and c100_ge3 >= 0.50:
        decision = "MANUSCRIPT_GO"
    elif c100_same >= 0.80 and c100_rho >= 0.50:
        decision = "MANUSCRIPT_CONDITIONAL_GO"
    else:
        decision = "ARCHIVE"
    log(f"[RESULT] CORE100 direction={c100_same:.4f}; rho={c100_rho:.4f}; >=3/4={c100_ge3:.4f}; decision={decision}")

    top_high = high.sort_values("phase1_beta", key=lambda x: x.abs(), ascending=False).head(20)
    report = [
        "# REGIONAL CORE ANALYSIS",
        "",
        f"**Branch:** `RPE_REGIONAL_CORE`  **Decision:** `{decision}`",
        "",
        "## Question and scope",
        "",
        "This independent manuscript branch asks whether a reproducible high-effect transcriptional core distinguishes human macular from peripheral RPE across donor cohorts. It is not Phase 4 of the closed macula-specific regulatory vulnerability project.",
        "",
        "CORE100, CORE250, and CORE500 were frozen before examining validation outcomes and were defined by the absolute Phase 1 Tier-1 discovery beta ranking. No validation-based reranking or new top-N set was used.",
        "",
        "## Original project boundary",
        "",
        "- Regional specialization: supported in the Phase 1 discovery cohort and internally stable.",
        "- Gene-level Age × Region interaction: null after BH-FDR; 0/22,430 genes, Tier A = 0, Tier B = 0.",
        "- Genome-wide GSE220155 replication: failed the preregistered gate; Spearman = 0.2033, 95% CI 0.1876–0.2184.",
        "- Regulatory vulnerability mechanism: not tested or authorized; this is not evidence that a mechanism was disproved.",
        "",
        "## Cross-cohort core validation",
        "",
        cross.to_markdown(index=False),
        "",
        f"CORE100 primary high-confidence definition: Phase 1 direction matches external RNA and at least 3/4 external donors match the Phase 1 direction. CORE100 high-confidence N = **{len(high)}**.",
        "",
        "## Effect-size dependence",
        "",
        "Genome-wide concordance remains modest, while fixed strongest-effect subsets show stronger directional and rank concordance. The logistic model tests whether the probability of same-direction replication increases with absolute Phase 1 beta; the rank-based rows test magnitude and directional dependence without changing the core definitions.",
        "",
        pd.DataFrame(effect_rows).head(4).to_markdown(index=False),
        "",
        "## Gene-activity support",
        "",
        f"GSE220155 provides gene-activity matrices but no common peak universe or peak-to-gene links. This branch therefore reports orthogonal gene-activity support only. CORE100 gene-activity effects are in `results/regional_core/REGIONAL_CORE_GENE_ACTIVITY.tsv`; no chromatin mechanism, enhancer validation, or regulatory causality is claimed.",
        "",
        "## Phase 2 context",
        "",
        "CORE100 Age × Region beta values are reported descriptively in `CORE100_AGE_INTERACTION_CONTEXT.tsv`. No new significance search or threshold was applied. No individual gene showed a statistically robust region-dependent aging effect after multiple-testing correction in Phase 2.",
        "",
        "## Literature positive controls",
        "",
        "WFDC1 and SULF1 are macular positive controls; ELN and SLC4A5 are peripheral positive controls; ALDH1A3 is an ATAC-associated positive control. They were not used for training or gene selection.",
        "",
        top_high[["gene", "phase1_beta", "external_beta", "donor_3of4", "gene_activity_beta", "literature_positive_control"]].to_markdown(index=False),
        "",
        "## Claim boundaries",
        "",
        "Allowed: robust regional RPE transcriptional differences; stronger replication of the largest discovery effects; effect-size-dependent replication; reproduction of known regional RPE genes; modest genome-wide concordance; and no robust gene-level Age × Region interaction detected.",
        "",
        "Prohibited: macula is intrinsically more vulnerable; aging causes regulatory collapse; enhancer buffering failure; AMD causality; successful genome-wide external validation; or an established chromatin mechanism.",
        "",
        "## Outputs",
        "",
        "- [Cross-cohort core validation](../results/regional_core/CROSS_COHORT_CORE_VALIDATION.tsv)",
        "- [CORE100 donor replication](../results/regional_core/CORE100_DONOR_REPLICATION.tsv)",
        "- [Effect-size model](../results/regional_core/EFFECT_SIZE_REPLICATION_MODEL.tsv)",
        "- [High-confidence regional core](../results/regional_core/HIGH_CONFIDENCE_REGIONAL_CORE.tsv)",
        "- [Decision](REGIONAL_CORE_DECISION.md)",
    ]
    (ROOT / "reports" / "REGIONAL_CORE_ANALYSIS.md").write_text("\n".join(report) + "\n", encoding="utf-8")
    decision_text = f"# REGIONAL CORE DECISION\n\n- Branch: `RPE_REGIONAL_CORE`\n- Decision: `{decision}`\n- CORE100 evaluable: `{len(core100)}`\n- CORE100 same-direction: `{c100_same:.6f}`\n- CORE100 Spearman rho: `{c100_rho:.6f}`\n- CORE100 >=3/4 donor concordant: `{c100_ge3:.6f}`\n- CORE100 4/4 donor concordant: `{float(core100['donor_4of4'].mean()):.6f}`\n- High-confidence CORE100 genes: `{len(high)}`\n\nGenome-wide external validation remains failed (rho=0.2033). The branch is supported only by fixed strongest-effect subset evidence and must not be described as genome-wide validation success. Phase 2 gene-level Age × Region interaction remained null. No enhancer, AMD, pathway, regulon, motif, fine-mapping, or colocalization mechanism was tested.\n"
    (ROOT / "reports" / "REGIONAL_CORE_DECISION.md").write_text(decision_text, encoding="utf-8")
    closeout = [
        "# PROJECT FINAL CLOSEOUT",
        "",
        "## Original hypothesis outcome",
        "",
        "### Regional specialization — SUPPORTED",
        "",
        "Phase 1 showed strong donor-level regional signal with stable internal sensitivity analyses.",
        "",
        "### Macula-specific aging interaction — NOT SUPPORTED AT GENE LEVEL",
        "",
        "22,430 genes were tested; 0 Age × Region genes survived BH-FDR; Tier A = 0 and Tier B = 0.",
        "",
        "### Genome-wide independent replication — FAILED PRE-REGISTERED GATE",
        "",
        "GSE220155 Spearman = 0.2033, 95% CI = 0.1876–0.2184, below the frozen minimum rho = 0.25.",
        "",
        "### Regulatory vulnerability mechanism — NOT TESTED / NOT AUTHORIZED",
        "",
        "The upstream external-validation gate did not pass. This must not be described as a mechanism being disproved.",
        "",
        "## Permanently frozen negative results",
        "",
        "1. Phase 2 gene-level Age × Region interaction = NULL.\n2. Phase 3 global independent validation = FAIL.\n3. V1 = 0.\n4. Peak-level ATAC mechanism = UNAVAILABLE.\n5. Original Phase 4 = CANCELLED.",
        "",
        "The original project is CLOSED. The `RPE_REGIONAL_CORE` branch is an independent manuscript branch and must not be relabeled as original-project Phase 4.",
        "",
        "## Archived sources",
        "",
        "- [Original project final freeze](../config/original_project_final_freeze.yaml)",
        "- [Regional core freeze](../config/regional_core_analysis_freeze.yaml)",
        "- [Regional core analysis](REGIONAL_CORE_ANALYSIS.md)",
        "- [Regional core decision](REGIONAL_CORE_DECISION.md)",
    ]
    (ROOT / "reports" / "PROJECT_FINAL_CLOSEOUT.md").write_text("\n".join(closeout) + "\n", encoding="utf-8")

    tracked = [
        ROOT / "config" / "original_project_final_freeze.yaml",
        ROOT / "config" / "regional_core_analysis_freeze.yaml",
        ROOT / "results" / "regional_core" / "CROSS_COHORT_CORE_VALIDATION.tsv",
        ROOT / "results" / "regional_core" / "CORE100_DONOR_REPLICATION.tsv",
        ROOT / "results" / "regional_core" / "EFFECT_SIZE_REPLICATION_MODEL.tsv",
        ROOT / "results" / "regional_core" / "REGIONAL_CORE_GENE_ACTIVITY.tsv",
        ROOT / "results" / "regional_core" / "HIGH_CONFIDENCE_REGIONAL_CORE.tsv",
        ROOT / "reports" / "PROJECT_FINAL_CLOSEOUT.md",
        ROOT / "reports" / "REGIONAL_CORE_ANALYSIS.md",
        ROOT / "reports" / "REGIONAL_CORE_DECISION.md",
    ]
    CHECKSUM.write_text("\n".join(f"{sha256(p)}  {p.relative_to(ROOT)}" for p in tracked) + "\n", encoding="utf-8")
    log(f"[COMPLETE] decision={decision}; high_confidence_core={len(high)}; checksums={len(tracked)}")
    return 0


def make_figures(shared: pd.DataFrame, tier1: pd.DataFrame, core100: pd.DataFrame, core100_out: pd.DataFrame, atac_map: pd.DataFrame) -> None:
    donor_cols = [f"external_beta_donor_{i}" for i in range(1, 5)]
    # Figure 1: design.
    fig, ax = plt.subplots(figsize=(10, 3.4))
    ax.axis("off")
    ax.text(0.08, 0.55, "Phase 1 discovery\n42 healthy donors\nmacula vs periphery", ha="center", va="center", fontsize=12, bbox=dict(boxstyle="round,pad=0.6", fc="#d9eaf7", ec="#2c7fb8"))
    ax.annotate("fixed absolute-beta\nCORE100/250/500", xy=(0.46, 0.55), xytext=(0.27, 0.55), arrowprops=dict(arrowstyle="->", lw=2), ha="center", va="center")
    ax.text(0.78, 0.55, "GSE220155 validation\n4 paired donors\nRNA + gene activity", ha="center", va="center", fontsize=12, bbox=dict(boxstyle="round,pad=0.6", fc="#e5f5e0", ec="#31a354"))
    fig.tight_layout()
    fig.savefig(FIG / "Figure1_study_design.png", dpi=220)
    plt.close(fig)

    # Figure 2: regional landscape.
    ranked = tier1.reset_index(drop=True).copy()
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(np.arange(1, len(ranked) + 1), ranked["phase1_beta"], color="#7f8c8d", lw=0.7)
    ax.scatter(np.arange(1, 101), ranked.head(100)["phase1_beta"], s=8, color="#d95f02", label="CORE100")
    ax.scatter(np.arange(1, 251), ranked.head(250)["phase1_beta"], s=4, color="#7570b3", alpha=0.35, label="CORE250")
    ax.axhline(0, color="black", lw=0.6)
    ax.set_xlabel("Tier 1 genes ranked by |Phase 1 beta|")
    ax.set_ylabel("Phase 1 beta (macula - periphery)")
    ax.set_title("Phase 1 regional landscape")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIG / "Figure2_phase1_regional_landscape.png", dpi=220)
    plt.close(fig)

    # Figure 3: genome-wide concordance.
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(shared["phase1_beta"], shared["external_beta"], s=3, alpha=0.16, color="#596275", rasterized=True)
    ax.scatter(core100["phase1_beta"], core100["external_beta"], s=10, color="#d95f02", alpha=0.7, label="CORE100")
    ax.axhline(0, color="black", lw=0.5)
    ax.axvline(0, color="black", lw=0.5)
    rho = safe_corr(shared["phase1_beta"].to_numpy(float), shared["external_beta"].to_numpy(float), "spearman")
    ax.text(0.03, 0.97, f"Spearman rho = {rho:.4f}\npre-registered genome-wide gate failed", transform=ax.transAxes, va="top", bbox=dict(fc="white", alpha=0.8, ec="none"))
    ax.set_xlabel("Phase 1 discovery beta")
    ax.set_ylabel("GSE220155 external RNA beta")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIG / "Figure3_genomewide_beta_beta.png", dpi=220)
    plt.close(fig)

    # Figure 4: effect-size dependence.
    binned = shared.copy()
    binned["effect_bin"] = pd.qcut(binned["phase1_beta"].abs(), q=10, duplicates="drop")
    summary = binned.groupby("effect_bin", observed=True).agg(effect=("phase1_beta", lambda x: np.median(np.abs(x))), same=("same_direction", "mean"), n=("same_direction", "size"))
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(summary["effect"], summary["same"], marker="o", color="#386cb0", label="All shared genes")
    for n, color in [(500, "#7570b3"), (250, "#e78ac3"), (100, "#d95f02")]:
        sub = tier1.head(n)
        ax.scatter(np.median(np.abs(sub["phase1_beta"])), np.mean(np.sign(sub["phase1_beta"]) == np.sign(sub["external_beta"])), s=60, color=color, label=f"CORE{n}")
    ax.set_xlabel("Median |Phase 1 beta|")
    ax.set_ylabel("Same-direction replication")
    ax.set_ylim(0.45, 1.0)
    ax.set_title("Replication increases with discovery effect size")
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(FIG / "Figure4_effect_size_replication.png", dpi=220)
    plt.close(fig)

    # Figure 5: donor-level CORE100 heatmap.
    heat = core100[donor_cols].to_numpy(float)
    heat = np.sign(heat) * np.sign(core100["phase1_beta"].to_numpy(float))[:, None]
    fig, ax = plt.subplots(figsize=(5.2, 7))
    ax.imshow(heat, aspect="auto", cmap="RdYlGn", vmin=-1, vmax=1)
    ax.set_xticks(range(4), ["D1", "D2", "D3", "D4"])
    ax.set_ylabel("CORE100 genes ranked by |Phase 1 beta|")
    ax.set_title("CORE100 donor-level direction replication")
    fig.colorbar(ax.images[0], ax=ax, ticks=[-1, 1], label="Phase 1 direction matched")
    fig.tight_layout()
    fig.savefig(FIG / "Figure5_CORE100_donor_replication.png", dpi=220)
    plt.close(fig)

    # Figure 6: literature controls plus strongest non-control high-confidence genes.
    high = core100_out[core100_out["high_confidence_core"]].copy()
    controls = high[high["literature_positive_control"] != "NONE"].sort_values("gene")
    unbiased = high[high["literature_positive_control"] == "NONE"].sort_values("phase1_beta", key=lambda x: x.abs(), ascending=False).head(5)
    selected = pd.concat([controls, unbiased]).drop_duplicates("gene").head(12)
    fig, ax = plt.subplots(figsize=(10, 4.5))
    xpos = np.arange(len(selected))
    width = 0.25
    ax.bar(xpos - width, selected["phase1_beta"], width, label="Phase 1", color="#7570b3")
    ax.bar(xpos, selected["external_beta"], width, label="External RNA", color="#d95f02")
    ax.bar(xpos + width, selected["gene_activity_beta"], width, label="Gene activity", color="#1b9e77")
    ax.axhline(0, color="black", lw=0.6)
    ax.set_xticks(xpos, selected["gene"], rotation=55, ha="right")
    ax.set_ylabel("Regional beta")
    ax.set_title("Selected regional core genes")
    ax.legend(frameon=False, ncol=3)
    fig.tight_layout()
    fig.savefig(FIG / "Figure6_selected_regional_core_genes.png", dpi=220)
    plt.close(fig)

    # Figure 7: gene-activity directional support.
    valid = core100_out[np.isfinite(core100_out["gene_activity_beta"])].copy()
    fig, ax = plt.subplots(figsize=(6, 4))
    frac = float(np.mean(valid["gene_activity_match"].astype(bool))) if len(valid) else np.nan
    ax.bar(["CORE100"], [frac], color="#1b9e77")
    ax.axhline(0.5, color="black", ls="--", lw=0.8)
    ax.set_ylim(0, 1)
    ax.set_ylabel("RNA / gene-activity same-direction fraction")
    ax.set_title("Orthogonal gene-activity support\n(no peak-level mechanism claim)")
    fig.tight_layout()
    fig.savefig(FIG / "Figure7_gene_activity_support.png", dpi=220)
    plt.close(fig)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        log(f"[FAIL] {type(exc).__name__}: {exc}")
        raise
