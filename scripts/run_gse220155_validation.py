#!/usr/bin/env python3
"""Frozen Phase 3 validation of Phase 1 regional RPE effects in GSE220155.

This script uses only author-provided processed integer count matrices.  It does
not discover genes, tune thresholds, or perform downstream mechanism/genetics.
The GSE220155 archive contains gene-activity matrices for ATAC, not a common
peak universe or peak-to-gene links; the script records that limitation
explicitly and never promotes gene activity to peak-level evidence.
"""

from __future__ import annotations

import hashlib
import json
import math
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "raw" / "GSE220155"
PROC = ROOT / "data" / "processed" / "pivot_phase3"
OUT = ROOT / "results" / "pivot_phase3"
FIG = ROOT / "figures" / "pivot_phase3"
LOG = ROOT / "logs" / "pivot_phase3" / "GSE220155_VALIDATION.log"
CHECKSUM = ROOT / "checksums" / "pivot_phase3.sha256"
PHASE1_MODEL = ROOT / "results" / "pivot_phase1" / "MODEL_A_ALL_ELIGIBLE.tsv.gz"
PHASE1_TIERS = ROOT / "results" / "pivot_phase1" / "REGIONAL_GENE_TIERS_ANNOTATED.tsv"
SEED = 20260922
BOOTSTRAP_REPS = 2000


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


def parse_external_name(path: Path) -> dict[str, str]:
    match = re.match(r"(GSM\d+)_donor_(\d+)_(macula|peripheral)(?:_geneactivity)?_counts\.csv\.gz$", path.name)
    if not match:
        raise ValueError(f"Unexpected GSE220155 processed file name: {path.name}")
    gsm, donor_number, region = match.groups()
    return {
        "gsm": gsm,
        "donor_id": f"donor_{donor_number}",
        "region": region.upper(),
        "sample_id": f"donor_{donor_number}_{region.upper()}",
    }


def external_files(modality: str) -> list[Path]:
    if modality == "RNA":
        files = sorted(RAW.glob("GSM*_donor_*_*_counts.csv.gz"))
        return [p for p in files if "geneactivity" not in p.name]
    if modality == "ATAC_GENE_ACTIVITY":
        return sorted(RAW.glob("GSM*_donor_*_*_geneactivity_counts.csv.gz"))
    raise ValueError(modality)


def load_count_units(modality: str) -> tuple[list[str], dict[str, dict], np.ndarray, dict, pd.DataFrame]:
    files = external_files(modality)
    if len(files) != 8:
        raise RuntimeError(f"Expected 8 {modality} files, found {len(files)}")
    first = pd.read_csv(files[0], nrows=0)
    meta_cols = ["barcode", "region", "new_cluster_number", "library"]
    if list(first.columns[:4]) != meta_cols:
        raise RuntimeError(f"Unexpected metadata columns in {files[0].name}: {first.columns[:4].tolist()}")
    genes = first.columns[4:].tolist()
    if len(genes) != len(set(genes)):
        raise RuntimeError(f"Duplicate gene/feature names in {files[0].name}")
    unit_order: list[str] = []
    unit_info: dict[str, dict] = {}
    totals: list[np.ndarray] = []
    cluster_totals: dict[str, np.ndarray] = {}
    manifest_rows: list[dict] = []
    for path in files:
        parsed = parse_external_name(path)
        df = pd.read_csv(path)
        if df.columns[4:].tolist() != genes:
            raise RuntimeError(f"Gene/feature columns differ in {path.name}")
        if not df["barcode"].is_unique:
            raise RuntimeError(f"Duplicate barcodes in {path.name}")
        if set(df["region"].astype(str).str.upper()) != {parsed["region"]}:
            raise RuntimeError(f"Region field mismatch in {path.name}")
        values = df.iloc[:, 4:].to_numpy(copy=False)
        if not np.isfinite(values).all() or (values < 0).any() or not np.equal(values, np.floor(values)).all():
            raise RuntimeError(f"Non-negative integer count audit failed for {path.name}")
        values = values.astype(np.int64, copy=False)
        total = values.sum(axis=0, dtype=np.int64)
        unit = parsed["sample_id"]
        unit_order.append(unit)
        unit_info[unit] = {
            **parsed,
            "modality": modality,
            "processed_file": str(path.relative_to(ROOT)),
            "processed_sha256": sha256(path),
            "cell_count": int(len(df)),
            "library_size": int(total.sum()),
            "cluster_labels": "|".join(sorted(df["new_cluster_number"].astype(str).unique(), key=str)),
        }
        totals.append(total)
        for cluster, idx in df.groupby("new_cluster_number", sort=True).groups.items():
            cluster_key = f"{unit}__cluster_{str(cluster)}"
            cluster_totals[cluster_key] = values[idx.to_numpy(), :].sum(axis=0, dtype=np.int64)
        manifest_rows.append(
            {
                "GSM": parsed["gsm"],
                "donor": parsed["donor_id"],
                "region": parsed["region"],
                "modality": "RNA" if modality == "RNA" else "ATAC",
                "processed_file": str(path.relative_to(ROOT)),
                "raw_file": "SRA linked from GEO; FASTQ not downloaded",
                "barcode_metadata": "barcode, region, new_cluster_number, library embedded in processed CSV",
                "cluster": unit_info[unit]["cluster_labels"],
                "cell_count": int(len(df)),
                "usable": True,
                "reason": "author processed integer count matrix; sample is RPE with auditable donor and region",
                "processed_sha256": unit_info[unit]["processed_sha256"],
                "normalized_companion": str(path).replace("_counts.csv.gz", "_normalized.csv.gz").replace(str(ROOT) + "/", ""),
            }
        )
    matrix = np.column_stack(totals)
    return unit_order, unit_info, matrix, cluster_totals, pd.DataFrame(manifest_rows)


def logcpm(counts: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    library_sizes = counts.sum(axis=0, dtype=np.float64)
    values = np.log2((counts.astype(np.float64) + 0.5) / library_sizes[None, :] * 1_000_000.0)
    return values, library_sizes


def paired_effect_table(feature_names: list[str], unit_order: list[str], counts: np.ndarray) -> pd.DataFrame:
    values, library_sizes = logcpm(counts)
    donors = sorted({u.split("_")[1] for u in unit_order}, key=int)
    rows = []
    donor_beta_cols = [f"external_beta_donor_{d}" for d in donors]
    for i, feature in enumerate(feature_names):
        betas = []
        for donor in donors:
            mac = unit_order.index(f"donor_{donor}_MACULA")
            per = unit_order.index(f"donor_{donor}_PERIPHERAL")
            betas.append(float(values[i, mac] - values[i, per]))
        b = np.asarray(betas, dtype=float)
        n = int(np.isfinite(b).sum())
        mean = float(np.nanmean(b)) if n else np.nan
        se = float(np.nanstd(b, ddof=1) / math.sqrt(n)) if n > 1 else np.nan
        df = n - 1
        tcrit = float(stats.t.ppf(0.975, df)) if df > 0 else np.nan
        tstat = mean / se if n > 1 and se > 0 else np.nan
        p = float(2 * stats.t.sf(abs(tstat), df)) if np.isfinite(tstat) else np.nan
        row = {
            "external_feature": feature,
            "n_donors": n,
            "external_beta": mean,
            "external_SE": se,
            "external_CI95_low": mean - tcrit * se if np.isfinite(tcrit) and np.isfinite(se) else np.nan,
            "external_CI95_high": mean + tcrit * se if np.isfinite(tcrit) and np.isfinite(se) else np.nan,
            "external_p_descriptive": p,
            "external_total_library_size": float(np.nanmean(library_sizes)),
        }
        row.update({col: val for col, val in zip(donor_beta_cols, b)})
        rows.append(row)
    return pd.DataFrame(rows)


def safe_corr(x: pd.Series, y: pd.Series, method: str) -> float:
    if len(x) < 3:
        return np.nan
    if method == "pearson":
        return float(stats.pearsonr(x.to_numpy(float), y.to_numpy(float)).statistic)
    return float(stats.spearmanr(x.to_numpy(float), y.to_numpy(float)).statistic)


def bootstrap_spearman(x: np.ndarray, y: np.ndarray, reps: int, seed: int) -> tuple[float, float]:
    rng = np.random.default_rng(seed)
    rx = stats.rankdata(x)
    ry = stats.rankdata(y)
    n = len(x)
    boot = np.empty(reps, dtype=float)
    for i in range(reps):
        idx = rng.integers(0, n, n)
        boot[i] = np.corrcoef(rx[idx], ry[idx])[0, 1]
    return float(np.quantile(boot, 0.025)), float(np.quantile(boot, 0.975))


def direction_summary(df: pd.DataFrame, label: str) -> dict:
    d = df.copy()
    d = d[np.isfinite(d["phase1_beta"]) & np.isfinite(d["external_beta"]) & (d["phase1_beta"] != 0)]
    same = np.sign(d["phase1_beta"]) == np.sign(d["external_beta"])
    n = len(d)
    k = int(same.sum())
    p = float(stats.binomtest(k, n, 0.5, alternative="greater").pvalue) if n else np.nan
    return {
        "set": label,
        "n_evaluable": n,
        "n_same_direction": k,
        "same_direction_fraction": k / n if n else np.nan,
        "binomial_p_greater_0.5": p,
        "median_external_beta": float(d["external_beta"].median()) if n else np.nan,
    }


def fisher_named_gene_table(shared: pd.DataFrame) -> pd.DataFrame:
    # The accessible article text names these regional genes; a complete author
    # supplement table was not included in the GEO processed archive.  We keep
    # this as a limited cross-check and do not compute a false enrichment test.
    genes = {
        "WFDC1": ("MACULA", "expression"),
        "SULF1": ("MACULA", "expression"),
        "ELN": ("PERIPHERAL", "expression"),
        "SLC4A5": ("PERIPHERAL", "expression"),
        "ALDH1A3": ("PERIPHERAL", "ATAC_peak"),
    }
    rows = []
    for gene, (direction, signal_type) in genes.items():
        match = shared[shared["gene_symbol"] == gene]
        row = {
            "gene_symbol": gene,
            "author_signal_type": signal_type,
            "author_published_direction": direction,
            "phase1_gene_id": match["phase1_gene_id"].iloc[0] if len(match) else "",
            "phase1_tier": match["evidence_tier"].iloc[0] if len(match) else "NOT_IN_PHASE1_UNIQUE_MAPPING",
            "phase1_beta": match["phase1_beta"].iloc[0] if len(match) else np.nan,
            "external_rna_beta": match["external_beta"].iloc[0] if len(match) else np.nan,
            "phase1_tier1_overlap": bool(len(match) and match["evidence_tier"].iloc[0] == "TIER_1_ROBUST_REGIONAL"),
            "expression_direction_agreement": np.nan,
            "enrichment_OR": np.nan,
            "enrichment_P": np.nan,
            "crosscheck_status": "LIMITED_MAIN_TEXT_NAMED_GENE_CROSSCHECK;_FULL_SUPPLEMENT_TABLE_NOT_RETRIEVED",
        }
        if len(match) and signal_type == "expression" and np.isfinite(row["external_rna_beta"]):
            expected = 1 if direction == "MACULA" else -1
            row["expression_direction_agreement"] = bool(np.sign(row["external_rna_beta"]) == expected)
        rows.append(row)
    return pd.DataFrame(rows)


def paired_cluster_composition(unit_info: dict, cluster_totals: dict, unit_order: list[str], total_counts: np.ndarray) -> pd.DataFrame:
    rows = []
    cluster_names = sorted({key.split("__cluster_", 1)[1] for key in cluster_totals}, key=str)
    total_by_unit = {u: float(total_counts[:, i].sum()) for i, u in enumerate(unit_order)}
    for cluster in cluster_names:
        deltas = []
        fractions = []
        for donor_number in range(1, 5):
            mac_u = f"donor_{donor_number}_MACULA"
            per_u = f"donor_{donor_number}_PERIPHERAL"
            mac_key = f"{mac_u}__cluster_{cluster}"
            per_key = f"{per_u}__cluster_{cluster}"
            mac_frac = cluster_totals.get(mac_key, np.zeros(total_counts.shape[0])).sum() / total_by_unit[mac_u]
            per_frac = cluster_totals.get(per_key, np.zeros(total_counts.shape[0])).sum() / total_by_unit[per_u]
            fractions.append((mac_frac, per_frac))
            deltas.append(mac_frac - per_frac)
        deltas = np.asarray(deltas)
        p = float(2 * stats.t.sf(abs(deltas.mean() / (deltas.std(ddof=1) / 2)), 3)) if deltas.std(ddof=1) > 0 else np.nan
        rows.append({
            "row_type": "composition",
            "cluster": cluster,
            "n_donors": 4,
            "macula_fraction_mean": float(np.mean([x[0] for x in fractions])),
            "peripheral_fraction_mean": float(np.mean([x[1] for x in fractions])),
            "paired_delta_mac_minus_peripheral": float(deltas.mean()),
            "paired_delta_P_descriptive": p,
            "note": "descriptive composition audit; not a cell-level inferential test",
        })
    return pd.DataFrame(rows)


def subtype_gene_sensitivity(
    top_genes: pd.DataFrame, phase1_effects: pd.DataFrame, unit_order: list[str], cluster_totals: dict, gene_index: dict[str, int]
) -> pd.DataFrame:
    rows = []
    clusters = sorted({key.split("__cluster_", 1)[1] for key in cluster_totals}, key=str)
    for _, gene_row in top_genes.iterrows():
        symbol = gene_row["gene_symbol"]
        if symbol not in gene_index:
            continue
        idx = gene_index[symbol]
        phase_dir = np.sign(gene_row["phase1_beta"])
        for cluster in clusters:
            betas = []
            for donor_number in range(1, 5):
                mac_key = f"donor_{donor_number}_MACULA__cluster_{cluster}"
                per_key = f"donor_{donor_number}_PERIPHERAL__cluster_{cluster}"
                mac = cluster_totals.get(mac_key)
                per = cluster_totals.get(per_key)
                if mac is None or per is None:
                    continue
                mac_lib = mac.sum()
                per_lib = per.sum()
                if mac_lib <= 0 or per_lib <= 0:
                    continue
                betas.append(float(np.log2((mac[idx] + 0.5) / mac_lib * 1e6) - np.log2((per[idx] + 0.5) / per_lib * 1e6)))
            b = np.asarray(betas, dtype=float)
            same = int(np.sum(np.sign(b) == phase_dir)) if len(b) else 0
            rows.append({
                "row_type": "gene_sensitivity",
                "gene_symbol": symbol,
                "phase1_beta": float(gene_row["phase1_beta"]),
                "cluster": cluster,
                "n_donors_complete": int(len(b)),
                "same_direction_fraction": same / len(b) if len(b) else np.nan,
                "mean_cluster_beta": float(np.mean(b)) if len(b) else np.nan,
                "direction_consistent": bool(len(b) and same / len(b) >= 0.75),
                "note": "fixed Phase1 top-100 sensitivity; no new subtype discovery",
            })
    return pd.DataFrame(rows)


def write_gz(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, sep="\t", index=False, compression="gzip", na_rep="NA")


def make_figures(shared: pd.DataFrame, tier1: pd.DataFrame, atac_support: pd.DataFrame, metrics: dict) -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(10, 4))
    donors = ["donor_1", "donor_2", "donor_3", "donor_4"]
    for i, donor in enumerate(donors):
        ax.plot([0, 1], [i, i], color="#b8c2cc", lw=2)
        ax.scatter([0, 1], [i, i], s=110, color=["#d95f02", "#1b9e77"], zorder=3)
        ax.text(-0.08, i, donor, ha="right", va="center")
    ax.set_xticks([0, 1], ["Macula", "Peripheral"])
    ax.set_yticks([])
    ax.set_title("GSE220155 design: four paired RPE donors\nRNA and ATAC processed modalities")
    ax.set_xlim(-0.3, 1.3)
    fig.tight_layout()
    fig.savefig(FIG / "P3-1_GSE220155_design.png", dpi=220)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6, 5))
    x = shared["phase1_beta"].to_numpy(float)
    y = shared["external_beta"].to_numpy(float)
    ax.scatter(x, y, s=4, alpha=0.18, color="#52616b", rasterized=True)
    t = tier1[np.isfinite(tier1["external_beta"])]
    ax.scatter(t["phase1_beta"], t["external_beta"], s=5, alpha=0.35, color="#d95f02", rasterized=True)
    ax.axhline(0, color="black", lw=0.6)
    ax.axvline(0, color="black", lw=0.6)
    ax.set_xlabel("Phase 1 Model A beta (macula - periphery)")
    ax.set_ylabel("GSE220155 paired RNA beta")
    ax.set_title(f"Global concordance; Spearman rho={metrics['spearman']:.3f}")
    fig.tight_layout()
    fig.savefig(FIG / "P3-2_phase1_vs_external_RNA.png", dpi=220)
    plt.close(fig)

    labels = ["Tier1 all", "Macula-up", "Periphery-up"]
    vals = [metrics["tier1_same"], metrics["macula_same"], metrics["periphery_same"]]
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.bar(labels, vals, color=["#7570b3", "#d95f02", "#1b9e77"])
    ax.axhline(0.5, color="black", ls="--", lw=0.8)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Same-direction fraction")
    ax.set_title("Tier 1 RNA direction concordance")
    fig.tight_layout()
    fig.savefig(FIG / "P3-3_Tier1_direction_concordance.png", dpi=220)
    plt.close(fig)

    top = tier1[np.isfinite(tier1["external_beta"])].sort_values("abs_phase1_beta", ascending=False).head(10)
    fig, axes = plt.subplots(max(1, len(top)), 1, figsize=(7, max(2.0, 1.6 * len(top))), squeeze=False)
    axes = axes[:, 0]
    donor_cols = [f"external_beta_donor_{i}" for i in range(1, 5)]
    for ax, (_, row) in zip(axes, top.iterrows()):
        b = [row.get(c, np.nan) for c in donor_cols]
        ax.axhline(0, color="black", lw=0.5)
        ax.plot(range(1, 5), b, marker="o", color="#d95f02")
        ax.set_xticks(range(1, 5), ["D1", "D2", "D3", "D4"])
        ax.set_ylabel(str(row["gene_symbol"]), rotation=0, ha="right", va="center", labelpad=34)
    if len(top):
        axes[-1].set_xlabel("External donor-specific beta")
    fig.suptitle("Top fixed Phase 1 Tier 1 genes: donor-paired external effects", y=1.0)
    fig.tight_layout()
    fig.savefig(FIG / "P3-4_top_validated_genes_donor_pairs.png", dpi=220, bbox_inches="tight")
    plt.close(fig)

    support = atac_support[np.isfinite(atac_support["gene_activity_beta"])]
    fig, ax = plt.subplots(figsize=(6, 5))
    ax.scatter(support["external_beta"], support["gene_activity_beta"], s=8, alpha=0.35, color="#386cb0", rasterized=True)
    ax.axhline(0, color="black", lw=0.6)
    ax.axvline(0, color="black", lw=0.6)
    ax.set_xlabel("External RNA beta")
    ax.set_ylabel("External ATAC gene-activity beta")
    ax.set_title("Orthogonal gene-activity comparison\n(no peak-level link asserted)")
    fig.tight_layout()
    fig.savefig(FIG / "P3-5_RNA_vs_ATAC_gene_activity.png", dpi=220)
    plt.close(fig)


def main() -> int:
    LOG.parent.mkdir(parents=True, exist_ok=True)
    LOG.write_text("", encoding="utf-8")
    for directory in [PROC, OUT, FIG, CHECKSUM.parent]:
        directory.mkdir(parents=True, exist_ok=True)
    log("[START] Phase 3 GSE220155 validation")

    rna_units, rna_info, rna_counts, rna_cluster_counts, rna_manifest = load_count_units("RNA")
    atac_units, atac_info, atac_counts, atac_cluster_counts, atac_manifest = load_count_units("ATAC_GENE_ACTIVITY")
    if rna_units != atac_units:
        raise RuntimeError("RNA and ATAC unit order differs")
    if sorted(rna_info) != sorted(atac_info):
        raise RuntimeError("RNA and ATAC donor-region unit identities differ")
    manifest = pd.concat([rna_manifest, atac_manifest], ignore_index=True)
    manifest.to_csv(OUT / "GSE220155_SAMPLE_MANIFEST.tsv", sep="\t", index=False)
    pairing = manifest.groupby(["donor", "region"]).size().unstack(fill_value=0)
    if set(pairing.index) != {"donor_1", "donor_2", "donor_3", "donor_4"} or not (pairing > 0).all().all():
        raise RuntimeError("4/4 donor-region pairing audit failed")
    log(f"[PASS] manifest: {len(manifest)} rows; 4 donors x 2 regions x 2 modalities")

    # Keep feature names explicit instead of relying on later column order.
    rna_feature_names = pd.read_csv(external_files("RNA")[0], nrows=0).columns[4:].tolist()
    atac_feature_names = pd.read_csv(external_files("ATAC_GENE_ACTIVITY")[0], nrows=0).columns[4:].tolist()
    rna_counts_df = pd.DataFrame(rna_counts, index=rna_feature_names, columns=rna_units)
    atac_counts_df = pd.DataFrame(atac_counts, index=atac_feature_names, columns=atac_units)
    rna_counts_df.reset_index(names="gene_symbol").to_csv(PROC / "GSE220155_RNA_PSEUDOBULK_COUNTS.tsv.gz", sep="\t", index=False, compression="gzip")
    rna_meta = pd.DataFrame([rna_info[u] for u in rna_units])
    rna_meta.insert(0, "sample_id", rna_meta.pop("sample_id"))
    rna_meta.to_csv(PROC / "GSE220155_RNA_PSEUDOBULK_METADATA.tsv", sep="\t", index=False)

    phase1_model = pd.read_csv(PHASE1_MODEL, sep="\t", compression="gzip")
    phase1_tiers = pd.read_csv(PHASE1_TIERS, sep="\t", usecols=["gene_id", "gene_symbol", "evidence_tier"])
    phase1 = phase1_model[["gene_id", "log2FC_MACULA_vs_PERIPHERY"]].rename(columns={"log2FC_MACULA_vs_PERIPHERY": "phase1_beta"}).merge(phase1_tiers, on="gene_id", how="left")
    phase1["gene_symbol"] = phase1["gene_symbol"].fillna(phase1["gene_id"])
    ext_symbol_counts = pd.Series(rna_feature_names).value_counts()
    phase1_symbol_counts = phase1["gene_symbol"].value_counts()
    map_rows = phase1[["gene_id", "gene_symbol", "evidence_tier"]].copy()
    map_rows["external_gene_symbol"] = map_rows["gene_symbol"]
    map_rows["external_matches"] = map_rows["external_gene_symbol"].map(ext_symbol_counts).fillna(0).astype(int)
    map_rows["phase1_symbol_multiplicity"] = map_rows["gene_symbol"].map(phase1_symbol_counts).fillna(0).astype(int)
    map_rows["shared"] = (map_rows["external_matches"] > 0) & (map_rows["phase1_symbol_multiplicity"] == 1)
    map_rows["mapping_status"] = np.select(
        [map_rows["phase1_symbol_multiplicity"] > 1, map_rows["external_matches"] == 0],
        ["AMBIGUOUS_PHASE1_SYMBOL", "NO_EXTERNAL_MATCH"],
        default="MATCH_UNIQUE",
    )
    map_rows.to_csv(OUT / "RNA_GENE_MAPPING.tsv", sep="\t", index=False)
    unique_map = map_rows[map_rows["shared"]].set_index("external_gene_symbol")
    log(f"[PASS] Phase1 Set A={len(phase1)}; unique shared genes={len(unique_map)}; Tier1={sum(phase1.evidence_tier == 'TIER_1_ROBUST_REGIONAL')}")

    rna_effects = paired_effect_table(rna_feature_names, rna_units, rna_counts)
    shared = phase1.merge(rna_effects, left_on="gene_symbol", right_on="external_feature", how="inner")
    shared = shared[shared["gene_symbol"].map(phase1_symbol_counts) == 1].copy()
    shared = shared.rename(columns={"external_feature": "external_gene_symbol"})
    shared["phase1_gene_id"] = shared["gene_id"]
    donor_cols = [f"external_beta_donor_{i}" for i in range(1, 5)]
    shared.to_csv(OUT / "GSE220155_RNA_PAIRED_EFFECTS.tsv.gz", sep="\t", index=False, compression="gzip", na_rep="NA")

    tier1 = shared[shared["evidence_tier"] == "TIER_1_ROBUST_REGIONAL"].copy()
    tier1["abs_phase1_beta"] = tier1["phase1_beta"].abs()
    phase1_dir = np.sign(tier1["phase1_beta"])
    donor_same = tier1[donor_cols].apply(lambda row: int(np.sum(np.sign(row.to_numpy(float)) == np.sign(tier1.loc[row.name, "phase1_beta"]))), axis=1)
    tier1["n_donors_same_direction"] = donor_same
    tier1["donor_direction_class"] = pd.cut(donor_same, bins=[-1, 1, 2, 3, 4], labels=["<=1/4", "2/4", "3/4", "4/4"])
    tier1.to_csv(OUT / "TIER1_DONOR_DIRECTION.tsv", sep="\t", index=False, na_rep="NA")

    global_x = shared[np.isfinite(shared["phase1_beta"]) & np.isfinite(shared["external_beta"])].copy()
    spearman = safe_corr(global_x["phase1_beta"], global_x["external_beta"], "spearman")
    pearson = safe_corr(global_x["phase1_beta"], global_x["external_beta"], "pearson")
    sp_low, sp_high = bootstrap_spearman(global_x["phase1_beta"].to_numpy(float), global_x["external_beta"].to_numpy(float), BOOTSTRAP_REPS, SEED)
    macula_up = tier1[tier1["phase1_beta"] > 0]
    periphery_up = tier1[tier1["phase1_beta"] < 0]
    summaries = [
        direction_summary(tier1, "TIER1_ALL"),
        direction_summary(macula_up, "TIER1_MACULA_UP"),
        direction_summary(periphery_up, "TIER1_PERIPHERY_UP"),
    ]
    tier1_ge3 = float((tier1["n_donors_same_direction"] >= 3).mean()) if len(tier1) else np.nan
    tier1_ge4 = float((tier1["n_donors_same_direction"] == 4).mean()) if len(tier1) else np.nan
    metrics_rows = [
        {"metric": "shared_genes_set_A", "value": len(global_x), "note": "unique exact symbol mapping"},
        {"metric": "spearman_beta_concordance", "value": spearman, "note": "primary metric"},
        {"metric": "spearman_bootstrap_CI_low", "value": sp_low, "note": "ordinary gene bootstrap; gene dependence limitation"},
        {"metric": "spearman_bootstrap_CI_high", "value": sp_high, "note": "ordinary gene bootstrap; gene dependence limitation"},
        {"metric": "pearson_beta_concordance", "value": pearson, "note": "secondary metric"},
        {"metric": "tier1_ge3_of_4_donor_fraction", "value": tier1_ge3, "note": "primary donor consistency metric"},
        {"metric": "tier1_4_of_4_donor_fraction", "value": tier1_ge4, "note": "descriptive"},
    ]
    metrics_rows.extend({"metric": f"{row['set']}_{key}", "value": value, "note": "direction summary"} for row in summaries for key, value in row.items() if key != "set")
    metrics_df = pd.DataFrame(metrics_rows)
    metrics_df.to_csv(OUT / "RNA_VALIDATION_METRICS.tsv", sep="\t", index=False, na_rep="NA")

    top_rows = []
    tier1_for_top = tier1.sort_values("abs_phase1_beta", ascending=False)
    for n in [100, 250, 500]:
        subset = tier1_for_top.head(n)
        s = direction_summary(subset, f"TOP_{n}")
        top_rows.append({"phase1_set": f"Top {n}", "n_selected": len(subset), "direction_concordance": s["same_direction_fraction"], "median_validation_beta": s["median_external_beta"], "rank_correlation": safe_corr(subset["phase1_beta"], subset["external_beta"], "spearman")})
    top_df = pd.DataFrame(top_rows)
    top_df.to_csv(OUT / "TOP_N_VALIDATION.tsv", sep="\t", index=False, na_rep="NA")

    composition = paired_cluster_composition(rna_info, rna_cluster_counts, rna_units, rna_counts)
    gene_index = {g: i for i, g in enumerate(rna_feature_names)}
    subtype = subtype_gene_sensitivity(tier1_for_top.head(100), tier1_for_top, rna_units, rna_cluster_counts, gene_index)
    pd.concat([composition, subtype], ignore_index=True, sort=False).to_csv(OUT / "RPE_SUBTYPE_VALIDATION.tsv", sep="\t", index=False, na_rep="NA")

    published = fisher_named_gene_table(shared)
    published.to_csv(OUT / "PUBLISHED_RESULT_CROSSCHECK.tsv", sep="\t", index=False, na_rep="NA")

    atac_effects = paired_effect_table(atac_feature_names, atac_units, atac_counts)
    atac_effects["feature_type"] = "gene_activity"
    atac_effects["peak_universe_available"] = False
    atac_effects["peak_to_gene_link_available"] = False
    atac_effects["phase1_gene_id"] = atac_effects["external_feature"].map(unique_map["gene_id"] if "gene_id" in unique_map else pd.Series(dtype=object))
    atac_effects["external_gene_symbol"] = atac_effects["external_feature"]
    atac_effects.to_csv(OUT / "GSE220155_ATAC_PAIRED_EFFECTS.tsv.gz", sep="\t", index=False, compression="gzip", na_rep="NA")
    atac_lookup = atac_effects.set_index("external_gene_symbol")
    support = tier1.copy()
    support["linked_peaks"] = 0
    support["best_peak_beta"] = np.nan
    support["RNA_ATAC_direction_consistent"] = "NA_NO_PEAK_LINK"
    support["gene_activity_beta"] = support["gene_symbol"].map(atac_lookup["external_beta"])
    support["gene_activity_direction_consistent"] = np.where(
        np.isfinite(support["gene_activity_beta"]),
        np.sign(support["external_beta"]) == np.sign(support["gene_activity_beta"]),
        np.nan,
    )
    support["evidence_level"] = "GENE_ACTIVITY_ONLY_NO_PEAK_LINK"
    support["note"] = "GEO processed archive provides gene activity, not a common peak universe or peak-to-gene links"
    support.to_csv(OUT / "TIER1_ATAC_SUPPORT.tsv", sep="\t", index=False, na_rep="NA")

    validation = support.copy()
    validation["v_tier"] = np.where((validation["n_donors_same_direction"] >= 3), "V2", "V3")
    validation["v_tier"] = np.where((validation["n_donors_same_direction"] >= 3) & (validation["RNA_ATAC_direction_consistent"] == True), "V1", validation["v_tier"])
    validation["v_tier_reason"] = np.where(validation["v_tier"] == "V2", "RNA same direction and >=3/4 donors; no usable peak link", "RNA direction concordant but donor consistency below 3/4")
    validation[["gene_id", "gene_symbol", "phase1_beta", "external_beta", "n_donors_same_direction", "v_tier", "v_tier_reason"]].rename(columns={"gene_id": "phase1_gene_id", "gene_symbol": "gene"}).to_csv(OUT / "EXTERNAL_VALIDATION_TIERS.tsv", sep="\t", index=False, na_rep="NA")

    lodo_rows = []
    for _, row in validation[validation["v_tier"].isin(["V1", "V2"])].iterrows():
        expected = np.sign(row["phase1_beta"])
        statuses = []
        for donor_number in range(1, 5):
            keep = [c for c in donor_cols if c != f"external_beta_donor_{donor_number}"]
            beta = float(row[keep].mean())
            stable = bool(np.sign(beta) == expected)
            statuses.append(stable)
            lodo_rows.append({"phase1_gene_id": row["gene_id"], "gene": row["gene_symbol"], "v_tier": row["v_tier"], "dropped_donor": f"donor_{donor_number}", "external_lodo_beta": beta, "phase1_direction": expected, "direction_stable": stable})
        for item in lodo_rows[-4:]:
            item["external_lodo_stable_all_donors"] = bool(all(statuses))
            item["external_lodo_status"] = "EXTERNAL_LODO_UNSTABLE" if not all(statuses) else "STABLE"
    lodo_df = pd.DataFrame(lodo_rows)
    lodo_df.to_csv(OUT / "EXTERNAL_LODO_STABILITY.tsv", sep="\t", index=False, na_rep="NA")

    gene_activity_evaluable = support[np.isfinite(support["gene_activity_beta"])]
    gene_activity_same = float(gene_activity_evaluable["gene_activity_direction_consistent"].mean()) if len(gene_activity_evaluable) else np.nan
    tier1_same = summaries[0]["same_direction_fraction"]
    if spearman >= 0.40 and tier1_same >= 0.70 and tier1_ge3 >= 0.60:
        validation_class = "STRONG"
    elif spearman >= 0.25 and tier1_same >= 0.65 and tier1_ge3 >= 0.50:
        validation_class = "MODERATE"
    elif spearman < 0.25 or tier1_same < 0.60:
        validation_class = "FAIL"
    else:
        validation_class = "BORDERLINE"
    if validation_class == "STRONG":
        decision = "GO_PHASE4"
    elif validation_class == "MODERATE" and np.isfinite(gene_activity_same) and gene_activity_same >= 0.60:
        decision = "CONDITIONAL_GO_PHASE4"
    elif validation_class in {"MODERATE", "BORDERLINE"} and spearman >= 0.25 and tier1_same >= 0.60 and tier1_ge3 >= 0.50:
        decision = "CONDITIONAL_GO_PHASE4"
    else:
        decision = "NO_GO_EXTERNAL_VALIDATION"
    log(f"[RESULT] Spearman={spearman:.4f}; Tier1 same={tier1_same:.4f}; Tier1 >=3/4={tier1_ge3:.4f}; class={validation_class}; decision={decision}")

    metrics_dict = {"spearman": spearman, "tier1_same": tier1_same, "macula_same": summaries[1]["same_direction_fraction"], "periphery_same": summaries[2]["same_direction_fraction"]}
    make_figures(shared, tier1, support, metrics_dict)

    # Final top externally validated table, capped at 20 genes.
    top_validated = validation[validation["v_tier"].isin(["V1", "V2"])].sort_values("abs_phase1_beta", ascending=False).head(20).copy()
    top_validated = top_validated[["gene_symbol", "phase1_beta", "external_beta", "n_donors_same_direction", "v_tier"]].rename(columns={"gene_symbol": "Gene", "phase1_beta": "Phase1_beta", "external_beta": "External_RNA_beta", "n_donors_same_direction": "Donor_agreement", "v_tier": "V_tier"})
    top_validated["ATAC_beta"] = top_validated["Gene"].map(atac_lookup["external_beta"])
    top_validated["ATAC_link"] = "NO_ATAC_LINK"
    top_validated.to_csv(OUT / "TOP_EXTERNALLY_VALIDATED_GENES.tsv", sep="\t", index=False, na_rep="NA")

    report_lines = [
        "# PIVOT PHASE 3 — GSE220155 VALIDATION",
        "",
        f"**Run date:** 2026-09-22  **Decision:** `{decision}`  **RNA class:** `{validation_class}`",
        "",
        "## Scope and freeze",
        "",
        "This is an independent RPE-only validation of the frozen Phase 1 macula-versus-periphery resource. GSE220155 was not used to select genes, retune thresholds, run Age × Region validation, or initiate pathway, regulatory-mechanism, or AMD-genetics analyses.",
        "",
        "The GEO processed archive contains sample-level integer RNA counts and ATAC gene-activity counts. It does not contain a common processed peak universe or peak-to-gene link table. ATAC results below are therefore labeled gene-activity-only; no gene is granted peak-level V1 evidence.",
        "",
        "## Cohort audit",
        "",
        f"- 4 donors; 4/4 have paired macula and peripheral regions in both modalities.",
        f"- RNA units: 8 donor-region pseudobulk units; ATAC units: 8 donor-region gene-activity units.",
        f"- RNA cells retained: {int(sum(rna_meta.cell_count))}; ATAC cells retained: {int(sum(atac_manifest.cell_count))}.",
        f"- Phase 1 Set A genes: {len(phase1)}; unique shared genes: {len(global_x)}; Tier 1 evaluable: {len(tier1)}.",
        "- Primary cell scope: author-processed high-quality RPE nuclei; cluster labels 1–4 were retained for descriptive composition sensitivity.",
        "",
        "## RNA validation",
        "",
        f"- Phase 1 vs external beta Spearman: **{spearman:.4f}**; ordinary gene-bootstrap 95% CI **[{sp_low:.4f}, {sp_high:.4f}]**.",
        f"- Pearson: **{pearson:.4f}**.",
        f"- Tier 1 same-direction: **{summaries[0]['same_direction_fraction']:.4f}** ({summaries[0]['n_same_direction']}/{summaries[0]['n_evaluable']}); macula-up **{summaries[1]['same_direction_fraction']:.4f}**; periphery-up **{summaries[2]['same_direction_fraction']:.4f}**.",
        f"- Tier 1 ≥3/4 donor concordant: **{tier1_ge3:.4f}**; 4/4 donor concordant: **{tier1_ge4:.4f}**.",
        f"- Validation class: **{validation_class}** under the pre-registered gate.",
        "",
        "## Top-N validation",
        "",
        top_df.to_markdown(index=False),
        "",
        "## ATAC support",
        "",
        f"- Peak-level features tested: **0** (processed peak universe unavailable); gene-activity features tested: **{len(atac_effects)}**.",
        f"- Tier 1 genes with matched gene-activity effect: **{len(gene_activity_evaluable)}**; RNA/gene-activity same-direction fraction: **{gene_activity_same:.4f}**." if np.isfinite(gene_activity_same) else "- Tier 1 gene-activity direction support: unavailable.",
        f"- V1 genes: **{int((validation.v_tier == 'V1').sum())}**; V2 genes: **{int((validation.v_tier == 'V2').sum())}**; V3 genes: **{int((validation.v_tier == 'V3').sum())}**.",
        "- `NO_ATAC_LINK` is not treated as an RNA validation failure; it prevents peak-level V1 assignment.",
        "",
        "## Author-published cross-check",
        "",
        "The accessible article text names WFDC1 and SULF1 as macula-enriched expression genes, ELN and SLC4A5 as periphery-enriched expression genes, and ALDH1A3 as a gene associated with the largest number of regional ATAC peaks. A complete author Supplement table was not included in the GEO processed archive, so full-table overlap, odds ratio, and Fisher P were left `NA` rather than inferred from five named genes.",
        "",
        published[["gene_symbol", "author_signal_type", "author_published_direction", "phase1_tier", "phase1_beta", "external_rna_beta", "phase1_tier1_overlap", "expression_direction_agreement"]].to_markdown(index=False),
        "",
        "## Leave-one-donor-out",
        "",
        f"- V1/V2 genes evaluated: {validation[validation.v_tier.isin(['V1','V2'])].shape[0]}; stable genes: {int(lodo_df.groupby('phase1_gene_id').external_lodo_stable_all_donors.first().sum()) if len(lodo_df) else 0}; unstable genes: {int(lodo_df.groupby('phase1_gene_id').external_lodo_stable_all_donors.first().eq(False).sum()) if len(lodo_df) else 0}.",
        "- Any unstable gene is labeled `EXTERNAL_LODO_UNSTABLE` in the machine-readable output.",
        "",
        "## Phase 2 context",
        "",
        "No gene-level Age × Region interaction survived BH-FDR in Phase 2; Phase 2 Tier A = 0 and Tier B = 0. GSE220155 is not an independent validation of aging interactions.",
        "",
        "## Decision and lock",
        "",
        f"**{decision}**. Phase 4, if unlocked by this decision, is restricted to the externally validated V1/V2 RPE genes and GSE262151 regulatory-mechanism triangulation. AMD genetics, pathway, regulon, motif/chromVAR, fine-mapping, and colocalization remain locked.",
        "",
        "## Outputs",
        "",
        "- [Sample manifest](../results/pivot_phase3/GSE220155_SAMPLE_MANIFEST.tsv)",
        "- [RNA validation metrics](../results/pivot_phase3/RNA_VALIDATION_METRICS.tsv)",
        "- [External validation tiers](../results/pivot_phase3/EXTERNAL_VALIDATION_TIERS.tsv)",
        "- [ATAC support audit](../results/pivot_phase3/TIER1_ATAC_SUPPORT.tsv)",
        "- [Decision report](PIVOT_PHASE3_DECISION.md)",
    ]
    (ROOT / "reports" / "PIVOT_PHASE3_GSE220155_VALIDATION.md").write_text("\n".join(report_lines) + "\n", encoding="utf-8")
    decision_text = f"# PIVOT PHASE 3 DECISION\n\n- Decision: `{decision}`\n- RNA validation class: `{validation_class}`\n- Spearman: `{spearman:.6f}`\n- Tier1 same direction: `{tier1_same:.6f}`\n- Tier1 >=3/4 donor fraction: `{tier1_ge3:.6f}`\n- Peak-level ATAC: `UNAVAILABLE`; gene-activity-only support recorded\n\nPhase 2 had no BH-FDR-significant gene-level Age × Region interaction; no aging-interaction validation claim is permitted.\n\nPhase 4 remains restricted to V1/V2 genes if the decision is GO or CONDITIONAL_GO. AMD genetics and pathway/regulon/motif analyses remain locked.\n"
    (ROOT / "reports" / "PIVOT_PHASE3_DECISION.md").write_text(decision_text, encoding="utf-8")

    checks = []
    tracked = [
        ROOT / "config" / "pivot_phase3_freeze.yaml",
        ROOT / "results" / "pivot_phase3" / "GSE220155_SAMPLE_MANIFEST.tsv",
        ROOT / "results" / "pivot_phase3" / "RNA_GENE_MAPPING.tsv",
        ROOT / "results" / "pivot_phase3" / "GSE220155_RNA_PAIRED_EFFECTS.tsv.gz",
        ROOT / "results" / "pivot_phase3" / "RNA_VALIDATION_METRICS.tsv",
        ROOT / "results" / "pivot_phase3" / "TIER1_DONOR_DIRECTION.tsv",
        ROOT / "results" / "pivot_phase3" / "PUBLISHED_RESULT_CROSSCHECK.tsv",
        ROOT / "results" / "pivot_phase3" / "RPE_SUBTYPE_VALIDATION.tsv",
        ROOT / "results" / "pivot_phase3" / "GSE220155_ATAC_PAIRED_EFFECTS.tsv.gz",
        ROOT / "results" / "pivot_phase3" / "TIER1_ATAC_SUPPORT.tsv",
        ROOT / "results" / "pivot_phase3" / "EXTERNAL_VALIDATION_TIERS.tsv",
        ROOT / "results" / "pivot_phase3" / "EXTERNAL_LODO_STABILITY.tsv",
        ROOT / "reports" / "PIVOT_PHASE3_GSE220155_VALIDATION.md",
        ROOT / "reports" / "PIVOT_PHASE3_DECISION.md",
    ]
    for path in tracked:
        if path.exists():
            checks.append(f"{sha256(path)}  {path.relative_to(ROOT)}")
    CHECKSUM.write_text("\n".join(checks) + "\n", encoding="utf-8")
    input_paths = [RAW / "GSE220155_RAW.tar"] + sorted(RAW.glob("GSM*.csv.gz"))
    input_checks = [f"{sha256(path)}  {path.relative_to(ROOT)}" for path in input_paths if path.exists()]
    (CHECKSUM.parent / "pivot_phase3_inputs.sha256").write_text("\n".join(input_checks) + "\n", encoding="utf-8")
    log(f"[COMPLETE] outputs={len(checks)}; decision={decision}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:  # keep a durable failure record
        log(f"[FAIL] {type(exc).__name__}: {exc}")
        raise
