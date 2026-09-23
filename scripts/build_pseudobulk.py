#!/usr/bin/env python3
"""Build donor-region broad-RPE pseudobulk counts and pre-model QC.

This script intentionally performs no differential testing.  It uses the
author-labelled broad RPE cells from the local CellxGene H5AD and the raw.X
count matrix, aggregates all cells within donor x region, and writes the
auditable sample-level inputs for the R-based model stage.
"""

from __future__ import annotations

import argparse
import gzip
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd


def norm_text(value) -> str:
    if pd.isna(value):
        return "NA"
    return str(value)


def parse_age(value) -> float:
    m = re.search(r"(\d+(?:\.\d+)?)", norm_text(value))
    return float(m.group(1)) if m else np.nan


def unique_join(values) -> str:
    vals = sorted({norm_text(x) for x in values if norm_text(x) not in {"NA", "nan"}})
    return ";".join(vals) if vals else "NA"


def safe_id(donor: str, region: str) -> str:
    donor = re.sub(r"[^A-Za-z0-9_.-]+", "_", donor)
    return f"RPE_{donor}_{region.upper()}"


def aggregate_rows(matrix, rows: np.ndarray, chunk_size: int = 512) -> np.ndarray:
    out = None
    for start in range(0, len(rows), chunk_size):
        block = matrix[rows[start : start + chunk_size], :]
        block_sum = np.asarray(block.sum(axis=0)).ravel().astype(np.float64)
        out = block_sum if out is None else out + block_sum
    if out is None:
        raise ValueError("Cannot aggregate an empty group")
    return out


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    # outdir is project/data/processed/pivot_phase1; keep QC in the project
    # level results tree, not under data/results.
    qcdir = outdir.parents[2] / "results" / "pivot_phase1"
    qcdir.mkdir(parents=True, exist_ok=True)

    atlas = ad.read_h5ad(args.h5ad, backed="r")
    obs = atlas.obs.copy()
    # Confirm the frozen local RPE metadata size and use the H5AD's matching
    # author/majorclass labels for count-row selection.
    local = pd.read_parquet(args.metadata)
    if len(local) != 69244:
        raise ValueError(f"Unexpected frozen RPE metadata rows: {len(local)}")

    rpe_mask = (obs["majorclass"].astype(str) == "RPE cells") & obs["author_cell_type"].astype(str).str.startswith("RPE_")
    if int(rpe_mask.sum()) != len(local):
        raise ValueError(f"H5AD RPE row count {int(rpe_mask.sum())} != local metadata {len(local)}")

    raw = atlas.raw.X if atlas.raw is not None else atlas.X
    genes = atlas.raw.var.copy() if atlas.raw is not None else atlas.var.copy()
    gene_names = pd.Index(genes.index.astype(str))
    feature_names = genes.get("feature_name", pd.Series(gene_names, index=genes.index)).astype(str).to_numpy()
    mito_mask = np.char.startswith(np.char.upper(feature_names.astype(str)), "MT-")
    ribo_mask = np.char.startswith(np.char.upper(feature_names.astype(str)), "RPS") | np.char.startswith(np.char.upper(feature_names.astype(str)), "RPL")

    # Keep every donor-region pseudobulk, including the one disease extension
    # group, but mark it explicitly. All formal models will use YES only.
    rpe_positions = np.flatnonzero(rpe_mask.to_numpy())
    rpe_obs = obs.iloc[rpe_positions].copy()
    rpe_obs["_position"] = rpe_positions
    rpe_obs["region_norm"] = rpe_obs["Region"].astype(str).map({"Macular": "MACULA", "Peripheral": "PERIPHERY"}).fillna(rpe_obs["Region"].astype(str).str.upper())
    rpe_obs["disease_status"] = np.where(rpe_obs["reported_diseases"].astype(str).str.lower().eq("none"), "HEALTHY", "DISEASE_EXTENSION")
    rpe_obs["primary_eligible"] = np.where(rpe_obs["disease_status"].eq("HEALTHY"), "YES", "NO")

    region_by_donor = rpe_obs[rpe_obs["primary_eligible"].eq("YES")].groupby("donor_id", observed=True)["region_norm"].unique().to_dict()
    rows = []
    meta_rows = []
    qc_rows = []
    count_arrays = []

    grouped = rpe_obs.groupby(["donor_id", "region_norm"], sort=True, observed=True)
    total_groups = len(grouped)
    for i, ((donor, region), grp) in enumerate(grouped, start=1):
        positions = grp["_position"].to_numpy(dtype=np.int64)
        counts = aggregate_rows(raw, positions)
        sample_id = safe_id(norm_text(donor), norm_text(region))
        age_values = pd.to_numeric(grp["donor_age"].map(parse_age), errors="coerce")
        age = float(age_values.dropna().iloc[0]) if not age_values.dropna().empty else np.nan
        healthy = bool(grp["primary_eligible"].eq("YES").all())
        paired = "PAIRED" if healthy and set(region_by_donor.get(donor, [])) == {"MACULA", "PERIPHERY"} else "UNPAIRED"
        assay = unique_join(grp["assay"])
        metadata = {
            "sample_id": sample_id,
            "donor_id": norm_text(donor),
            "region": norm_text(region),
            "age": age,
            "age_precision": "exact_years" if np.isfinite(age) else "missing",
            "sex": unique_join(grp["sex"]),
            "ancestry_ethnicity": unique_join(grp["self_reported_ethnicity"]),
            "assay": assay,
            "batch_if_available": "NA",
            "n_RPE_cells": int(len(grp)),
            "n_libraries": int(grp["library_id"].nunique()),
            "n_samples": int(grp["sample_id"].nunique()),
            "library_ids": unique_join(grp["library_id"]),
            "sample_ids": unique_join(grp["sample_id"]),
            "disease_status": "HEALTHY" if healthy else "DISEASE_EXTENSION",
            "paired_status": paired,
            "primary_eligible": "YES" if healthy else "NO",
            "author_cell_type": unique_join(grp["author_cell_type"]),
            "source_h5ad": str(Path(args.h5ad)),
            "count_matrix": "raw.X",
        }
        meta_rows.append(metadata)
        count_arrays.append(counts)
        total = float(counts.sum())
        detected = int(np.count_nonzero(counts))
        mito = float(counts[mito_mask].sum() / total) if total else np.nan
        ribo = float(counts[ribo_mask].sum() / total) if total else np.nan
        qc_rows.append({
            **metadata,
            "total_counts": int(round(total)),
            "detected_genes": detected,
            "zero_genes": int(len(counts) - detected),
            "mitochondrial_fraction": mito,
            "ribosomal_fraction": ribo,
            "detected_genes_per_1000_counts": (detected / total * 1000) if total else np.nan,
            "mean_counts_per_cell": total / len(grp),
            "median_cell_counts": float(pd.to_numeric(grp["nCount_RNA"], errors="coerce").median()),
            "median_cell_detected_genes": float(pd.to_numeric(grp["nFeature_RNA"], errors="coerce").median()),
            "assay_source": assay,
        })
        if i == 1 or i % 10 == 0 or i == total_groups:
            print(f"aggregated {i}/{total_groups} donor-region groups", flush=True)

    metadata_df = pd.DataFrame(meta_rows).sort_values(["primary_eligible", "region", "donor_id"], ascending=[False, True, True])
    qc_df = pd.DataFrame(qc_rows).merge(metadata_df[["sample_id", "primary_eligible"]], on=["sample_id", "primary_eligible"], how="left")
    # Align count columns to the metadata order and write genes x samples.
    order = metadata_df["sample_id"].tolist()
    original_by_sample = {m["sample_id"]: arr for m, arr in zip(meta_rows, count_arrays)}
    matrix = np.column_stack([original_by_sample[s] for s in order])
    count_path = outdir / "RPE_PSEUDOBULK_COUNTS.tsv.gz"
    with gzip.open(count_path, "wt") as handle:
        handle.write("gene_id\t" + "\t".join(order) + "\n")
        for gene, row in zip(gene_names, matrix):
            vals = np.rint(row).astype(np.int64)
            handle.write(str(gene) + "\t" + "\t".join(map(str, vals)) + "\n")
    metadata_df.to_csv(outdir / "RPE_PSEUDOBULK_METADATA.tsv", sep="\t", index=False)
    qc_df.drop_duplicates(subset=["sample_id"]).sort_values(["primary_eligible", "region", "donor_id"], ascending=[False, True, True]).to_csv(qcdir / "PSEUDOBULK_SAMPLE_QC.tsv", sep="\t", index=False)

    primary = metadata_df[metadata_df["primary_eligible"].eq("YES")].copy()
    thresholds = [10, 20, 30, 50, 100]
    dist = []
    for threshold in thresholds:
        eligible = primary[primary["n_RPE_cells"] >= threshold]
        dist.append({
            "threshold_cells": threshold,
            "donor_region_total": int(len(eligible)),
            "macula_donor_region": int((eligible["region"] == "MACULA").sum()),
            "periphery_donor_region": int((eligible["region"] == "PERIPHERY").sum()),
            "paired_donors_retained": int(eligible.groupby("donor_id")["region"].nunique().ge(2).sum()),
            "min_cells_retained": int(eligible["n_RPE_cells"].min()) if len(eligible) else 0,
            "median_cells_retained": float(eligible["n_RPE_cells"].median()) if len(eligible) else np.nan,
            "median_total_counts": float(qc_df.merge(eligible[["sample_id"]], on="sample_id")["total_counts"].median()) if len(eligible) else np.nan,
            "median_detected_genes": float(qc_df.merge(eligible[["sample_id"]], on="sample_id")["detected_genes"].median()) if len(eligible) else np.nan,
        })
    pd.DataFrame(dist).to_csv(qcdir / "RPE_CELLCOUNT_DISTRIBUTION.tsv", sep="\t", index=False)
    print(f"wrote {count_path}", flush=True)
    print(f"wrote {len(metadata_df)} donor-region samples; primary={int((metadata_df.primary_eligible == 'YES').sum())}", flush=True)


if __name__ == "__main__":
    main()
