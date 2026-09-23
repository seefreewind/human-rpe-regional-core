#!/usr/bin/env python3
from pathlib import Path
import anndata as ad
import pandas as pd

project = Path(__file__).resolve().parents[2]
atlas = ad.read_h5ad(project / "data/raw/reference/cellxgene_00ad64f9-ff32-445d-a51b-05e499fa3491.h5ad", backed="r")
var = atlas.raw.var.copy() if atlas.raw is not None else atlas.var.copy()
mapping = pd.DataFrame({"gene_id": var.index.astype(str), "gene_symbol": var["feature_name"].astype(str).values})
mapping.to_csv(project / "results/pivot_phase1/GENE_ID_SYMBOL_MAP.tsv", sep="\t", index=False)
tiers = pd.read_csv(project / "results/pivot_phase1/REGIONAL_GENE_TIERS.tsv", sep="\t")
tiers = tiers.merge(mapping, on="gene_id", how="left")
cols = ["gene_id", "gene_symbol"] + [c for c in tiers.columns if c not in {"gene_id", "gene_symbol"}]
tiers[cols].to_csv(project / "results/pivot_phase1/REGIONAL_GENE_TIERS_ANNOTATED.tsv", sep="\t", index=False)
print("mapped", len(mapping), "genes; annotated", len(tiers), "tier rows")
atlas.file.close()

