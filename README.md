# Human RPE Regional Core

Reproducible donor-level analysis of regional transcriptional specialization in human retinal pigment epithelium.

## Scientific scope

This repository contains the analysis scripts, frozen configuration and selected derived tables used to evaluate regional RPE transcriptional effects and their transfer to an independent public cohort. The package supports donor-region pseudobulk construction, regional and sensitivity models, Age × Region testing, external RNA/gene-activity validation, fixed CORE100/250/500 comparisons, effect-size-dependent replication, and the 71-gene high-confidence regional core.

The global external concordance gate was not met (Spearman rho = 0.2033; pre-specified gate rho >= 0.25). The package therefore supports reproducibility and bounded interpretation; it does not claim peak-level regulation, causal AMD inference or clinical validation.

## Contents

- `config/`: frozen branch configuration.
- `scripts/`: analysis scripts relevant to pseudobulk, regional modelling, external validation and core definition.
- `results/`: selected frozen derived tables.
- `CODE_MANIFEST.tsv`: script-to-analysis mapping.
- `DATA_MANIFEST.tsv`: source and derived-data mapping.
- `REPRODUCIBILITY_README.md`: reproduction order and runtime notes.

## Source data

The discovery source is the public CELLxGENE collection identified in `DATA_MANIFEST.tsv`. The external validation source is public GEO accession GSE220155. Raw source matrices are not mirrored here; use the accessions and access routes in `DATA_MANIFEST.tsv`.

## Reproducibility

Validate package versions and input checksums before rerunning. The scripts and derived tables in this repository represent the frozen analysis branch; no new downstream analyses are part of this release.

## Citation

Version 1.0 is archived on Zenodo: [10.5281/zenodo.22915573](https://doi.org/10.5281/zenodo.22915573). The concept DOI for all versions is [10.5281/zenodo.22915572](https://doi.org/10.5281/zenodo.22915572). The archived release corresponds to the GitHub tag [`1.0`](https://github.com/seefreewind/human-rpe-regional-core/tree/1.0). Citation metadata are provided in `CITATION.cff`.

## Data and code status

This repository contains derived tables and code intended for public reuse. Reused public datasets remain governed by their original access and licensing terms.
