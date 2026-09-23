# Reproducibility README

## Intended order

1. Obtain the public CELLxGENE discovery source and the GEO GSE220155 source using their accession records.
2. Review the frozen metadata and donor-region eligibility rules.
3. Run the pseudobulk and regional model scripts using the frozen configuration.
4. Run the external validation script against the processed GSE220155 source.
5. Run the regional-core script using the fixed nested CORE sets and frozen branch configuration.
6. Compare outputs with the derived tables in `results/` and the source-data mappings in `DATA_MANIFEST.tsv`.

## Current reproducibility limitation

This release does not yet contain a versioned environment lock or session-information file. Validate the runtime and package versions before rerunning. Exact bitwise reproduction is not claimed without that information.

## Data boundary

The package does not redistribute raw public datasets that are already available through CELLxGENE and GEO. It contains selected derived tables only. Any restricted, third-party or duplicated raw data remain governed by the original source terms.
