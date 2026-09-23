#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(pheatmap)
})

project <- normalizePath(".", mustWork = TRUE)
counts_file <- file.path(project, "data/processed/pivot_phase1/RPE_PSEUDOBULK_COUNTS.tsv.gz")
meta_file <- file.path(project, "data/processed/pivot_phase1/RPE_PSEUDOBULK_METADATA.tsv")
out_dir <- file.path(project, "results/pivot_phase1")
fig_dir <- file.path(project, "figures/pivot_phase1")
log_dir <- file.path(project, "logs/pivot_phase1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

write_tsv_gz <- function(x, path) {
  con <- gzfile(path, "wt")
  on.exit(close(con), add = TRUE)
  write.table(x, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

message("Reading donor-region counts")
counts_df <- read.delim(gzfile(counts_file), check.names = FALSE, stringsAsFactors = FALSE)
gene_ids <- make.unique(as.character(counts_df$gene_id))
counts_df$gene_id <- NULL
counts <- as.matrix(counts_df)
storage.mode(counts) <- "numeric"
rownames(counts) <- gene_ids
meta <- read.delim(meta_file, check.names = FALSE, stringsAsFactors = FALSE)
meta <- meta[match(colnames(counts), meta$sample_id), , drop = FALSE]
if (anyNA(meta$sample_id) || any(meta$sample_id != colnames(counts))) stop("Counts and metadata are not aligned")

meta$region <- factor(meta$region, levels = c("PERIPHERY", "MACULA"))
meta$sex <- factor(meta$sex)
meta$ancestry_ethnicity <- factor(meta$ancestry_ethnicity)
meta$assay_covariate <- factor(ifelse(grepl(";", meta$assay), "MIXED_ASSAY", "10x_3p_v3"))
meta$age <- as.numeric(meta$age)
meta$donor_id <- factor(meta$donor_id)

make_subset <- function(d, threshold = 30, adult_only = FALSE, age_overlap = FALSE, exclude_mixed = FALSE) {
  keep <- d$primary_eligible == "YES" & d$n_RPE_cells >= threshold
  if (adult_only) keep <- keep & d$age >= 18
  if (age_overlap) {
    rmin <- tapply(d$age[keep], d$region[keep], min, na.rm = TRUE)
    rmax <- tapply(d$age[keep], d$region[keep], max, na.rm = TRUE)
    lo <- max(rmin[c("MACULA", "PERIPHERY")])
    hi <- min(rmax[c("MACULA", "PERIPHERY")])
    keep <- keep & d$age >= lo & d$age <= hi
  }
  if (exclude_mixed) keep <- keep & d$assay_covariate != "MIXED_ASSAY"
  ans <- d[keep, , drop = FALSE]
  rownames(ans) <- ans$sample_id
  ans
}

make_design_a <- function(d) {
  d$Region <- relevel(droplevels(factor(d$region)), ref = "PERIPHERY")
  d$Age_c <- d$age - mean(d$age, na.rm = TRUE)
  d$Sex <- droplevels(factor(d$sex))
  d$Ancestry <- droplevels(factor(d$ancestry_ethnicity))
  d$Assay <- droplevels(factor(d$assay_covariate))
  terms <- c("Region", "Age_c")
  for (v in c("Sex", "Ancestry", "Assay")) {
    if (nlevels(d[[v]]) > 1) terms <- c(terms, v)
  }
  design <- model.matrix(as.formula(paste("~", paste(terms, collapse = " + "))), d)
  if (qr(design)$rank < ncol(design)) {
    for (v in rev(c("Assay", "Ancestry", "Sex"))) {
      if (v %in% terms) {
        trial_terms <- setdiff(terms, v)
        trial <- model.matrix(as.formula(paste("~", paste(trial_terms, collapse = " + "))), d)
        if (qr(trial)$rank == ncol(trial)) {
          terms <- trial_terms
          design <- trial
          break
        }
      }
    }
  }
  coef_name <- grep("^RegionMACULA$", colnames(design), value = TRUE)
  if (length(coef_name) != 1 || qr(design)$rank < ncol(design)) stop("Model A design is not estimable")
  list(data = d, design = design, coef = coef_name, terms = terms)
}

make_design_b <- function(d) {
  d$Region <- relevel(droplevels(factor(d$region)), ref = "PERIPHERY")
  d$Donor <- droplevels(factor(d$donor_id))
  design <- model.matrix(~ Donor + Region, d)
  coef_name <- grep("^RegionMACULA$", colnames(design), value = TRUE)
  if (length(coef_name) != 1 || qr(design)$rank < ncol(design)) stop("Model B design is not estimable")
  list(data = d, design = design, coef = coef_name, terms = colnames(design))
}

extract_fit_table <- function(fit, coef_name, genes) {
  tt <- topTable(fit, coef = coef_name, number = Inf, sort.by = "none")
  se <- fit$stdev.unscaled[, coef_name] * sqrt(fit$s2.post)
  df <- fit$df.total
  beta <- as.numeric(fit$coefficients[, coef_name])
  data.frame(
    gene_id = genes,
    log2FC_MACULA_vs_PERIPHERY = beta,
    SE = as.numeric(se),
    CI95_low = as.numeric(beta - qt(0.975, df = df) * se),
    CI95_high = as.numeric(beta + qt(0.975, df = df) * se),
    raw_P = as.numeric(tt$P.Value),
    FDR = as.numeric(tt$adj.P.Val),
    stringsAsFactors = FALSE
  )
}

fit_model_a <- function(d, gene_universe = NULL, counts_all = counts, primary_filter = FALSE) {
  md <- make_design_a(d)
  ids <- md$data$sample_id
  y <- DGEList(counts = counts_all[, ids, drop = FALSE])
  if (primary_filter) {
    keep <- filterByExpr(y, design = md$design)
    genes <- rownames(y)[keep]
  } else {
    genes <- intersect(rownames(y), gene_universe)
    y <- y[genes, , keep.lib.sizes = FALSE]
    keep2 <- filterByExpr(y, design = md$design)
    genes <- genes[keep2]
    y <- y[genes, , keep.lib.sizes = FALSE]
    keep <- rep(FALSE, nrow(counts_all))
    names(keep) <- rownames(counts_all)
    keep[genes] <- TRUE
  }
  if (length(genes) < 100) stop("Too few genes after donor-aware filtering")
  y <- calcNormFactors(y)
  v0 <- voom(y, md$design, plot = FALSE)
  corfit <- duplicateCorrelation(v0, md$design, block = md$data$donor_id)
  cor_value <- as.numeric(corfit$consensus)
  if (!is.finite(cor_value)) cor_value <- 0
  v <- voom(y, md$design, plot = FALSE, block = md$data$donor_id, correlation = cor_value)
  fit <- lmFit(v, md$design, block = md$data$donor_id, correlation = cor_value)
  fit <- eBayes(fit, robust = TRUE)
  tab <- extract_fit_table(fit, md$coef, genes)
  list(table = tab, fit = fit, y = y, voom = v, design = md$design, data = md$data, coef = md$coef, terms = md$terms, correlation = cor_value, keep = keep)
}

fit_model_b <- function(d, gene_universe, counts_all = counts) {
  md <- make_design_b(d)
  ids <- md$data$sample_id
  y0 <- DGEList(counts = counts_all[, ids, drop = FALSE])
  genes <- intersect(rownames(y0), gene_universe)
  y <- y0[genes, , keep.lib.sizes = FALSE]
  keep2 <- filterByExpr(y, design = md$design)
  genes <- genes[keep2]
  y <- y[genes, , keep.lib.sizes = FALSE]
  y <- calcNormFactors(y)
  v <- voom(y, md$design, plot = FALSE)
  fit <- lmFit(v, md$design)
  fit <- eBayes(fit, robust = TRUE)
  list(table = extract_fit_table(fit, md$coef, genes), fit = fit, y = y, voom = v, design = md$design, data = md$data, coef = md$coef, terms = md$terms, correlation = NA_real_)
}

complete_with_universe <- function(tab, universe) {
  out <- data.frame(gene_id = universe, stringsAsFactors = FALSE)
  out <- merge(out, tab, by = "gene_id", all.x = TRUE, sort = FALSE)
  out[match(universe, out$gene_id), , drop = FALSE]
}

write_model <- function(tab, path, model_name, d, method, terms) {
  tab$model <- model_name
  tab$method <- method
  tab$n_samples <- nrow(d)
  tab$n_donors <- length(unique(d$donor_id))
  tab$n_macula <- sum(d$region == "MACULA")
  tab$n_periphery <- sum(d$region == "PERIPHERY")
  tab$design_terms <- paste(terms, collapse = " + ")
  tab <- tab[, c("gene_id", "log2FC_MACULA_vs_PERIPHERY", "SE", "CI95_low", "CI95_high", "raw_P", "FDR", "model", "method", "n_samples", "n_donors", "n_macula", "n_periphery", "design_terms")]
  write_tsv_gz(tab, path)
  tab
}

primary <- make_subset(meta, 30)
if (nrow(primary) < 40) stop("Primary threshold leaves too few samples")
design_primary <- make_design_a(primary)
y_primary <- DGEList(counts = counts[, primary$sample_id, drop = FALSE])
gene_keep <- filterByExpr(y_primary, design = design_primary$design)
gene_universe <- rownames(y_primary)[gene_keep]
gene_filter <- data.frame(
  metric = c("genes_before", "genes_after", "fraction_retained", "primary_min_cells", "lower_sensitivity_min_cells", "higher_sensitivity_min_cells", "filter_rule", "primary_samples", "primary_donors", "primary_macula", "primary_periphery"),
  value = c(nrow(counts), length(gene_universe), length(gene_universe) / nrow(counts), 30, 20, 50, "edgeR::filterByExpr(counts, design=primary_Model_A_design); one shared universe", nrow(primary), length(unique(primary$donor_id)), sum(primary$region == "MACULA"), sum(primary$region == "PERIPHERY")),
  stringsAsFactors = FALSE
)
write.table(gene_filter, file.path(out_dir, "GENE_FILTER_FREEZE.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

message("Fitting Model A")
a_fit <- fit_model_a(primary, gene_universe = gene_universe, counts_all = counts)
a_tab <- write_model(complete_with_universe(a_fit$table, gene_universe), file.path(out_dir, "MODEL_A_ALL_ELIGIBLE.tsv.gz"), "MODEL_A_ALL_ELIGIBLE", primary, "limma_voom_duplicateCorrelation_donor_block", a_fit$terms)

paired <- primary[primary$donor_id %in% names(which(table(primary$donor_id) == 2)), , drop = FALSE]
paired <- paired[order(paired$donor_id, paired$region), , drop = FALSE]
message("Fitting paired Model B on ", length(unique(paired$donor_id)), " donors")
b_fit <- fit_model_b(paired, gene_universe, counts)
b_tab <- write_model(complete_with_universe(b_fit$table, gene_universe), file.path(out_dir, "MODEL_B_PAIRED.tsv.gz"), "MODEL_B_PAIRED", paired, "limma_voom_fixed_donor_block", b_fit$terms)

adult <- make_subset(meta, 30, adult_only = TRUE)
adult_fit <- fit_model_a(adult, gene_universe = gene_universe, counts_all = counts)
adult_tab <- write_model(complete_with_universe(adult_fit$table, gene_universe), file.path(out_dir, "MODEL_A_ADULT_ONLY.tsv.gz"), "MODEL_A_ADULT_ONLY", adult, "limma_voom_duplicateCorrelation_donor_block", adult_fit$terms)

adult_paired <- adult[adult$donor_id %in% names(which(table(adult$donor_id) == 2)), , drop = FALSE]
adult_paired <- adult_paired[order(adult_paired$donor_id, adult_paired$region), , drop = FALSE]
adult_b_tab <- NULL
if (length(unique(adult_paired$donor_id)) >= 5) {
  adult_b_fit <- fit_model_b(adult_paired, gene_universe, counts)
  adult_b_tab <- write_model(complete_with_universe(adult_b_fit$table, gene_universe), file.path(out_dir, "MODEL_B_ADULT_PAIRED.tsv.gz"), "MODEL_B_ADULT_PAIRED", adult_paired, "limma_voom_fixed_donor_block", adult_b_fit$terms)
}

age_overlap <- make_subset(meta, 30, age_overlap = TRUE)
age_fit <- fit_model_a(age_overlap, gene_universe = gene_universe, counts_all = counts)
age_tab <- write_model(complete_with_universe(age_fit$table, gene_universe), file.path(out_dir, "MODEL_A_AGE_OVERLAP.tsv.gz"), "MODEL_A_AGE_OVERLAP", age_overlap, "limma_voom_duplicateCorrelation_donor_block", age_fit$terms)

lower <- make_subset(meta, 20)
lower_fit <- fit_model_a(lower, gene_universe = gene_universe, counts_all = counts)
lower_tab <- write_model(complete_with_universe(lower_fit$table, gene_universe), file.path(out_dir, "MODEL_A_CELLCOUNT_20.tsv.gz"), "MODEL_A_CELLCOUNT_20", lower, "limma_voom_duplicateCorrelation_donor_block", lower_fit$terms)
higher <- make_subset(meta, 50)
higher_fit <- fit_model_a(higher, gene_universe = gene_universe, counts_all = counts)
higher_tab <- write_model(complete_with_universe(higher_fit$table, gene_universe), file.path(out_dir, "MODEL_A_CELLCOUNT_50.tsv.gz"), "MODEL_A_CELLCOUNT_50", higher, "limma_voom_duplicateCorrelation_donor_block", higher_fit$terms)

mixed_excl <- make_subset(meta, 30, exclude_mixed = TRUE)
mixed_fit <- fit_model_a(mixed_excl, gene_universe = gene_universe, counts_all = counts)
mixed_tab <- write_model(complete_with_universe(mixed_fit$table, gene_universe), file.path(out_dir, "MODEL_A_MIXED_ASSAY_EXCLUDED.tsv.gz"), "MODEL_A_MIXED_ASSAY_EXCLUDED", mixed_excl, "limma_voom_duplicateCorrelation_donor_block", mixed_fit$terms)

cor_summary <- function(x, y, label) {
  z <- merge(x[, c("gene_id", "log2FC_MACULA_vs_PERIPHERY", "FDR")], y[, c("gene_id", "log2FC_MACULA_vs_PERIPHERY")], by = "gene_id", suffixes = c("_x", "_y"))
  z <- z[is.finite(z$log2FC_MACULA_vs_PERIPHERY_x) & is.finite(z$log2FC_MACULA_vs_PERIPHERY_y), , drop = FALSE]
  hit <- z$FDR < 0.05
  rho <- if (nrow(z) >= 3) suppressWarnings(cor(z$log2FC_MACULA_vs_PERIPHERY_x, z$log2FC_MACULA_vs_PERIPHERY_y, method = "spearman")) else NA_real_
  dir_all <- if (nrow(z) > 0) mean(sign(z$log2FC_MACULA_vs_PERIPHERY_x) == sign(z$log2FC_MACULA_vs_PERIPHERY_y)) else NA_real_
  dir_hit <- if (sum(hit) > 0) mean(sign(z$log2FC_MACULA_vs_PERIPHERY_x[hit]) == sign(z$log2FC_MACULA_vs_PERIPHERY_y[hit])) else NA_real_
  data.frame(comparison = label, n_genes = nrow(z), spearman_rho = rho, direction_concordance_all = dir_all, n_model_a_fdr05 = sum(hit), same_direction_among_model_a_hits = dir_hit, stringsAsFactors = FALSE)
}
conc <- rbind(
  cor_summary(a_tab, b_tab, "MODEL_A_vs_MODEL_B_PAIRED"),
  cor_summary(a_tab, adult_tab, "MODEL_A_vs_ADULT_ONLY"),
  cor_summary(a_tab, age_tab, "MODEL_A_vs_AGE_OVERLAP"),
  cor_summary(a_tab, lower_tab, "MODEL_A_vs_CELLCOUNT_20"),
  cor_summary(a_tab, higher_tab, "MODEL_A_vs_CELLCOUNT_50"),
  cor_summary(a_tab, mixed_tab, "MODEL_A_vs_MIXED_ASSAY_EXCLUDED")
)
write.table(conc, file.path(out_dir, "MODEL_CONCORDANCE.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

merge_effects <- function(x, y, suffix) {
  out <- y[, c("gene_id", "log2FC_MACULA_vs_PERIPHERY", "SE", "CI95_low", "CI95_high", "raw_P", "FDR")]
  names(out)[-1] <- paste0(names(out)[-1], suffix)
  merge(x, out, by = "gene_id", all.x = TRUE, sort = FALSE)
}
tiers <- a_tab[, c("gene_id", "log2FC_MACULA_vs_PERIPHERY", "SE", "CI95_low", "CI95_high", "raw_P", "FDR")]
names(tiers) <- c("gene_id", "Model_A_log2FC", "Model_A_SE", "Model_A_CI95_low", "Model_A_CI95_high", "Model_A_raw_P", "Model_A_FDR")
tiers <- merge_effects(tiers, b_tab, "_paired")
tiers <- merge_effects(tiers, adult_tab, "_adult")
tiers <- merge_effects(tiers, age_tab, "_age_overlap")
same_pair <- is.finite(tiers$log2FC_MACULA_vs_PERIPHERY_paired) & sign(tiers$Model_A_log2FC) == sign(tiers$log2FC_MACULA_vs_PERIPHERY_paired)
same_adult <- is.finite(tiers$log2FC_MACULA_vs_PERIPHERY_adult) & sign(tiers$Model_A_log2FC) == sign(tiers$log2FC_MACULA_vs_PERIPHERY_adult)
tiers$evidence_tier <- "NOT_SIGNIFICANT"
sig <- is.finite(tiers$Model_A_FDR) & tiers$Model_A_FDR < 0.05
tiers$evidence_tier[sig & same_pair & same_adult & abs(tiers$Model_A_log2FC) >= 0.25] <- "TIER_1_ROBUST_REGIONAL"
tiers$evidence_tier[sig & same_pair & tiers$evidence_tier == "NOT_SIGNIFICANT"] <- "TIER_2_SUPPORTED_REGIONAL"
tiers$evidence_tier[sig & tiers$evidence_tier == "NOT_SIGNIFICANT"] <- "TIER_3_MODEL_A_ONLY"
tiers$paired_same_direction <- same_pair
tiers$adult_same_direction <- same_adult
tiers$primary_effect_threshold_abs_log2FC <- 0.25
tiers <- tiers[order(tiers$evidence_tier, tiers$Model_A_FDR, -abs(tiers$Model_A_log2FC)), , drop = FALSE]
write.table(tiers, file.path(out_dir, "REGIONAL_GENE_TIERS.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

tier1_genes <- tiers$gene_id[tiers$evidence_tier == "TIER_1_ROBUST_REGIONAL"]
lodo_rows <- list()
if (length(tier1_genes) > 0) {
  full_sign <- sign(tiers$Model_A_log2FC[match(tier1_genes, tiers$gene_id)])
  donors <- unique(as.character(primary$donor_id))
  # Refit all Tier 1 genes in one limma call per left-out donor. This is an
  # exact donor-level coefficient re-estimation on the fixed primary voom
  # scale and avoids thousands of redundant one-gene voom fits.
  for (donor in donors) {
    keep_sample <- as.character(primary$donor_id) != donor
    d2 <- primary[keep_sample, , drop = FALSE]
    md2 <- make_design_a(d2)
    gene_idx_lodo <- match(tier1_genes, rownames(a_fit$voom$E))
    sample_idx_lodo <- match(d2$sample_id, colnames(a_fit$voom$E))
    E2 <- a_fit$voom$E[gene_idx_lodo, sample_idx_lodo, drop = FALSE]
    W2 <- a_fit$voom$weights[gene_idx_lodo, sample_idx_lodo, drop = FALSE]
    f2 <- tryCatch(lmFit(E2, md2$design, weights = W2, block = md2$data$donor_id, correlation = a_fit$correlation), error = function(e) NULL)
    beta <- if (is.null(f2)) rep(NA_real_, length(tier1_genes)) else as.numeric(f2$coefficients[, md2$coef])
    lodo_rows[[length(lodo_rows) + 1]] <- data.frame(
      gene_id = tier1_genes, donor_left_out = donor, lodo_log2FC = beta,
      full_model_log2FC = tiers$Model_A_log2FC[match(tier1_genes, tiers$gene_id)],
      direction_flip = is.finite(beta) & sign(beta) != full_sign,
      stringsAsFactors = FALSE
    )
    message("LODO donor ", donor)
  }
}
if (length(lodo_rows) > 0) {
  lodo <- do.call(rbind, lodo_rows)
  stability <- do.call(rbind, lapply(tier1_genes, function(g) {
    z <- lodo[lodo$gene_id == g, , drop = FALSE]
    ok <- is.finite(z$lodo_log2FC)
    data.frame(gene_id = g, full_model_log2FC = z$full_model_log2FC[1], min_lodo_log2FC = if (any(ok)) min(z$lodo_log2FC[ok]) else NA_real_, max_lodo_log2FC = if (any(ok)) max(z$lodo_log2FC[ok]) else NA_real_, n_direction_flips = sum(z$direction_flip, na.rm = TRUE), LODO_status = if (any(ok) && sum(z$direction_flip, na.rm = TRUE) == 0) "STABLE" else "UNSTABLE_OR_UNESTIMABLE", stringsAsFactors = FALSE)
  }))
  write.table(stability, file.path(out_dir, "TIER1_LODO_STABILITY.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
} else {
  write.table(data.frame(gene_id = character(), full_model_log2FC = numeric(), min_lodo_log2FC = numeric(), max_lodo_log2FC = numeric(), n_direction_flips = integer(), LODO_status = character()), file.path(out_dir, "TIER1_LODO_STABILITY.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}

qc <- read.delim(file.path(out_dir, "PSEUDOBULK_SAMPLE_QC.tsv"), check.names = FALSE, stringsAsFactors = FALSE)
qc <- qc[qc$primary_eligible == "YES" & qc$n_RPE_cells >= 30, , drop = FALSE]
qc$region <- factor(qc$region, levels = c("PERIPHERY", "MACULA"))
p_lib <- ggplot(qc, aes(x = reorder(sample_id, total_counts), y = total_counts, colour = region)) + geom_point(size = 2.2) + scale_y_log10() + coord_flip() + theme_bw(base_size = 9) + labs(x = "Donor-region pseudobulk", y = "Total raw counts (log10)", colour = "Region")
ggsave(file.path(fig_dir, "P1_QC_library_size.png"), p_lib, width = 8, height = 8, dpi = 220)
p_det <- ggplot(qc, aes(x = reorder(sample_id, detected_genes), y = detected_genes, colour = region)) + geom_point(size = 2.2) + coord_flip() + theme_bw(base_size = 9) + labs(x = "Donor-region pseudobulk", y = "Detected genes", colour = "Region")
ggsave(file.path(fig_dir, "P1_QC_detected_genes.png"), p_det, width = 8, height = 8, dpi = 220)

logcpm <- cpm(a_fit$y, log = TRUE, prior.count = 2)
var_genes <- head(order(apply(logcpm, 1, var, na.rm = TRUE), decreasing = TRUE), min(5000, nrow(logcpm)))
pca <- prcomp(t(logcpm[var_genes, , drop = FALSE]), scale. = FALSE)
pca_df <- data.frame(sample_id = rownames(pca$x), PC1 = pca$x[, 1], PC2 = pca$x[, 2], region = primary[rownames(pca$x), "region"], age = primary[rownames(pca$x), "age"])
p_pca <- ggplot(pca_df, aes(PC1, PC2, colour = region)) + geom_point(size = 3) + theme_bw() + labs(title = "Donor-level RPE pseudobulk PCA", subtitle = paste0("Top variable genes: ", length(var_genes)))
ggsave(file.path(fig_dir, "P1_PCA.png"), p_pca, width = 7, height = 5.5, dpi = 220)
ggsave(file.path(fig_dir, "P1-2_pseudobulk_PCA.png"), p_pca, width = 7, height = 5.5, dpi = 220)
mds <- plotMDS(logcpm[var_genes, , drop = FALSE], plot = FALSE)
mds_df <- data.frame(sample_id = colnames(logcpm), MDS1 = mds$x, MDS2 = mds$y, region = primary[colnames(logcpm), "region"])
p_mds <- ggplot(mds_df, aes(MDS1, MDS2, colour = region)) + geom_point(size = 3) + theme_bw() + labs(title = "Donor-level RPE pseudobulk MDS")
ggsave(file.path(fig_dir, "P1_MDS.png"), p_mds, width = 7, height = 5.5, dpi = 220)
corr <- cor(logcpm[var_genes, , drop = FALSE], method = "pearson")
png(file.path(fig_dir, "P1_sample_correlation_heatmap.png"), width = 1800, height = 1600, res = 220)
pheatmap(corr, show_rownames = FALSE, show_colnames = FALSE, main = "Pseudobulk sample correlation")
dev.off()

cohort_df <- primary
cohort_df$donor_label <- reorder(cohort_df$donor_id, cohort_df$age)
p_cohort <- ggplot(cohort_df, aes(x = age, y = donor_label, colour = region, size = n_RPE_cells)) + geom_point(alpha = 0.85) + theme_bw(base_size = 9) + labs(title = "Primary donor-region structure", x = "Age (years)", y = "Donor", size = "RPE cells", colour = "Region")
ggsave(file.path(fig_dir, "P1-1_cohort_donor_region_structure.png"), p_cohort, width = 8, height = 9, dpi = 220)

pdat <- merge(a_tab, tiers[, c("gene_id", "evidence_tier")], by = "gene_id", all.x = TRUE)
p_volcano <- ggplot(pdat, aes(x = log2FC_MACULA_vs_PERIPHERY, y = -log10(pmax(FDR, .Machine$double.xmin)), colour = evidence_tier)) + geom_point(alpha = 0.65, size = 1.2) + geom_vline(xintercept = c(-0.25, 0.25), linetype = 2, colour = "grey50") + geom_hline(yintercept = -log10(0.05), linetype = 2, colour = "grey50") + theme_bw() + labs(title = "Model A regional effect", x = "log2FC (macula vs periphery)", y = "-log10 FDR", colour = "Evidence")
ggsave(file.path(fig_dir, "P1-3_model_A_volcano.png"), p_volcano, width = 7, height = 5.5, dpi = 220)

bb <- merge(a_tab[, c("gene_id", "log2FC_MACULA_vs_PERIPHERY")], b_tab[, c("gene_id", "log2FC_MACULA_vs_PERIPHERY")], by = "gene_id", suffixes = c("_A", "_B"))
p_beta <- ggplot(bb, aes(log2FC_MACULA_vs_PERIPHERY_A, log2FC_MACULA_vs_PERIPHERY_B)) + geom_point(alpha = 0.5, size = 1) + geom_abline(slope = 1, intercept = 0, linetype = 2) + theme_bw() + labs(title = "Model A versus paired Model B", x = "Model A log2FC", y = "Paired Model B log2FC")
ggsave(file.path(fig_dir, "P1-4_model_A_vs_paired_beta_beta.png"), p_beta, width = 6.5, height = 5.5, dpi = 220)

top1 <- tiers$gene_id[tiers$evidence_tier == "TIER_1_ROBUST_REGIONAL"]
if (length(top1) == 0) {
  top1 <- head(tiers$gene_id[order(tiers$Model_A_FDR, -abs(tiers$Model_A_log2FC), na.last = NA)], 20)
  heat_title <- "Top regional genes; no Tier 1 robust genes"
} else {
  top1 <- head(top1, 30)
  heat_title <- "Tier 1 robust regional genes"
}
heat <- logcpm[top1, , drop = FALSE]
heat <- t(scale(t(heat)))
png(file.path(fig_dir, "P1-5_top_regional_genes_heatmap.png"), width = 2200, height = 1700, res = 220)
pheatmap(heat, show_rownames = TRUE, show_colnames = FALSE, annotation_col = data.frame(Region = primary[colnames(heat), "region"], row.names = colnames(heat)), main = heat_title)
dev.off()

writeLines(c(
  paste0("primary_samples=", nrow(primary)),
  paste0("primary_donors=", length(unique(primary$donor_id))),
  paste0("paired_donors=", length(unique(paired$donor_id))),
  paste0("adult_paired_donors=", length(unique(adult_paired$donor_id))),
  paste0("age_overlap_min=", min(age_overlap$age), ";age_overlap_max=", max(age_overlap$age)),
  paste0("primary_gene_universe=", length(gene_universe)),
  paste0("model_a_design=", paste(a_fit$terms, collapse = "+")),
  paste0("model_a_duplicateCorrelation=", a_fit$correlation)
), file.path(log_dir, "MODEL_RUN_SUMMARY.txt"))
message("Phase 1 model outputs complete")
