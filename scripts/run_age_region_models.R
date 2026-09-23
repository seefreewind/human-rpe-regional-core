#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
  library(ggplot2)
  library(pheatmap)
  library(splines)
})

project <- normalizePath(".", mustWork = TRUE)
counts_file <- file.path(project, "data/processed/pivot_phase1/RPE_PSEUDOBULK_COUNTS.tsv.gz")
meta_file <- file.path(project, "data/processed/pivot_phase1/RPE_PSEUDOBULK_METADATA.tsv")
phase1_tiers_file <- file.path(project, "results/pivot_phase1/REGIONAL_GENE_TIERS.tsv")
phase1_annot_file <- file.path(project, "results/pivot_phase1/REGIONAL_GENE_TIERS_ANNOTATED.tsv")
out_dir <- file.path(project, "results/pivot_phase2")
fig_dir <- file.path(project, "figures/pivot_phase2")
log_dir <- file.path(project, "logs/pivot_phase2")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

AGE_CENTER <- 53.0

write_tsv_gz <- function(x, path) {
  con <- gzfile(path, "wt")
  on.exit(close(con), add = TRUE)
  write.table(x, con, sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}

message("Reading Phase 1 pseudobulk and inherited gene universe")
counts_df <- read.delim(gzfile(counts_file), check.names = FALSE, stringsAsFactors = FALSE)
gene_ids <- make.unique(as.character(counts_df$gene_id))
counts_df$gene_id <- NULL
counts <- as.matrix(counts_df)
storage.mode(counts) <- "numeric"
rownames(counts) <- gene_ids
meta <- read.delim(meta_file, check.names = FALSE, stringsAsFactors = FALSE)
meta <- meta[match(colnames(counts), meta$sample_id), , drop = FALSE]
if (anyNA(meta$sample_id) || any(meta$sample_id != colnames(counts))) stop("Counts and metadata are not aligned")
phase1_tiers <- read.delim(phase1_tiers_file, check.names = FALSE, stringsAsFactors = FALSE)
gene_universe <- intersect(rownames(counts), unique(phase1_tiers$gene_id))
if (length(gene_universe) != 22430) stop("Inherited Phase 1 gene universe is not 22430 genes")
counts_u <- counts[gene_universe, , drop = FALSE]

meta$region <- factor(meta$region, levels = c("PERIPHERY", "MACULA"))
meta$sex <- droplevels(factor(meta$sex))
meta$ancestry_ethnicity <- droplevels(factor(meta$ancestry_ethnicity))
meta$assay_covariate <- droplevels(factor(ifelse(grepl(";", meta$assay), "MIXED_ASSAY", "10x_3p_v3")))
meta$age <- as.numeric(meta$age)
meta$Age_c <- meta$age - AGE_CENTER
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
  ans$region <- droplevels(ans$region)
  ans$sex <- droplevels(ans$sex)
  ans$ancestry_ethnicity <- droplevels(ans$ancestry_ethnicity)
  ans$assay_covariate <- droplevels(ans$assay_covariate)
  ans$donor_id <- droplevels(ans$donor_id)
  rownames(ans) <- ans$sample_id
  ans
}

build_linear_design <- function(d) {
  d$Region <- relevel(droplevels(factor(d$region)), ref = "PERIPHERY")
  d$Age_c <- d$age - AGE_CENTER
  d$Sex <- droplevels(factor(d$sex))
  d$Ancestry <- droplevels(factor(d$ancestry_ethnicity))
  d$Assay <- droplevels(factor(d$assay_covariate))
  optional <- c("Sex", "Ancestry", "Assay")
  retained <- optional[nlevels(d$Sex) > 1 | optional != "Sex"]
  retained <- optional[sapply(optional, function(x) nlevels(d[[x]]) > 1)]
  make <- function(opts) {
    formula_text <- paste("~ Region * Age_c", if (length(opts)) paste("+", paste(opts, collapse = " + ")) else "")
    model.matrix(as.formula(formula_text), d)
  }
  design <- make(retained)
  if (qr(design)$rank < ncol(design)) {
    for (v in rev(retained)) {
      trial_opts <- setdiff(retained, v)
      trial <- make(trial_opts)
      if (qr(trial)$rank == ncol(trial)) {
        retained <- trial_opts
        design <- trial
        break
      }
    }
  }
  int_name <- grep("^RegionMACULA:Age_c$", colnames(design), value = TRUE)
  if (length(int_name) != 1 || !all(c("RegionMACULA", "Age_c") %in% colnames(design)) || qr(design)$rank < ncol(design)) {
    stop("Mandatory linear Age x Region design is not estimable")
  }
  list(data = d, design = design, interaction = int_name, retained_optional = retained, terms = colnames(design))
}

build_linear_contrasts <- function(design, interaction) {
  L <- matrix(0, nrow = ncol(design), ncol = 3, dimnames = list(colnames(design), c("beta_age_periphery", "beta_age_macula", "beta_interaction")))
  L["Age_c", c("beta_age_periphery", "beta_age_macula")] <- 1
  L[interaction, c("beta_age_macula", "beta_interaction")] <- 1
  L
}

extract_linear_table <- function(fit_raw, design_info, genes) {
  L <- build_linear_contrasts(design_info$design, design_info$interaction)
  fit_c <- contrasts.fit(fit_raw, L)
  fit_c <- eBayes(fit_c, robust = TRUE)
  tt <- topTable(fit_c, coef = "beta_interaction", number = Inf, sort.by = "none")
  df <- fit_c$df.total
  beta_i <- as.numeric(fit_c$coefficients[, "beta_interaction"])
  beta_p <- as.numeric(fit_c$coefficients[, "beta_age_periphery"])
  beta_m <- as.numeric(fit_c$coefficients[, "beta_age_macula"])
  se_i <- as.numeric(fit_c$stdev.unscaled[, "beta_interaction"] * sqrt(fit_c$s2.post))
  se_p <- as.numeric(fit_c$stdev.unscaled[, "beta_age_periphery"] * sqrt(fit_c$s2.post))
  se_m <- as.numeric(fit_c$stdev.unscaled[, "beta_age_macula"] * sqrt(fit_c$s2.post))
  data.frame(
    gene_id = genes,
    beta_age_periphery = beta_p,
    SE_age_periphery = se_p,
    CI95_low_age_periphery = beta_p - qt(0.975, df = df) * se_p,
    CI95_high_age_periphery = beta_p + qt(0.975, df = df) * se_p,
    beta_age_macula = beta_m,
    SE_age_macula = se_m,
    CI95_low_age_macula = beta_m - qt(0.975, df = df) * se_m,
    CI95_high_age_macula = beta_m + qt(0.975, df = df) * se_m,
    beta_interaction = beta_i,
    SE_interaction = se_i,
    CI95_low_interaction = beta_i - qt(0.975, df = df) * se_i,
    CI95_high_interaction = beta_i + qt(0.975, df = df) * se_i,
    P_interaction = as.numeric(tt$P.Value),
    BH_FDR_interaction = as.numeric(tt$adj.P.Val),
    stringsAsFactors = FALSE
  )
}

fit_linear_model <- function(d, label, counts_matrix = counts_u, genes = gene_universe, correlation = NULL) {
  di <- build_linear_design(d)
  ids <- di$data$sample_id
  y <- DGEList(counts = counts_matrix[, ids, drop = FALSE])
  y <- calcNormFactors(y)
  v0 <- voom(y, di$design, plot = FALSE)
  corfit <- duplicateCorrelation(v0, di$design, block = di$data$donor_id)
  cor_value <- if (is.null(correlation)) as.numeric(corfit$consensus) else correlation
  if (!is.finite(cor_value)) cor_value <- 0
  v <- voom(y, di$design, plot = FALSE, block = di$data$donor_id, correlation = cor_value)
  fit_raw <- lmFit(v, di$design, block = di$data$donor_id, correlation = cor_value)
  tab <- extract_linear_table(fit_raw, di, genes)
  tab$model <- label
  tab$n_samples <- nrow(di$data)
  tab$n_donors <- length(unique(di$data$donor_id))
  tab$n_macula <- sum(di$data$region == "MACULA")
  tab$n_periphery <- sum(di$data$region == "PERIPHERY")
  tab$age_center_years <- AGE_CENTER
  tab$design_terms <- paste(di$terms, collapse = " + ")
  list(table = tab, data = di$data, design_info = di, y = y, voom = v, fit_raw = fit_raw, correlation = cor_value)
}

fit_paired_delta <- function(d, label = "MODEL_PAIRED_DELTA_AGE") {
  paired_donors <- names(which(table(d$donor_id) == 2))
  pd <- d[as.character(d$donor_id) %in% paired_donors, , drop = FALSE]
  pd <- pd[order(pd$donor_id, pd$region), , drop = FALSE]
  if (length(unique(pd$donor_id)) < 5) stop("Too few paired donors for paired-delta model")
  ids <- pd$sample_id
  y <- DGEList(counts = counts_u[, ids, drop = FALSE])
  y <- calcNormFactors(y)
  logcpm <- cpm(y, log = TRUE, prior.count = 2)
  donors <- unique(as.character(pd$donor_id))
  delta <- sapply(donors, function(x) {
    mac <- pd$sample_id[as.character(pd$donor_id) == x & pd$region == "MACULA"]
    per <- pd$sample_id[as.character(pd$donor_id) == x & pd$region == "PERIPHERY"]
    if (length(mac) != 1 || length(per) != 1) stop("Paired donor has unexpected regions")
    logcpm[, mac] - logcpm[, per]
  })
  colnames(delta) <- donors
  dm <- pd[match(donors, as.character(pd$donor_id)), , drop = FALSE]
  dm$Donor <- droplevels(factor(dm$donor_id))
  dm$Sex <- droplevels(factor(dm$sex))
  dm$Ancestry <- droplevels(factor(dm$ancestry_ethnicity))
  dm$Age_c <- dm$age - AGE_CENTER
  opts <- c("Sex", "Ancestry")[sapply(c("Sex", "Ancestry"), function(x) nlevels(dm[[x]]) > 1)]
  make <- function(o) model.matrix(as.formula(paste("~ Age_c", if (length(o)) paste("+", paste(o, collapse = " + ")) else "")), dm)
  design <- make(opts)
  if (qr(design)$rank < ncol(design)) {
    for (v in rev(opts)) {
      trial <- make(setdiff(opts, v))
      if (qr(trial)$rank == ncol(trial)) {
        opts <- setdiff(opts, v)
        design <- trial
        break
      }
    }
  }
  if (!("Age_c" %in% colnames(design)) || qr(design)$rank < ncol(design)) stop("Paired delta design is not estimable")
  fit <- eBayes(lmFit(delta, design), robust = TRUE)
  tt <- topTable(fit, coef = "Age_c", number = Inf, sort.by = "none")
  beta <- as.numeric(fit$coefficients[, "Age_c"])
  se <- as.numeric(fit$stdev.unscaled[, "Age_c"] * sqrt(fit$s2.post))
  df <- fit$df.total
  tab <- data.frame(
    gene_id = rownames(delta),
    beta_paired_delta_age = beta,
    beta_interaction = beta,
    SE_paired_delta_age = se,
    CI95_low_paired_delta_age = beta - qt(0.975, df = df) * se,
    CI95_high_paired_delta_age = beta + qt(0.975, df = df) * se,
    P_paired_delta_age = as.numeric(tt$P.Value),
    FDR_paired_delta_age = as.numeric(tt$adj.P.Val),
    model = label,
    n_samples = ncol(delta),
    n_donors = length(donors),
    n_macula = length(donors),
    n_periphery = length(donors),
    age_center_years = AGE_CENTER,
    design_terms = paste(colnames(design), collapse = " + "),
    stringsAsFactors = FALSE
  )
  list(table = tab, data = dm, delta = delta, design = design)
}

fit_spline_diagnostic <- function(d, counts_matrix = counts_u, genes = gene_universe) {
  di <- build_linear_design(d)
  dd <- di$data
  dd$Region <- relevel(droplevels(factor(dd$region)), ref = "PERIPHERY")
  dd$Age_c <- dd$age - AGE_CENTER
  dd$Sex <- droplevels(factor(dd$sex))
  dd$Ancestry <- droplevels(factor(dd$ancestry_ethnicity))
  dd$Assay <- droplevels(factor(dd$assay_covariate))
  optional <- c("Sex", "Ancestry", "Assay")
  optional <- optional[sapply(optional, function(x) nlevels(dd[[x]]) > 1)]
  make <- function(opts) {
    formula_text <- paste("~ Region * splines::ns(Age_c, df=3)", if (length(opts)) paste("+", paste(opts, collapse = " + ")) else "")
    model.matrix(as.formula(formula_text), dd)
  }
  design <- make(optional)
  if (qr(design)$rank < ncol(design)) {
    for (v in rev(optional)) {
      trial <- make(setdiff(optional, v))
      if (qr(trial)$rank == ncol(trial)) {
        optional <- setdiff(optional, v)
        design <- trial
        break
      }
    }
  }
  int_cols <- grep("^RegionMACULA:.*ns", colnames(design), value = TRUE)
  if (length(int_cols) != 3 || qr(design)$rank < ncol(design)) stop("Spline diagnostic design is not estimable")
  y <- DGEList(counts = counts_matrix[, dd$sample_id, drop = FALSE])
  y <- calcNormFactors(y)
  v0 <- voom(y, design, plot = FALSE)
  corfit <- duplicateCorrelation(v0, design, block = dd$donor_id)
  cor_value <- as.numeric(corfit$consensus)
  if (!is.finite(cor_value)) cor_value <- 0
  v <- voom(y, design, plot = FALSE, block = dd$donor_id, correlation = cor_value)
  fit_raw <- lmFit(v, design, block = dd$donor_id, correlation = cor_value)
  L <- matrix(0, nrow = ncol(design), ncol = length(int_cols), dimnames = list(colnames(design), int_cols))
  for (i in seq_along(int_cols)) L[int_cols[i], i] <- 1
  fit_c <- eBayes(contrasts.fit(fit_raw, L), robust = TRUE)
  tt <- topTableF(fit_c, number = Inf, sort.by = "none")
  data.frame(gene_id = genes, spline_F = tt$F, spline_P = tt$P.Value, spline_FDR = tt$adj.P.Val, spline_design_terms = paste(colnames(design), collapse = " + "), stringsAsFactors = FALSE)
}

add_direction <- function(primary_beta, x) {
  is.finite(primary_beta) & is.finite(x) & sign(primary_beta) == sign(x) & sign(primary_beta) != 0
}

write_model <- function(tab, path) {
  fixed <- c("gene_id", "beta_age_periphery", "SE_age_periphery", "CI95_low_age_periphery", "CI95_high_age_periphery", "beta_age_macula", "SE_age_macula", "CI95_low_age_macula", "CI95_high_age_macula", "beta_interaction", "SE_interaction", "CI95_low_interaction", "CI95_high_interaction", "P_interaction", "BH_FDR_interaction", "model", "n_samples", "n_donors", "n_macula", "n_periphery", "age_center_years", "design_terms")
  write_tsv_gz(tab[, intersect(fixed, names(tab)), drop = FALSE], path)
}

primary <- make_subset(meta, 30)
message("Fitting primary Age x Region model")
primary_fit <- fit_linear_model(primary, "MODEL_AGE_REGION_PRIMARY")
primary_tab <- primary_fit$table
write_model(primary_tab, file.path(out_dir, "MODEL_AGE_REGION_PRIMARY.tsv.gz"))

message("Fitting adult-only, age-overlap, cell-count and mixed-assay models")
adult_fit <- fit_linear_model(make_subset(meta, 30, adult_only = TRUE), "MODEL_AGE_REGION_ADULT")
adult_tab <- adult_fit$table
write_model(adult_tab, file.path(out_dir, "MODEL_AGE_REGION_ADULT.tsv.gz"))

overlap_fit <- fit_linear_model(make_subset(meta, 30, age_overlap = TRUE), "MODEL_AGE_REGION_OVERLAP")
overlap_tab <- overlap_fit$table
write_model(overlap_tab, file.path(out_dir, "MODEL_AGE_REGION_OVERLAP.tsv.gz"))

cell20_fit <- fit_linear_model(make_subset(meta, 20), "MODEL_AGE_REGION_CELLCOUNT_20")
cell20_tab <- cell20_fit$table
write_model(cell20_tab, file.path(out_dir, "MODEL_AGE_REGION_CELLCOUNT_20.tsv.gz"))

cell50_fit <- fit_linear_model(make_subset(meta, 50), "MODEL_AGE_REGION_CELLCOUNT_50")
cell50_tab <- cell50_fit$table
write_model(cell50_tab, file.path(out_dir, "MODEL_AGE_REGION_CELLCOUNT_50.tsv.gz"))

mixed_fit <- fit_linear_model(make_subset(meta, 30, exclude_mixed = TRUE), "MODEL_AGE_REGION_MIXED_ASSAY_EXCLUDED")
mixed_tab <- mixed_fit$table
write_model(mixed_tab, file.path(out_dir, "MODEL_AGE_REGION_MIXED_ASSAY_EXCLUDED.tsv.gz"))

paired_fit <- fit_paired_delta(primary)
paired_tab <- paired_fit$table
write_tsv_gz(paired_tab, file.path(out_dir, "MODEL_PAIRED_DELTA_AGE.tsv.gz"))

slopes <- primary_tab[, c("gene_id", "beta_age_periphery", "SE_age_periphery", "CI95_low_age_periphery", "CI95_high_age_periphery", "beta_age_macula", "SE_age_macula", "CI95_low_age_macula", "CI95_high_age_macula", "beta_interaction", "SE_interaction", "CI95_low_interaction", "CI95_high_interaction", "P_interaction", "BH_FDR_interaction")]
write_tsv_gz(slopes, file.path(out_dir, "REGION_SPECIFIC_AGE_SLOPES.tsv.gz"))

message("Running natural-spline nonlinearity diagnostic")
spline_tab <- fit_spline_diagnostic(primary)
diag_tab <- merge(primary_tab[, c("gene_id", "P_interaction", "BH_FDR_interaction", "beta_interaction")], spline_tab, by = "gene_id", sort = FALSE)
diag_tab$diagnostic_class <- "NO_PRIMARY_INTERACTION"
diag_tab$diagnostic_class[diag_tab$BH_FDR_interaction < 0.05 & diag_tab$spline_FDR < 0.05] <- "LINEAR_AND_SPLINE_SUPPORTED"
diag_tab$diagnostic_class[diag_tab$BH_FDR_interaction < 0.05 & diag_tab$spline_FDR >= 0.05] <- "NONLINEAR_AGE_DEPENDENCE"
diag_tab$diagnostic_class[diag_tab$BH_FDR_interaction >= 0.05 & diag_tab$spline_FDR < 0.05] <- "SPLINE_ONLY_DIAGNOSTIC"
write.table(diag_tab, file.path(out_dir, "NONLINEAR_AGE_DIAGNOSTIC.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

concordance_row <- function(x, y, label, beta_x = "beta_interaction", beta_y = "beta_interaction", fdr_x = "BH_FDR_interaction") {
  z <- merge(x[, c("gene_id", beta_x, fdr_x)], y[, c("gene_id", beta_y)], by = "gene_id", sort = FALSE)
  names(z) <- c("gene_id", "beta_primary", "fdr_primary", "beta_sensitivity")
  z <- z[is.finite(z$beta_primary) & is.finite(z$beta_sensitivity), , drop = FALSE]
  hit <- z$fdr_primary < 0.05
  data.frame(
    comparison = label,
    n_genes = nrow(z),
    spearman_rho = if (nrow(z) >= 3) suppressWarnings(cor(z$beta_primary, z$beta_sensitivity, method = "spearman")) else NA_real_,
    direction_concordance_all = if (nrow(z)) mean(sign(z$beta_primary) == sign(z$beta_sensitivity) & sign(z$beta_primary) != 0) else NA_real_,
    n_primary_fdr05 = sum(hit),
    same_direction_among_primary_hits = if (sum(hit)) mean(sign(z$beta_primary[hit]) == sign(z$beta_sensitivity[hit]) & sign(z$beta_primary[hit]) != 0) else NA_real_,
    stringsAsFactors = FALSE
  )
}
conc <- rbind(
  concordance_row(primary_tab, adult_tab, "PRIMARY_vs_ADULT"),
  concordance_row(primary_tab, overlap_tab, "PRIMARY_vs_AGE_OVERLAP"),
  concordance_row(primary_tab, paired_tab, "PRIMARY_vs_PAIRED_DELTA", beta_y = "beta_interaction"),
  concordance_row(primary_tab, cell20_tab, "PRIMARY_vs_CELLCOUNT_20"),
  concordance_row(primary_tab, cell50_tab, "PRIMARY_vs_CELLCOUNT_50"),
  concordance_row(primary_tab, mixed_tab, "PRIMARY_vs_MIXED_ASSAY_EXCLUDED")
)
write.table(conc, file.path(out_dir, "AGE_INTERACTION_SENSITIVITY_CONCORDANCE.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

message("Building interaction tiers and cross-classification")
merge_beta <- function(x, y, suffix) {
  out <- y[, c("gene_id", "beta_interaction", "P_interaction", "BH_FDR_interaction")]
  names(out)[-1] <- paste0(names(out)[-1], suffix)
  merge(x, out, by = "gene_id", all.x = TRUE, sort = FALSE)
}
tiers <- primary_tab
names(tiers)[names(tiers) == "P_interaction"] <- "primary_P_interaction"
names(tiers)[names(tiers) == "BH_FDR_interaction"] <- "primary_FDR_interaction"
tiers <- merge_beta(tiers, adult_tab, "_adult")
tiers <- merge_beta(tiers, overlap_tab, "_overlap")
tiers <- merge_beta(tiers, cell20_tab, "_cell20")
tiers <- merge_beta(tiers, cell50_tab, "_cell50")
tiers <- merge_beta(tiers, mixed_tab, "_mixed")
paired_merge <- paired_tab[, c("gene_id", "beta_interaction", "P_paired_delta_age", "FDR_paired_delta_age")]
names(paired_merge) <- c("gene_id", "beta_interaction_paired_delta", "P_paired_delta", "FDR_paired_delta")
tiers <- merge(tiers, paired_merge, by = "gene_id", all.x = TRUE, sort = FALSE)

sens_cols <- c("beta_interaction_adult", "beta_interaction_overlap", "beta_interaction_paired_delta", "beta_interaction_cell20", "beta_interaction_cell50", "beta_interaction_mixed")
same_cols <- paste0("same_direction_", c("adult", "overlap", "paired_delta", "cell20", "cell50", "mixed"))
for (i in seq_along(sens_cols)) tiers[[same_cols[i]]] <- add_direction(tiers$beta_interaction, tiers[[sens_cols[i]]])
tiers$n_sensitivity_same <- rowSums(tiers[, same_cols, drop = FALSE], na.rm = TRUE)
tiers$n_sensitivity_estimable <- rowSums(is.finite(as.matrix(tiers[, sens_cols, drop = FALSE])))
tiers$major_sensitivity_opposite <- rowSums(!tiers[, paste0("same_direction_", c("adult", "overlap", "paired_delta")), drop = FALSE], na.rm = TRUE)
tiers$interaction_significant <- is.finite(tiers$primary_FDR_interaction) & tiers$primary_FDR_interaction < 0.05
tiers$interaction_class <- "NOT_PRIMARY_SIGNIFICANT"
sig <- tiers$interaction_significant
divergent <- sig & is.finite(tiers$beta_age_periphery) & is.finite(tiers$beta_age_macula) & (tiers$beta_age_periphery * tiers$beta_age_macula < 0)
tiers$interaction_class[divergent] <- "DIVERGENT"
tiers$interaction_class[sig & !divergent & tiers$beta_interaction > 0 & tiers$beta_age_periphery >= 0 & tiers$beta_age_macula >= 0] <- "MACULA_ACCELERATED_INCREASE"
tiers$interaction_class[sig & !divergent & tiers$beta_interaction < 0 & tiers$beta_age_periphery <= 0 & tiers$beta_age_macula <= 0] <- "MACULA_ACCELERATED_DECLINE"
tiers$interaction_class[sig & tiers$interaction_class == "NOT_PRIMARY_SIGNIFICANT"] <- "DIFFERENTIAL_MAGNITUDE"

lodo_candidates <- tiers$gene_id[tiers$interaction_significant]
lodo_rows <- list()
if (length(lodo_candidates) > 0) {
  donors <- unique(as.character(primary$donor_id))
  full_sign <- sign(tiers$beta_interaction[match(lodo_candidates, tiers$gene_id)])
  for (donor in donors) {
    keep_sample <- as.character(primary$donor_id) != donor
    d2 <- primary[keep_sample, , drop = FALSE]
    di2 <- build_linear_design(d2)
    gene_idx <- match(lodo_candidates, rownames(primary_fit$voom$E))
    sample_idx <- match(d2$sample_id, colnames(primary_fit$voom$E))
    E2 <- primary_fit$voom$E[gene_idx, sample_idx, drop = FALSE]
    W2 <- primary_fit$voom$weights[gene_idx, sample_idx, drop = FALSE]
    f2 <- tryCatch(lmFit(E2, di2$design, weights = W2, block = di2$data$donor_id, correlation = primary_fit$correlation), error = function(e) NULL)
    beta <- rep(NA_real_, length(lodo_candidates))
    if (!is.null(f2)) {
      L2 <- build_linear_contrasts(di2$design, di2$interaction)
      beta <- as.numeric(contrasts.fit(f2, L2)$coefficients[, "beta_interaction"])
    }
    lodo_rows[[length(lodo_rows) + 1]] <- data.frame(gene_id = lodo_candidates, donor_left_out = donor, lodo_beta_interaction = beta, full_model_beta_interaction = tiers$beta_interaction[match(lodo_candidates, tiers$gene_id)], direction_flip = is.finite(beta) & sign(beta) != full_sign, stringsAsFactors = FALSE)
    message("LODO donor ", donor)
  }
}
if (length(lodo_rows) > 0) {
  lodo <- do.call(rbind, lodo_rows)
  lodo_summary <- do.call(rbind, lapply(lodo_candidates, function(g) {
    z <- lodo[lodo$gene_id == g, , drop = FALSE]
    ok <- is.finite(z$lodo_beta_interaction)
    data.frame(gene_id = g, full_model_beta_interaction = z$full_model_beta_interaction[1], min_lodo_beta_interaction = if (any(ok)) min(z$lodo_beta_interaction[ok]) else NA_real_, max_lodo_beta_interaction = if (any(ok)) max(z$lodo_beta_interaction[ok]) else NA_real_, n_direction_flips = sum(z$direction_flip, na.rm = TRUE), n_lodo_fits = sum(ok), LODO_status = if (any(ok) && sum(z$direction_flip, na.rm = TRUE) == 0) "STABLE" else "LODO_UNSTABLE_OR_UNESTIMABLE", stringsAsFactors = FALSE)
  }))
} else {
  lodo_summary <- data.frame(gene_id = character(), full_model_beta_interaction = numeric(), min_lodo_beta_interaction = numeric(), max_lodo_beta_interaction = numeric(), n_direction_flips = integer(), n_lodo_fits = integer(), LODO_status = character())
}
write.table(lodo_summary, file.path(out_dir, "TIER_AB_LODO_STABILITY.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
lodo_status <- lodo_summary$LODO_status[match(tiers$gene_id, lodo_summary$gene_id)]
tiers$LODO_status <- lodo_status
tiers$Tier <- "NOT_PRIMARY_SIGNIFICANT"
tiers$Tier[sig & tiers$n_sensitivity_same >= 4 & tiers$major_sensitivity_opposite <= 1 & (is.na(tiers$LODO_status) | tiers$LODO_status %in% c("STABLE", "LODO_UNSTABLE_OR_UNESTIMABLE"))] <- "TIER_B_SUPPORTED_AGE_LABILE"
all_six <- sig & tiers$n_sensitivity_estimable == 6 & tiers$n_sensitivity_same == 6 & tiers$LODO_status == "STABLE"
tiers$Tier[all_six] <- "TIER_A_ROBUST_MACULAR_AGE_LABILE"

annot <- if (file.exists(phase1_annot_file)) read.delim(phase1_annot_file, check.names = FALSE, stringsAsFactors = FALSE)[, c("gene_id", "gene_symbol")] else data.frame(gene_id = gene_universe, gene_symbol = gene_universe)
tiers <- merge(tiers, annot, by = "gene_id", all.x = TRUE, sort = FALSE)
tiers <- tiers[order(tiers$Tier, tiers$primary_FDR_interaction, na.last = TRUE), , drop = FALSE]
write.table(tiers, file.path(out_dir, "AGE_LABILE_GENE_TIERS.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

phase1_class <- phase1_tiers[, c("gene_id", "evidence_tier")]
phase1_class$phase1_regional <- phase1_class$evidence_tier %in% c("TIER_1_ROBUST_REGIONAL", "TIER_2_SUPPORTED_REGIONAL")
cross <- merge(tiers[, c("gene_id", "gene_symbol", "primary_FDR_interaction", "Tier", "interaction_class")], phase1_class, by = "gene_id", all.x = TRUE, sort = FALSE)
cross$phase2_interaction <- cross$primary_FDR_interaction < 0.05
cross$cross_classification <- "NEITHER"
cross$cross_classification[cross$phase1_regional & !cross$phase2_interaction] <- "REGION_ONLY"
cross$cross_classification[!cross$phase1_regional & cross$phase2_interaction] <- "AGE_INTERACTION_ONLY"
cross$cross_classification[cross$phase1_regional & cross$phase2_interaction] <- "REGION_AND_AGE_INTERACTION"
write.table(cross, file.path(out_dir, "REGION_AGE_CROSSCLASSIFICATION.tsv"), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")

message("Writing model summary")
writeLines(c(
  paste0("primary_samples=", nrow(primary)),
  paste0("primary_donors=", length(unique(primary$donor_id))),
  paste0("primary_macula=", sum(primary$region == "MACULA")),
  paste0("primary_periphery=", sum(primary$region == "PERIPHERY")),
  paste0("paired_donors=", length(unique(as.character(paired_fit$data$donor_id)))),
  paste0("adult_samples=", nrow(adult_fit$data)),
  paste0("adult_paired_donors=", length(unique(as.character(adult_fit$data$donor_id[adult_fit$data$donor_id %in% adult_fit$data$donor_id[duplicated(adult_fit$data$donor_id)]])))),
  paste0("age_center_years=", AGE_CENTER),
  paste0("genes_tested=", length(gene_universe)),
  paste0("primary_design_terms=", primary_tab$design_terms[1]),
  paste0("primary_duplicateCorrelation=", primary_fit$correlation),
  paste0("primary_fdr05=", sum(primary_tab$BH_FDR_interaction < 0.05)),
  paste0("tier_a=", sum(tiers$Tier == "TIER_A_ROBUST_MACULAR_AGE_LABILE")),
  paste0("tier_b=", sum(tiers$Tier == "TIER_B_SUPPORTED_AGE_LABILE")),
  paste0("tier_c=", sum(tiers$Tier == "NOT_PRIMARY_SIGNIFICANT" & tiers$interaction_significant))
), file.path(log_dir, "MODEL_RUN_SUMMARY.txt"))
message("Phase 2 interaction outputs complete")

