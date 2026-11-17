#!/usr/bin/env Rscript
# =============================================================================
# Analysis 2: L vs non-L haplogroup differential expression (PE cohort)
#
# Runs TWO models, per the agreed design:
#   (a) FULL N (~410 samples): design ~ L_status, no sex/age covariate
#   (b) RESTRICTED N (~238 samples with complete clinical data):
#       design ~ sex + age + L_status
#
# L vs non-L: any haplogroup starting with "L" (L0-L6, the African-lineage
# macrohaplogroup in PhyloTree nomenclature) is "L"; everything else is
# "non-L". No other haplogroup family uses this letter, so this split is
# unambiguous directly from the haplogroup string.
#
# Inputs:
#   mtdna_expression_matrix.tsv  — 37 genes x samples (from 27_aggregate_mtdna_counts.py)
#   pe_metadata_unified.tsv      — sample_id, cohort, condition, sex, age,
#                                   haplogroup, haplogroup_quality
#
# Usage:
#   Rscript deseq2_PE_L_vs_nonL.R
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
})

# ---- EDIT THESE PATHS ----
EXPR_MATRIX   <- "mtdna_expression_matrix.tsv"
METADATA      <- "pe_metadata_unified.tsv"
OUTDIR        <- "analysis2_L_vs_nonL"
# ---------------------------

dir.create(OUTDIR, showWarnings = FALSE)

# ---- Load data ----
expr <- read.delim(EXPR_MATRIX, stringsAsFactors = FALSE, check.names = FALSE)
meta <- read.delim(METADATA, stringsAsFactors = FALSE)

cat("Expression matrix:", nrow(expr), "genes x", ncol(expr) - 3, "samples\n")
cat("Metadata:", nrow(meta), "PE samples\n")

# ---- Derive L vs non-L grouping ----
meta$L_status <- ifelse(startsWith(meta$haplogroup, "L"), "L", "non_L")
cat("\nL vs non-L counts (all PE samples):\n")
print(table(meta$L_status))

# ---- Build count matrix: genes as rows, samples as columns ----
gene_info <- expr[, c("gene_id", "gene_name", "gene_type")]
count_cols <- setdiff(colnames(expr), c("gene_id", "gene_name", "gene_type"))
count_mat_full <- as.matrix(expr[, count_cols])
rownames(count_mat_full) <- expr$gene_id
mode(count_mat_full) <- "numeric"  # coerces "NA" strings to NA

run_deseq <- function(sample_ids, design_formula, colData_extra_cols, label) {
  cat("\n=====", label, "=====\n")

  # Restrict to samples present in both the count matrix and this sample list,
  # and with no NA counts (drops the handful of samples missing personalized
  # remap output)
  sample_ids <- intersect(sample_ids, colnames(count_mat_full))
  sub_mat <- count_mat_full[, sample_ids, drop = FALSE]
  complete <- colSums(is.na(sub_mat)) == 0
  n_dropped_na <- sum(!complete)
  sub_mat <- sub_mat[, complete, drop = FALSE]
  sample_ids <- colnames(sub_mat)

  cat("Samples after NA filtering:", length(sample_ids),
      "(", n_dropped_na, "dropped for missing counts )\n")

  storage.mode(sub_mat) <- "integer"

  col_data <- meta[match(sample_ids, meta$sample_id), c("sample_id", "L_status", colData_extra_cols)]
  rownames(col_data) <- col_data$sample_id
  col_data$L_status <- factor(col_data$L_status, levels = c("non_L", "L"))  # non_L = reference

  stopifnot(all(colnames(sub_mat) == rownames(col_data)))

  cat("Design:", deparse(design_formula), "\n")
  cat("Group counts in this model:\n")
  print(table(col_data$L_status))

  dds <- DESeqDataSetFromMatrix(countData = sub_mat, colData = col_data, design = design_formula)
  dds <- DESeq(dds)
  res <- results(dds, name = "L_status_L_vs_non_L")
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  res_df <- merge(gene_info, res_df, by = "gene_id")
  res_df <- res_df[order(res_df$pvalue), ]

  out_path <- file.path(OUTDIR, paste0(label, "_results.tsv"))
  write.table(res_df, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("Wrote:", out_path, "\n")

  # Also export DESeq2-normalized counts (library-size-corrected), needed
  # for plotting expression levels — raw counts aren't comparable across
  # samples with different sequencing depth.
  norm_counts <- counts(dds, normalized = TRUE)
  norm_counts_df <- as.data.frame(norm_counts)
  norm_counts_df$gene_id <- rownames(norm_counts_df)
  norm_counts_df <- merge(gene_info, norm_counts_df, by = "gene_id")
  norm_out_path <- file.path(OUTDIR, paste0(label, "_normalized_counts.tsv"))
  write.table(norm_counts_df, norm_out_path, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("Wrote:", norm_out_path, "\n")

  n_sig <- sum(res_df$padj < 0.05, na.rm = TRUE)
  cat("Genes with padj < 0.05:", n_sig, "/", nrow(res_df), "\n")

  res_df
}

# ---- Model (a): FULL N, no covariates ----
full_ids <- meta$sample_id
res_full <- run_deseq(full_ids, ~ L_status, character(0), "full_N")

# ---- Model (b): RESTRICTED N, with sex + age covariates ----
restricted_ids <- meta$sample_id[meta$condition != "NA" & !is.na(meta$condition)]
meta$age <- suppressWarnings(as.numeric(meta$age))
meta$age_centered <- meta$age - mean(meta$age, na.rm = TRUE)  # centering silences
                                                                 # DESeq2's GLM
                                                                 # convergence warning
res_restricted <- run_deseq(restricted_ids, ~ sex + age_centered + L_status,
                             c("sex", "age_centered"), "restricted_N")

# ---- Sensitivity check: exclude the 19 low-HaploGrep3-quality samples ----
# (haplogroup_quality < 0.90 — see conversation notes on the 19 flagged
# samples from the original HaploGrep3 QC step)
low_quality_ids <- meta$sample_id[!is.na(meta$haplogroup_quality) &
                                     as.numeric(meta$haplogroup_quality) < 0.90]
cat("\nExcluding", length(low_quality_ids), "low-haplogroup-quality samples for sensitivity check\n")

full_ids_hq <- setdiff(full_ids, low_quality_ids)
res_full_hq <- run_deseq(full_ids_hq, ~ L_status, character(0), "full_N_highQualityOnly")

restricted_ids_hq <- setdiff(restricted_ids, low_quality_ids)
res_restricted_hq <- run_deseq(restricted_ids_hq, ~ sex + age_centered + L_status,
                                c("sex", "age_centered"), "restricted_N_highQualityOnly")

# ---- Compare full_N vs restricted_N ----
cat("\n===== Comparing full_N vs restricted_N results =====\n")
merged <- merge(res_full[, c("gene_name", "log2FoldChange", "padj")],
                 res_restricted[, c("gene_name", "log2FoldChange", "padj")],
                 by = "gene_name", suffixes = c("_full", "_restricted"))
merged <- merged[order(merged$padj_full), ]
print(merged)
write.table(merged, file.path(OUTDIR, "full_vs_restricted_comparison.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# ---- Compare full_N vs full_N with low-quality samples excluded ----
cat("\n===== Sensitivity check: full_N vs full_N (high quality only) =====\n")
merged_sens <- merge(res_full[, c("gene_name", "log2FoldChange", "padj")],
                      res_full_hq[, c("gene_name", "log2FoldChange", "padj")],
                      by = "gene_name", suffixes = c("_all", "_highQualityOnly"))
merged_sens <- merged_sens[order(merged_sens$padj_all), ]
print(merged_sens)
write.table(merged_sens, file.path(OUTDIR, "full_N_quality_sensitivity_comparison.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nDone. Check", OUTDIR, "for results.\n")
