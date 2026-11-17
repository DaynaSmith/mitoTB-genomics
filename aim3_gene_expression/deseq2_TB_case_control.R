#!/usr/bin/env Rscript
# =============================================================================
# Analysis 3: mtDNA gene DE between TB Case and Control, PE discovery / SE
# replication. Covariates: sex + age (paper's own TB-analysis methods
# weren't available to us -- this is a project-specific extension, so we
# use the same BH/FDR standard as Analysis 2, not the paper's Bonferroni
# convention which was specific to their eQTL/co-expression methods).
#
# L_status is included as an additional covariate ONLY IF the balance check
# (61_check_TB_covariate_balance.R) showed a meaningful Case/Control
# imbalance -- see INCLUDE_L_STATUS_COVARIATE below, set based on that
# script's output.
#
# Usage:
#   Rscript deseq2_TB_case_control.R
# =============================================================================

suppressPackageStartupMessages(library(DESeq2))

# ---- EDIT THESE PATHS ----
EXPR_MATRIX   <- "mtdna_expression_matrix.tsv"
PE_METADATA   <- "pe_metadata_unified.tsv"
SE_METADATA   <- "se_metadata_unified.tsv"
OUTDIR        <- "analysis3_TB_case_control"
INCLUDE_L_STATUS_COVARIATE <- FALSE  # EDIT based on 61's balance check output
# ---------------------------

dir.create(OUTDIR, showWarnings = FALSE)

expr <- read.delim(EXPR_MATRIX, stringsAsFactors = FALSE, check.names = FALSE)
gene_info <- expr[, c("gene_id", "gene_name", "gene_type")]
count_cols <- setdiff(colnames(expr), c("gene_id", "gene_name", "gene_type"))
count_mat_full <- as.matrix(expr[, count_cols])
rownames(count_mat_full) <- expr$gene_id
mode(count_mat_full) <- "numeric"

run_TB_deseq <- function(meta, label) {
  cat("\n=====", label, "=====\n")

  meta$L_status <- ifelse(startsWith(meta$haplogroup, "L"), "L", "non_L")
  meta$age <- suppressWarnings(as.numeric(meta$age))
  meta_cc <- meta[meta$condition %in% c("Case", "Control"), ]
  meta_cc$age_centered <- meta_cc$age - mean(meta_cc$age, na.rm = TRUE)

  sample_ids <- intersect(meta_cc$sample_id, colnames(count_mat_full))
  sub_mat <- count_mat_full[, sample_ids, drop = FALSE]
  complete <- colSums(is.na(sub_mat)) == 0
  n_dropped <- sum(!complete)
  sub_mat <- sub_mat[, complete, drop = FALSE]
  sample_ids <- colnames(sub_mat)
  storage.mode(sub_mat) <- "integer"

  extra_cols <- c("sex", "age_centered")
  if (INCLUDE_L_STATUS_COVARIATE) extra_cols <- c(extra_cols, "L_status")

  col_data <- meta_cc[match(sample_ids, meta_cc$sample_id),
                       c("sample_id", "condition", extra_cols)]
  rownames(col_data) <- col_data$sample_id
  col_data$condition <- factor(col_data$condition, levels = c("Control", "Case"))
  stopifnot(all(colnames(sub_mat) == rownames(col_data)))

  cat("Samples (", n_dropped, " dropped for missing counts ):", length(sample_ids), "\n", sep = "")
  print(table(col_data$condition))

  design_formula <- if (INCLUDE_L_STATUS_COVARIATE) {
    ~ sex + age_centered + L_status + condition
  } else {
    ~ sex + age_centered + condition
  }
  cat("Design:", deparse(design_formula), "\n")

  dds <- DESeqDataSetFromMatrix(countData = sub_mat, colData = col_data, design = design_formula)
  dds <- DESeq(dds)
  res <- results(dds, name = "condition_Case_vs_Control")
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  res_df <- merge(gene_info, res_df, by = "gene_id")
  res_df <- res_df[order(res_df$pvalue), ]

  out_path <- file.path(OUTDIR, paste0(label, "_results.tsv"))
  write.table(res_df, out_path, sep = "\t", quote = FALSE, row.names = FALSE)

  norm_counts <- counts(dds, normalized = TRUE)
  norm_counts_df <- as.data.frame(norm_counts)
  norm_counts_df$gene_id <- rownames(norm_counts_df)
  norm_counts_df <- merge(gene_info, norm_counts_df, by = "gene_id")
  write.table(norm_counts_df, file.path(OUTDIR, paste0(label, "_normalized_counts.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  n_sig <- sum(res_df$padj < 0.05, na.rm = TRUE)
  cat("Genes with padj < 0.05:", n_sig, "/", nrow(res_df), "\n")
  cat("Wrote:", out_path, "\n")
  res_df
}

pe_meta <- read.delim(PE_METADATA, stringsAsFactors = FALSE)
se_meta <- read.delim(SE_METADATA, stringsAsFactors = FALSE)

pe_res <- run_TB_deseq(pe_meta, "PE")
se_res <- run_TB_deseq(se_meta, "SE")

# ---- Replication check ----
cat("\n===== PE vs SE replication =====\n")
merged <- merge(pe_res[, c("gene_name", "log2FoldChange", "padj")],
                 se_res[, c("gene_name", "log2FoldChange", "padj")],
                 by = "gene_name", suffixes = c("_PE", "_SE"))
merged$same_direction <- sign(merged$log2FoldChange_PE) == sign(merged$log2FoldChange_SE)
merged <- merged[order(merged$padj_PE), ]
print(merged)
write.table(merged, file.path(OUTDIR, "PE_vs_SE_replication.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nDone. Check", OUTDIR, "for results.\n")