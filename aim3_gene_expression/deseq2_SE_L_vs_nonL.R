#!/usr/bin/env Rscript
# =============================================================================
# SE replication of Analysis 2 (L vs non-L DE) -- mirrors 30_deseq2_L_vs_nonL.R's
# PE design exactly, applied to the SE cohort. This was flagged complete
# earlier but never actually run.
#
# Usage:
#   Rscript deseq2_SE_L_vs_nonL.R
# =============================================================================

suppressPackageStartupMessages(library(DESeq2))

# ---- EDIT THESE PATHS ----
EXPR_MATRIX   <- "mtdna_expression_matrix.tsv"
SE_METADATA   <- "se_metadata_unified.tsv"
OUTDIR        <- "analysis2_L_vs_nonL"
# ---------------------------

dir.create(OUTDIR, showWarnings = FALSE)

expr <- read.delim(EXPR_MATRIX, stringsAsFactors = FALSE, check.names = FALSE)
se_meta <- read.delim(SE_METADATA, stringsAsFactors = FALSE)

se_meta$L_status <- ifelse(startsWith(se_meta$haplogroup, "L"), "L", "non_L")
cat("SE L vs non-L counts:\n")
print(table(se_meta$L_status))

gene_info <- expr[, c("gene_id", "gene_name", "gene_type")]
count_cols <- setdiff(colnames(expr), c("gene_id", "gene_name", "gene_type"))
count_mat_full <- as.matrix(expr[, count_cols])
rownames(count_mat_full) <- expr$gene_id
mode(count_mat_full) <- "numeric"

run_deseq <- function(sample_ids, design_formula, colData_extra_cols, label) {
  cat("\n=====", label, "=====\n")
  sample_ids <- intersect(sample_ids, colnames(count_mat_full))
  sub_mat <- count_mat_full[, sample_ids, drop = FALSE]
  complete <- colSums(is.na(sub_mat)) == 0
  n_dropped_na <- sum(!complete)
  sub_mat <- sub_mat[, complete, drop = FALSE]
  sample_ids <- colnames(sub_mat)
  cat("Samples after NA filtering:", length(sample_ids), "(", n_dropped_na, "dropped )\n")
  storage.mode(sub_mat) <- "integer"

  col_data <- se_meta[match(sample_ids, se_meta$sample_id), c("sample_id", "L_status", colData_extra_cols)]
  rownames(col_data) <- col_data$sample_id
  col_data$L_status <- factor(col_data$L_status, levels = c("non_L", "L"))
  stopifnot(all(colnames(sub_mat) == rownames(col_data)))

  cat("Design:", deparse(design_formula), "\n")
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

# ---- Full N (all SE samples, no covariates) ----
res_full <- run_deseq(se_meta$sample_id, ~ L_status, character(0), "SE_full_N")

# ---- Restricted N (samples with complete clinical data, sex+age covariates) ----
se_meta$age <- suppressWarnings(as.numeric(se_meta$age))
se_meta$age_centered <- se_meta$age - mean(se_meta$age, na.rm = TRUE)
restricted_ids <- se_meta$sample_id[se_meta$condition != "NA" & !is.na(se_meta$condition)]
res_restricted <- run_deseq(restricted_ids, ~ sex + age_centered + L_status,
                             c("sex", "age_centered"), "SE_restricted_N")

# ---- Compare against PE's original results, for the full PE vs SE replication check ----
pe_full <- read.delim(file.path(OUTDIR, "full_N_results.tsv"), stringsAsFactors = FALSE)
merged <- merge(pe_full[, c("gene_name", "log2FoldChange", "padj")],
                 res_full[, c("gene_name", "log2FoldChange", "padj")],
                 by = "gene_name", suffixes = c("_PE", "_SE"))
merged$same_direction <- sign(merged$log2FoldChange_PE) == sign(merged$log2FoldChange_SE)
merged <- merged[order(merged$padj_PE), ]

cat("\n===== PE vs SE replication (full_N) =====\n")
print(merged)
write.table(merged, file.path(OUTDIR, "PE_vs_SE_full_N_replication.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nDone. Check", OUTDIR, "for results.\n")