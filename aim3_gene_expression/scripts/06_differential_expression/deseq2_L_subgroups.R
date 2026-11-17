#!/usr/bin/env Rscript
# =============================================================================
# L-subhaplogroup-resolved DE (L0, L1, L2, L3, L4 vs non-L), extending the
# binary L vs non-L analysis with finer resolution.
#
# PE: all 5 subgroups tested (each has reasonable N: L0=140, L2=109, L3=75,
#     L1=28, L4=6 -- L4 is thin, treated as exploratory).
# SE: ONLY L0 vs non-L is tested (SE's other subgroups have too few samples
#     -- L2=8, L1=2, L3=2, L4=1 -- for any meaningful test or replication
#     claim; we don't pretend otherwise).
# L5/L6 (1 sample each in PE, none in SE) excluded from formal testing --
# cannot estimate within-group variance from n=1.
#
# Usage:
#   Rscript deseq2_L_subgroups.R
# =============================================================================

suppressPackageStartupMessages(library(DESeq2))

# ---- EDIT THESE PATHS ----
EXPR_MATRIX   <- "mtdna_expression_matrix.tsv"
PE_METADATA   <- "pe_metadata_unified.tsv"
SE_METADATA   <- "se_metadata_unified.tsv"
OUTDIR        <- "analysis2_L_subgroups"
# ---------------------------

dir.create(OUTDIR, showWarnings = FALSE)

expr <- read.delim(EXPR_MATRIX, stringsAsFactors = FALSE, check.names = FALSE)
pe_meta <- read.delim(PE_METADATA, stringsAsFactors = FALSE)
se_meta <- read.delim(SE_METADATA, stringsAsFactors = FALSE)

get_major_haplogroup <- function(hap) {
  ifelse(startsWith(hap, "L"), substr(hap, 1, 2), "non_L")
}
pe_meta$L_subgroup <- get_major_haplogroup(pe_meta$haplogroup)
se_meta$L_subgroup <- get_major_haplogroup(se_meta$haplogroup)

gene_info <- expr[, c("gene_id", "gene_name", "gene_type")]
count_cols <- setdiff(colnames(expr), c("gene_id", "gene_name", "gene_type"))
count_mat_full <- as.matrix(expr[, count_cols])
rownames(count_mat_full) <- expr$gene_id
mode(count_mat_full) <- "numeric"

run_subgroup_deseq <- function(meta, sample_ids, subgroups_to_test, label) {
  cat("\n=====", label, "=====\n")

  # Exclude L5/L6 (or anything not in subgroups_to_test + non_L) -- can't
  # test groups with n=1
  keep_levels <- c("non_L", subgroups_to_test)
  meta_use <- meta[meta$L_subgroup %in% keep_levels, ]
  sample_ids <- intersect(sample_ids, meta_use$sample_id)

  sample_ids <- intersect(sample_ids, colnames(count_mat_full))
  sub_mat <- count_mat_full[, sample_ids, drop = FALSE]
  complete <- colSums(is.na(sub_mat)) == 0
  n_dropped <- sum(!complete)
  sub_mat <- sub_mat[, complete, drop = FALSE]
  sample_ids <- colnames(sub_mat)
  storage.mode(sub_mat) <- "integer"

  col_data <- meta_use[match(sample_ids, meta_use$sample_id), c("sample_id", "L_subgroup")]
  rownames(col_data) <- col_data$sample_id
  col_data$L_subgroup <- factor(col_data$L_subgroup, levels = keep_levels)  # non_L = reference
  stopifnot(all(colnames(sub_mat) == rownames(col_data)))

  cat("Samples (NA-filtered, ", n_dropped, " dropped):", length(sample_ids), "\n")
  print(table(col_data$L_subgroup))

  dds <- DESeqDataSetFromMatrix(countData = sub_mat, colData = col_data, design = ~ L_subgroup)
  dds <- DESeq(dds)

  # Omnibus LRT: does ANY subgroup differ from non-L / from each other?
  dds_lrt <- DESeq(dds, test = "LRT", reduced = ~1)
  lrt_res <- as.data.frame(results(dds_lrt))
  lrt_res$gene_id <- rownames(lrt_res)
  lrt_res <- merge(gene_info, lrt_res, by = "gene_id")
  lrt_res <- lrt_res[order(lrt_res$pvalue), ]
  write.table(lrt_res, file.path(OUTDIR, paste0(label, "_omnibus_LRT.tsv")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("Omnibus LRT (any subgroup differs): ", sum(lrt_res$padj < 0.05, na.rm = TRUE),
      "/", nrow(lrt_res), "genes significant\n")

  # Per-subgroup contrasts vs non-L
  all_contrasts <- list()
  for (sg in subgroups_to_test) {
    if (!(sg %in% col_data$L_subgroup)) next
    res <- results(dds, contrast = c("L_subgroup", sg, "non_L"))
    res_df <- as.data.frame(res)
    res_df$gene_id <- rownames(res_df)
    res_df <- merge(gene_info, res_df, by = "gene_id")
    res_df$subgroup <- sg
    res_df <- res_df[order(res_df$pvalue), ]
    n_sig <- sum(res_df$padj < 0.05, na.rm = TRUE)
    cat(sg, "vs non_L: n=", sum(col_data$L_subgroup == sg), " | padj<0.05:", n_sig, "/", nrow(res_df), "\n")
    write.table(res_df, file.path(OUTDIR, paste0(label, "_", sg, "_vs_nonL_results.tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    all_contrasts[[sg]] <- res_df
  }

  list(dds = dds, contrasts = all_contrasts, lrt = lrt_res)
}

# ---- PE: all subgroups with reasonable N ----
pe_results <- run_subgroup_deseq(pe_meta, pe_meta$sample_id,
                                  subgroups_to_test = c("L0", "L1", "L2", "L3", "L4"),
                                  label = "PE")

# ---- SE: L0 only (only subgroup with viable N) ----
se_results <- run_subgroup_deseq(se_meta, se_meta$sample_id,
                                  subgroups_to_test = c("L0"),
                                  label = "SE")

# ---- Replication check for L0 (the only subgroup testable in both) ----
cat("\n===== L0 vs non-L: PE/SE replication =====\n")
pe_L0 <- pe_results$contrasts[["L0"]]
se_L0 <- se_results$contrasts[["L0"]]
merged <- merge(pe_L0[, c("gene_name", "log2FoldChange", "padj")],
                 se_L0[, c("gene_name", "log2FoldChange", "padj")],
                 by = "gene_name", suffixes = c("_PE", "_SE"))
merged$same_direction <- sign(merged$log2FoldChange_PE) == sign(merged$log2FoldChange_SE)
merged <- merged[order(merged$padj_PE), ]
print(merged)
write.table(merged, file.path(OUTDIR, "L0_vs_nonL_PE_SE_replication.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\nDone. Check", OUTDIR, "for results.\n")
cat("\nNOTE: L2/L3/L1/L4 results are PE-only (exploratory/discovery) -- SE\n")
cat("sample sizes for these subgroups (L2=8, L1=2, L3=2, L4=1) are too small\n")
cat("for a meaningful replication test, so none is claimed.\n")
