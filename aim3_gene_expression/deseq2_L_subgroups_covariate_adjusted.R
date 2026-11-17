#!/usr/bin/env Rscript
# =============================================================================
# Covariate-adjusted (sex + age) version of the L-subhaplogroup DE analysis,
# to check whether L3's weak signal (1/37 genes unadjusted) holds up once
# the sex imbalance flagged in the diagnostic (L3=71% female, vs non_L=37%)
# is controlled for. Restricted to samples with complete condition/sex/age,
# matching the "restricted_N" pattern used throughout this project.
#
# Usage:
#   Rscript deseq2_L_subgroups_covariate_adjusted.R
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

pe_meta$age <- suppressWarnings(as.numeric(pe_meta$age))
pe_meta$age_centered <- pe_meta$age - mean(pe_meta$age, na.rm = TRUE)
se_meta$age <- suppressWarnings(as.numeric(se_meta$age))
se_meta$age_centered <- se_meta$age - mean(se_meta$age, na.rm = TRUE)

gene_info <- expr[, c("gene_id", "gene_name", "gene_type")]
count_cols <- setdiff(colnames(expr), c("gene_id", "gene_name", "gene_type"))
count_mat_full <- as.matrix(expr[, count_cols])
rownames(count_mat_full) <- expr$gene_id
mode(count_mat_full) <- "numeric"

run_subgroup_deseq_adjusted <- function(meta, subgroups_to_test, label) {
  cat("\n=====", label, "=====\n")

  keep_levels <- c("non_L", subgroups_to_test)
  meta_use <- meta[meta$L_subgroup %in% keep_levels & meta$condition != "NA" & !is.na(meta$condition), ]
  sample_ids <- intersect(meta_use$sample_id, colnames(count_mat_full))

  sub_mat <- count_mat_full[, sample_ids, drop = FALSE]
  complete <- colSums(is.na(sub_mat)) == 0
  n_dropped <- sum(!complete)
  sub_mat <- sub_mat[, complete, drop = FALSE]
  sample_ids <- colnames(sub_mat)
  storage.mode(sub_mat) <- "integer"

  col_data <- meta_use[match(sample_ids, meta_use$sample_id),
                        c("sample_id", "L_subgroup", "sex", "age_centered")]
  rownames(col_data) <- col_data$sample_id
  col_data$L_subgroup <- factor(col_data$L_subgroup, levels = keep_levels)
  stopifnot(all(colnames(sub_mat) == rownames(col_data)))

  cat("Samples (restricted to complete clinical data, NA-filtered, ", n_dropped, " dropped):",
      length(sample_ids), "\n", sep = "")
  print(table(col_data$L_subgroup))
  cat("Sex balance by subgroup:\n")
  print(table(col_data$L_subgroup, col_data$sex))

  if (any(table(col_data$L_subgroup) < 5)) {
    cat("WARNING: at least one group has <5 samples in this restricted subset -- treat with caution\n")
  }

  dds <- DESeqDataSetFromMatrix(countData = sub_mat, colData = col_data,
                                 design = ~ sex + age_centered + L_subgroup)
  dds <- DESeq(dds)

  all_contrasts <- list()
  for (sg in subgroups_to_test) {
    if (!(sg %in% col_data$L_subgroup)) next
    n_sg <- sum(col_data$L_subgroup == sg)
    if (n_sg < 5) {
      cat(sg, "vs non_L: n=", n_sg, " -- SKIPPED (too few samples)\n")
      next
    }
    res <- results(dds, contrast = c("L_subgroup", sg, "non_L"))
    res_df <- as.data.frame(res)
    res_df$gene_id <- rownames(res_df)
    res_df <- merge(gene_info, res_df, by = "gene_id")
    res_df <- res_df[order(res_df$pvalue), ]
    n_sig <- sum(res_df$padj < 0.05, na.rm = TRUE)
    cat(sg, "vs non_L (adjusted): n=", n_sg, " | padj<0.05:", n_sig, "/", nrow(res_df), "\n")
    write.table(res_df, file.path(OUTDIR, paste0(label, "_", sg, "_vs_nonL_ADJUSTED_results.tsv")),
                sep = "\t", quote = FALSE, row.names = FALSE)
    all_contrasts[[sg]] <- res_df
  }
  all_contrasts
}

# ---- PE: all subgroups, sex+age adjusted ----
pe_adjusted <- run_subgroup_deseq_adjusted(pe_meta, c("L0", "L1", "L2", "L3", "L4"), "PE")

# ---- SE: L0 only, sex+age adjusted ----
se_adjusted <- run_subgroup_deseq_adjusted(se_meta, c("L0"), "SE")

# ---- Compare L3's unadjusted vs adjusted results directly ----
cat("\n===== L3 vs non-L: unadjusted vs sex/age-adjusted comparison =====\n")
l3_unadjusted <- read.delim(file.path(OUTDIR, "PE_L3_vs_nonL_results.tsv"), stringsAsFactors = FALSE)
l3_adjusted <- pe_adjusted[["L3"]]
if (!is.null(l3_adjusted)) {
  merged <- merge(l3_unadjusted[, c("gene_name", "log2FoldChange", "padj")],
                   l3_adjusted[, c("gene_name", "log2FoldChange", "padj")],
                   by = "gene_name", suffixes = c("_unadjusted", "_adjusted"))
  merged <- merged[order(merged$padj_unadjusted), ]
  print(merged)
  write.table(merged, file.path(OUTDIR, "L3_unadjusted_vs_adjusted_comparison.tsv"),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat("\nSignificant (padj<0.05) unadjusted:", sum(merged$padj_unadjusted < 0.05, na.rm=TRUE),
      " | adjusted:", sum(merged$padj_adjusted < 0.05, na.rm=TRUE), "\n")
}

cat("\nDone. Check", OUTDIR, "for results.\n")