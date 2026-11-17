#!/usr/bin/env Rscript
# =============================================================================
# Genome-wide nuclear co-expression analysis: correlate ALL nuclear genes
# against mtDNA gene expression, matching the paper's actual design (they
# ran co-expression FIRST, genome-wide, then tested the significantly
# co-expressed genes for genetic association -- the reverse order of what
# we did with the MitoCarta-restricted eQTL).
#
# Run separately in PE and SE (discovery/replication), matching the
# framework used throughout this project.
#
# Usage:
#   Rscript nuclear_coexpression.R
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
})

# ---- EDIT THESE PATHS ----
NUCLEAR_EXPR   <- "nuclear_expression_matrix.tsv"
MTDNA_EXPR     <- "mtdna_expression_matrix.tsv"
PE_METADATA    <- "pe_metadata_unified.tsv"
SE_METADATA    <- "se_metadata_unified.tsv"
OUTDIR         <- "analysis1_eqtl"
FDR_THRESHOLD  <- 0.05
# ---------------------------

dir.create(OUTDIR, showWarnings = FALSE)

normalize_matrix <- function(path, extra_cols = c("gene_name", "gene_type")) {
  raw <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  count_cols <- setdiff(colnames(raw), c("gene_id", extra_cols))
  mat <- as.matrix(raw[, count_cols])
  rownames(mat) <- if ("gene_name" %in% colnames(raw)) raw$gene_name else raw$gene_id
  mode(mat) <- "numeric"
  complete <- colSums(is.na(mat)) == 0
  mat <- mat[, complete]
  mat <- mat[rowSums(mat) > 0, ]
  storage.mode(mat) <- "integer"
  dds <- DESeqDataSetFromMatrix(mat, colData = data.frame(row.names = colnames(mat)), design = ~1)
  dds <- estimateSizeFactors(dds)
  log2(counts(dds, normalized = TRUE) + 1)
}

cat("Normalizing nuclear expression matrix (genome-wide -- this will take a while)...\n")
nuclear_norm <- normalize_matrix(NUCLEAR_EXPR)
cat("Nuclear matrix:", nrow(nuclear_norm), "genes x", ncol(nuclear_norm), "samples\n")

cat("Normalizing mtDNA expression matrix...\n")
mtdna_norm <- normalize_matrix(MTDNA_EXPR)
cat("mtDNA matrix:", nrow(mtdna_norm), "genes x", ncol(mtdna_norm), "samples\n")

pe_meta <- read.delim(PE_METADATA, stringsAsFactors = FALSE)
se_meta <- read.delim(SE_METADATA, stringsAsFactors = FALSE)

run_coexpression_genomewide <- function(sample_ids, cohort_label) {
  samples <- Reduce(intersect, list(sample_ids, colnames(nuclear_norm), colnames(mtdna_norm)))
  cat("\n=====", cohort_label, "=====\n")
  cat("Samples:", length(samples), "\n")

  nuc_sub <- t(nuclear_norm[, samples])
  mt_sub  <- t(mtdna_norm[, samples])
  n <- length(samples)

  cat("Computing", ncol(nuc_sub), "x", ncol(mt_sub), "=",
      format(ncol(nuc_sub) * ncol(mt_sub), big.mark = ","), "correlations...\n")
  r_mat <- cor(nuc_sub, mt_sub, method = "pearson")

  t_mat <- r_mat * sqrt((n - 2) / (1 - r_mat^2))
  p_mat <- 2 * pt(-abs(t_mat), df = n - 2)

  res <- data.frame(
    nuclear_gene = rep(rownames(r_mat), times = ncol(r_mat)),
    mtdna_gene   = rep(colnames(r_mat), each = nrow(r_mat)),
    r = as.vector(r_mat),
    pvalue = as.vector(p_mat),
    n = n,
    cohort = cohort_label
  )
  res <- res[!is.na(res$pvalue), ]
  res$padj <- p.adjust(res$pvalue, method = "BH")
  res
}

pe_res <- run_coexpression_genomewide(pe_meta$sample_id, "PE")
se_res <- run_coexpression_genomewide(se_meta$sample_id, "SE")

pe_sig <- pe_res[pe_res$padj < FDR_THRESHOLD, ]
cat("\nPE: significant nuclear-mtDNA co-expression pairs (padj <", FDR_THRESHOLD, "):", nrow(pe_sig), "\n")
cat("PE: unique nuclear genes involved:", length(unique(pe_sig$nuclear_gene)), "\n")

write.table(pe_res, file.path(OUTDIR, "PE_nuclear_coexpression_full.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(se_res, file.path(OUTDIR, "SE_nuclear_coexpression_full.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

merged <- merge(pe_sig, se_res, by = c("nuclear_gene", "mtdna_gene"), suffixes = c("_PE", "_SE"))
merged$same_direction <- sign(merged$r_PE) == sign(merged$r_SE)
merged$sig_in_SE <- merged$padj_SE < FDR_THRESHOLD

cat("\nReplication check:\n")
cat("  PE significant pairs:", nrow(pe_sig), "\n")
cat("  Present in SE:", nrow(merged), "\n")
cat("  Same direction in SE:", sum(merged$same_direction, na.rm = TRUE), "\n")
cat("  Also significant in SE:", sum(merged$sig_in_SE, na.rm = TRUE), "\n")

replicating_genes <- unique(merged$nuclear_gene[merged$sig_in_SE])
cat("  Unique replicating nuclear genes:", length(replicating_genes), "\n")

write.table(merged, file.path(OUTDIR, "nuclear_coexpression_replication_check.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(replicating_genes, file.path(OUTDIR, "replicating_coexpressed_nuclear_genes.txt"))

cat("\nDone. Replicating gene list written to:",
    file.path(OUTDIR, "replicating_coexpressed_nuclear_genes.txt"), "\n")
