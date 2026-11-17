#!/usr/bin/env Rscript
# =============================================================================
# Phase F: nuclear co-expression + L/non-L DE + GO enrichment -- CLEAN REBUILD
# matching the paper's actual two-stage design:
#
#   "To assess co-expression of nDNA-encoded genes with the mtDNA-encoded
#    ones, Pearson correlation was employed (p < 5.23e-8, Bonferroni-
#    corrected)... co-expression was sought between the 15 mtDNA-encoded
#    mRNA and rRNA genes and 63,662 nDNA-encoded genes. Then, Matrix eQTL
#    was employed to identify differential expression between L and non-L
#    genetic backgrounds among the significantly co-expressed genes...
#    Finally, only nDNA-encoded genes that were significantly co-expressed,
#    AND differentially expressed [L vs non-L]... were subjected to GO
#    analysis."
#
# STAGE 1: Pearson co-expression, mtDNA mRNA+rRNA genes (15, matching the
#          paper's count exactly) vs all presence-filtered nuclear genes.
#          Bonferroni threshold computed from OUR actual gene counts.
#          Required significant in BOTH PE and SE (this project's
#          discovery/replication framework, consistent with every other
#          analysis here -- paper's original design was single-cohort).
# STAGE 2: Among STAGE 1 genes, Matrix eQTL: L vs non-L status (binary,
#          treated as genotype) vs nuclear gene expression. Bonferroni
#          corrected. Covariates: sex (paper also used mtDNA copy number +
#          lab-of-origin, not available here -- documented gap).
# STAGE 3: GO enrichment (gprofiler2, local -- see conversation notes on why
#          not the HPC/RCurl route) on genes passing BOTH stages, in BOTH
#          cohorts.
#
# Usage:
#   Rscript coexpression_clean.R
# =============================================================================

if (!requireNamespace("MatrixEQTL", quietly = TRUE)) install.packages("MatrixEQTL")
suppressPackageStartupMessages({
  library(DESeq2)
  library(MatrixEQTL)
})

# ============================== CONFIG ======================================
MTDNA_EXPR    <- "mtdna_expression_matrix.tsv"
NUCLEAR_EXPR  <- "nuclear_expression_matrix.tsv"
PE_METADATA   <- "pe_metadata_unified.tsv"
SE_METADATA   <- "se_metadata_unified.tsv"
OUTDIR        <- "phase_f_coexpression"
MIN_PRESENCE  <- 0.90
# =============================================================================

dir.create(OUTDIR, showWarnings = FALSE)
pe_meta <- read.delim(PE_METADATA, stringsAsFactors = FALSE)
se_meta <- read.delim(SE_METADATA, stringsAsFactors = FALSE)
pe_meta$L_status <- ifelse(startsWith(pe_meta$haplogroup, "L"), 1, 0)
se_meta$L_status <- ifelse(startsWith(se_meta$haplogroup, "L"), 1, 0)

# -----------------------------------------------------------------------------
# Load + normalize both expression matrices, with gene-presence filtering
# -----------------------------------------------------------------------------
cat("\n############ Loading and normalizing expression data ############\n")

load_and_normalize <- function(path, keep_types = NULL) {
  raw <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  extra_cols <- intersect(c("gene_id", "gene_name", "gene_type"), colnames(raw))
  if (!is.null(keep_types) && "gene_type" %in% colnames(raw)) {
    raw <- raw[raw$gene_type %in% keep_types, ]
  }
  count_cols <- setdiff(colnames(raw), extra_cols)
  mat <- as.matrix(raw[, count_cols])
  rownames(mat) <- raw$gene_id
  mode(mat) <- "numeric"
  complete <- colSums(is.na(mat)) == 0
  mat <- mat[, complete, drop = FALSE]
  mat <- mat[rowSums(mat) > 0, , drop = FALSE]
  storage.mode(mat) <- "integer"
  dds <- DESeqDataSetFromMatrix(mat, colData = data.frame(row.names = colnames(mat)), design = ~1)
  dds <- estimateSizeFactors(dds)
  list(norm = log2(counts(dds, normalized = TRUE) + 1), gene_info = raw[, extra_cols])
}

# mtDNA side: mRNA + rRNA ONLY (matches paper's "15 mtDNA-encoded mRNA and
# rRNA genes" exactly -- 13 protein-coding + 2 rRNA in our data)
mtdna_loaded <- load_and_normalize(MTDNA_EXPR, keep_types = c("protein_coding", "Mt_rRNA"))
mtdna_norm <- mtdna_loaded$norm
cat("mtDNA mRNA+rRNA genes:", nrow(mtdna_norm), "(paper: 15)\n")

nuclear_loaded <- load_and_normalize(NUCLEAR_EXPR)
nuclear_norm <- nuclear_loaded$norm
cat("Nuclear genes (all):", nrow(nuclear_norm), "\n")

# ---- Exclude NUMT pseudogenes (nuclear mtDNA-homologous loci) -- reads
# from real mtDNA transcripts can misalign here, creating a spurious
# "co-expressed with mtDNA + differs by L status" signal that's actually
# just spillover, not real nuclear regulatory biology. See conversation
# notes: this is exactly what happened in the first run of this pipeline.
NUMT_BLOCKLIST <- "numt_pseudogene_blocklist.txt"
if (file.exists(NUMT_BLOCKLIST)) {
  numt_ids <- readLines(NUMT_BLOCKLIST)
  numt_ids_in_matrix <- intersect(numt_ids, rownames(nuclear_norm))
  cat("Excluding", length(numt_ids_in_matrix), "NUMT pseudogenes from the nuclear candidate list\n")
  nuclear_norm <- nuclear_norm[!(rownames(nuclear_norm) %in% numt_ids), , drop = FALSE]
  cat("Nuclear genes after NUMT exclusion:", nrow(nuclear_norm), "\n")
} else {
  cat("WARNING:", NUMT_BLOCKLIST, "not found -- NUMT pseudogenes NOT excluded.\n")
  cat("Run 53b_build_numt_blocklist.sh first to generate it.\n")
}

presence_filter <- function(norm_mat, raw_path, sample_ids, cohort_label) {
  raw <- read.delim(raw_path, stringsAsFactors = FALSE, check.names = FALSE)
  extra_cols <- intersect(c("gene_id", "gene_name", "gene_type"), colnames(raw))
  count_cols <- setdiff(colnames(raw), extra_cols)
  samples <- intersect(sample_ids, count_cols)
  sub <- as.matrix(raw[, samples, drop = FALSE]); mode(sub) <- "numeric"
  rownames(sub) <- raw$gene_id
  presence_frac <- rowMeans(!is.na(sub) & sub > 0)
  passing <- rownames(sub)[presence_frac >= MIN_PRESENCE]
  passing <- intersect(passing, rownames(norm_mat))
  cat(cohort_label, ": ", length(passing), "/", nrow(norm_mat), " genes pass presence filter\n", sep = "")
  passing
}

pe_nuclear_genes <- presence_filter(nuclear_norm, NUCLEAR_EXPR, pe_meta$sample_id, "PE nuclear")
se_nuclear_genes <- presence_filter(nuclear_norm, NUCLEAR_EXPR, se_meta$sample_id, "SE nuclear")
pe_mtdna_genes   <- presence_filter(mtdna_norm, MTDNA_EXPR, pe_meta$sample_id, "PE mtDNA mRNA/rRNA")
se_mtdna_genes   <- presence_filter(mtdna_norm, MTDNA_EXPR, se_meta$sample_id, "SE mtDNA mRNA/rRNA")

# -----------------------------------------------------------------------------
# STAGE 1: Co-expression (Pearson, Bonferroni), required significant in BOTH
# -----------------------------------------------------------------------------
cat("\n############ STAGE 1: Nuclear-mtDNA co-expression ############\n")

run_coexpression <- function(sample_ids, cohort_label, nuclear_genes, mtdna_genes) {
  samples <- Reduce(intersect, list(sample_ids, colnames(nuclear_norm), colnames(mtdna_norm)))
  nuc_sub <- t(nuclear_norm[nuclear_genes, samples, drop = FALSE])
  mt_sub  <- t(mtdna_norm[mtdna_genes, samples, drop = FALSE])
  n <- length(samples)
  n_tests <- length(nuclear_genes) * length(mtdna_genes)
  bonferroni_alpha <- 0.05 / n_tests

  cat(cohort_label, "- samples:", n, " | nuclear genes:", length(nuclear_genes),
      " | mtDNA genes:", length(mtdna_genes), " | tests:", n_tests,
      " | Bonferroni alpha:", format(bonferroni_alpha, scientific = TRUE), "\n")

  r_mat <- cor(nuc_sub, mt_sub, method = "pearson")
  t_mat <- r_mat * sqrt((n - 2) / (1 - r_mat^2))
  p_mat <- 2 * pt(-abs(t_mat), df = n - 2)

  res <- data.frame(
    nuclear_gene = rep(rownames(r_mat), times = ncol(r_mat)),
    mtdna_gene   = rep(colnames(r_mat), each = nrow(r_mat)),
    r = as.vector(r_mat), pvalue = as.vector(p_mat)
  )
  res <- res[!is.na(res$pvalue), ]
  res$bonferroni_sig <- res$pvalue < bonferroni_alpha
  res
}

pe_coexp <- run_coexpression(pe_meta$sample_id, "PE", pe_nuclear_genes, pe_mtdna_genes)
se_coexp <- run_coexpression(se_meta$sample_id, "SE", se_nuclear_genes, se_mtdna_genes)

pe_coexp_sig_genes <- unique(pe_coexp$nuclear_gene[pe_coexp$bonferroni_sig])
se_coexp_sig_genes <- unique(se_coexp$nuclear_gene[se_coexp$bonferroni_sig])
coexp_both <- intersect(pe_coexp_sig_genes, se_coexp_sig_genes)
cat("\nCo-expressed in PE:", length(pe_coexp_sig_genes),
    " | in SE:", length(se_coexp_sig_genes),
    " | in BOTH:", length(coexp_both), "\n")

write.table(pe_coexp[pe_coexp$bonferroni_sig, ], file.path(OUTDIR, "PE_coexpression_significant.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(se_coexp[se_coexp$bonferroni_sig, ], file.path(OUTDIR, "SE_coexpression_significant.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(coexp_both, file.path(OUTDIR, "coexpressed_in_both_cohorts.txt"))

# -----------------------------------------------------------------------------
# STAGE 2: L vs non-L DE among co-expressed genes (Matrix eQTL, matching
# paper's method of using Matrix eQTL for this step specifically)
# -----------------------------------------------------------------------------
cat("\n############ STAGE 2: L vs non-L DE among co-expressed genes ############\n")

if (length(coexp_both) == 0) {
  cat("No genes co-expressed in both cohorts -- stopping here.\n")
  quit(save = "no")
}

run_L_vs_nonL <- function(sample_ids, meta, cohort_label, candidate_genes) {
  samples <- Reduce(intersect, list(sample_ids, colnames(nuclear_norm), meta$sample_id))
  meta_sub <- meta[match(samples, meta$sample_id), ]

  expr_sub <- nuclear_norm[candidate_genes, samples, drop = FALSE]

  sex_numeric <- as.numeric(factor(meta_sub$sex)) - 1
  cvrt <- SlicedData$new()
  if (all(!is.na(sex_numeric)) && length(unique(sex_numeric)) > 1) {
    m <- matrix(sex_numeric, nrow = 1, dimnames = list("sex", samples))
    cvrt$CreateFromMatrix(m)
  }

  l_status <- matrix(meta_sub$L_status, nrow = 1, dimnames = list("L_status", samples))
  snps <- SlicedData$new(); snps$CreateFromMatrix(l_status)
  gene <- SlicedData$new(); gene$CreateFromMatrix(expr_sub)

  n_tests <- length(candidate_genes)
  bonferroni_alpha <- 0.05 / n_tests
  cat(cohort_label, "- samples:", length(samples), " | genes tested:", n_tests,
      " | Bonferroni alpha:", format(bonferroni_alpha, scientific = TRUE), "\n")

  out_file <- file.path(OUTDIR, paste0(cohort_label, "_L_vs_nonL_DE_results.tsv"))
  me <- Matrix_eQTL_engine(snps = snps, gene = gene, cvrt = cvrt, output_file_name = out_file,
                            pvOutputThreshold = 1, useModel = modelLINEAR,
                            errorCovariance = numeric(), verbose = FALSE, pvalue.hist = FALSE)
  res <- me$all$eqtls
  res$bonferroni_sig <- res$pvalue < bonferroni_alpha
  write.table(res, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("Significant:", sum(res$bonferroni_sig), "/", nrow(res), "\n")
  res
}

pe_de <- run_L_vs_nonL(pe_meta$sample_id, pe_meta, "PE", coexp_both)
se_de <- run_L_vs_nonL(se_meta$sample_id, se_meta, "SE", coexp_both)

pe_de_sig <- unique(pe_de$gene[pe_de$bonferroni_sig])
se_de_sig <- unique(se_de$gene[se_de$bonferroni_sig])
final_genes <- intersect(pe_de_sig, se_de_sig)
cat("\nDE-significant in PE:", length(pe_de_sig), " | in SE:", length(se_de_sig),
    " | in BOTH (final gene set):", length(final_genes), "\n")

writeLines(final_genes, file.path(OUTDIR, "final_genes_coexpressed_AND_DE.txt"))

# -----------------------------------------------------------------------------
# STAGE 3: GO enrichment on the final gene set
# -----------------------------------------------------------------------------
cat("\n############ STAGE 3: GO enrichment ############\n")

if (length(final_genes) == 0) {
  cat("No genes passed both stages -- nothing to run GO enrichment on.\n")
} else {
  clean_ids <- unique(sub("\\.[0-9]+$", "", final_genes))
  cat("Genes passing both stages:", length(clean_ids), "\n")
  cat("STOPPING HERE if running on the HPC -- GO enrichment (gprofiler2) needs\n")
  cat("RCurl, which this cluster can't build (missing curl-config). Transfer\n")
  cat(file.path(OUTDIR, "final_genes_coexpressed_AND_DE.txt"),
      "to your Mac and run 58_go_enrichment_final_local.R there instead.\n")
}

cat("\n\nDone. All output in", OUTDIR, "\n")
