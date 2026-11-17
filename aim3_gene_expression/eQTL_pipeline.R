#!/usr/bin/env Rscript
# =============================================================================
# Analysis 1 (Matrix eQTL) -- rewritten to match the original paper's stated
# methodology as closely as our available data allows.
#
# From the paper's Methods:
#   "Only genes presenting with reads in more than 90% of the samples were
#    used" -- gene-presence filter, applied per cohort.
#   "eQTL mapping was calculated according to the allele assignment... A
#    Bonferroni correction was employed... we focused on SNPs shared by at
#    least 10 individuals... major allele frequency will not exceed 95%"
#    -- MAC>=10, major allele freq <=95%, Bonferroni (not FDR), and
#    genotype treated as a categorical GROUP variable (modelANOVA) for the
#    mtDNA-SNP model, matching "divided into groups according to allele
#    assignment" rather than a continuous linear dosage.
#
# DOCUMENTED DEVIATIONS from the paper (data constraints, not oversights):
#   - Covariates: paper used gender + mtDNA copy number + lab-of-origin.
#     We only have gender (sex). No matched DNA-seq means no true mtDNA
#     copy number; genotype PCs are used as an ancestry-structure proxy for
#     the nDNA-SNP model instead of lab-of-origin (NOT applied to the
#     mtDNA-SNP model, where mtDNA-linked ancestry IS the signal of
#     interest, not a confound to remove).
#   - nDNA candidate genes: paper used a specific "mitochondrial RNA-binding
#     activity" gene list + 4 named TFs (c-Jun, JunD, CEBPb, Mef2D), not
#     MitoCarta3.0. We used MitoCarta3.0 (broader) since the paper's exact
#     source list wasn't readily available to us.
#   - nDNA SNPs: paper restricted to known 1000 Genomes variants. We did not
#     cross-reference against 1000 Genomes; our SNPs are whatever GATK
#     called directly from RNA-seq.
#
# Usage:
#   Rscript analysis1_eqtl_clean.R
# =============================================================================

if (!requireNamespace("MatrixEQTL", quietly = TRUE)) install.packages("MatrixEQTL")
suppressPackageStartupMessages({
  library(MatrixEQTL)
  library(DESeq2)
})

# modelANOVA defaults to a max of 3 genotype categories -- our mtDNA calls
# (now correctly capturing multiallelic sites, matching the paper's note on
# up to 3 alleles per position) can exceed that in edge cases. Set generously
# high to avoid this failing again on data specifics we haven't fully audited.
options(MatrixEQTL.ANOVA.categories = 10)

# ============================== CONFIG ======================================
MTDNA_EXPR     <- "mtdna_expression_matrix.tsv"
MTDNA_GENO     <- "mtdna_genotype_matrix.tsv"
NDNA_GENO      <- "ndna_genotype_matrix.tsv"
PE_METADATA    <- "pe_metadata_unified.tsv"
SE_METADATA    <- "se_metadata_unified.tsv"
OUTDIR         <- "analysis1_eqtl_v3"
MIN_PRESENCE   <- 0.90   # paper's gene-presence QC criterion
MIN_MAC        <- 10     # paper: "SNPs shared by at least 10 individuals"
MAX_MAJOR_FREQ <- 0.95   # paper: major allele frequency <= 95% (mtDNA only)
N_GENO_PCS     <- 5
# =============================================================================

dir.create(OUTDIR, showWarnings = FALSE)
pe_meta <- read.delim(PE_METADATA, stringsAsFactors = FALSE)
se_meta <- read.delim(SE_METADATA, stringsAsFactors = FALSE)

# -----------------------------------------------------------------------------
# STEP 1: Gene-presence filter on mtDNA expression (per cohort)
# Paper's benchmark: this left 16/22 tRNA genes in their data -- worth
# comparing our own count against, as a sanity check on filter behavior.
# -----------------------------------------------------------------------------
cat("\n############ STEP 1: Gene-presence filter (mtDNA expression) ############\n")

expr_raw <- read.delim(MTDNA_EXPR, stringsAsFactors = FALSE, check.names = FALSE)
gene_info <- expr_raw[, c("gene_id", "gene_name", "gene_type")]
expr_count_cols <- setdiff(colnames(expr_raw), c("gene_id", "gene_name", "gene_type"))
expr_mat_raw <- as.matrix(expr_raw[, expr_count_cols])
rownames(expr_mat_raw) <- expr_raw$gene_id
mode(expr_mat_raw) <- "numeric"

get_presence_pass <- function(sample_ids, cohort_label) {
  samples <- intersect(sample_ids, colnames(expr_mat_raw))
  sub <- expr_mat_raw[, samples, drop = FALSE]
  presence_frac <- rowMeans(!is.na(sub) & sub > 0)
  passing <- rownames(sub)[presence_frac >= MIN_PRESENCE]

  passing_info <- gene_info[gene_info$gene_id %in% passing, ]
  n_trna_pass <- sum(passing_info$gene_type == "Mt_tRNA")
  cat(cohort_label, ": ", length(passing), "/", nrow(sub), " genes pass (",
      n_trna_pass, " tRNA genes -- paper's benchmark was 16/22)\n", sep = "")
  passing
}

pe_genes_pass <- get_presence_pass(pe_meta$sample_id, "PE")
se_genes_pass <- get_presence_pass(se_meta$sample_id, "SE")

# -----------------------------------------------------------------------------
# STEP 2: DESeq2-normalize mtDNA expression
# -----------------------------------------------------------------------------
cat("\n############ STEP 2: Normalize mtDNA expression ############\n")

complete <- colSums(is.na(expr_mat_raw)) == 0
expr_mat_clean <- expr_mat_raw[, complete]
storage.mode(expr_mat_clean) <- "integer"
dds <- DESeqDataSetFromMatrix(expr_mat_clean, colData = data.frame(row.names = colnames(expr_mat_clean)), design = ~1)
dds <- estimateSizeFactors(dds)
norm_expr <- log2(counts(dds, normalized = TRUE) + 1)
cat("Normalized:", nrow(norm_expr), "genes x", ncol(norm_expr), "samples\n")

# -----------------------------------------------------------------------------
# STEP 3: Genotype QC filters
#   - MAC >= 10 (both genotype types)
#   - Major allele frequency <= 95% (mtDNA only -- paper's multiallelic-site
#     criterion; not described for nDNA)
# -----------------------------------------------------------------------------
cat("\n############ STEP 3: Genotype QC filters ############\n")

load_geno <- function(path) {
  g <- read.delim(path, stringsAsFactors = FALSE, check.names = FALSE, row.names = 1)
  as.matrix(g)
}

coerce_geno <- function(sub) {
  sub_char <- sub
  mode(sub_char) <- "character"
  is_na_flag <- sub_char == "NA" | is.na(sub_char)
  sub_num <- suppressWarnings(matrix(as.numeric(sub_char), nrow = nrow(sub_char), dimnames = dimnames(sub_char)))
  sub_num[is_na_flag] <- NA
  sub_num
}

geno_qc_filter <- function(geno_mat_raw, sample_ids, cohort_label, ploidy_max,
                            apply_major_freq_filter = FALSE) {
  samples <- intersect(sample_ids, colnames(geno_mat_raw))
  sub <- geno_mat_raw[, samples, drop = FALSE]
  sub_num <- coerce_geno(sub)

  cat("\n--", cohort_label, "--\n")
  cat("Samples:", length(samples), " | SNPs (raw):", nrow(sub_num), "\n")

  stats <- t(apply(sub_num, 1, function(x) {
    x <- x[!is.na(x)]
    if (length(x) == 0) return(c(mac = NA, major_freq = NA))
    allele_count <- sum(x)
    total_possible <- length(x) * ploidy_max
    mac <- min(allele_count, total_possible - allele_count)
    major_freq <- max(allele_count, total_possible - allele_count) / total_possible
    c(mac = mac, major_freq = major_freq)
  }))

  keep <- !is.na(stats[, "mac"]) & stats[, "mac"] >= MIN_MAC
  cat("Passing MAC >=", MIN_MAC, ":", sum(keep), "/", nrow(sub_num), "\n")

  if (apply_major_freq_filter) {
    keep_freq <- !is.na(stats[, "major_freq"]) & stats[, "major_freq"] <= MAX_MAJOR_FREQ
    cat("Also passing major allele freq <=", MAX_MAJOR_FREQ, ":", sum(keep & keep_freq), "\n")
    keep <- keep & keep_freq
  }

  cat("Final passing SNPs:", sum(keep), "\n")
  rownames(sub_num)[keep]
}

mtdna_geno_raw <- load_geno(MTDNA_GENO)
ndna_geno_raw  <- load_geno(NDNA_GENO)

pe_mtdna_snp_pass <- geno_qc_filter(mtdna_geno_raw, pe_meta$sample_id, "PE mtDNA SNPs", ploidy_max = 1, apply_major_freq_filter = TRUE)
se_mtdna_snp_pass <- geno_qc_filter(mtdna_geno_raw, se_meta$sample_id, "SE mtDNA SNPs", ploidy_max = 1, apply_major_freq_filter = TRUE)
pe_ndna_snp_pass  <- geno_qc_filter(ndna_geno_raw,  pe_meta$sample_id, "PE nDNA SNPs",  ploidy_max = 2, apply_major_freq_filter = FALSE)
se_ndna_snp_pass  <- geno_qc_filter(ndna_geno_raw,  se_meta$sample_id, "SE nDNA SNPs",  ploidy_max = 2, apply_major_freq_filter = FALSE)

# -----------------------------------------------------------------------------
# STEP 4: Genotype PCs (ancestry proxy, nDNA-SNP model only -- our
# substitute for the paper's "lab of origin" covariate, which we don't have)
# -----------------------------------------------------------------------------
cat("\n############ STEP 4: Genotype PCs (nDNA, ancestry proxy) ############\n")

compute_geno_pcs <- function(geno_mat_raw, snp_pass_list, sample_ids, cohort_label) {
  samples <- intersect(sample_ids, colnames(geno_mat_raw))
  sub <- geno_mat_raw[snp_pass_list, samples, drop = FALSE]
  sub <- coerce_geno(sub)
  complete <- rowSums(is.na(sub)) == 0
  sub <- sub[complete, , drop = FALSE]
  snp_var <- apply(sub, 1, var)
  sub <- sub[snp_var > 0, , drop = FALSE]

  max_snps <- 50000
  if (nrow(sub) > max_snps) {
    set.seed(42)
    sub <- sub[sample(nrow(sub), max_snps), , drop = FALSE]
  }
  cat(cohort_label, "- PCA on", nrow(sub), "SNPs x", ncol(sub), "samples\n")

  pca <- prcomp(t(sub), center = TRUE, scale. = TRUE)
  pcs <- pca$x[, 1:min(N_GENO_PCS, ncol(pca$x)), drop = FALSE]
  colnames(pcs) <- paste0("genoPC", 1:ncol(pcs))
  var_exp <- round(100 * pca$sdev[1:ncol(pcs)]^2 / sum(pca$sdev^2), 1)
  cat("  Variance explained:", paste(var_exp, collapse = "%, "), "%\n")
  data.frame(sample_id = rownames(pcs), pcs, check.names = FALSE)
}

pe_geno_pcs <- compute_geno_pcs(ndna_geno_raw, pe_ndna_snp_pass, pe_meta$sample_id, "PE")
se_geno_pcs <- compute_geno_pcs(ndna_geno_raw, se_ndna_snp_pass, se_meta$sample_id, "SE")

# Save these to disk -- needed by downstream robustness/sensitivity checks
# (e.g. 64_pcbd2_robustness_checks.R) so they use the EXACT SAME PCs as
# this headline analysis, not an independently-recomputed (and potentially
# slightly different) version.
write.table(pe_geno_pcs, file.path(OUTDIR, "PE_genotype_pcs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(se_geno_pcs, file.path(OUTDIR, "SE_genotype_pcs.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("Wrote genotype PCs to", OUTDIR, "\n")

# -----------------------------------------------------------------------------
# STEP 5: Matrix eQTL
#   mtDNA-SNP model: modelANOVA (categorical allele-assignment groups,
#     matching "divided into groups according to their allele assignment"
#     -- appropriate for haploid, potentially triallelic mtDNA genotypes)
#   nDNA-SNP model: modelLINEAR (standard diploid dosage)
#   Both: Bonferroni correction (not FDR), matching the paper
# -----------------------------------------------------------------------------
cat("\n############ STEP 5: Matrix eQTL ############\n")

run_eqtl <- function(geno_mat_raw, snp_pass, cohort_samples, meta, label, geno_pcs, gene_pass,
                      use_model = modelLINEAR) {
  cat("\n=====", label, "=====\n")
  samples <- Reduce(intersect, list(cohort_samples, colnames(geno_mat_raw), colnames(norm_expr), meta$sample_id))
  if (!is.null(geno_pcs)) samples <- intersect(samples, geno_pcs$sample_id)
  cat("Samples:", length(samples), "\n")
  if (length(samples) < 10) { cat("Too few samples -- skipping\n"); return(NULL) }

  geno_sub <- geno_mat_raw[rownames(geno_mat_raw) %in% snp_pass, samples, drop = FALSE]
  geno_sub <- coerce_geno(geno_sub)
  cat("SNPs (post QC filter):", nrow(geno_sub), "\n")

  expr_sub <- norm_expr[rownames(norm_expr) %in% gene_pass, samples, drop = FALSE]
  cat("Genes (post presence filter):", nrow(expr_sub), "\n")

  if (nrow(geno_sub) == 0 || nrow(expr_sub) == 0) {
    cat("ERROR: empty genotype or expression matrix after filtering -- skipping this model\n")
    return(NULL)
  }

  n_tests <- nrow(geno_sub) * nrow(expr_sub)
  bonferroni_alpha <- 0.05 / n_tests
  cat("Total tests:", n_tests, " | Bonferroni alpha:", format(bonferroni_alpha, scientific = TRUE), "\n")

  meta_sub <- meta[match(samples, meta$sample_id), ]
  sex_numeric <- as.numeric(factor(meta_sub$sex)) - 1

  cvrt_rows <- list(); cvrt_names <- c()
  if (all(!is.na(sex_numeric)) && length(unique(sex_numeric)) > 1) {
    cvrt_rows[[1]] <- sex_numeric; cvrt_names <- "sex"
  }
  if (!is.null(geno_pcs)) {
    pcs_sub <- geno_pcs[match(samples, geno_pcs$sample_id), ]
    for (pc in grep("^genoPC", colnames(pcs_sub), value = TRUE)) {
      cvrt_rows[[length(cvrt_rows) + 1]] <- pcs_sub[[pc]]
      cvrt_names <- c(cvrt_names, pc)
    }
  }

  snps <- SlicedData$new(); snps$CreateFromMatrix(geno_sub); snps$ResliceCombined(sliceSize = 5000)
  gene <- SlicedData$new(); gene$CreateFromMatrix(expr_sub); gene$ResliceCombined(sliceSize = 500)
  cvrt <- SlicedData$new()
  if (length(cvrt_rows) > 0) {
    m <- do.call(rbind, cvrt_rows); rownames(m) <- cvrt_names; colnames(m) <- samples
    cvrt$CreateFromMatrix(m)
  }

  out_file <- file.path(OUTDIR, paste0(label, "_results.tsv"))
  me <- Matrix_eQTL_engine(snps = snps, gene = gene, cvrt = cvrt, output_file_name = out_file,
                            pvOutputThreshold = 1, useModel = use_model,
                            errorCovariance = numeric(), verbose = FALSE, pvalue.hist = FALSE)
  res <- me$all$eqtls
  res$bonferroni_sig <- res$pvalue < bonferroni_alpha
  write.table(res, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("Significant (Bonferroni p<", format(bonferroni_alpha, scientific = TRUE), "):",
      sum(res$bonferroni_sig), "\n")
  res
}

pe_mtdna <- run_eqtl(mtdna_geno_raw, pe_mtdna_snp_pass, pe_meta$sample_id, pe_meta, "PE_mtDNA_SNP_eQTL",
                      NULL, pe_genes_pass, use_model = modelANOVA)
pe_ndna  <- run_eqtl(ndna_geno_raw,  pe_ndna_snp_pass,  pe_meta$sample_id, pe_meta, "PE_nDNA_SNP_eQTL",
                      pe_geno_pcs, pe_genes_pass, use_model = modelLINEAR)
se_mtdna <- run_eqtl(mtdna_geno_raw, se_mtdna_snp_pass, se_meta$sample_id, se_meta, "SE_mtDNA_SNP_eQTL",
                      NULL, se_genes_pass, use_model = modelANOVA)
se_ndna  <- run_eqtl(ndna_geno_raw,  se_ndna_snp_pass,  se_meta$sample_id, se_meta, "SE_nDNA_SNP_eQTL",
                      se_geno_pcs, se_genes_pass, use_model = modelLINEAR)

# -----------------------------------------------------------------------------
# STEP 6: Replication check
# -----------------------------------------------------------------------------
cat("\n############ STEP 6: Replication check ############\n")

check_replication <- function(pe_res, se_res, label) {
  if (is.null(pe_res) || is.null(se_res)) return(invisible(NULL))
  pe_sig <- pe_res[pe_res$bonferroni_sig, ]
  if (nrow(pe_sig) == 0) { cat("\n", label, ": no PE discovery hits\n"); return(invisible(NULL)) }
  merged <- merge(pe_sig, se_res, by = c("snps", "gene"), suffixes = c("_PE", "_SE"))
  cat("\n", label, ":\n")
  cat("  PE Bonferroni-sig hits:", nrow(pe_sig), " | Present in SE:", nrow(merged),
      " | Also Bonferroni-sig in SE:", sum(merged$bonferroni_sig_SE, na.rm = TRUE), "\n")
  write.table(merged, file.path(OUTDIR, paste0(label, "_replication.tsv")), sep = "\t", quote = FALSE, row.names = FALSE)
}

check_replication(pe_mtdna, se_mtdna, "mtDNA_SNP_eQTL")
check_replication(pe_ndna, se_ndna, "nDNA_SNP_eQTL")

cat("\n\nDone. All output in", OUTDIR, "\n")
cat("\nReminder of documented deviations from the paper's exact methods --\n")
cat("see comments at the top of this script.\n")