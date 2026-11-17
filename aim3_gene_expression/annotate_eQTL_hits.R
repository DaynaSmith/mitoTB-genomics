#!/usr/bin/env Rscript
# =============================================================================
# Annotate the TRULY replicating hits (Bonferroni-sig in BOTH PE and SE --
# the *_replication.tsv files contain all PE-sig hits merely PRESENT in SE,
# not pre-filtered to double-significance, so we filter explicitly here).
#
# Usage:
#   Rscript annotate_v3_hits.R
# =============================================================================

# ---- EDIT THESE PATHS ----
OUTDIR         <- "analysis1_eqtl_v3"
MTDNA_GTF      <- "chrM_genes.gtf"
NDNA_BED       <- "mito_nuclear_genes.bed"
GENE_INFO_FILE <- "mtdna_expression_matrix.tsv"
# ---------------------------

gene_info <- read.delim(GENE_INFO_FILE, stringsAsFactors = FALSE)[, c("gene_id", "gene_name")]

# Length-safe attribute extraction -- regmatches()/regexpr() silently drops
# non-matching elements rather than returning NA, which shortens the vector
# and breaks column assignment. sapply here guarantees one output per input.
# NOTE: this GTF's attributes are UNQUOTED (e.g. "gene_id ENSG...;" not
# 'gene_id "ENSG...";'), unlike standard GENCODE format -- pattern matches
# everything between the field name and the next semicolon, no quotes.
extract_gtf_field <- function(attr_col, field) {
  sapply(attr_col, function(s) {
    m <- regmatches(s, regexpr(paste0(field, " [^;]+;"), s))
    if (length(m) == 0) return(NA_character_)
    val <- sub(paste0("^", field, " "), "", m)
    val <- sub(";$", "", val)
    trimws(val)
  }, USE.NAMES = FALSE)
}

annotate_snp_gene <- function(rep_file, bed_or_gtf_path, is_gtf, label) {
  rep <- read.delim(file.path(OUTDIR, rep_file), stringsAsFactors = FALSE)
  cat("\n=====", label, "=====\n")
  cat("Rows in replication file (PE-sig, present in SE):", nrow(rep), "\n")

  # Filter to TRUE double-replication: Bonferroni-significant in BOTH cohorts
  sig_col <- grep("^bonferroni_sig_SE$|^fdr_sig_SE$", colnames(rep), value = TRUE)[1]
  if (is.na(sig_col)) stop("Could not find an SE-significance column -- check actual column names: ",
                            paste(colnames(rep), collapse = ", "))
  rep <- rep[!is.na(rep[[sig_col]]) & rep[[sig_col]], ]
  cat("Rows after filtering to double-significant (", sig_col, "):", nrow(rep), "\n", sep = "")

  if (nrow(rep) == 0) {
    cat("No doubly-significant hits -- nothing to annotate.\n")
    return(invisible(NULL))
  }

  if (is_gtf) {
    # mtDNA-specific: SNP IDs are "NC_012920.1_POS_REF_ALT" -- the chromosome
    # name itself contains an underscore, so naive split-by-underscore
    # incorrectly grabs "NC" as chrom and "012920.1" as pos. Also,
    # chrM_genes.gtf (from GENCODE) uses UCSC-style "chrM" naming, while our
    # genotype matrix uses RefSeq-style "NC_012920.1" -- same genome,
    # different label, so we normalize to "chrM" to match the GTF.
    rep$chrom <- "chrM"
    rep$pos <- as.numeric(sub("^NC_012920\\.1_([0-9]+)_.*$", "\\1", rep$snps))
  } else {
    # nDNA: chrom names (chr1, chr5, etc.) have no internal underscore, so
    # simple split-by-underscore correctly separates chrom from pos.
    snp_parts <- strsplit(rep$snps, "_")
    rep$chrom <- sapply(snp_parts, `[`, 1)
    rep$pos <- as.numeric(sapply(snp_parts, `[`, 2))
  }

  if (is_gtf) {
    bed <- read.delim(bed_or_gtf_path, header = FALSE, comment.char = "#", stringsAsFactors = FALSE)
    bed <- bed[bed$V3 == "gene", ]
    cat("GTF gene rows found:", nrow(bed), "\n")
    cat("Example raw attribute string:", substr(bed$V9[1], 1, 150), "\n")

    # Extract gene_id (simpler/more standardized than gene_name) and join
    # against the already-clean gene_id -> gene_name mapping from
    # mtdna_expression_matrix.tsv, rather than re-parsing gene_name
    # directly from GTF attributes (which was silently failing).
    bed$gene_id_gtf <- extract_gtf_field(bed$V9, "gene_id")
    cat("Gene IDs successfully extracted:", sum(!is.na(bed$gene_id_gtf)), "/", nrow(bed), "\n")
    bed <- merge(bed, gene_info, by.x = "gene_id_gtf", by.y = "gene_id", all.x = TRUE)
    cat("Gene names resolved via join:", sum(!is.na(bed$gene_name)), "/", nrow(bed), "\n")
    bed <- data.frame(chrom = bed$V1, start = bed$V4, end = bed$V5, gene_name = bed$gene_name)
  } else {
    bed <- read.delim(bed_or_gtf_path, header = FALSE, stringsAsFactors = FALSE,
                       col.names = c("chrom", "start", "end", "gene_name"))
  }

  match_gene <- function(chrom, pos) {
    hits <- bed[bed$chrom == chrom & bed$start <= pos & bed$end >= pos, "gene_name"]
    if (length(hits) == 0) return(NA)
    paste(unique(hits), collapse = ";")
  }
  rep$snp_gene_name <- mapply(match_gene, rep$chrom, rep$pos)

  rep <- merge(rep, gene_info, by.x = "gene", by.y = "gene_id", all.x = TRUE)
  colnames(rep)[colnames(rep) == "gene_name"] <- "mtdna_outcome_gene"

  n_unique_snp_genes <- length(unique(na.omit(rep$snp_gene_name)))
  n_unique_outcome_genes <- length(unique(na.omit(rep$mtdna_outcome_gene)))
  cat("Unique SNP-side genes:", n_unique_snp_genes, "\n")
  cat("Unique mtDNA outcome genes:", n_unique_outcome_genes, "\n")

  out_cols <- intersect(c("snps", "chrom", "pos", "snp_gene_name", "mtdna_outcome_gene",
                           "beta_PE", "pvalue_PE", "beta_SE", "pvalue_SE"), colnames(rep))
  rep <- rep[order(rep$pvalue_PE), out_cols]

  cat("\nAll doubly-significant replicating pairs:\n")
  print(rep)

  out_path <- file.path(OUTDIR, paste0(label, "_annotated.tsv"))
  write.table(rep, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("\nWrote:", out_path, "\n")
}

annotate_snp_gene("mtDNA_SNP_eQTL_replication.tsv", MTDNA_GTF, is_gtf = TRUE, label = "mtDNA_SNP")
annotate_snp_gene("nDNA_SNP_eQTL_replication.tsv", NDNA_BED, is_gtf = FALSE, label = "nDNA_SNP")