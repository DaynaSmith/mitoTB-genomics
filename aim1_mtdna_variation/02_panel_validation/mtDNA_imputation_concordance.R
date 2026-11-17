#!/usr/bin/env Rscript

# ==============================================================================
# mtDNA IMPUTATION GENOTYPE CONCORDANCE
#
# Purpose:
#   Evaluate genotype accuracy specifically at mtSNPs that were NOT directly
#   represented on the original genotyping array.
#
# Comparison:
#   Newly imputed mtSNPs vs complete mitochondrial sequence truth set
#
# Key outputs:
#   1. Number of sequence-derived untyped mtSNPs
#   2. Number/proportion recovered by imputation
#   3. Site-specific TP, TN, FP, FN
#   4. Site-specific Matthews correlation coefficient (MCC)
#   5. Mean and median site-specific MCC across estimable sites
#   6. Overall genotype concordance rate
#
# Notes:
#   - Directly genotyped array POSITIONS are excluded from imputation accuracy.
#   - Only biallelic SNP records are included in MCC calculations.
#   - mtDNA is treated as haploid.
#   - Homozygous diploid encodings (0/0, 1/1) are converted to haploid calls.
#   - Heterozygous/heteroplasmic calls (0/1, 1/0) are treated as missing.
# ==============================================================================


# ------------------------------------------------------------------------------
# 0. USER SETTINGS
# ------------------------------------------------------------------------------

array_vcf_path <- "MEGA_MT_geno0.05_mind0.1_SNPsOnly.vcf"
imputed_vcf_path <- "MEGA_imputed_geno0.05_mind0.1_SNPsOnly_INFO0.7.recode.vcf.gz"
sequence_vcf_path <- "mtDNA_allSamples_PASS_ExcludeBlacklistSNPs_SamplesForFinalAnalysis_RenameSamples.vcf"

out_dir <- "mtdna_imputation_concordance"

# Set to NULL to auto-detect MT / chrM / chrMT / M
mito_chr <- NULL


# ------------------------------------------------------------------------------
# 1. READ VCF
# ------------------------------------------------------------------------------

read_vcf <- function(path) {
  
  message("Reading: ", path)
  
  con <- if (grepl("\\.gz$", path, ignore.case = TRUE)) {
    gzfile(path, "rt")
  } else {
    file(path, "rt")
  }
  
  on.exit(close(con))
  
  lines <- readLines(con)
  
  header_idx <- which(startsWith(lines, "#CHROM"))
  
  if (length(header_idx) != 1) {
    stop("Could not identify a unique #CHROM header in: ", path)
  }
  
  col_names <- strsplit(
    sub("^#", "", lines[header_idx]),
    "\t"
  )[[1]]
  
  data_lines <- lines[(header_idx + 1):length(lines)]
  data_lines <- data_lines[nzchar(data_lines)]
  
  if (length(data_lines) == 0) {
    stop("No variant records found in: ", path)
  }
  
  fields <- strsplit(data_lines, "\t")
  
  if (any(lengths(fields) != length(col_names))) {
    stop("Malformed VCF records detected in: ", path)
  }
  
  mat <- matrix(
    unlist(fields),
    ncol = length(col_names),
    byrow = TRUE
  )
  
  colnames(mat) <- col_names
  
  fixed_cols <- c(
    "CHROM", "POS", "ID", "REF",
    "ALT", "QUAL", "FILTER", "INFO"
  )
  
  fix <- as.data.frame(
    mat[, fixed_cols, drop = FALSE],
    stringsAsFactors = FALSE
  )
  
  fix$POS <- as.integer(fix$POS)
  
  sample_names <- setdiff(
    col_names,
    c(fixed_cols, "FORMAT")
  )
  
  if (length(sample_names) == 0) {
    stop("No sample columns found in: ", path)
  }
  
  format_col <- mat[, "FORMAT"]
  
  geno_block <- mat[
    ,
    sample_names,
    drop = FALSE
  ]
  
  # Extract GT from FORMAT
  format_split <- strsplit(format_col, ":")
  
  gt_index <- vapply(
    format_split,
    function(x) {
      idx <- which(x == "GT")
      if (length(idx) == 0) NA_integer_ else idx[1]
    },
    integer(1)
  )
  
  gt <- matrix(
    NA_character_,
    nrow = nrow(geno_block),
    ncol = ncol(geno_block),
    dimnames = dimnames(geno_block)
  )
  
  for (i in seq_len(nrow(geno_block))) {
    
    idx <- gt_index[i]
    
    if (is.na(idx)) next
    
    vals <- strsplit(geno_block[i, ], ":")
    
    gt[i, ] <- vapply(
      vals,
      function(v) {
        if (length(v) >= idx) v[idx] else NA_character_
      },
      character(1)
    )
  }
  
  list(
    fix = fix,
    gt = gt,
    samples = sample_names
  )
}


# ------------------------------------------------------------------------------
# 2. RESTRICT TO MITOCHONDRIAL CONTIG
# ------------------------------------------------------------------------------

filter_mito <- function(vcf, mito_chr = NULL) {
  
  chrom <- vcf$fix$CHROM
  
  if (is.null(mito_chr)) {
    
    accepted <- c("MT", "chrM", "chrMT", "M")
    
    keep <- chrom %in% accepted
    
    if (!any(keep)) {
      stop(
        "Could not identify mitochondrial contig. Observed contigs include: ",
        paste(head(unique(chrom), 20), collapse = ", ")
      )
    }
    
  } else {
    
    keep <- chrom == mito_chr
    
    if (!any(keep)) {
      stop("No records found on mitochondrial contig: ", mito_chr)
    }
  }
  
  vcf$fix <- vcf$fix[keep, , drop = FALSE]
  vcf$gt <- vcf$gt[keep, , drop = FALSE]
  
  vcf
}


# ------------------------------------------------------------------------------
# 3. RETAIN BIALLELIC SNPs
# ------------------------------------------------------------------------------

filter_biallelic_snps <- function(vcf) {
  
  # Single REF nucleotide
  ref_snp <- nchar(vcf$fix$REF) == 1
  
  # Single ALT nucleotide
  alt_snp <- nchar(vcf$fix$ALT) == 1
  
  # Exclude multiallelic records
  single_alt <- !grepl(",", vcf$fix$ALT, fixed = TRUE)
  
  # Exclude symbolic/missing ALT
  valid_alt <- !(vcf$fix$ALT %in% c(".", "*"))
  
  keep <- ref_snp & alt_snp & single_alt & valid_alt
  
  message(
    "Retaining ",
    sum(keep),
    " of ",
    length(keep),
    " mitochondrial records as biallelic SNPs."
  )
  
  vcf$fix <- vcf$fix[keep, , drop = FALSE]
  vcf$gt <- vcf$gt[keep, , drop = FALSE]
  
  vcf
}


# ------------------------------------------------------------------------------
# 4. STANDARDISE mtDNA GENOTYPES
# ------------------------------------------------------------------------------

standardize_gt <- function(gt_mat) {
  
  out <- matrix(
    NA_integer_,
    nrow = nrow(gt_mat),
    ncol = ncol(gt_mat),
    dimnames = dimnames(gt_mat)
  )
  
  for (i in seq_len(nrow(gt_mat))) {
    
    for (j in seq_len(ncol(gt_mat))) {
      
      g <- gt_mat[i, j]
      
      if (
        is.na(g) ||
        g == "." ||
        g == "./." ||
        g == ".|."
      ) {
        next
      }
      
      alleles <- unlist(strsplit(g, "[/|]"))
      
      alleles <- alleles[
        alleles != "." &
          !is.na(alleles)
      ]
      
      if (length(alleles) == 0) next
      
      unique_alleles <- unique(alleles)
      
      # Heterozygous/heteroplasmic call -> missing
      if (length(unique_alleles) != 1) next
      
      allele <- suppressWarnings(
        as.integer(unique_alleles)
      )
      
      # Biallelic SNPs only: 0 = REF, 1 = ALT
      if (!is.na(allele) && allele %in% c(0, 1)) {
        out[i, j] <- allele
      }
    }
  }
  
  out
}


# ------------------------------------------------------------------------------
# 5. LOAD AND PREPARE DATA
# ------------------------------------------------------------------------------

prepare_vcf <- function(path, mito_chr = NULL) {
  
  v <- read_vcf(path)
  
  v <- filter_mito(
    v,
    mito_chr
  )
  
  v <- filter_biallelic_snps(v)
  
  calls <- standardize_gt(v$gt)
  
  # Harmonised key for exact variant comparison
  variant_id <- paste(
    v$fix$POS,
    toupper(v$fix$REF),
    toupper(v$fix$ALT),
    sep = "_"
  )
  
  # Remove duplicate variant records if present
  duplicated_id <- duplicated(variant_id)
  
  if (any(duplicated_id)) {
    
    warning(
      sum(duplicated_id),
      " duplicate variant IDs detected; keeping first occurrence."
    )
    
    keep <- !duplicated_id
    
    v$fix <- v$fix[keep, , drop = FALSE]
    calls <- calls[keep, , drop = FALSE]
    variant_id <- variant_id[keep]
  }
  
  rownames(calls) <- variant_id
  
  v$fix$variant_id <- variant_id
  
  list(
    fix = v$fix,
    calls = calls,
    samples = colnames(calls),
    variant_id = variant_id
  )
}


# ------------------------------------------------------------------------------
# 6. MCC FUNCTION
# ------------------------------------------------------------------------------

calculate_binary_mcc <- function(truth, predicted) {
  
  complete <- !is.na(truth) & !is.na(predicted)
  
  truth <- truth[complete]
  predicted <- predicted[complete]
  
  n <- length(truth)
  
  if (n == 0) {
    return(
      list(
        TP = 0,
        TN = 0,
        FP = 0,
        FN = 0,
        n = 0,
        concordance = NA_real_,
        mcc = NA_real_
      )
    )
  }
  
  # Force counts to numeric/double rather than integer
  TP <- as.numeric(sum(truth == 1 & predicted == 1))
  TN <- as.numeric(sum(truth == 0 & predicted == 0))
  FP <- as.numeric(sum(truth == 0 & predicted == 1))
  FN <- as.numeric(sum(truth == 1 & predicted == 0))
  
  concordance <- (TP + TN) / n
  
  # MCC denominator calculated using doubles to avoid integer overflow
  denominator <- sqrt(
    as.numeric(TP + FP) *
      as.numeric(TP + FN) *
      as.numeric(TN + FP) *
      as.numeric(TN + FN)
  )
  
  # MCC is undefined when one of the required classes is absent
  if (
    is.na(denominator) ||
    !is.finite(denominator) ||
    denominator == 0
  ) {
    mcc <- NA_real_
  } else {
    
    numerator <- TP * TN - FP * FN
    
    mcc <- numerator / denominator
  }
  
  list(
    TP = TP,
    TN = TN,
    FP = FP,
    FN = FN,
    n = n,
    concordance = concordance,
    mcc = mcc
  )
}

# ------------------------------------------------------------------------------
# 7. LOAD DATASETS
# ------------------------------------------------------------------------------

message("\nLoading sequence truth set ...")

sequence <- prepare_vcf(
  sequence_vcf_path,
  mito_chr
)

message("\nLoading original array data ...")

array <- prepare_vcf(
  array_vcf_path,
  mito_chr
)

message("\nLoading post-imputation data ...")

imputed <- prepare_vcf(
  imputed_vcf_path,
  mito_chr
)


# ------------------------------------------------------------------------------
# 8. DEFINE DIRECTLY GENOTYPED ARRAY POSITIONS
# ------------------------------------------------------------------------------

# IMPORTANT:
# We exclude by POSITION rather than exact REF/ALT representation.
#
# If an mtDNA position was interrogated by the original array, it is not treated
# as a newly imputed site even if REF/ALT representation differs after
# normalisation or decomposition.

array_positions <- unique(array$fix$POS)

message(
  "\nOriginal array contained ",
  length(array_positions),
  " unique mitochondrial SNP positions."
)


# ------------------------------------------------------------------------------
# 9. IDENTIFY SEQUENCE-DERIVED UNTYPED SNPs
# ------------------------------------------------------------------------------

sequence_untyped_ids <- sequence$fix$variant_id[
  !(sequence$fix$POS %in% array_positions)
]

sequence_untyped_positions <- unique(
  sequence$fix$POS[
    !(sequence$fix$POS %in% array_positions)
  ]
)

message(
  "Sequence truth set contained ",
  length(sequence_untyped_ids),
  " biallelic mtSNP records at positions not represented on the array."
)


# ------------------------------------------------------------------------------
# 10. IDENTIFY NEWLY IMPUTED SNPs
# ------------------------------------------------------------------------------

imputed_only_ids <- imputed$fix$variant_id[
  !(imputed$fix$POS %in% array_positions)
]

message(
  "Post-imputation dataset contained ",
  length(imputed_only_ids),
  " biallelic mtSNP records at positions not represented on the array."
)


# ------------------------------------------------------------------------------
# 11. IDENTIFY RECOVERED UNTYPED SEQUENCE SNPs
# ------------------------------------------------------------------------------

recovered_ids <- intersect(
  sequence_untyped_ids,
  imputed_only_ids
)

n_sequence_untyped <- length(sequence_untyped_ids)
n_recovered <- length(recovered_ids)

recovery_rate <- if (
  n_sequence_untyped > 0
) {
  n_recovered / n_sequence_untyped
} else {
  NA_real_
}

message(
  "\nRecovered ",
  n_recovered,
  " of ",
  n_sequence_untyped,
  " untyped sequence-derived mtSNPs (",
  round(100 * recovery_rate, 2),
  "%)."
)


# ------------------------------------------------------------------------------
# 12. MATCH SAMPLES
# ------------------------------------------------------------------------------

common_samples <- intersect(
  sequence$samples,
  imputed$samples
)

if (length(common_samples) == 0) {
  stop("No overlapping samples between sequence and imputed datasets.")
}

message(
  "Number of samples available for genotype comparison: ",
  length(common_samples)
)


# ------------------------------------------------------------------------------
# 13. SITE-SPECIFIC MCC
# ------------------------------------------------------------------------------

if (length(recovered_ids) == 0) {
  stop(
    "No newly imputed mtSNPs overlap the sequence truth set."
  )
}

site_results <- do.call(
  rbind,
  lapply(
    recovered_ids,
    function(v) {
      
      truth <- sequence$calls[
        v,
        common_samples
      ]
      
      predicted <- imputed$calls[
        v,
        common_samples
      ]
      
      res <- calculate_binary_mcc(
        truth,
        predicted
      )
      
      seq_row <- sequence$fix[
        sequence$fix$variant_id == v,
        ,
        drop = FALSE
      ][1, ]
      
      data.frame(
        variant_id = v,
        POS = seq_row$POS,
        REF = seq_row$REF,
        ALT = seq_row$ALT,
        n_samples_compared = res$n,
        TP = res$TP,
        TN = res$TN,
        FP = res$FP,
        FN = res$FN,
        genotype_concordance = res$concordance,
        MCC = res$mcc,
        stringsAsFactors = FALSE
      )
    }
  )
)


# ------------------------------------------------------------------------------
# 14. SUMMARISE SITE-SPECIFIC PERFORMANCE
# ------------------------------------------------------------------------------

estimable_mcc <- site_results$MCC[
  !is.na(site_results$MCC)
]

mean_mcc <- if (
  length(estimable_mcc) > 0
) {
  mean(estimable_mcc)
} else {
  NA_real_
}

median_mcc <- if (
  length(estimable_mcc) > 0
) {
  median(estimable_mcc)
} else {
  NA_real_
}

min_mcc <- if (
  length(estimable_mcc) > 0
) {
  min(estimable_mcc)
} else {
  NA_real_
}

max_mcc <- if (
  length(estimable_mcc) > 0
) {
  max(estimable_mcc)
} else {
  NA_real_
}


# Overall genotype concordance across all comparable genotype calls
total_TP <- sum(site_results$TP)
total_TN <- sum(site_results$TN)
total_FP <- sum(site_results$FP)
total_FN <- sum(site_results$FN)

total_pairs <- (
  total_TP +
    total_TN +
    total_FP +
    total_FN
)

overall_concordance <- if (
  total_pairs > 0
) {
  (total_TP + total_TN) / total_pairs
} else {
  NA_real_
}


summary_results <- data.frame(
  
  n_array_positions = length(array_positions),
  
  n_sequence_untyped_mtSNPs = n_sequence_untyped,
  
  n_recovered_imputed_mtSNPs = n_recovered,
  
  recovery_rate = recovery_rate,
  
  n_common_samples = length(common_samples),
  
  n_sites_with_estimable_MCC = length(estimable_mcc),
  
  mean_site_MCC = mean_mcc,
  
  median_site_MCC = median_mcc,
  
  minimum_site_MCC = min_mcc,
  
  maximum_site_MCC = max_mcc,
  
  overall_genotype_concordance = overall_concordance,
  
  stringsAsFactors = FALSE
)


# ------------------------------------------------------------------------------
# 15. WRITE OUTPUTS
# ------------------------------------------------------------------------------

if (!dir.exists(out_dir)) {
  dir.create(
    out_dir,
    recursive = TRUE
  )
}

write.csv(
  site_results,
  file.path(
    out_dir,
    "imputed_only_site_specific_MCC.csv"
  ),
  row.names = FALSE
)

write.csv(
  summary_results,
  file.path(
    out_dir,
    "imputed_only_MCC_summary.csv"
  ),
  row.names = FALSE
)

writeLines(
  recovered_ids,
  file.path(
    out_dir,
    "recovered_imputed_variant_ids.txt"
  )
)


# ------------------------------------------------------------------------------
# 16. REPORT
# ------------------------------------------------------------------------------

cat(
  "\n============================================================\n",
  " mtDNA IMPUTATION GENOTYPE CONCORDANCE\n",
  "============================================================\n\n",
  sep = ""
)

cat(
  "Directly genotyped array positions: ",
  length(array_positions),
  "\n",
  sep = ""
)

cat(
  "Sequence-derived mtSNPs at untyped positions: ",
  n_sequence_untyped,
  "\n",
  sep = ""
)

cat(
  "Untyped sequence-derived mtSNPs recovered by imputation: ",
  n_recovered,
  "\n",
  sep = ""
)

cat(
  "Recovery rate: ",
  round(100 * recovery_rate, 2),
  "%\n",
  sep = ""
)

cat(
  "Samples compared: ",
  length(common_samples),
  "\n",
  sep = ""
)

cat(
  "Recovered sites with estimable MCC: ",
  length(estimable_mcc),
  " / ",
  n_recovered,
  "\n",
  sep = ""
)

cat(
  "Mean site-specific MCC: ",
  round(mean_mcc, 4),
  "\n",
  sep = ""
)

cat(
  "Median site-specific MCC: ",
  round(median_mcc, 4),
  "\n",
  sep = ""
)

cat(
  "Site-specific MCC range: ",
  round(min_mcc, 4),
  " to ",
  round(max_mcc, 4),
  "\n",
  sep = ""
)

cat(
  "Overall genotype concordance: ",
  round(overall_concordance, 4),
  "\n",
  sep = ""
)

cat(
  "\nResults written to: ",
  normalizePath(out_dir),
  "\n",
  sep = ""
)

cat(
  "\n============================================================\n"
)
