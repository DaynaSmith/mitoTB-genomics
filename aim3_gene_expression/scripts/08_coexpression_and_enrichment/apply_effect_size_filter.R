#!/usr/bin/env Rscript
# =============================================================================
# Apply a minimum effect-size threshold (|r|) on top of statistical
# significance for the nuclear co-expression replication results.
#
# With N=407 (PE) and 2.28M tests, statistical significance alone is reached
# by trivially small correlations -- this is standard practice in
# co-expression analysis (e.g. WGCNA-style methods always gate on
# correlation magnitude, not p-value alone) to separate meaningful
# co-regulation from noise/technical confounds.
#
# Re-filters the ALREADY-COMPUTED results -- no need to rerun the expensive
# genome-wide correlation step.
#
# Usage:
#   Rscript apply_effect_size_filter.R
# =============================================================================

# ---- EDIT THESE PATHS ----
REPLICATION_FILE <- "analysis1_eqtl/nuclear_coexpression_replication_check.tsv"
OUT_GENE_LIST     <- "analysis1_eqtl/replicating_coexpressed_nuclear_genes_filtered.txt"
OUT_FULL          <- "analysis1_eqtl/nuclear_coexpression_replication_filtered.tsv"
MIN_ABS_R         <- 0.5   # strict threshold -- see conversation notes on
                             # why 0.3 still left too many genes (6,270) to
                             # confidently distinguish real coordinated
                             # mitochondrial biogenesis biology from a
                             # persistent technical confound
# ---------------------------

merged <- read.delim(REPLICATION_FILE, stringsAsFactors = FALSE)
cat("Total replicating pairs before effect-size filter:", nrow(merged), "\n")
cat("Unique genes before filter:", length(unique(merged$nuclear_gene[merged$sig_in_SE])), "\n")

filtered <- merged[merged$sig_in_SE &
                      abs(merged$r_PE) >= MIN_ABS_R &
                      abs(merged$r_SE) >= MIN_ABS_R, ]

cat("\nAfter requiring |r| >=", MIN_ABS_R, "in BOTH cohorts:\n")
cat("  Pairs:", nrow(filtered), "\n")
replicating_genes <- unique(filtered$nuclear_gene)
cat("  Unique genes:", length(replicating_genes), "\n")

write.table(filtered, OUT_FULL, sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(sort(replicating_genes), OUT_GENE_LIST)

cat("\nWrote:", OUT_FULL, "\n")
cat("Wrote:", OUT_GENE_LIST, "\n")

# Show the sensitivity to different thresholds, for context
cat("\nSensitivity to different |r| thresholds (both cohorts):\n")
for (thresh in c(0.1, 0.2, 0.3, 0.4, 0.5)) {
  n <- sum(merged$sig_in_SE & abs(merged$r_PE) >= thresh & abs(merged$r_SE) >= thresh)
  ng <- length(unique(merged$nuclear_gene[merged$sig_in_SE & abs(merged$r_PE) >= thresh & abs(merged$r_SE) >= thresh]))
  cat("  |r| >=", thresh, ":", n, "pairs,", ng, "unique genes\n")
}
