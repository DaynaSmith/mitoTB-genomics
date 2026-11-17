#!/usr/bin/env Rscript
# =============================================================================
# GO enrichment on the replicating co-expressed nuclear genes, calling
# g:Profiler's REST API directly via wget (bypassing the gprofiler2/RCurl
# package, which requires curl-config -- a system library not available on
# this HPC and not installable without root).
#
# Needs outbound internet access (wget -> biit.cs.ut.ee) from wherever this
# runs -- same consideration as everything else we've called out over the
# network on this cluster.
#
# Usage:
#   Rscript 46_go_enrichment.R
# =============================================================================

if (!requireNamespace("jsonlite", quietly = TRUE)) {
  install.packages("jsonlite", repos = "https://cloud.r-project.org",
                    lib = Sys.getenv("R_LIBS_USER"))
}
library(jsonlite)

# ---- EDIT THESE PATHS ----
REPLICATING_GENES <- "analysis1_eqtl/replicating_coexpressed_nuclear_genes_filtered.txt"
OUT_FILE           <- "analysis1_eqtl/go_enrichment_results.tsv"
REQUEST_JSON        <- "analysis1_eqtl/gprofiler_request.json"
RESPONSE_JSON       <- "analysis1_eqtl/gprofiler_response.json"
API_URL             <- "https://biit.cs.ut.ee/gprofiler/api/gost/profile/"
# ---------------------------

target_genes <- sort(unique(readLines(REPLICATING_GENES)))
cat("Genes for GO enrichment:", length(target_genes), "genes\n")

payload <- list(
  organism = "hsapiens",
  query = target_genes,
  sources = list("GO:BP", "GO:MF", "GO:CC"),
  user_threshold = 0.05,
  all_results = TRUE
)

json_body <- toJSON(payload, auto_unbox = TRUE)
writeLines(json_body, REQUEST_JSON)

cat("Sending request to g:Profiler API via wget...\n")
exit_code <- system2("wget",
  args = c(
    "--header=Content-Type: application/json",
    paste0("--post-file=", REQUEST_JSON),
    "-O", RESPONSE_JSON,
    API_URL
  )
)

if (exit_code != 0 || !file.exists(RESPONSE_JSON)) {
  stop("wget request failed (exit code ", exit_code, ") -- check network access from this node")
}

resp <- fromJSON(RESPONSE_JSON, flatten = TRUE)

if (is.null(resp$result) || (is.data.frame(resp$result) && nrow(resp$result) == 0)) {
  cat("\nNo GO terms returned. Full response for debugging:\n")
  str(resp)
} else {
  res <- resp$result
  res <- res[order(res$p_value), ]

  keep_cols <- intersect(c("source", "term_id", "term_name", "p_value",
                            "term_size", "intersection_size"), colnames(res))
  write.table(res[, keep_cols], OUT_FILE, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("\nWrote:", OUT_FILE, "\n")

  n_sig <- sum(res$p_value < 0.05, na.rm = TRUE)
  cat("Terms with p <0.05:", n_sig, "/", nrow(res), "total terms returned\n")

  cat("\nTop 20 terms by p-value:\n")
  print(head(res[, keep_cols], 20))
}
