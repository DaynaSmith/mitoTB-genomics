#!/usr/bin/env Rscript
# Plot corrected counts; this script does NOT filter or match VCF records.
# Input: platform,n_imputed,n_seq,n_common (one row per platform).
# Both sets must contain biallelic SNP records at positions absent from that
# platform's array, using the same allele-aware matching rule as the analysis.
# Keep recovered records with undefined MCC in these counts.
# Usage: Rscript venn_diagrams_updated.R corrected_counts.csv venn_L0d

for (pkg in c("VennDiagram", "gridExtra")) {
  if (!requireNamespace(pkg, quietly = TRUE))
    stop("Install the required package: install.packages('", pkg, "')")
}
library(grid)

match_column <- function(df, aliases, patterns) {
  headers <- tolower(names(df))
  hits <- which(headers %in% aliases)
  if (!length(hits)) hits <- which(vapply(headers, function(x)
    any(vapply(patterns, function(p) grepl(p, x), logical(1))), logical(1)))
  if (length(hits) != 1L)
    stop("Expected exactly one column for ", aliases[1],
         "; found ", length(hits), ". Use the canonical headers: ",
         "platform,n_imputed,n_seq,n_common")
  names(df)[hits]
}

load_counts <- function(path) {
  df <- read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!nrow(df)) stop("The input CSV contains no rows.")
  keys <- c(
    match_column(df, c("platform"), c("^array$", "^chip$", "^name$")),
    match_column(df, c("n_imputed", "total_imputed"), c("imput")),
    match_column(df, c("n_seq", "total_sequenced"), c("sequenc", "^seq")),
    match_column(df, c("n_common", "total_common"), c("common", "overlap", "shared", "intersect")))
  if (anyDuplicated(keys)) stop("A column was matched to multiple quantities.")
  out <- data.frame(platform = trimws(as.character(df[[keys[1]]])),
                    n_imputed = suppressWarnings(as.numeric(df[[keys[2]]])),
                    n_seq = suppressWarnings(as.numeric(df[[keys[3]]])),
                    n_common = suppressWarnings(as.numeric(df[[keys[4]]])))
  if (anyNA(out$platform) || any(!nzchar(out$platform)) || anyDuplicated(out$platform))
    stop("Platform names must be nonempty and unique.")
  vals <- as.matrix(out[, c("n_imputed", "n_seq", "n_common")])
  if (any(!is.finite(vals)) || any(vals < 0) || any(vals != floor(vals)))
    stop("All counts must be finite, nonnegative integers.")
  if (any(out$n_seq == 0)) stop("Sequenced totals must be positive to calculate recovery.")
  if (any(out$n_common > pmin(out$n_imputed, out$n_seq)))
    stop("Shared counts cannot exceed either set total.")
  out$n_seq_only <- out$n_seq - out$n_common
  out$n_imputed_only <- out$n_imputed - out$n_common
  out$recovery_percent <- 100 * out$n_common / out$n_seq
  out$unrecovered_percent <- 100 - out$recovery_percent
  out
}

platform_label <- function(x) {
  labels <- c("Affymetrix6.0" = "Affymetrix 6.0", "H3Africa" = "H3Africa",
              "InfiniumGDA" = "Infinium GDA", "InfiniumOmni5-4" = "Infinium Omni5-4",
              "MultiEthnicAMR-AFR" = "Multi-Ethnic AMR/AFR",
              "MultiEthnicGlobal" = "Multi-Ethnic Global")
  if (x %in% names(labels)) unname(labels[x]) else x
}

make_venn_grob <- function(row) {
  # Sequenced reference = grey; imputed = teal. Fixed geometry makes panels
  # comparable and accommodates disjoint, contained, identical or empty sets.
  # Counts, rather than circle areas, communicate set sizes.
  venn <- VennDiagram::draw.pairwise.venn(
    area1 = row$n_seq, area2 = row$n_imputed, cross.area = row$n_common,
    category = c("Sequenced", "Imputed"),
    fill = c("#BDBDBD", "#35978F"), alpha = 0.8,
    col = c("#666666", "#23665F"), lwd = 1,
    fontfamily = "sans", cex = 1.15, cat.cex = 1,
    cat.fontfamily = "sans", cat.pos = c(-20, 20), cat.dist = c(0.04, 0.04),
    scaled = FALSE, euler.d = FALSE, ind = FALSE)
  body <- gTree(children = venn, vp = viewport(width = 0.9, height = 0.82))
  gridExtra::arrangeGrob(body,
    top = textGrob(platform_label(row$platform), gp = gpar(fontsize = 13, fontface = "bold")),
    bottom = textGrob(sprintf("Recovery: %.1f%%", row$recovery_percent),
                      gp = gpar(fontsize = 10)))
}

save_plot <- function(g, stem, width, height) {
  pdf(paste0(stem, ".pdf"), width = width, height = height,
      family = "sans", useDingbats = FALSE, bg = "white")
  tryCatch({grid.newpage(); grid.draw(g)}, finally = dev.off())
  png(paste0(stem, ".png"), width = width, height = height,
      units = "in", res = 600, bg = "white")
  tryCatch({grid.newpage(); grid.draw(g)}, finally = dev.off())
}

run_venn_diagrams <- function(input_csv, out_dir = "venn_diagrams") {
  counts <- load_counts(input_csv)
  order_keys <- c("Affymetrix6.0", "H3Africa", "InfiniumGDA", "InfiniumOmni5-4",
                  "MultiEthnicAMR-AFR", "MultiEthnicGlobal")
  ranks <- match(counts$platform, order_keys)
  counts <- counts[order(is.na(ranks), ranks, counts$platform), , drop = FALSE]
  safe <- gsub("[^A-Za-z0-9_-]+", "_", counts$platform)
  if (anyDuplicated(safe)) stop("Platform names produce duplicate output filenames.")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  grobs <- lapply(seq_len(nrow(counts)), function(i) make_venn_grob(counts[i, ]))
  for (i in seq_along(grobs))
    save_plot(grobs[[i]], file.path(out_dir, paste0("venn_", safe[i])), 4.2, 4.2)
  nc <- min(3L, length(grobs)); nr <- ceiling(length(grobs) / nc)
  combined <- gridExtra::arrangeGrob(grobs = grobs, ncol = nc,
    bottom = textGrob("Biallelic mtSNP records at untyped positions; circle areas are not proportional to counts.",
                      gp = gpar(fontsize = 9)))
  save_plot(combined, file.path(out_dir, "venn_all_platforms_combined"),
            4.2 * nc, 4.2 * nr + 0.3)
  write.csv(counts, file.path(out_dir, "variant_recovery_summary.csv"), row.names = FALSE)
  print(counts)
  message("Saved PDF, 600-dpi PNG and recovery summary to: ", normalizePath(out_dir))
  invisible(counts)
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly = TRUE)
  run_venn_diagrams(if (length(args) >= 1) args[1] else "corrected_counts.csv",
                    if (length(args) >= 2) args[2] else "venn_diagrams")
}
