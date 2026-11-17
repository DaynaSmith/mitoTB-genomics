## =============================================================================
## poly_C_length_heteroplasmy_composition.R
##
## Constructs per-sample length-heteroplasmy compositional profiles at ALL
## FOUR confirmed common-heteroplasmy positions (m.302, m.567, m.955,
## m.16182) in one run, following the method described in Chen et al. 2023,
## Nature (Methods: "Determination of chrM:302 length heteroplasmy
## composition"):
##
##   1. For each QC-passed sample, gather ALL calls (pre- and post-call-QC)
##      at the target position.
##   2. "QC-fail" estimate = sum of heteroplasmy fractions (AF) across calls
##      that failed call-level QC (call_pass_qc == FALSE).
##   3. "Reference" fraction = 1 - (sum of QC-pass allele AFs + QC-fail
##      estimate). Samples with NO calls at the position get reference = 1.
##   4. The named tracked alleles are reported individually; every other
##      allele (including the QC-fail estimate) is pooled into "Other".
##   5. Final step: any of {reference, tracked alleles} with a value <0.05
##      is relabelled into "Other" (matches the pipeline's own 5% HF floor).
##
## TRACKED ALLELES PER POSITION -- determined empirically from the final
## pipeline's QC-passed call set (variant_calls_annotated.csv), NOT assumed
## to be the same three alleles at every position. Total calls summed
## across all three cohorts:
##   m.302:   A>AC (200), A>ACC (172), A>ACCC (29)      -- matches Chen et
##            al.'s own choice for this position exactly.
##   m.955:   A>AC (18), A>ACC (16), A>ACCC (15)          -- same short-ladder
##            pattern as m.302.
##   m.567:   A>ACCCC (16), A>ACCCCC (15), A>ACCCCCC (15) -- NOTE: the
##            *longer* insertions dominate at this position (A>AC and
##            A>ACC are rare here, 3 and 2 calls respectively) -- using
##            the m.302 allele set here would have tracked the wrong,
##            least-common alleles.
##   m.16182: A>AC (24), A>ACC (14), AACC>A (6)           -- only two
##            insertion alleles were observed at this position, so the
##            third tracked category is a 4bp DELETION (AACC>A), not a
##            third insertion length. This is a genuinely different kind
##            of variant from the other three positions' pure insertion
##            ladders -- reported as its own category rather than folded
##            into "Other", since at 6 calls it's too common to treat as
##            residual noise, but its deletion/insertion status should be
##            called out explicitly wherever this position's composition
##            is discussed.
##
## KNOWN DATA LIMITATION: the depth<100 sample-exclusion criterion described
## in the source paper can only be directly applied to samples that already
## have >=1 call at the target position -- this VCF format does not record
## per-sample depth at positions/records where that sample has no variant
## call (FORMAT fields are ".:.:." for uncalled genotypes, confirmed by
## inspection of the raw VCF). Samples with zero calls at a position are
## therefore INCLUDED by default without a verified depth check. If a
## separate per-position depth file (e.g. samtools depth/mosdepth output)
## exists from the original alignment/calling pipeline, use it to apply the
## depth filter properly to all samples, not just called ones.
##
## Run after 01_mtDNA_cohort_QC.R (only needs variant_calls_all_flagged.rds
## and sample_metadata_QC_passed.rds; does not require scripts 02/03).
## =============================================================================

library(dplyr); library(tidyr); library(purrr); library(readr); library(tibble); library(ggplot2)

qc_dir <- "mtDNA_QC_outputs"
out_subdir <- file.path(qc_dir, "length_heteroplasmy_composition")
dir.create(out_subdir, showWarnings = FALSE)

vcf_long_all <- readRDS(file.path(qc_dir, "variant_calls_all_flagged.rds"))
sample_meta_pass <- readRDS(file.path(qc_dir, "sample_metadata_QC_passed.rds"))

## -----------------------------------------------------------------------
## 1. CONFIG
## -----------------------------------------------------------------------
min_position_depth <- 100   ## per Chen et al.; stricter than the pipeline-wide min_site_depth = 30
hf_relabel_cutoff <- 0.05

position_configs <- list(
  list(pos = 302,   tracked_alleles = c("AC", "ACC", "ACCC")),
  list(pos = 567,   tracked_alleles = c("ACCCC", "ACCCCC", "ACCCCCC")),
  list(pos = 955,   tracked_alleles = c("AC", "ACC", "ACCC")),
  list(pos = 16182, tracked_alleles = c("AC", "ACC", "A"))  ## "A" here is the AACC>A deletion allele
)

## cohort colour palette (soft/natural tones), keyed to short cohort names --
## consistent with the other pipeline scripts. Not applied to the stacked
## composition plot, since its fill maps to composition CATEGORY (reference/
## tracked alleles/other), not cohort -- no palette has been established for
## that variable, so it's left at ggplot's default.
cohort_colors_raw <- c(CAAPA = "#3B6FB6", BantuSpeakingAfrican = "#4F8A3F", SouthAfricanColoured = "#D77A2D")
cohort_name_map <- c(Bantu = "BantuSpeakingAfrican", SAC = "SouthAfricanColoured", CAAPA = "CAAPA")
cohort_colors <- stats::setNames(cohort_colors_raw[cohort_name_map], names(cohort_name_map))

theme_set(theme_bw(base_size = 11))

## -----------------------------------------------------------------------
## Relabel-small-to-other helper (shared across all positions)
## -----------------------------------------------------------------------
relabel_small_to_other <- function(df, cols, cutoff = hf_relabel_cutoff) {
  df <- as.data.frame(df)
  moved <- rep(0, nrow(df))
  for (cl in cols) {
    small <- df[[cl]] > 0 & df[[cl]] < cutoff
    moved[small] <- moved[small] + df[[cl]][small]
    df[[cl]][small] <- 0
  }
  df$other <- df$other + moved
  tibble::as_tibble(df)
}

## -----------------------------------------------------------------------
## Main per-position routine
## -----------------------------------------------------------------------
profile_position <- function(target_position, tracked_alleles) {
  
  message(sprintf("\n##### chrM:%d (tracked alleles: %s) #####", target_position, paste(tracked_alleles, collapse = ", ")))
  
  ## 2. Pull all calls at the target position, QC-pass samples only
  calls_at_pos <- vcf_long_all %>%
    dplyr::filter(POS == target_position, sample_id %in% sample_meta_pass$sample_id) %>%
    dplyr::mutate(AF = suppressWarnings(as.numeric(AF)))
  message(sprintf("Found %d call rows across %d distinct QC-passed samples.",
                  nrow(calls_at_pos), dplyr::n_distinct(calls_at_pos$sample_id)))
  
  ## 3. Depth filter (see header caveat)
  depth_at_pos <- calls_at_pos %>% dplyr::filter(!is.na(site_DP)) %>%
    dplyr::group_by(cohort, sample_id) %>% dplyr::summarise(depth_at_pos = max(site_DP, na.rm = TRUE), .groups = "drop")
  
  all_samples <- sample_meta_pass %>% dplyr::select(cohort, sample_id) %>%
    dplyr::left_join(depth_at_pos, by = c("cohort", "sample_id")) %>%
    dplyr::mutate(has_call = sample_id %in% unique(calls_at_pos$sample_id[!is.na(calls_at_pos$AF)]),
                  depth_known = !is.na(depth_at_pos))
  
  n_excluded_low_depth <- sum(all_samples$has_call & all_samples$depth_at_pos < min_position_depth, na.rm = TRUE)
  n_depth_unknown_no_call <- sum(!all_samples$has_call)
  message(sprintf("%d sample(s) excluded for depth <%d (verified, has a call there). %d sample(s) have no call -- depth unverifiable, included by default.",
                  n_excluded_low_depth, min_position_depth, n_depth_unknown_no_call))
  
  included_samples <- all_samples %>%
    dplyr::filter(!(has_call & depth_at_pos < min_position_depth)) %>%
    dplyr::select(cohort, sample_id)
  
  ## 4. Per-sample composition: QC-pass alleles (wide), QC-fail sum
  qc_pass_wide <- calls_at_pos %>%
    dplyr::filter(sample_id %in% included_samples$sample_id, call_pass_qc %in% TRUE, !is.na(AF)) %>%
    dplyr::group_by(cohort, sample_id, ALT) %>% dplyr::summarise(AF = sum(AF), .groups = "drop") %>%
    tidyr::pivot_wider(names_from = ALT, values_from = AF, values_fill = 0)
  
  qc_fail_sum <- calls_at_pos %>%
    dplyr::filter(sample_id %in% included_samples$sample_id, !(call_pass_qc %in% TRUE), !is.na(AF)) %>%
    dplyr::group_by(cohort, sample_id) %>% dplyr::summarise(qc_fail = sum(AF, na.rm = TRUE), .groups = "drop")
  
  composition <- included_samples %>%
    dplyr::left_join(qc_pass_wide, by = c("cohort", "sample_id")) %>%
    dplyr::left_join(qc_fail_sum, by = c("cohort", "sample_id")) %>%
    dplyr::mutate(dplyr::across(-c(cohort, sample_id), ~ tidyr::replace_na(.x, 0)))
  
  for (al in tracked_alleles) if (!al %in% names(composition)) composition[[al]] <- 0
  if (!"qc_fail" %in% names(composition)) composition$qc_fail <- 0
  
  all_qc_pass_allele_cols <- setdiff(names(qc_pass_wide), c("cohort", "sample_id"))
  other_allele_cols <- setdiff(all_qc_pass_allele_cols, tracked_alleles)
  
  composition <- composition %>%
    dplyr::rowwise() %>%
    dplyr::mutate(
      other_qc_pass_sum = sum(dplyr::c_across(dplyr::all_of(other_allele_cols)), na.rm = TRUE),
      tracked_sum = sum(dplyr::c_across(dplyr::all_of(tracked_alleles)), na.rm = TRUE),
      total_alt = tracked_sum + other_qc_pass_sum + qc_fail,
      reference = pmax(0, 1 - total_alt),
      other = other_qc_pass_sum + qc_fail
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(cohort, sample_id, reference, dplyr::all_of(tracked_alleles), other)
  
  ## 5. Final relabelling: any tracked category (incl. reference) <0.05 -> Other
  composition_final <- relabel_small_to_other(composition, c("reference", tracked_alleles))
  
  ## sanity check -- composition should sum to ~1 per sample
  composition_final <- composition_final %>%
    dplyr::mutate(row_sum = reference + rowSums(dplyr::across(dplyr::all_of(tracked_alleles))) + other)
  n_bad_sum <- sum(abs(composition_final$row_sum - 1) > 1e-6)
  if (n_bad_sum > 0) warning(sprintf("chrM:%d -- %d sample(s) have a composition not summing to 1.", target_position, n_bad_sum))
  
  readr::write_csv(composition_final %>% dplyr::select(-row_sum),
                   file.path(out_subdir, sprintf("composition_chrM_%d.csv", target_position)))
  
  ## 6. Summary + plots
  mean_comp <- composition_final %>% dplyr::group_by(cohort) %>%
    dplyr::summarise(n = dplyr::n(), dplyr::across(c(reference, dplyr::all_of(tracked_alleles), other), mean), .groups = "drop")
  message("Mean composition by cohort:")
  print(mean_comp, n = Inf)
  
  carrier_n <- composition_final %>% dplyr::group_by(cohort) %>%
    dplyr::summarise(dplyr::across(c(dplyr::all_of(tracked_alleles), other), ~ sum(.x > 0)), .groups = "drop")
  message("N carriers (composition value > 0) by cohort:")
  print(carrier_n, n = Inf)
  
  plot_long <- composition_final %>%
    dplyr::select(cohort, sample_id, reference, dplyr::all_of(tracked_alleles), other) %>%
    tidyr::pivot_longer(-c(cohort, sample_id), names_to = "category", values_to = "hf") %>%
    dplyr::mutate(category = factor(category, levels = c("reference", tracked_alleles, "other")))
  
  p_stack <- plot_long %>%
    dplyr::group_by(cohort, sample_id) %>%
    dplyr::arrange(sample_id, .by_group = TRUE) %>%
    ggplot(aes(x = reorder(sample_id, hf * (category == "reference")), y = hf, fill = category)) +
    geom_col(width = 1) + facet_wrap(~ cohort, scales = "free_x", ncol = 1) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
    labs(title = sprintf("Length-heteroplasmy composition at chrM:%d", target_position),
         x = "Sample (sorted by reference fraction)", y = "Fraction of composition", fill = NULL)
  ggsave(file.path(out_subdir, sprintf("composition_stacked_chrM_%d.png", target_position)), p_stack, width = 10, height = 8, dpi = 300)
  
  p_dist <- plot_long %>% dplyr::filter(category != "reference") %>%
    ggplot(aes(cohort, hf, fill = cohort)) + geom_boxplot(outlier.shape = NA, alpha = .7) +
    geom_jitter(width = .15, size = .5, alpha = .3) + facet_wrap(~ category, scales = "free_y") +
    scale_fill_manual(values = cohort_colors) +
    labs(title = sprintf("Distribution of non-reference composition at chrM:%d", target_position),
         x = NULL, y = "Fraction") + theme(legend.position = "none")
  ggsave(file.path(out_subdir, sprintf("composition_distribution_chrM_%d.png", target_position)), p_dist, width = 9, height = 6, dpi = 300)
  
  mean_comp %>% dplyr::mutate(POS = target_position, .before = 1)
}

## -----------------------------------------------------------------------
## Run for all four confirmed positions, collect a combined summary
## -----------------------------------------------------------------------
combined_summary <- purrr::map_dfr(position_configs, ~ profile_position(.x$pos, .x$tracked_alleles))
readr::write_csv(combined_summary, file.path(out_subdir, "composition_mean_summary_all_positions.csv"))

message(sprintf("\nDone. Per-position outputs and combined summary in: %s", normalizePath(out_subdir)))
