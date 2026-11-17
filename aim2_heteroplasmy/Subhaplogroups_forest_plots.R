## =============================================================================
## irr_forest_subhaplogroup_tidy.R
##
## Rebuilds the subhaplogroup IRR forest plot (irr_forest_subhaplogroup_by_cohort.png)
## from the saved regression results CSV, with tidied labels throughout:
##   - cohort facet labels: "Bantu" -> "Bantu-speaking African", "SAC" ->
##     "South African Coloured" (CAAPA unchanged)
##   - outcome facet labels: n_intermediate_het / n_low_het / n_total_het ->
##     full readable phrases, no underscores
##   - term (y-axis) labels: broad_haplogroup_grpL0a1 etc -> "Subhaplogroup: L0a1"
##     etc, no underscores
##
## Same free x AND y scales as the original (scales = "free") -- each
## cohort genuinely has a different set of subhaplogroups, so the y-axis
## is not shared across rows; this is intentional, not a labelling
## oversight, and is preserved here.
##
## No colour/significance mapping, matching the original: these are
## exploratory, NOT Bonferroni-corrected coefficients, so they're shown as
## plain black points/intervals rather than flagged red/grey like the
## macro-haplogroup plot.
##
## Does NOT change any statistics -- purely a labelling pass.
##
## Reads: mtDNA_QC_outputs/association_analysis/negbin_regression_results_by_cohort.csv
## (written by 03_heteroplasmy_haplogroup_association.R)
## =============================================================================

library(dplyr); library(stringr); library(ggplot2); library(forcats)

stats_dir <- "mtDNA_QC_outputs/association_analysis"
out_path <- file.path(stats_dir, "plots", "irr_forest_subhaplogroup_by_cohort_tidy.png")

tidy_all <- read.csv(file.path(stats_dir, "negbin_regression_results_by_cohort.csv"), stringsAsFactors = FALSE)

## -----------------------------------------------------------------------
## Label lookups
## -----------------------------------------------------------------------
cohort_full <- c(
  Bantu = "Bantu-speaking African",
  SAC   = "South African Coloured",
  CAAPA = "CAAPA"
)

outcome_full <- c(
  n_low_het          = "Low-level Heteroplasmy",
  n_intermediate_het = "Intermediate-level Heteroplasmy",
  n_total_het        = "Total Heteroplasmy"
)

## -----------------------------------------------------------------------
## Prepare plotting data
## -----------------------------------------------------------------------
sub_df <- tidy_all %>%
  dplyr::filter(haplogroup_resolution == "sub", stringr::str_starts(term, "broad_haplogroup_grp")) %>%
  dplyr::mutate(
    cohort_label = factor(cohort_full[cohort], levels = cohort_full[c("Bantu", "SAC", "CAAPA")]),
    outcome_label = factor(outcome_full[outcome],
                           levels = outcome_full[c("n_low_het", "n_intermediate_het", "n_total_het")]),
    subhaplogroup = stringr::str_remove(term, "^broad_haplogroup_grp"),
    term_label = paste0("Subhaplogroup: ", subhaplogroup)
  )

## Deliberate y-axis order (ascending, so the first level here renders at
## the BOTTOM of each panel and the last at the TOP): alphabetical by
## subhaplogroup code, with "Other" forced to the very end/top regardless
## of where it would otherwise alphabetise to, matching the same
## "Other last" convention used in the tidied macro-haplogroup plot
subhaplo_order_raw <- sub_df %>% dplyr::distinct(subhaplogroup) %>%
  dplyr::arrange(subhaplogroup == "Other", subhaplogroup) %>% dplyr::pull(subhaplogroup)
term_order <- paste0("Subhaplogroup: ", subhaplo_order_raw)

sub_df <- sub_df %>% dplyr::mutate(term_label = factor(term_label, levels = term_order))

## -----------------------------------------------------------------------
## Plot (same structure/statistics as the original, tidied labels only;
## no colour mapping -- these coefficients are exploratory/uncorrected)
## -----------------------------------------------------------------------
p_forest_sub_tidy <- ggplot(sub_df, aes(x = IRR, y = term_label)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_pointrange(aes(xmin = IRR_lower, xmax = IRR_upper)) +
  scale_x_log10() +
  facet_grid(cohort_label ~ outcome_label, scales = "free") +
  labs(
    title = "Subhaplogroup Coefficients, per Cohort (Exploratory -- Not Bonferroni-corrected)",
    x = "IRR (log scale)", y = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

dir.create(file.path(stats_dir, "plots"), showWarnings = FALSE)
ggsave(out_path, p_forest_sub_tidy, width = 13, height = 10, dpi = 300, bg = "white")

message(sprintf("Saved: %s", normalizePath(out_path)))