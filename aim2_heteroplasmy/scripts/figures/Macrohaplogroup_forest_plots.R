## =============================================================================
## irr_forest_macro_tidy.R
##
## Rebuilds the macro-haplogroup IRR forest plot (irr_forest_macro_by_cohort.png)
## from the saved regression results CSV, with tidied labels throughout:
##   - cohort facet labels: "Bantu" -> "Bantu-speaking African", "SAC" ->
##     "South African Coloured" (CAAPA unchanged, already a proper name)
##   - outcome facet labels: n_intermediate_het / n_low_het / n_total_het ->
##     full readable phrases, no underscores
##   - term (y-axis) labels: sexMale/mtcn/mean_depth/macro_haplogroup_grp*
##     -> full readable phrases, no underscores, in a deliberate display
##     order (demographic/biological covariates first, then macro-haplogroup
##     levels) rather than the original's incidental alphabetical order
##
## Does NOT change any statistics -- same IRRs, CIs, and Bonferroni
## significance flags as the original plot, purely a labelling pass.
##
## Reads: mtDNA_QC_outputs/association_analysis/negbin_regression_results_by_cohort.csv
## (written by 03_heteroplasmy_haplogroup_association.R)
## =============================================================================

library(dplyr); library(stringr); library(ggplot2); library(forcats)

stats_dir <- "mtDNA_QC_outputs/association_analysis"
out_path <- file.path(stats_dir, "plots", "irr_forest_macro_by_cohort_tidy.png")

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

## Term relabeling: exact matches for the fixed covariates, pattern-based
## for the macro-haplogroup levels (level suffix varies: Other/L1/L2/L3/...)
clean_term_label <- function(term) {
  dplyr::case_when(
    term == "age"              ~ "Age",
    term == "sexMale"          ~ "Sex (Male)",
    term == "mtcn"              ~ "mtDNA Copy Number",
    term == "mean_depth"       ~ "Mean Sequencing Depth",
    stringr::str_starts(term, "macro_haplogroup_grp") ~
      paste0("Macro-haplogroup: ", stringr::str_remove(term, "^macro_haplogroup_grp")),
    TRUE ~ term
  )
}

## -----------------------------------------------------------------------
## Prepare plotting data: clean labels, deliberate display order
## -----------------------------------------------------------------------
macro_df <- tidy_all %>%
  dplyr::filter(haplogroup_resolution == "macro") %>%
  dplyr::mutate(
    cohort_label = factor(cohort_full[cohort], levels = cohort_full[c("Bantu", "SAC", "CAAPA")]),
    outcome_label = factor(outcome_full[outcome],
                           levels = outcome_full[c("n_low_het", "n_intermediate_het", "n_total_het")]),
    term_label = clean_term_label(term)
  )

## Deliberate y-axis order: fixed covariates first (in a natural reading
## order), then macro-haplogroup levels sorted (numeric L-levels first,
## then any single-letter macro codes, then "Other" last)
fixed_order <- c("Age", "Sex (Male)", "mtDNA Copy Number", "Mean Sequencing Depth")
haplo_levels_present <- macro_df %>%
  dplyr::filter(stringr::str_starts(term_label, "Macro-haplogroup: ")) %>%
  dplyr::pull(term_label) %>% unique()
haplo_suffix <- stringr::str_remove(haplo_levels_present, "^Macro-haplogroup: ")
haplo_order <- haplo_levels_present[order(haplo_suffix == "Other", haplo_suffix)]  ## "Other" sorted last

term_order <- rev(c(fixed_order, haplo_order))  ## rev() so first-listed ends up at top of plot
macro_df <- macro_df %>% dplyr::mutate(term_label = factor(term_label, levels = term_order))

## -----------------------------------------------------------------------
## Plot (same structure/statistics as the original, tidied labels only)
## -----------------------------------------------------------------------
p_forest_macro_tidy <- ggplot(macro_df, aes(x = IRR, y = term_label, colour = sig_bonferroni)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_pointrange(aes(xmin = IRR_lower, xmax = IRR_upper)) +
  scale_x_log10() +
  scale_colour_manual(values = c(`TRUE` = "firebrick", `FALSE` = "grey40"), name = "Bonferroni p < 0.05") +
  facet_grid(cohort_label ~ outcome_label, scales = "free_x") +
  labs(
    title = "Predictors of Heteroplasmic Burden, per Cohort (Macro-haplogroup Model)",
    x = "IRR (log scale)", y = NULL
  ) +
  theme_bw(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

dir.create(file.path(stats_dir, "plots"), showWarnings = FALSE)
ggsave(out_path, p_forest_macro_tidy, width = 13, height = 9, dpi = 300, bg = "white")

message(sprintf("Saved: %s", normalizePath(out_path)))
