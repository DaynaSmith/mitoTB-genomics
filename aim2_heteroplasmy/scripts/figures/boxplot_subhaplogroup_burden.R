## =============================================================================
## boxplot_subhaplogroup_burden_tidy.R
##
## Rebuilds the three per-cohort subhaplogroup burden boxplots
## (boxplot_subhaplogroup_burden_Bantu/CAAPA/SAC.png) with tidied labels:
##   - title: cohort code -> full name (Bantu -> Bantu-speaking African,
##     SAC -> South African Coloured; CAAPA unchanged), Title Case
##   - facet strip labels: n_low_het/n_intermediate_het/n_total_het -> full
##     readable phrases, no underscores
##   - both axis titles now explicit ("Subhaplogroup" / "Number of
##     Heteroplasmic Variants") -- the originals left the categorical axis
##     completely unlabelled (x = NULL); added here for a fully
##     self-explanatory chart
##
## Same underlying data, colours (haplo_colors_sub), boxplot/jitter layers,
## and per-cohort ordering (subhaplogroups ordered by median total burden)
## as the original -- purely a labelling pass, no statistics or ordering
## logic changed beyond the ordering already present in cohort_data_prepared.rds.
##
## Reads: mtDNA_QC_outputs/association_analysis/cohort_data_prepared.rds
## (written by 03_heteroplasmy_haplogroup_association.R)
## =============================================================================

library(dplyr); library(tidyr); library(stringr); library(ggplot2); library(forcats)

stats_dir <- "mtDNA_QC_outputs/association_analysis"
plot_dir <- file.path(stats_dir, "plots")
dir.create(plot_dir, showWarnings = FALSE)

cohort_data <- readRDS(file.path(stats_dir, "cohort_data_prepared.rds"))

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
outcomes <- names(outcome_full)  ## matches script 03's outcome ordering

## -----------------------------------------------------------------------
## Colour palette (same as script 03's haplo_colors_sub)
## -----------------------------------------------------------------------
haplo_colors <- c(
  "L0a1" = "#DFC27D", "L0a2" = "#BF812D",
  "L0d1" = "#E9A03B", "L0d2" = "#F5C46B", "L0d3" = "#FEE090",
  "L0f1" = "#BFA26A",
  "L1b1" = "#4D9221", "L1c1" = "#2B7A0B", "L1c2" = "#2E7D32", "L1c3" = "#1B5E20",
  "L5a2" = "#18786F", "L5b1" = "#36BBA7",
  "L2a1" = "#35978F", "L2a5" = "#80CDC1", "L2b" = "#1C9099", "L2b2" = "#7BC8C3",
  "L2b3" = "#9EDBD7", "L2c" = "#4EB3A5", "L2c2" = "#36BBA7", "L2c3" = "#14BDA7", "L2e1" = "#109482",
  "L3b1" = "#3BB8DB", "L3d1" = "#42D3F2", "L3d2" = "#21BCFF", "L3d3" = "#34A6F4", "L3d6" = "#2984D1",
  "L3e1" = "#1C69A8", "L3e2" = "#51A2FF", "L3e3" = "#2B7FFF", "L3f1" = "#155DFC",
  "L2'3'4'6" = "#615FFF",
  "L4b1" = "#8E51FF", "L4b2" = "#7F22FE",
  "M" = "#EC407A", "N" = "#E67597", "R" = "#C51B7D"
)
haplo_colors_sub <- c(haplo_colors, "Other" = "#BDBDBD")

## -----------------------------------------------------------------------
## Rebuild each cohort's plot
## -----------------------------------------------------------------------
purrr::iwalk(cohort_data, function(d, cohort_name) {
  plot_df <- d %>%
    dplyr::mutate(broad_haplogroup_grp = forcats::fct_reorder(broad_haplogroup_grp, n_total_het, .fun = median)) %>%
    tidyr::pivot_longer(c(n_low_het, n_intermediate_het, n_total_het), names_to = "outcome", values_to = "n") %>%
    dplyr::mutate(outcome_label = factor(outcome_full[outcome], levels = outcome_full[outcomes]))
  
  p <- ggplot(plot_df, aes(broad_haplogroup_grp, n, fill = broad_haplogroup_grp)) +
    geom_boxplot(outlier.shape = NA, alpha = .7) +
    geom_jitter(width = .15, size = .7, alpha = .4) +
    scale_fill_manual(values = haplo_colors_sub) +
    facet_wrap(~ outcome_label, scales = "free_x", ncol = 3) +
    coord_flip() +
    labs(
      title = sprintf("%s: Heteroplasmic Burden by Subhaplogroup", cohort_full[cohort_name]),
      x = "Subhaplogroup", y = "Number of Heteroplasmic Variants"
    ) +
    theme_bw(base_size = 12) +
    theme(
      strip.text = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      legend.position = "none"
    )
  
  out_path <- file.path(plot_dir, sprintf("boxplot_subhaplogroup_burden_%s_tidy.png", cohort_name))
  ggsave(out_path, p, width = 11, height = max(4, 0.5 * dplyr::n_distinct(d$broad_haplogroup_grp)), dpi = 300, bg = "white")
  message(sprintf("Saved: %s", normalizePath(out_path)))
})
