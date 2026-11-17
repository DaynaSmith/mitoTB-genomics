## =============================================================================
## 02_heteroplasmy_diagnostic_plots_by_cohort.R
##
## Reproduces a 6-panel diagnostic figure (heteroplasmy count vs. sex, age,
## mean depth, mtDNA-CN, macro-haplogroup, subhaplogroup), generated
## SEPARATELY for each cohort. Adapted from a reference figure with panels
## for Cohort/Sex/Age/Depth/Ancestry/Haplogroup:
##   - the "Cohort" panel is dropped -- meaningless once each figure is
##     already restricted to one cohort
##   - "Ancestry" (binary) is replaced with "Macrohaplogroup" (multi-level,
##     so tested with Kruskal-Wallis rather than Wilcoxon)
##   - "Haplogroup" is kept as-is, now representing SUBhaplogroup
##     (broad_haplogroup_grp) to sit alongside the new macro-haplogroup panel
##   - mtDNA-CN is added as a new continuous covariate panel (Spearman)
##   - no smoother line on the continuous panels: a loess curve implies a
##     shape the Spearman statistic isn't testing (rank-monotonic, not
##     curve-fit), and risks overfitting noise at n=105-159 per cohort;
##     these panels show the raw points and let the Spearman rho/p speak
##     for themselves
##
## Statistical tests: Wilcoxon rank-sum for sex, Spearman correlation for
## continuous covariates (age, depth, mtDNA-CN), Kruskal-Wallis for
## multi-category comparisons (macro-haplogroup, subhaplogroup).
##
## Requires cohort_data_prepared.rds, written by
## 03_heteroplasmy_haplogroup_association.R -- run that script first.
## =============================================================================

## install.packages("patchwork") if not already available
library(dplyr); library(ggplot2); library(patchwork)

qc_dir <- "mtDNA_QC_outputs"
stats_dir <- file.path(qc_dir, "association_analysis")
plot_dir <- file.path(stats_dir, "plots", "diagnostics")
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

cohort_data <- readRDS(file.path(stats_dir, "cohort_data_prepared.rds"))

## -----------------------------------------------------------------------
## COLOUR PALETTES (user-provided, soft/natural tones)
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
## "Other" is not a real haplogroup and has no entry above -- add a neutral
## grey so the collapsed rare-category bucket doesn't fall back to
## ggplot's default palette and clash with the rest.
haplo_colors_sub <- c(haplo_colors, "Other" = "#BDBDBD")

## Your palette has no entries for the bare trunk codes used at MACRO
## resolution ("L0"/"L1"/.../"L4" -- only subclade-level keys like L0a1
## exist). This is a constructed fallback, not something you specified:
## one representative shade per branch, picked FROM your own subclade
## colours since each branch already occupies a distinct hue family
## (gold/tan = L0, green = L1, teal = L2, blue = L3, purple = L4).
## Override this vector directly if these picks aren't what you'd choose.
macro_colors <- c(
  "L0" = unname(haplo_colors["L0d1"]),
  "L1" = unname(haplo_colors["L1b1"]),
  "L2" = unname(haplo_colors["L2a1"]),
  "L3" = unname(haplo_colors["L3e1"]),
  "L4" = unname(haplo_colors["L4b1"]),
  "L5" = unname(haplo_colors["L5a2"]),
  "M" = unname(haplo_colors["M"]), "N" = unname(haplo_colors["N"]), "R" = unname(haplo_colors["R"]),
  "Other" = "#BDBDBD"
)

cohort_colors <- c(
  CAAPA = "#3B6FB6",
  BantuSpeakingAfrican = "#4F8A3F",
  SouthAfricanColoured = "#D77A2D"
)
## bridge from this script's short internal cohort names to the long-form
## keys used in cohort_colors above
cohort_name_map <- c(Bantu = "BantuSpeakingAfrican", SAC = "SouthAfricanColoured", CAAPA = "CAAPA")

## no sex-specific palette was provided -- two soft, muted, natural tones
## chosen to sit comfortably alongside the palettes above; swap these for
## an established sex palette if you have one
sex_colors <- c(Female = "#D9A5A0", Male = "#8FA6B2")

## -----------------------------------------------------------------------
## CONFIG -- which outcome(s) to plot. Default is total heteroplasmic
## burden only, matching the reference figure's single "Number of
## Heteroplasmies" outcome; add "n_low_het"/"n_intermediate_het" here to
## also generate figures for those.
## -----------------------------------------------------------------------
outcomes_to_plot <- c("n_total_het")
outcome_labels <- c(n_low_het = "N low-level heteroplasmies",
                    n_intermediate_het = "N intermediate-level heteroplasmies",
                    n_total_het = "N heteroplasmies")

fmt_p <- function(p) if (is.na(p)) "NA" else if (p < 0.001) formatC(p, format = "e", digits = 2) else sprintf("%.4g", p)

## -----------------------------------------------------------------------
## Panel builder
## -----------------------------------------------------------------------
make_diagnostic_panels <- function(df, outcome_col, cohort_name) {
  y_lab <- outcome_labels[[outcome_col]]
  theme_diag <- theme_bw(base_size = 10)
  cohort_accent <- unname(cohort_colors[cohort_name_map[cohort_name]])
  
  ## Sex (Wilcoxon rank-sum)
  sex_test <- tryCatch(wilcox.test(df[[outcome_col]] ~ df$sex), error = function(e) NULL)
  p_sex <- ggplot(df, aes(sex, .data[[outcome_col]], fill = sex)) +
    geom_boxplot(outlier.shape = NA, alpha = .7) +
    geom_jitter(width = .15, size = .8, alpha = .4) +
    scale_fill_manual(values = sex_colors) +
    labs(title = if (!is.null(sex_test)) sprintf("Wilcoxon rank-sum p = %s", fmt_p(sex_test$p.value)) else "Wilcoxon: insufficient groups",
         x = NULL, y = y_lab) + theme_diag + theme(legend.position = "none")
  
  ## Age (Spearman) -- points tinted with this cohort's own colour, no smoother
  age_test <- suppressWarnings(cor.test(df$age, df[[outcome_col]], method = "spearman"))
  p_age <- ggplot(df, aes(age, .data[[outcome_col]])) +
    geom_jitter(height = .1, alpha = .5, colour = cohort_accent) +
    labs(title = sprintf("Spearman rho = %.3g, p = %s", unname(age_test$estimate), fmt_p(age_test$p.value)),
         x = "Age (years)", y = y_lab) + theme_diag
  
  ## Mean read depth (Spearman)
  depth_test <- suppressWarnings(cor.test(df$mean_depth, df[[outcome_col]], method = "spearman"))
  p_depth <- ggplot(df, aes(mean_depth, .data[[outcome_col]])) +
    geom_jitter(height = .1, alpha = .5, colour = cohort_accent) +
    labs(title = sprintf("Spearman rho = %.3g, p = %s", unname(depth_test$estimate), fmt_p(depth_test$p.value)),
         x = "Mean read depth", y = y_lab) + theme_diag
  
  ## mtDNA-CN (Spearman)
  cn_test <- suppressWarnings(cor.test(df$mtcn, df[[outcome_col]], method = "spearman"))
  p_cn <- ggplot(df, aes(mtcn, .data[[outcome_col]])) +
    geom_jitter(height = .1, alpha = .5, colour = cohort_accent) +
    labs(title = sprintf("Spearman rho = %.3g, p = %s", unname(cn_test$estimate), fmt_p(cn_test$p.value)),
         x = "mtDNA copy number", y = y_lab) + theme_diag
  
  ## Macrohaplogroup (Kruskal-Wallis)
  macro_test <- tryCatch(kruskal.test(df[[outcome_col]], df$macro_haplogroup_grp), error = function(e) NULL)
  p_macro <- ggplot(df, aes(macro_haplogroup_grp, .data[[outcome_col]], fill = macro_haplogroup_grp)) +
    geom_boxplot(outlier.shape = NA, alpha = .7) + geom_jitter(width = .15, size = .6, alpha = .4) +
    scale_fill_manual(values = macro_colors) +
    labs(title = if (!is.null(macro_test)) sprintf("Kruskal-Wallis p = %s", fmt_p(macro_test$p.value)) else "Kruskal-Wallis: insufficient groups",
         x = "Macrohaplogroup", y = y_lab) + theme_diag +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
  
  ## Subhaplogroup (Kruskal-Wallis)
  sub_test <- tryCatch(kruskal.test(df[[outcome_col]], df$broad_haplogroup_grp), error = function(e) NULL)
  p_sub <- ggplot(df, aes(broad_haplogroup_grp, .data[[outcome_col]], fill = broad_haplogroup_grp)) +
    geom_boxplot(outlier.shape = NA, alpha = .7) + geom_jitter(width = .15, size = .6, alpha = .4) +
    scale_fill_manual(values = haplo_colors_sub) +
    labs(title = if (!is.null(sub_test)) sprintf("Kruskal-Wallis p = %s", fmt_p(sub_test$p.value)) else "Kruskal-Wallis: insufficient groups",
         x = "Subhaplogroup", y = y_lab) + theme_diag +
    theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1))
  
  list(sex = p_sex, age = p_age, depth = p_depth, mtcn = p_cn, macro = p_macro, sub = p_sub)
}

## -----------------------------------------------------------------------
## Generate one combined figure per cohort x outcome
## -----------------------------------------------------------------------
for (cohort_name in names(cohort_data)) {
  df <- cohort_data[[cohort_name]]
  for (outcome_col in outcomes_to_plot) {
    panels <- make_diagnostic_panels(df, outcome_col, cohort_name)
    
    combined <- (panels$sex | panels$age) /
      (panels$depth | panels$mtcn) /
      (panels$macro | panels$sub) +
      patchwork::plot_annotation(title = sprintf("%s: %s", cohort_name, outcome_labels[[outcome_col]]))
    
    out_path <- file.path(plot_dir, sprintf("diagnostics_%s_%s.png", cohort_name, outcome_col))
    ggsave(out_path, combined, width = 10, height = 12, dpi = 300)
    message(sprintf("Saved %s", out_path))
  }
}

message(sprintf("\nDone. Outputs in: %s", normalizePath(plot_dir)))
