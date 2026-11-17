## =============================================================================
## 04_heteroplasmy_haplogroup_association.R
##
## Heteroplasmic burden association analysis, run SEPARATELY per cohort
## (Bantu, SAC, CAAPA) rather than pooled with a cohort covariate.
##
## WHY STRATIFY RATHER THAN POOL-AND-ADJUST: mean sequencing depth and
## cohort were confounded in this dataset (median depth 1849x/1184x/845x
## for Bantu/CAAPA/SAC; non-overlapping 5th-95th percentile depth ranges
## for some cohort pairs -- see section 4). Earlier pipeline versions tried
## to handle this within a single pooled model, first with a spline-depth
## term, then with linear depth once the spline was found to add nothing
## once other data-quality issues were fixed. Analysing each cohort
## separately sidesteps the confound structurally rather than statistically
## -- there is no cross-cohort depth comparison happening inside any one
## cohort's model, so a plain linear mean_depth covariate is sufficient
## within each cohort without needing a cohort term or a flexible
## depth-response shape.
##
## Each cohort gets FULL models at both macro-haplogroup and subhaplogroup
## resolution (6 negative binomial models per cohort: 2 resolutions x 3
## outcomes; 18 total), each with mean_depth as a plain linear covariate:
##   outcome ~ age + sex + mtcn + mean_depth + haplogroup_term
##
## Other consolidated decisions, unchanged from prior versions:
##   - burden excludes haplogroup-defining homoplasmic markers
##   - macro-haplogroup and subhaplogroup categories are now collapsed
##     ("Other" if rare) and releveled to the modal category WITHIN each
##     cohort separately, not pooled across cohorts as before -- consistent
##     with treating each cohort as its own independent analysis population
##   - Wald CIs computed manually from glm.nb's coefficient table (not
##     broom's profile-likelihood confint(), which fails on sparse
##     factor-level combinations)
##   - Bonferroni correction applied within each cohort's own planned
##     comparison family (age, sex, mtcn, mean_depth, haplogroup x 3
##     outcomes = 15 tests per cohort; NOT further multiplied across the 3
##     cohorts, each treated as a self-contained analysis)
##   - position-diversity check counts DISTINCT SAMPLES per position, not
##     raw call rows (see prior version's notes on the m.955 counting bug)
##   - burden reflects the finalized artifact-filter decision from script
##     01: m.302, m.567, m.955, m.16182 are genuine poly-C
##     length-heteroplasmy (Chen et al. 2023 Nature), fully reinstated,
##     each length-allele counted separately
##
## Run after 03_variant_annotation_mitomap.R.
## =============================================================================

library(dplyr); library(tidyr); library(stringr); library(purrr); library(readr)
library(tibble); library(forcats); library(ggplot2); library(MASS); library(broom)

qc_dir <- "mtDNA_QC_outputs"
stats_dir <- file.path(qc_dir, "association_analysis")
dir.create(stats_dir, showWarnings = FALSE); dir.create(file.path(stats_dir, "plots"), showWarnings = FALSE)

variant_calls <- readRDS(file.path(qc_dir, "variant_calls_annotated_dedup_for_burden.rds"))
sample_meta   <- readRDS(file.path(qc_dir, "sample_metadata_QC_passed.rds"))

## -----------------------------------------------------------------------
## 1. CONFIG
## -----------------------------------------------------------------------
min_macro_haplogroup_n <- 10   ## applied WITHIN each cohort now
min_subhaplogroup_n <- 8       ## applied WITHIN each cohort
n_planned_tests <- 5 * 3       ## age, sex, mtcn, mean_depth, haplogroup x 3 outcomes -- per cohort
cohort_meta_paths <- list(Bantu = "../Outputs/Bantu_metadata.csv", SAC = "../Outputs/SAC_metadata.csv", CAAPA = "../Outputs/CAAPA_metadata.csv")
sample_id_col <- "sample_id"
broad_haplogroup_col <- "broad_haplogroup"
cohort_names <- names(cohort_meta_paths)

outcomes <- c("n_low_het", "n_intermediate_het", "n_total_het")
hg_terms <- c(macro = "macro_haplogroup_grp", sub = "broad_haplogroup_grp")

## -----------------------------------------------------------------------
## 2. Per-individual heteroplasmic burden (excludes haplogroup markers)
## -----------------------------------------------------------------------
het_calls <- variant_calls %>%
  dplyr::filter(hf_class %in% c("low_het", "intermediate_het")) %>%
  dplyr::filter(is.na(haplogroup_defining) | !haplogroup_defining)

burden <- het_calls %>%
  dplyr::group_by(sample_id, cohort, hf_class) %>% dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = hf_class, values_from = n, values_fill = 0) %>%
  dplyr::rename_with(~ paste0("n_", .x), c(dplyr::any_of(c("low_het", "intermediate_het"))))

burden <- sample_meta %>% dplyr::select(sample_id, cohort) %>%
  dplyr::left_join(burden, by = c("sample_id", "cohort")) %>%
  dplyr::mutate(n_low_het = tidyr::replace_na(n_low_het, 0), n_intermediate_het = tidyr::replace_na(n_intermediate_het, 0),
                n_total_het = n_low_het + n_intermediate_het)

## -----------------------------------------------------------------------
## 3. Merge covariates + derive macro-haplogroup + merge subhaplogroup
## (collapsing/releveling deferred to per-cohort preparation, section 5)
## -----------------------------------------------------------------------
derive_macro_haplogroup <- function(hg) {
  hg <- as.character(hg)
  dplyr::case_when(is.na(hg) ~ NA_character_,
                   stringr::str_starts(hg, "L") ~ stringr::str_extract(hg, "^L[0-6]?"),
                   TRUE ~ stringr::str_extract(hg, "^[A-Z]{1,2}"))
}

read_broad_haplogroup <- function(path, cohort_name) {
  df <- readr::read_csv(path, show_col_types = FALSE)
  if (!broad_haplogroup_col %in% names(df)) stop(sprintf("[%s] no column '%s' -- columns present: %s",
                                                         cohort_name, broad_haplogroup_col, paste(names(df), collapse = ", ")))
  df %>% dplyr::select(sample_id = !!sample_id_col, broad_haplogroup = !!broad_haplogroup_col) %>% dplyr::mutate(cohort = cohort_name)
}
subhg_all <- purrr::imap(cohort_meta_paths, read_broad_haplogroup) %>% dplyr::bind_rows()

analysis_df_base <- burden %>%
  dplyr::left_join(sample_meta %>% dplyr::select(sample_id, age, sex, haplogroup, mtcn, mean_depth), by = "sample_id") %>%
  dplyr::left_join(subhg_all %>% dplyr::select(sample_id, broad_haplogroup), by = "sample_id") %>%
  dplyr::mutate(macro_haplogroup = derive_macro_haplogroup(haplogroup))

n_missing_subhg <- sum(is.na(analysis_df_base$broad_haplogroup))
if (n_missing_subhg > 0) warning(sprintf("%d samples unmatched on broad_haplogroup after join -- check sample_id consistency.", n_missing_subhg))

readr::write_csv(analysis_df_base, file.path(stats_dir, "per_sample_heteroplasmic_burden_with_haplogroups.csv"))

## -----------------------------------------------------------------------
## 4. Depth / cohort confounding diagnostic -- retained as the DOCUMENTED
## JUSTIFICATION for stratifying by cohort (see header), not as an input
## to choosing a depth functional form within a pooled model.
## -----------------------------------------------------------------------
depth_by_cohort <- analysis_df_base %>% dplyr::group_by(cohort) %>%
  dplyr::summarise(median_depth = median(mean_depth), p5 = quantile(mean_depth, .05),
                   p95 = quantile(mean_depth, .95), n = dplyr::n(), .groups = "drop")
message("\n=== Depth distribution by cohort (motivates per-cohort stratification) ==="); print(depth_by_cohort)

cohort_pairs <- combn(cohort_names, 2, simplify = FALSE)
overlap_check <- purrr::map_dfr(cohort_pairs, function(pair) {
  d1 <- depth_by_cohort %>% dplyr::filter(cohort == pair[1]); d2 <- depth_by_cohort %>% dplyr::filter(cohort == pair[2])
  lo <- max(d1$p5, d2$p5); hi <- min(d1$p95, d2$p95)
  tibble(cohort_A = pair[1], cohort_B = pair[2], overlapping_range = hi > lo, overlap_lo = lo, overlap_hi = hi)
})
readr::write_csv(overlap_check, file.path(stats_dir, "depth_cohort_overlap_check.csv"))
print(overlap_check)

## -----------------------------------------------------------------------
## 5. Per-cohort data preparation: collapse rare haplogroup categories and
## relevel to the modal category WITHIN each cohort
## -----------------------------------------------------------------------
relevel_to_mode <- function(x) {
  ref <- names(sort(table(x), decreasing = TRUE))[1]
  forcats::fct_relevel(factor(x), ref)
}

prepare_cohort_data <- function(df, this_cohort) {
  d <- df %>% dplyr::filter(cohort == this_cohort)
  
  macro_counts <- d %>% dplyr::count(macro_haplogroup, name = "n_samples")
  rare_macro <- macro_counts$macro_haplogroup[macro_counts$n_samples < min_macro_haplogroup_n]
  d <- d %>% dplyr::mutate(macro_haplogroup_grp = dplyr::if_else(macro_haplogroup %in% rare_macro | is.na(macro_haplogroup), "Other", macro_haplogroup))
  
  sub_counts <- d %>% dplyr::count(broad_haplogroup, name = "n_samples")
  rare_sub <- sub_counts$broad_haplogroup[sub_counts$n_samples < min_subhaplogroup_n]
  d <- d %>% dplyr::mutate(broad_haplogroup_grp = dplyr::if_else(broad_haplogroup %in% rare_sub | is.na(broad_haplogroup), "Other", broad_haplogroup))
  
  d %>% dplyr::mutate(
    macro_haplogroup_grp = relevel_to_mode(macro_haplogroup_grp),
    broad_haplogroup_grp = relevel_to_mode(broad_haplogroup_grp),
    sex = forcats::fct_relevel(factor(sex), sort(unique(as.character(sex)))[1])
  )
}

cohort_data <- purrr::map(cohort_names, ~ prepare_cohort_data(analysis_df_base, .x)) %>% purrr::set_names(cohort_names)
saveRDS(cohort_data, file.path(stats_dir, "cohort_data_prepared.rds"))

purrr::iwalk(cohort_data, function(d, cohort_name) {
  message(sprintf("[%s] n = %d. Reference macro-haplogroup = '%s'; reference subhaplogroup = '%s'.",
                  cohort_name, nrow(d), levels(d$macro_haplogroup_grp)[1], levels(d$broad_haplogroup_grp)[1]))
})

## -----------------------------------------------------------------------
## 6. Model fitting -- per cohort, per resolution, per outcome. mean_depth
## is a plain linear covariate; no cohort term (each model is already
## within one cohort).
## -----------------------------------------------------------------------
fit_nb_model <- function(outcome, hg_term, df) {
  f <- as.formula(paste(outcome, "~ age + sex + mtcn + mean_depth +", hg_term))
  model_df <- df %>% tidyr::drop_na(age, sex, mtcn, mean_depth, macro_haplogroup_grp, broad_haplogroup_grp)
  fit <- tryCatch(MASS::glm.nb(f, data = model_df), error = function(e) {
    message(sprintf("glm.nb failed (%s, %s) -- quasi-Poisson fallback.", outcome, hg_term))
    stats::glm(f, data = model_df, family = quasipoisson())
  })
  list(fit = fit, n = nrow(model_df))
}

tidy_nb_safe <- function(fit) {
  s <- summary(fit)$coefficients
  tibble::tibble(
    term = rownames(s), estimate = exp(s[, 1]), std.error = s[, 2],
    statistic = s[, 3], p.value = s[, 4],
    conf.low = exp(s[, 1] - 1.96 * s[, 2]), conf.high = exp(s[, 1] + 1.96 * s[, 2])
  )
}

combo_grid <- tidyr::expand_grid(cohort = cohort_names, outcome = outcomes, hg = names(hg_terms))
all_fits <- purrr::pmap(combo_grid, function(cohort, outcome, hg) fit_nb_model(outcome, hg_terms[[hg]], cohort_data[[cohort]]))
names(all_fits) <- paste(combo_grid$cohort, combo_grid$outcome, combo_grid$hg, sep = "__")
saveRDS(all_fits, file.path(stats_dir, "all_models_by_cohort.rds"))

## -----------------------------------------------------------------------
## 7. Model comparison: macro vs sub-haplogroup resolution, WITHIN each
## cohort (AIC)
## -----------------------------------------------------------------------
model_comparison <- purrr::imap_dfr(all_fits, function(m, key) {
  parts <- stringr::str_split(key, "__")[[1]]
  tibble(cohort = parts[1], outcome = parts[2], haplogroup_resolution = parts[3], n = m$n,
         AIC = tryCatch(AIC(m$fit), error = function(e) NA))
})
message("\n=== Model comparison: macro vs sub-haplogroup resolution, per cohort ===")
print(model_comparison %>% dplyr::arrange(cohort, outcome, AIC), n = Inf)
readr::write_csv(model_comparison, file.path(stats_dir, "model_comparison_by_cohort.csv"))

## -----------------------------------------------------------------------
## 8. Tidy coefficient tables, both resolutions, per cohort
## -----------------------------------------------------------------------
tidy_all <- purrr::imap_dfr(all_fits, function(m, key) {
  parts <- stringr::str_split(key, "__")[[1]]
  tidy_nb_safe(m$fit) %>%
    dplyr::mutate(cohort = parts[1], outcome = parts[2], haplogroup_resolution = parts[3], n_samples = m$n, .before = 1)
}) %>%
  dplyr::filter(!stringr::str_detect(term, "Intercept")) %>%
  dplyr::mutate(p_bonferroni = pmin(p.value * n_planned_tests, 1), sig_bonferroni = p_bonferroni < 0.05) %>%
  dplyr::rename(IRR = estimate, IRR_lower = conf.low, IRR_upper = conf.high)

readr::write_csv(tidy_all, file.path(stats_dir, "negbin_regression_results_by_cohort.csv"))

message("\n=== Macro-haplogroup results, per cohort ===")
tidy_all %>% dplyr::filter(haplogroup_resolution == "macro") %>%
  dplyr::select(cohort, outcome, term, IRR, IRR_lower, IRR_upper, p.value, p_bonferroni, sig_bonferroni) %>%
  dplyr::arrange(cohort, outcome) %>% print(n = Inf)

message("\n=== Subhaplogroup results, per cohort ===")
tidy_all %>% dplyr::filter(haplogroup_resolution == "sub", stringr::str_starts(term, "broad_haplogroup_grp")) %>%
  dplyr::select(cohort, outcome, term, IRR, IRR_lower, IRR_upper, p.value) %>%
  dplyr::arrange(cohort, outcome) %>% print(n = Inf)
message("NOTE: subhaplogroup p-values above are NOT Bonferroni-corrected -- treat as exploratory,")
message("especially for categories with small n. See Kruskal-Wallis results below for a")
message("covariate-free, more robust subhaplogroup test.")

## -----------------------------------------------------------------------
## 9. Kruskal-Wallis tests, subhaplogroup, per cohort (covariate-free, most
## robust given small/uneven subhaplogroup group sizes)
## -----------------------------------------------------------------------
kw_results <- purrr::imap_dfr(cohort_data, function(d, cohort_name) {
  if (dplyr::n_distinct(d$broad_haplogroup_grp) < 2) {
    return(tibble(cohort = cohort_name, outcome = outcomes, kruskal_p = NA_real_, n_groups = 1))
  }
  purrr::map_dfr(outcomes, function(outcome) {
    test <- kruskal.test(d[[outcome]], d$broad_haplogroup_grp)
    tibble(cohort = cohort_name, outcome = outcome, kruskal_p = test$p.value, n_groups = dplyr::n_distinct(d$broad_haplogroup_grp))
  })
})
message("\n=== Kruskal-Wallis: burden ~ subhaplogroup, per cohort ===")
print(kw_results, n = Inf)
readr::write_csv(kw_results, file.path(stats_dir, "kruskal_wallis_subhaplogroup_by_cohort.csv"))

## -----------------------------------------------------------------------
## 10. Position-diversity check (counts distinct samples, not raw call
## rows). Run manually for any haplogroup/subhaplogroup category with a
## notably elevated coefficient above.
## -----------------------------------------------------------------------
check_position_diversity <- function(sample_ids, calls_df = variant_calls) {
  calls_df %>%
    dplyr::filter(sample_id %in% sample_ids, hf_class %in% c("low_het", "intermediate_het")) %>%
    dplyr::group_by(POS) %>%
    dplyr::summarise(n_distinct_samples = dplyr::n_distinct(sample_id), n_calls = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(n_distinct_samples))
}
message("\nReminder: for any haplogroup/subhaplogroup category with a notably elevated coefficient,")
message('run: ids <- cohort_data[["CAAPA"]] %>% filter(broad_haplogroup_grp == "GROUP_NAME") %>% pull(sample_id)')
message("      check_position_diversity(ids)")
message("Check n_distinct_samples, not n_calls.")

## -----------------------------------------------------------------------
## 11. Plots
## -----------------------------------------------------------------------
theme_set(theme_bw(base_size = 11))

## subhaplogroup colour palette (soft/natural tones) for the box-and-whisker
## plot below -- NOT applied to the forest plots (p_forest_macro,
## p_subhg_coefs), which keep their significance-driven colouring.
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

## macro-haplogroup forest plot, faceted by cohort x outcome
p_forest_macro <- tidy_all %>% dplyr::filter(haplogroup_resolution == "macro") %>%
  ggplot(aes(x = IRR, y = term, colour = sig_bonferroni)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_pointrange(aes(xmin = IRR_lower, xmax = IRR_upper)) + scale_x_log10() +
  scale_colour_manual(values = c(`TRUE` = "firebrick", `FALSE` = "grey40"), name = "Bonferroni p<0.05") +
  facet_grid(cohort ~ outcome, scales = "free_x") +
  labs(title = "Predictors of heteroplasmic burden, per cohort (macro-haplogroup model)", x = "IRR (log scale)", y = NULL)
ggsave(file.path(stats_dir, "plots", "irr_forest_macro_by_cohort.png"), p_forest_macro, width = 13, height = 9, dpi = 300)

## subhaplogroup forest plot, faceted by cohort x outcome
p_subhg_coefs <- tidy_all %>% dplyr::filter(haplogroup_resolution == "sub", stringr::str_starts(term, "broad_haplogroup_grp")) %>%
  ggplot(aes(x = IRR, y = term)) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_pointrange(aes(xmin = IRR_lower, xmax = IRR_upper)) + scale_x_log10() +
  facet_grid(cohort ~ outcome, scales = "free") +
  labs(title = "Subhaplogroup coefficients, per cohort (exploratory -- not Bonferroni-corrected)", x = "IRR (log scale)", y = NULL)
ggsave(file.path(stats_dir, "plots", "irr_forest_subhaplogroup_by_cohort.png"), p_subhg_coefs, width = 13, height = 10, dpi = 300)

## box-and-whisker plots of burden by subhaplogroup -- ONE PER COHORT, since
## haplogroup collapsing/releveling ("Other") is now cohort-specific and the
## categories differ by cohort; separate files keep each cohort's own set
## of subhaplogroups clearly readable rather than forcing a shared axis.
purrr::iwalk(cohort_data, function(d, cohort_name) {
  plot_df <- d %>%
    dplyr::mutate(broad_haplogroup_grp = forcats::fct_reorder(broad_haplogroup_grp, n_total_het, .fun = median)) %>%
    tidyr::pivot_longer(c(n_low_het, n_intermediate_het, n_total_het), names_to = "outcome", values_to = "n") %>%
    dplyr::mutate(outcome = factor(outcome, levels = outcomes))
  
  p <- ggplot(plot_df, aes(broad_haplogroup_grp, n, fill = broad_haplogroup_grp)) +
    geom_boxplot(outlier.shape = NA, alpha = .7) +
    geom_jitter(width = .15, size = .7, alpha = .4) +
    scale_fill_manual(values = haplo_colors_sub) +
    facet_wrap(~ outcome, scales = "free_x", ncol = 3) + coord_flip() +
    labs(title = sprintf("%s: heteroplasmic burden by subhaplogroup", cohort_name),
         x = NULL, y = "N heteroplasmic variants") +
    theme(legend.position = "none")
  ggsave(file.path(stats_dir, "plots", sprintf("boxplot_subhaplogroup_burden_%s.png", cohort_name)),
         p, width = 11, height = max(4, 0.4 * dplyr::n_distinct(d$broad_haplogroup_grp)), dpi = 300)
})

message(sprintf("\nDone. Outputs in: %s", normalizePath(stats_dir)))