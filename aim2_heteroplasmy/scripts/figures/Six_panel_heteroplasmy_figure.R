# ============================================================
# FUNCTION: SIX-PANEL HETEROPLASMY FIGURE
# ============================================================

make_burden_figure <- function(
    data,
    cohort_name,
    outcome = "total_heteroplasmy") {
  
  d <- data %>%
    filter(cohort_label == cohort_name)
  
  # ----------------------------------------------------------
  # A. SEX
  # ----------------------------------------------------------
  
  p_sex <- ggplot(
    d,
    aes(
      x = sex,
      y = .data[[outcome]],
      fill = sex
    )
  ) +
    geom_boxplot(
      width = 0.65,
      outlier.shape = NA,
      alpha = 0.70
    ) +
    geom_jitter(
      aes(colour = sex),
      width = 0.15,
      height = 0,
      size = 1.2,
      alpha = 0.65
    ) +
    scale_fill_manual(
      values = sex_colors,
      drop = FALSE
    ) +
    scale_colour_manual(
      values = sex_colors,
      drop = FALSE
    ) +
    labs(
      title = get_wilcox_label(d, outcome),
      x = NULL,
      y = "N heteroplasmies"
    ) +
    theme_thesis
  
  
  # ----------------------------------------------------------
  # B. AGE
  # ----------------------------------------------------------
  
  p_age <- ggplot(
    d,
    aes(
      x = age,
      y = .data[[outcome]]
    )
  ) +
    geom_point(
      size = 1.5,
      alpha = 0.7,
      colour = "#5F5F5F"
    ) +
    labs(
      title = get_spearman_label(
        d,
        "age",
        outcome
      ),
      x = "Age (years)",
      y = "N heteroplasmies"
    ) +
    theme_thesis
  
  
  # ----------------------------------------------------------
  # C. MEAN READ DEPTH
  # ----------------------------------------------------------
  
  p_depth <- ggplot(
    d,
    aes(
      x = mean_depth,
      y = .data[[outcome]]
    )
  ) +
    geom_point(
      size = 1.5,
      alpha = 0.7,
      colour = "#5F5F5F"
    ) +
    labs(
      title = get_spearman_label(
        d,
        "mean_depth",
        outcome
      ),
      x = "Mean read depth",
      y = "N heteroplasmies"
    ) +
    theme_thesis
  
  
  # ----------------------------------------------------------
  # D. mtDNA COPY NUMBER
  # ----------------------------------------------------------
  
  p_mtcn <- ggplot(
    d,
    aes(
      x = mtcn,
      y = .data[[outcome]]
    )
  ) +
    geom_point(
      size = 1.5,
      alpha = 0.7,
      colour = "#5F5F5F"
    ) +
    labs(
      title = get_spearman_label(
        d,
        "mtcn",
        outcome
      ),
      x = "mtDNA copy number",
      y = "N heteroplasmies"
    ) +
    theme_thesis
  
  
  # ----------------------------------------------------------
  # E. MACROHAPLOGROUP
  # ----------------------------------------------------------
  
  p_macro <- ggplot(
    d %>% filter(!is.na(macrohaplogroup)),
    aes(
      x = macrohaplogroup,
      y = .data[[outcome]],
      fill = macrohaplogroup
    )
  ) +
    geom_boxplot(
      width = 0.7,
      outlier.shape = NA,
      alpha = 0.75
    ) +
    geom_jitter(
      aes(colour = macrohaplogroup),
      width = 0.15,
      height = 0,
      size = 1.1,
      alpha = 0.65
    ) +
    scale_fill_manual(
      values = macro_colors,
      drop = FALSE
    ) +
    scale_colour_manual(
      values = macro_colors,
      drop = FALSE
    ) +
    labs(
      title = get_kw_label(
        d,
        "macrohaplogroup",
        outcome
      ),
      x = "Macrohaplogroup",
      y = "N heteroplasmies"
    ) +
    theme_thesis
  
  
  # ----------------------------------------------------------
  # F. SUBHAPLOGROUP / BROAD HAPLOGROUP
  # ----------------------------------------------------------
  
  p_sub <- ggplot(
    d %>% filter(!is.na(broad_haplogroup)),
    aes(
      x = broad_haplogroup,
      y = .data[[outcome]],
      fill = broad_haplogroup
    )
  ) +
    geom_boxplot(
      width = 0.7,
      outlier.shape = NA,
      alpha = 0.75
    ) +
    geom_jitter(
      aes(colour = broad_haplogroup),
      width = 0.15,
      height = 0,
      size = 1.1,
      alpha = 0.65
    ) +
    scale_fill_manual(
      values = haplo_colors,
      drop = FALSE,
      na.value = "#BDBDBD"
    ) +
    scale_colour_manual(
      values = haplo_colors,
      drop = FALSE,
      na.value = "#BDBDBD"
    ) +
    labs(
      title = get_kw_label(
        d,
        "broad_haplogroup",
        outcome
      ),
      x = "Sub-haplogroup",
      y = "N heteroplasmies"
    ) +
    theme_thesis +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      )
    )
  
  # ----------------------------------------------------------
  # COMBINE
  # ----------------------------------------------------------
  
  final_plot <-
    (p_sex | p_age) /
    (p_depth | p_mtcn) /
    (p_macro | p_sub) +
    plot_annotation(
      title = cohort_name,
      tag_levels = "A",
      theme = theme(
        plot.title = element_text(
          face = "bold",
          size = 14
        )
      )
    )
  
  return(final_plot)
}

fig_LTBI <- make_burden_figure(
  plot_data,
  cohort_name = "LTBI"
)

fig_ResisTB <- make_burden_figure(
  plot_data,
  cohort_name = "ResisTB"
)

fig_CAAPA <- make_burden_figure(
  plot_data,
  cohort_name = "CAAPA"
)

fig_LTBI
fig_ResisTB
fig_CAAPA

ggsave(
  "LTBI_heteroplasmy_burden.png",
  fig_LTBI,
  width = 8,
  height = 10,
  dpi = 600
)

ggsave(
  "ResisTB_heteroplasmy_burden.png",
  fig_ResisTB,
  width = 8,
  height = 10,
  dpi = 600
)

ggsave(
  "CAAPA_heteroplasmy_burden.png",
  fig_CAAPA,
  width = 8,
  height = 10,
  dpi = 600
)

ggsave(
  "LTBI_heteroplasmy_burden.pdf",
  fig_LTBI,
  width = 8,
  height = 10
)

ggsave(
  "ResisTB_heteroplasmy_burden.pdf",
  fig_ResisTB,
  width = 8,
  height = 10
)

ggsave(
  "CAAPA_heteroplasmy_burden.pdf",
  fig_CAAPA,
  width = 8,
  height = 10
)

