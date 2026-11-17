p48a <- ggplot(
  plot_data,
  aes(
    x = cohort_label,
    y = total_heteroplasmy,
    fill = cohort_label
  )
) +
  geom_boxplot(
    width = 0.65,
    outlier.shape = NA,
    alpha = 0.75
  ) +
  geom_jitter(
    width = 0.15,
    size = 1.3,
    alpha = 0.55
  ) +
  scale_fill_manual(values = c("lightgrey", "lightgrey", "lightgrey")) +
  labs(
    x = NULL,
    y = "Number of heteroplasmies"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "none"
  )

burden_long <- plot_data %>%
  select(
    sample_id,
    cohort_label,
    low,
    intermediate
  ) %>%
  pivot_longer(
    cols = c(low, intermediate),
    names_to = "heteroplasmy_level",
    values_to = "burden"
  ) %>%
  mutate(
    heteroplasmy_level = factor(
      heteroplasmy_level,
      levels = c("low", "intermediate"),
      labels = c("Low-level", "Intermediate-level")
    )
  )

p48b <- ggplot(
  burden_long,
  aes(
    x = cohort_label,
    y = burden,
    fill = heteroplasmy_level
  )
) +
  geom_boxplot(
    position = position_dodge(width = 0.75),
    width = 0.65,
    outlier.shape = NA,
    alpha = 0.75
  ) +
  scale_fill_manual(values = c("#ec8e75", "#82a18f")) +
  labs(
    x = NULL,
    y = "Number of heteroplasmies",
    fill = "Heteroplasmy level"
  ) +
  theme_classic(base_size = 11) +
  theme(
    legend.position = "right"
  )

fig_4_8 <- p48a | p48b

fig_4_8
