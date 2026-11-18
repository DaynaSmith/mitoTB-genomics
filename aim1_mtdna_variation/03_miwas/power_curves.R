
library(ggplot2)
library(dplyr)
library(scales)

# 1. Load simulation results
power_data <- read.csv("miWAS_Firth_power_curves.csv")

# 2. Prepare labels
power_data <- power_data %>%
  mutate(
    Carrier_group = paste0(carriers, " carriers"),
    Carrier_group = factor(
      Carrier_group,
      levels = paste0(sort(unique(carriers)), " carriers")
    ),
    Power_percent = power * 100
  )

# 3. Plot statistical power curves
power_plot <- ggplot(
  power_data,
  aes(
    x = OR,
    y = Power_percent,
    colour = Carrier_group,
    group = Carrier_group
  )
) +
  geom_hline(
    yintercept = 80,
    linetype = "dashed",
    colour = "grey40",
    linewidth = 0.6
  ) +
  geom_line(linewidth = 1.1) +
  geom_point(size = 2.2) +
  annotate(
    "text",
    x = max(power_data$OR),
    y = 83,
    label = "80% power",
    hjust = 1,
    size = 3.5,
    colour = "grey30"
  ) +
  scale_x_continuous(
    breaks = sort(unique(power_data$OR))
  ) +
  scale_y_continuous(
    breaks = seq(0, 100, 20),
    limits = c(0, 105),
    labels = function(x) paste0(x, "%")
  ) +
  scale_colour_viridis_d(
    option = "plasma",
    end = 0.9
  ) +
  labs(
    title = "Statistical power of the mitochondrial-wide association study",
    subtitle = paste(
      "Firth logistic regression |",
      "n = 815 (412 cases, 403 controls)"
    ),
    x = "Odds ratio (OR)",
    y = "Statistical power",
    colour = "Variant frequency",
    caption = paste(
      "Bonferroni threshold: p = 2.15 × 10⁻⁴;",
      "300 simulations per scenario."
    )
  ) +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13),
    plot.subtitle = element_text(size = 10),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold"),
    plot.caption = element_text(
      size = 9,
      colour = "grey40",
      hjust = 0
    )
  )

print(power_plot)

# 4. Export publication-quality figures
ggsave(
  "miWAS_power_curves.pdf",
  plot = power_plot,
  width = 9,
  height = 6
)

ggsave(
  "miWAS_power_curves.png",
  plot = power_plot,
  width = 9,
  height = 6,
  dpi = 600,
  bg = "white"
)
