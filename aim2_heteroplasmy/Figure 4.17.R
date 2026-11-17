# ============================================================
# Figure 4.17
# Read-level distribution of m.955 poly-C tract lengths
# among L0a1b1a1 individuals in the ResisTB cohort
# ============================================================

library(tidyverse)

# -----------------------------
# 1. Read input data
# -----------------------------

dat <- read.delim(
  "m955_consensus_vs_FASTQ.tsv",
  header = TRUE,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

# -----------------------------
# 2. Convert C-length fractions
#    from wide to long format
# -----------------------------

plot_dat <- dat %>%
  select(
    sample,
    modal_C_length,
    mean_C_length,
    consensus_C_length,
    starts_with("C")
  ) %>%
  select(
    sample,
    modal_C_length,
    mean_C_length,
    consensus_C_length,
    matches("^C[0-9]+_fraction$")
  ) %>%
  pivot_longer(
    cols = matches("^C[0-9]+_fraction$"),
    names_to = "tract_length",
    values_to = "fraction"
  ) %>%
  mutate(
    # Extract numeric tract length from names such as C8_fraction
    C_length = as.numeric(
      str_extract(tract_length, "[0-9]+")
    ),
    percent = fraction * 100
  )

# -----------------------------
# 3. Order samples by modal tract
#    length, then mean tract length
# -----------------------------

sample_order <- dat %>%
  arrange(modal_C_length, mean_C_length) %>%
  pull(sample)

plot_dat <- plot_dat %>%
  mutate(
    sample = factor(sample, levels = sample_order)
  )

consensus_dat <- dat %>%
  mutate(
    sample = factor(sample, levels = sample_order)
  )

# -----------------------------
# 4. Keep states represented by
#    >=5% of pure-C reads
# -----------------------------

plot_dat_5pct <- plot_dat %>%
  filter(percent >= 5)

# -----------------------------
# 5. Build Figure 4.17
# -----------------------------

p <- ggplot() +
  
  # Bubble layer:
  # size = percentage of pure-C reads supporting each length
  geom_point(
    data = plot_dat_5pct,
    aes(
      x = sample,
      y = C_length,
      size = percent
    ),
    alpha = 0.75
  ) +
  
  # Consensus tract length:
  # open diamond so it remains visible over bubbles
  geom_point(
    data = consensus_dat,
    aes(
      x = sample,
      y = consensus_C_length
    ),
    shape = 23,
    size = 3.4,
    fill = "white",
    stroke = 0.8
  ) +
  
  # Integer C-tract lengths only
  scale_y_continuous(
    breaks = 8:20,
    limits = c(7.5, 20.5)
  ) +
  
  # Bubble-size legend
  scale_size_continuous(
    name = "Reads supporting\ntract length (%)",
    range = c(2.5, 10),
    breaks = c(5, 10, 25, 50, 75)
  ) +
  
  labs(
    x = "Individual",
    y = "Uninterrupted poly-C tract length",
    caption = "\u25C7 Consensus mitochondrial tract length"
  ) +
  
  theme_classic(base_size = 12) +
  
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      vjust = 1
    ),
    axis.title = element_text(
      face = "bold"
    ),
    legend.title = element_text(
      face = "bold"
    ),
    plot.caption = element_text(
      hjust = 0,
      size = 10
    )
  )

p

# -----------------------------
# 6. Save high-resolution files
# -----------------------------

ggsave(
  "Figure_4.17_m955_read_level_bubbleplot.pdf",
  p,
  width = 9,
  height = 5.8
)

ggsave(
  "Figure_4.17_m955_read_level_bubbleplot.png",
  p,
  width = 9,
  height = 5.8,
  dpi = 600
)
