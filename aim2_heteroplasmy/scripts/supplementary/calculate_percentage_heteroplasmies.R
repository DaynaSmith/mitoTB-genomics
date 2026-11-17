# ============================================================
# Calculate percentage of individuals heteroplasmic
# at each mitochondrial position within each cohort
# ============================================================

library(dplyr)

# 1. Get number of QC-passed individuals in each cohort
cohort_n <- sample_metadata_QC_passed %>%
  filter(sample_pass_qc == TRUE) %>%
  count(cohort, name = "n_total")

cohort_n

position_frequency <- variant_calls_annotated_dedup_for_burden %>%
  
  # Keep QC-passed heteroplasmic calls
  filter(
    sample_pass_qc == TRUE,
    call_pass_qc == TRUE,
    hf_class %in% c("low_het", "intermediate_het")
  ) %>%
  
  # One individual should only contribute once per position
  distinct(
    cohort,
    sample_id,
    POS
  ) %>%
  
  # Count individuals heteroplasmic at each position
  count(
    cohort,
    POS,
    name = "n_heteroplasmic"
  ) %>%
  
  # Add cohort denominator
  left_join(
    cohort_n,
    by = "cohort"
  ) %>%
  
  # Calculate percentage
  mutate(
    pct_samples = 100 * n_heteroplasmic / n_total,
    
    cohort_label = dplyr::recode(
      cohort,
      "Bantu" = "ResisTB",
      "SAC"   = "LTBI",
      "CAAPA" = "CAAPA"
    ),
    
    cohort_label = factor(
      cohort_label,
      levels = c("LTBI", "ResisTB", "CAAPA")
    )
  )

head(position_frequency)

position_frequency %>%
  arrange(desc(pct_samples)) %>%
  head(20)


# ============================================================
# FIGURE 4.9
# Distribution of recurrent mitochondrial heteroplasmy
# across LTBI, ResisTB and CAAPA
# ============================================================

library(dplyr)
library(ggplot2)
library(ggrepel)

# ------------------------------------------------------------
# 1. Dissertation cohort colours
# ------------------------------------------------------------

cohort_colors <- c(
  "CAAPA"   = "#86ADD1",
  "LTBI"    = "#D9A18F",
  "ResisTB" = "#E6A15A"
)

# ------------------------------------------------------------
# 2. Positions to highlight
# ------------------------------------------------------------

highlight_positions <- c(302, 567, 955, 16182)

# ------------------------------------------------------------
# 3. Prepare plotting dataset
#
# Assumes your position-level dataframe contains:
#   cohort
#   POS
#   pct_samples
#
# If you already have cohort_label, remove the recode step.
# ------------------------------------------------------------

plot_4_9_data <- position_frequency %>%
  mutate(
    cohort_label = dplyr::recode(
      cohort,
      "SAC"   = "LTBI",
      "Bantu" = "ResisTB",
      "CAAPA" = "CAAPA"
    ),
    
    cohort_label = factor(
      cohort_label,
      levels = c("LTBI", "ResisTB", "CAAPA")
    ),
    
    highlight = POS %in% highlight_positions
  )

# ------------------------------------------------------------
# 4. Separate background and highlighted positions
# ------------------------------------------------------------

background_sites <- plot_4_9_data %>%
  filter(!highlight)

highlight_sites <- plot_4_9_data %>%
  filter(highlight)

# Only label highlighted sites occurring in >=5% of individuals
label_sites <- highlight_sites %>%
  filter(pct_samples >= 5)

# ------------------------------------------------------------
# 5. Generate figure
# ------------------------------------------------------------

fig_4_9 <- ggplot() +
  
  # Background heteroplasmic positions
  geom_point(
    data = background_sites,
    aes(
      x = POS,
      y = pct_samples
    ),
    colour = "grey75",
    size = 1.5,
    alpha = 0.55
  ) +
  
  # 5% recurrence threshold
  geom_hline(
    yintercept = 5,
    linetype = "dashed",
    linewidth = 0.45,
    colour = "grey45"
  ) +
  
  # Highlight recurrent positions
  geom_point(
    data = highlight_sites,
    aes(
      x = POS,
      y = pct_samples,
      colour = cohort_label
    ),
    size = 3,
    alpha = 0.95
  ) +
  
  # Labels for important recurrent sites
  geom_text_repel(
    data = label_sites,
    aes(
      x = POS,
      y = pct_samples,
      label = paste0("m.", POS),
      colour = cohort_label
    ),
    size = 3.4,
    fontface = "italic",
    show.legend = FALSE,
    box.padding = 0.4,
    point.padding = 0.3,
    min.segment.length = 0,
    max.overlaps = Inf,
    seed = 42
  ) +
  
  # One panel per cohort
  facet_wrap(
    ~ cohort_label,
    ncol = 1
  ) +
  
  # Dissertation colours
  scale_colour_manual(
    values = cohort_colors
  ) +
  
  # Square-root scale makes low-frequency sites easier to see
  scale_y_sqrt(
    breaks = c(
      0, 1, 2, 5, 10,
      20, 40, 60
    ),
    labels = c(
      "0", "1", "2", "5", "10",
      "20", "40", "60"
    ),
    expand = expansion(
      mult = c(0.02, 0.08)
    )
  ) +
  
  # mtDNA genome coordinates
  scale_x_continuous(
    limits = c(1, 16569),
    breaks = c(
      1,
      2500,
      5000,
      7500,
      10000,
      12500,
      15000,
      16569
    ),
    labels = c(
      "1",
      "2,500",
      "5,000",
      "7,500",
      "10,000",
      "12,500",
      "15,000",
      "16,569"
    ),
    expand = expansion(
      mult = c(0.01, 0.01)
    )
  ) +
  
  labs(
    x = "Mitochondrial genome position (rCRS)",
    y = "Individuals with heteroplasmy (%)"
  ) +
  
  theme_classic(
    base_size = 12
  ) +
  
  theme(
    legend.position = "none",
    
    strip.background = element_blank(),
    
    strip.text = element_text(
      face = "bold",
      size = 12,
      hjust = 0
    ),
    
    axis.title = element_text(
      size = 12
    ),
    
    axis.text = element_text(
      size = 10
    ),
    
    axis.line = element_line(
      linewidth = 0.5
    ),
    
    axis.ticks = element_line(
      linewidth = 0.4
    ),
    
    panel.spacing = unit(
      0.8,
      "lines"
    ),
    
    plot.margin = margin(
      10, 15, 10, 10
    )
  )

# View figure
fig_4_9
