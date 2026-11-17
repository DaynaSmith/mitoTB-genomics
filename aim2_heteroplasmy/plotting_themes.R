# ============================================================
# PLOTTING THEME
# ============================================================

theme_thesis <- theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(
      size = 10,
      hjust = 0
    ),
    axis.title = element_text(size = 10),
    axis.text = element_text(size = 9),
    axis.line = element_line(linewidth = 0.4),
    axis.ticks = element_line(linewidth = 0.4),
    legend.position = "none",
    plot.margin = margin(5, 8, 5, 5)
  )

# ============================================================
# DISSERTATION COLOUR PALETTES
# ============================================================

# Cohorts
cohort_colors <- c(
  "LTBI"    = "#D8A48F",
  "ResisTB" = "#E6A15A",
  "CAAPA"   = "#6B9AC4"
)

# Sex
sex_colors <- c(
  "Female" = "#D8A48F",
  "Male"   = "#8FA8B3"
)

# Macrohaplogroups
macro_colors <- c(
  "L0" = "#D9A441",
  "L1" = "#4D9221",
  "L2" = "#66C2A4",
  "L3" = "#35978F",
  "L4" = "#80CDC1",
  "L5" = "#5AB4AC",
  "M"  = "#8073AC",
  "N"  = "#C2A5CF",
  "R"  = "#D6604D"
)

# Fine-scale haplogroups
haplo_colors <- c(
  # L0a
  "L0a'g" = "#F6D89A",
  "L0a"   = "#DFC27D",
  "L0a1"  = "#BF812D",
  "L0a2"  = "#E6AB52",
  
  # L0d
  "L0d1" = "#E9A03B",
  "L0d2" = "#F5C46B",
  "L0d3" = "#FEE090",
  
  # Other L0
  "L0f" = "#BFA26A",
  "L0k" = "#AEBF7A",
  
  # L1
  "L1b" = "#4D9221",
  "L1c" = "#2B7A0B",
  
  # L2
  "L2a" = "#80CFA9",
  "L2b" = "#66C2A4",
  "L2c" = "#31A354",
  
  # L3
  "L3"  = "#35978F",
  "L3d" = "#2A9D8F",
  "L3e" = "#5AB4AC",
  
  # non-L
  "M" = "#8073AC",
  "N" = "#C2A5CF",
  "R" = "#D6604D",
  
  "Other" = "#BDBDBD"
)

