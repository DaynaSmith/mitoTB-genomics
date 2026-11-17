# ============================================================
# Figure 4.15
# Multi-allelic length heteroplasmy at m.955 among L0a1
# individuals in the ResisTB cohort
# ============================================================

library(tidyverse)
library(patchwork)
library(scales)
library(forcats)

# ------------------------------------------------------------
# 1. FILE PATHS
# ------------------------------------------------------------

# Tract-composition output
composition_file <- "length_heteroplasmy_composition/composition_chrM_955.csv"

# QC-passed call-level variant file.
# This is the file containing columns such as:
# POS, REF, ALT, sample_id, AF, cohort, call_pass_qc, hf_class
variant_file <- "length_heteroplasmy_composition/het_calls__Bantu-speakingL0a1_individuals.csv"


# ------------------------------------------------------------
# 2. LOAD DATA
# ------------------------------------------------------------

composition <- read_csv(
  composition_file,
  show_col_types = FALSE
)

variants <- read_csv(
  variant_file,
  show_col_types = FALSE
)


# ------------------------------------------------------------
# 3. RESTRICT TO ResisTB / BANTU-SPEAKING COHORT
# ------------------------------------------------------------

composition_bantu <- composition %>%
  filter(cohort == "Bantu")

variants_bantu <- variants %>%
  filter(
    cohort == "Bantu",
    call_pass_qc == TRUE
  )


# ------------------------------------------------------------
# 4. IDENTIFY THE 15 m.955-POSITIVE INDIVIDUALS
# ------------------------------------------------------------

# In the ResisTB cohort, these correspond to the 15 L0a1
# individuals identified in the position-level analysis.

l0a1_ids <- composition_bantu %>%
  mutate(
    nonref_fraction = 1 - reference
  ) %>%
  filter(nonref_fraction > 0) %>%
  pull(sample_id)

cat(
  "Number of m.955-positive/L0a1 individuals:",
  length(l0a1_ids),
  "\n"
)

print(l0a1_ids)


# ------------------------------------------------------------
# 5. EXTRACT m.955 CALLS
# ------------------------------------------------------------

m955_calls <- variants_bantu %>%
  filter(
    sample_id %in% l0a1_ids,
    POS == 955
  )


# ------------------------------------------------------------
# 6. CLASSIFY ALLELES BY LENGTH
# ------------------------------------------------------------

# REF at m.955 is A.
#
# ALT states:
# AC       = +1 C
# ACC      = +2 C
# ACCC     = +3 C
# ACCCC    = +4 C
# ACCCCC   = +5 C
# ACCCCCCC = +6 C
#
# A>C is retained separately because it is a substitution
# representation rather than an insertion-length state.

m955_calls <- m955_calls %>%
  mutate(
    allele = case_when(
      ALT == "AC"      ~ "+1 C",
      ALT == "ACC"     ~ "+2 C",
      ALT == "ACCC"    ~ "+3 C",
      ALT == "ACCCC"   ~ "+4 C",
      ALT == "ACCCCC"  ~ "+5 C",
      ALT == "ACCCCCC" ~ "+6 C",
      ALT == "C"       ~ "A>C",
      TRUE             ~ paste0("Other: ", ALT)
    )
  )


# ------------------------------------------------------------
# 7. SET ALLELE ORDER
# ------------------------------------------------------------

allele_order <- c(
  "+1 C",
  "+2 C",
  "+3 C",
  "+4 C",
  "+5 C",
  "+6 C",
  "A>C"
)

m955_calls <- m955_calls %>%
  mutate(
    allele = factor(
      allele,
      levels = allele_order
    )
  )


# ------------------------------------------------------------
# 8. CREATE PER-SAMPLE SUMMARY FOR ORDERING
# ------------------------------------------------------------

# Order individuals according to:
# 1. total m.955 alternative allele fraction
# 2. longest detected insertion
#
# This produces a biologically meaningful left-to-right order
# rather than arbitrary sample-ID order.

allele_length_key <- c(
  "+1 C" = 1,
  "+2 C" = 2,
  "+3 C" = 3,
  "+4 C" = 4,
  "+5 C" = 5,
  "+6 C" = 6,
  "A>C"  = 0
)

sample_order <- m955_calls %>%
  mutate(
    allele_length = unname(
      allele_length_key[as.character(allele)]
    )
  ) %>%
  group_by(sample_id) %>%
  summarise(
    total_alt_fraction = sum(AF, na.rm = TRUE),
    max_length = max(allele_length, na.rm = TRUE),
    n_alleles = n(),
    .groups = "drop"
  ) %>%
  arrange(
    desc(max_length),
    desc(total_alt_fraction)
  ) %>%
  pull(sample_id)


# ------------------------------------------------------------
# 9. APPLY SAMPLE ORDER
# ------------------------------------------------------------

m955_calls <- m955_calls %>%
  mutate(
    sample_id = factor(
      sample_id,
      levels = sample_order
    )
  )

composition_l0a1 <- composition_bantu %>%
  filter(sample_id %in% l0a1_ids) %>%
  mutate(
    sample_id = factor(
      sample_id,
      levels = sample_order
    )
  )


# ------------------------------------------------------------
# 10. SUMMARY COUNTS
# ------------------------------------------------------------

allele_summary <- m955_calls %>%
  group_by(allele) %>%
  summarise(
    n_individuals = n_distinct(sample_id),
    .groups = "drop"
  )

print(allele_summary)

# Expected approximate pattern from your final call-level data:
#
# +1 C : 15 individuals
# +2 C : 14 individuals
# +3 C : 13 individuals
# +4 C :  7 individuals based on retained calls in supplied file
# +5 C :  4 individuals
# +6 C :  2 individuals
# A>C  :  3 individuals
#
# IMPORTANT:
# These counts are calculated directly from the input file,
# so use the printed output rather than hard-coding manuscript
# values if the final QC dataset changes.


# ------------------------------------------------------------
# 11. LABEL ALLELES WITH CARRIER COUNTS
# ------------------------------------------------------------

allele_labels <- allele_summary %>%
  mutate(
    label = paste0(
      as.character(allele),
      "  (n=", n_individuals, ")"
    )
  ) %>%
  dplyr::select(allele, label) %>%
  deframe()


# ------------------------------------------------------------
# 12. PANEL A — BUBBLE PLOT
# ------------------------------------------------------------

p_bubble <- ggplot(
  m955_calls,
  aes(
    x = sample_id,
    y = allele,
    size = AF
  )
) +
  geom_point(
    alpha = 0.78
  ) +
  scale_size_continuous(
    name = "Variant allele\nfraction",
    range = c(2.5, 10),
    limits = c(0.05, max(m955_calls$AF, na.rm = TRUE)),
    breaks = c(0.05, 0.10, 0.20, 0.30)
  ) +
  scale_y_discrete(
    labels = allele_labels,
    drop = FALSE
  ) +
  scale_x_discrete(
    drop = FALSE
  ) +
  labs(
    title = "A   Multiple m.955 length alleles co-occur within L0a1 individuals",
    x = NULL,
    y = "Alternative m.955 allele"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 12,
      hjust = 0
    ),
    
    axis.title.y = element_text(
      face = "bold"
    ),
    
    axis.text.y = element_text(
      size = 9
    ),
    
    axis.text.x = element_text(
      angle = 50,
      hjust = 1,
      vjust = 1,
      size = 8
    ),
    
    legend.position = "right",
    
    legend.title = element_text(
      face = "bold",
      size = 9
    ),
    
    legend.text = element_text(
      size = 8
    ),
    
    plot.margin = margin(
      t = 8,
      r = 8,
      b = 8,
      l = 8
    )
  )


# ------------------------------------------------------------
# 13. PREPARE COMPOSITION DATA FOR PANEL B
# ------------------------------------------------------------

# composition_chrM_955.csv contains:
# reference, AC, ACC, ACCC, other
#
# The "other" category includes the longer and/or less common
# sequence states not separately represented in this table.

composition_long <- composition_l0a1 %>%
  dplyr::select(
    sample_id,
    reference,
    AC,
    ACC,
    ACCC,
    other
  ) %>%
  pivot_longer(
    cols = c(
      reference,
      AC,
      ACC,
      ACCC,
      other
    ),
    names_to = "tract_state",
    values_to = "fraction"
  ) %>%
  mutate(
    tract_state = dplyr::recode(
      tract_state,
      reference = "Reference",
      AC        = "+1 C",
      ACC       = "+2 C",
      ACCC      = "+3 C",
      other     = "Other states"
    ),
    tract_state = factor(
      tract_state,
      levels = c(
        "Reference",
        "+1 C",
        "+2 C",
        "+3 C",
        "Other states"
      )
    )
  )


# ------------------------------------------------------------
# 14. COLOURS FOR TRACT COMPOSITION
# ------------------------------------------------------------

# Natural, muted progression consistent with your thesis palette.

tract_colors <- c(
  "Reference"    = "#D9D9D9",
  "+1 C"         = "#DFC27D",
  "+2 C"         = "#BF812D",
  "+3 C"         = "#80CFA9",
  "Other states" = "#35978F"
)


# ------------------------------------------------------------
# 15. PANEL B — STACKED COMPOSITION
# ------------------------------------------------------------

p_composition <- ggplot(
  composition_long,
  aes(
    x = sample_id,
    y = fraction,
    fill = tract_state
  )
) +
  geom_col(
    width = 0.82
  ) +
  scale_fill_manual(
    values = tract_colors,
    drop = FALSE
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),
    labels = number_format(accuracy = 0.01),
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    title = "B   Within-individual sequence composition of the m.955 tract",
    x = "L0a1 individual",
    y = "Read fraction",
    fill = "Sequence state"
  ) +
  theme_classic(base_size = 11) +
  theme(
    plot.title = element_text(
      face = "bold",
      size = 12,
      hjust = 0
    ),
    
    axis.title = element_text(
      face = "bold"
    ),
    
    axis.text.x = element_text(
      angle = 50,
      hjust = 1,
      vjust = 1,
      size = 8
    ),
    
    legend.position = "right",
    
    legend.title = element_text(
      face = "bold",
      size = 9
    ),
    
    legend.text = element_text(
      size = 8
    ),
    
    plot.margin = margin(
      t = 8,
      r = 8,
      b = 8,
      l = 8
    )
  )


# ------------------------------------------------------------
# 16. COMBINE PANELS
# ------------------------------------------------------------

figure_415 <-
  p_bubble /
  p_composition +
  plot_layout(
    heights = c(1.05, 1)
  ) +
  plot_annotation(
    title =
      "Multi-allelic length heteroplasmy at m.955 among L0a1 individuals",
    
    subtitle =
      "ResisTB cohort; n = 15 L0a1 individuals",
    
    theme = theme(
      plot.title = element_text(
        face = "bold",
        size = 16,
        hjust = 0
      ),
      
      plot.subtitle = element_text(
        size = 11,
        hjust = 0,
        margin = margin(b = 10)
      )
    )
  )


# ------------------------------------------------------------
# 17. DISPLAY
# ------------------------------------------------------------

figure_415


# ------------------------------------------------------------
# 18. EXPORT
# ------------------------------------------------------------

ggsave(
  filename = "Figure_4_15_m955_length_heteroplasmy.pdf",
  plot = figure_415,
  width = 11,
  height = 8.5,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = "Figure_4_15_m955_length_heteroplasmy.png",
  plot = figure_415,
  width = 11,
  height = 8.5,
  units = "in",
  dpi = 600,
  bg = "white"
)


# ------------------------------------------------------------
# 19. OPTIONAL TABLE OF ALLELE FREQUENCIES
# ------------------------------------------------------------

allele_table <- m955_calls %>%
  mutate(
    sample_id = as.character(sample_id)
  ) %>%
  group_by(allele) %>%
  summarise(
    n_individuals = n_distinct(sample_id),
    percent_individuals =
      100 * n_distinct(sample_id) / length(l0a1_ids),
    
    median_VAF =
      median(AF, na.rm = TRUE),
    
    min_VAF =
      min(AF, na.rm = TRUE),
    
    max_VAF =
      max(AF, na.rm = TRUE),
    
    .groups = "drop"
  )

print(allele_table)


# ------------------------------------------------------------
# 20. OPTIONAL: NUMBER OF ALLELES PER INDIVIDUAL
# ------------------------------------------------------------

per_individual <- m955_calls %>%
  group_by(sample_id) %>%
  summarise(
    n_retained_alleles = n(),
    total_alt_fraction = sum(AF, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(n_retained_alleles))

print(per_individual)
