## =============================================================================
## region_type_and_consequence_barplots.R
##
## Two grouped bar charts summarising the functional annotation results:
##   1. Distinct variant positions by mtDNA region type, per cohort
##   2. Synonymous vs. missense substitution counts (protein-coding), per cohort
##
## Data are hardcoded directly from the reported summary counts (not reread
## from pipeline output) since this script's only job is to visualise those
## already-finalised numbers. If you rerun the pipeline and the numbers
## change, update the two data frames below to match.
## =============================================================================

library(dplyr); library(tidyr); library(ggplot2); library(forcats)

out_dir <- "."  ## change to write plots elsewhere, e.g. "mtDNA_QC_outputs/plots"

## -----------------------------------------------------------------------
## Cohort colour palette (soft/natural tones), consistent with the rest of
## the pipeline's figures
## -----------------------------------------------------------------------
cohort_colors <- c(
  "Bantu-speaking African" = "#4F8A3F",
  "South African Coloured" = "#D77A2D",
  "CAAPA" = "#3B6FB6"
)
cohort_levels <- names(cohort_colors)

theme_pub <- theme_bw(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", size = 17, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, colour = "grey30"),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.title = element_blank(),
    legend.position = "top",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

## =============================================================================
## Plot 1: distinct variant positions by mtDNA region type, per cohort
## =============================================================================
region_data <- tribble(
  ~region_type,       ~cohort,                    ~n_positions,
  "Protein-coding",   "Bantu-speaking African",    488,
  "Protein-coding",   "CAAPA",                      528,
  "Protein-coding",   "South African Coloured",     494,
  "Control region",   "Bantu-speaking African",    141,
  "Control region",   "CAAPA",                      146,
  "Control region",   "South African Coloured",     134,
  "rRNA",             "Bantu-speaking African",     71,
  "rRNA",             "CAAPA",                       68,
  "rRNA",             "South African Coloured",      64,
  "tRNA",             "Bantu-speaking African",     45,
  "tRNA",             "CAAPA",                       43,
  "tRNA",             "South African Coloured",      41,
  "Intergenic",       "Bantu-speaking African",     12,
  "Intergenic",       "CAAPA",                        7,
  "Intergenic",       "South African Coloured",      10,
  "Origin of light-strand\nreplication (OL)", "Bantu-speaking African", 1,
  "Origin of light-strand\nreplication (OL)", "CAAPA",                   2,
  "Origin of light-strand\nreplication (OL)", "South African Coloured",  1
) %>%
  mutate(
    cohort = factor(cohort, levels = cohort_levels),
    region_type = fct_reorder(region_type, n_positions, .fun = max, .desc = TRUE)
  )

p1 <- ggplot(region_data, aes(x = region_type, y = n_positions, fill = cohort)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7, colour = "grey30", linewidth = 0.3, alpha=0.8) +
  geom_text(aes(label = n_positions), position = position_dodge(width = 0.75),
            vjust = -0.4, size = 3.4, colour = "grey20") +
  scale_fill_manual(values = cohort_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "Distinct variant positions by mtDNA region type",
    x = NULL, y = "Number of distinct positions"
  ) +
  theme_pub

ggsave(file.path(out_dir, "barplot_region_type_by_cohort.png"), p1, width = 9, height = 6, dpi = 300, bg = "white")

## =============================================================================
## Plot 2: synonymous vs. missense counts (protein-coding), per cohort
## =============================================================================
consequence_data <- tribble(
  ~consequence,  ~cohort,                    ~n_variants,
  "Synonymous",  "Bantu-speaking African",    333,
  "Synonymous",  "CAAPA",                      381,
  "Synonymous",  "South African Coloured",     355,
  "Missense",    "Bantu-speaking African",     158,
  "Missense",    "CAAPA",                      152,
  "Missense",    "South African Coloured",     143
) %>%
  mutate(
    cohort = factor(cohort, levels = cohort_levels),
    consequence = factor(consequence, levels = c("Synonymous", "Missense"))
  )

p2 <- ggplot(consequence_data, aes(x = consequence, y = n_variants, fill = cohort)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7, colour = "grey30", linewidth = 0.3, alpha=0.8) +
  geom_text(aes(label = n_variants), position = position_dodge(width = 0.75),
            vjust = -0.4, size = 3.8, colour = "grey20") +
  scale_fill_manual(values = cohort_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "Synonymous vs. missense protein-coding variants",
    subtitle = "Synonymous substitutions occur roughly twice as often as missense across all cohorts",
    x = NULL, y = "Number of distinct variants"
  ) +
  theme_pub

ggsave(file.path(out_dir, "barplot_synonymous_vs_missense.png"), p2, width = 7, height = 6, dpi = 300, bg = "white")

## -----------------------------------------------------------------------
## Optional: combined two-panel figure (requires patchwork)
## -----------------------------------------------------------------------

combined <- p1 + p2
combined
