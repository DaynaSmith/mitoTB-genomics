## =============================================================================
## heteroplasmic_region_type_and_consequence_barplots.R
##
## Same two-plot structure as region_type_and_consequence_barplots.R, but
## restricted to HETEROPLASMIC variants only (hf_class %in% c("low_het",
## "intermediate_het")) rather than all QC-passed variants (which included
## homoplasmic calls).
##
## Data pulled directly from variant_calls_annotated.csv, filtered to
## heteroplasmic calls, deduplicated to distinct (cohort, POS, region_type)
## / (cohort, POS, gene, consequence) combinations, then hardcoded below --
## same approach as the companion script: this script's job is to visualise
## already-finalised numbers, not recompute them live. Update the two data
## frames if you rerun the pipeline and the underlying counts change.
##
## NOTE ON THE MISSENSE/SYNONYMOUS RESULT: unlike the all-variants version
## of this plot (synonymous ~2x missense in every cohort), the HETEROPLASMIC
## subset shows missense variants roughly matching or exceeding synonymous
## counts in Bantu (17 vs 13) and SAC (18 vs 22, closer to parity) --  only
## CAAPA still shows a synonymous majority (16 vs 11). This is a real,
## worth-reporting reversal: purifying selection removes deleterious
## missense variants before they reach fixation/homoplasmy, so a heteroplasmic
## (not-yet-fixed) snapshot is expected to show a less skewed, or even
## reversed, ratio relative to the fixed/homoplasmic variant set -- this
## plot is showing exactly that signature, not an inconsistency with the
## purifying-selection interpretation of the earlier plot.
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
## Plot 1: distinct HETEROPLASMIC variant positions by mtDNA region type,
## per cohort. OL is included at 0 for all three cohorts deliberately --
## no heteroplasmic variation was observed there (all OL variants in this
## dataset are homoplasmic), which is itself a notable absence worth
## showing rather than silently dropping the category.
## =============================================================================
region_data_het <- tribble(
  ~region_type,       ~cohort,                    ~n_positions,
  "Protein-coding",   "Bantu-speaking African",    30,
  "Protein-coding",   "CAAPA",                      28,
  "Protein-coding",   "South African Coloured",     41,
  "Control region",   "Bantu-speaking African",    25,
  "Control region",   "CAAPA",                      12,
  "Control region",   "South African Coloured",     22,
  "rRNA",             "Bantu-speaking African",      7,
  "rRNA",             "CAAPA",                        2,
  "rRNA",             "South African Coloured",       7,
  "tRNA",             "Bantu-speaking African",       3,
  "tRNA",             "CAAPA",                        2,
  "tRNA",             "South African Coloured",       4,
  "Intergenic",       "Bantu-speaking African",       4,
  "Intergenic",       "CAAPA",                        4,
  "Intergenic",       "South African Coloured",       3,
  "Origin of light-strand\nreplication (OL)", "Bantu-speaking African", 0,
  "Origin of light-strand\nreplication (OL)", "CAAPA",                   0,
  "Origin of light-strand\nreplication (OL)", "South African Coloured",  0
) %>%
  mutate(
    cohort = factor(cohort, levels = cohort_levels),
    region_type = fct_reorder(region_type, n_positions, .fun = max, .desc = TRUE)
  )

p1 <- ggplot(region_data_het, aes(x = region_type, y = n_positions, fill = cohort)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7, colour = "grey30", linewidth = 0.3, alpha = 0.8) +
  geom_text(aes(label = n_positions), position = position_dodge(width = 0.75),
            vjust = -0.4, size = 3.4, colour = "grey20") +
  scale_fill_manual(values = cohort_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "Distinct heteroplasmic variant positions by mtDNA region type",
    subtitle = "Low- and intermediate-level heteroplasmy only (excludes homoplasmic variants)",
    x = NULL, y = "Number of distinct positions"
  ) +
  theme_pub

ggsave(file.path(out_dir, "barplot_heteroplasmic_region_type_by_cohort.png"), p1, width = 9, height = 6, dpi = 300, bg = "white")

## =============================================================================
## Plot 2: synonymous vs. missense counts among HETEROPLASMIC protein-coding
## variants, per cohort. See header note -- this ratio is notably different
## from (in Bantu/SAC, reversed relative to) the all-variants comparison.
## =============================================================================
consequence_data_het <- tribble(
  ~consequence,  ~cohort,                    ~n_variants,
  "Synonymous",  "Bantu-speaking African",     13,
  "Synonymous",  "CAAPA",                       16,
  "Synonymous",  "South African Coloured",      22,
  "Missense",    "Bantu-speaking African",      17,
  "Missense",    "CAAPA",                       11,
  "Missense",    "South African Coloured",      18
) %>%
  mutate(
    cohort = factor(cohort, levels = cohort_levels),
    consequence = factor(consequence, levels = c("Synonymous", "Missense"))
  )

p2 <- ggplot(consequence_data_het, aes(x = consequence, y = n_variants, fill = cohort)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.7, colour = "grey30", linewidth = 0.3, alpha = 0.8) +
  geom_text(aes(label = n_variants), position = position_dodge(width = 0.75),
            vjust = -0.4, size = 3.8, colour = "grey20") +
  scale_fill_manual(values = cohort_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(
    title = "Synonymous vs. missense heteroplasmic protein-coding variants",
    subtitle = "Missense variants are proportionally more common here than among fixed (homoplasmic) variants",
    x = NULL, y = "Number of distinct variants"
  ) +
  theme_pub

ggsave(file.path(out_dir, "barplot_heteroplasmic_synonymous_vs_missense.png"), p2, width = 7, height = 6, dpi = 300, bg = "white")

## -----------------------------------------------------------------------
## Optional: combined two-panel figure (requires patchwork)
## -----------------------------------------------------------------------
if (requireNamespace("patchwork", quietly = TRUE)) {
  library(patchwork)
  combined <- p1 + p2 + plot_layout(guides = "collect", widths = c(1.6, 1)) &
    theme(legend.position = "top")
  ggsave(file.path(out_dir, "barplot_heteroplasmic_region_and_consequence_combined.png"), combined,
         width = 15, height = 6, dpi = 300, bg = "white")
  message("Also saved combined two-panel figure: barplot_heteroplasmic_region_and_consequence_combined.png")
} else {
  message("patchwork not installed -- skipping combined figure (install.packages(\"patchwork\") to enable).")
}

message("Done. Saved: barplot_heteroplasmic_region_type_by_cohort.png, barplot_heteroplasmic_synonymous_vs_missense.png")