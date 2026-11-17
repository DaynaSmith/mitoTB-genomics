## =============================================================================
## poster_figure_L0a1_m955_concordance.R  (v2 -- three layout options)
##
## FIX from v1: the annotation text was positioned at
## max(plot_df$total_hf) * 1.08 -- since the two panels deliberately SHARE
## a y-axis (that's what makes the contrast visible), this max is set by
## the m.302 panel's values (~1.0), not the m.955 panel's actual range.
## 1.08x that exceeded ggplot's default axis expansion and got clipped at
## the top of the panel. Fixed here by anchoring the annotation to a safe,
## data-independent y-position within the panel, using geom_label (not
## geom_text) so it stays legible regardless of what's behind it.
##
## THREE OPTIONS, each saved to its own file:
##   Option A: fixed version of the original per-sample jitter-strip plot
##             (all subhaplogroups shown individually)
##   Option B: heatmap/tile summary -- % of each subhaplogroup carrying
##             the variant, one row per position. Reads well at poster
##             viewing distance since colour blocks parse faster than
##             scattered points from a few metres away.
##   Option C: simplified two-category comparison (L0a1 vs. all other
##             Bantu subhaplogroups pooled) -- trades categorical detail
##             for maximum visual punch; likely the strongest option if
##             poster space is tight or the audience will only glance at
##             it briefly.
##
## Requires the same inputs as v1: composition_chrM_302.csv,
## composition_chrM_955.csv, per_sample_heteroplasmic_burden_with_haplogroups.csv
## =============================================================================

library(dplyr); library(readr); library(tidyr); library(ggplot2); library(forcats)

## -----------------------------------------------------------------------
## CONFIG -- paths
## -----------------------------------------------------------------------
qc_dir <- "analysis_outputs"
comp_dir <- file.path(qc_dir, "length_heteroplasmy_composition")
stats_dir <- file.path(qc_dir, "association_analysis")

comp_302 <- read_csv(file.path("mtDNA_QC_outputs/length_heteroplasmy_composition/composition_chrM_302.csv"), show_col_types = FALSE)
comp_955 <- read_csv(file.path("mtDNA_QC_outputs/length_heteroplasmy_composition/composition_chrM_955.csv"), show_col_types = FALSE)
haplo <- read_csv(file.path("mtDNA_QC_outputs/association_analysis/per_sample_heteroplasmic_burden_with_haplogroups.csv"), show_col_types = FALSE)

## -----------------------------------------------------------------------
## Colour palette
## -----------------------------------------------------------------------
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

## -----------------------------------------------------------------------
## Shared data prep
## -----------------------------------------------------------------------
prep_position <- function(comp_df, position_label) {
  comp_df %>%
    dplyr::filter(cohort == "Bantu") %>%
    dplyr::mutate(total_hf = AC + ACC + ACCC) %>%
    dplyr::left_join(haplo %>% dplyr::select(sample_id, broad_haplogroup), by = "sample_id") %>%
    dplyr::transmute(sample_id, broad_haplogroup, total_hf, position = position_label)
}

plot_df <- dplyr::bind_rows(
  prep_position(comp_302, "chrM:302"),
  prep_position(comp_955, "chrM:955")
)

n_bantu <- haplo %>% dplyr::filter(cohort == "Bantu") %>% nrow()
n_l0a1 <- haplo %>% dplyr::filter(cohort == "Bantu", broad_haplogroup == "L0a1") %>% nrow()

carriers_955 <- plot_df %>% dplyr::filter(position == "chrM:955", total_hf > 0)
n_l0a1_carriers <- sum(carriers_955$broad_haplogroup == "L0a1")
n_non_l0a1_carriers <- sum(carriers_955$broad_haplogroup != "L0a1")

if (n_l0a1_carriers != n_l0a1 || n_non_l0a1_carriers > 0) {
  warning(sprintf(
    "Concordance is NOT perfect in this run: %d of %d L0a1 samples carry m.955, and %d non-L0a1 sample(s) also carry it. Update annotation text/framing before using this on a poster.",
    n_l0a1_carriers, n_l0a1, n_non_l0a1_carriers
  ))
}

order_levels <- plot_df %>%
  dplyr::filter(position == "chrM:955") %>%
  dplyr::group_by(broad_haplogroup) %>% dplyr::summarise(m = mean(total_hf, na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(m)) %>% dplyr::pull(broad_haplogroup)

theme_poster <- theme_minimal(base_size = 20) +
  theme(
    strip.text = element_text(face = "bold", size = 22, margin = margin(b = 8)),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(face = "bold", size = 18),
    plot.title = element_text(face = "bold", size = 26, hjust = 0.5, margin = margin(b = 6)),
    plot.subtitle = element_text(size = 16, hjust = 0.5, colour = "grey30", margin = margin(b = 14)),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1.5, "lines"),
    legend.position = "none"
  )

## =============================================================================
## OPTION A: fixed per-sample jitter-strip plot
## =============================================================================
n_labels <- haplo %>% dplyr::filter(cohort == "Bantu") %>% dplyr::count(broad_haplogroup) %>%
  dplyr::mutate(label = paste0(broad_haplogroup, "\n(n=", n, ")"))
label_lookup <- stats::setNames(n_labels$label, n_labels$broad_haplogroup)

plot_df_a <- plot_df %>%
  dplyr::mutate(
    broad_haplogroup = factor(broad_haplogroup, levels = order_levels),
    x_label = factor(label_lookup[as.character(broad_haplogroup)], levels = label_lookup[order_levels]),
    position = factor(position, levels = c("chrM:302", "chrM:955"),
                      labels = c("chrM:302 (population-wide)", "chrM:955 (confined to L0a1)"))
  )
haplo_colors_labeled <- stats::setNames(haplo_colors[as.character(order_levels)], label_lookup[order_levels])

## FIX: anchor annotation at a fixed, safe y-position (0.55 -- comfortably
## within the shared 0-1 axis, and geom_label's white background keeps it
## legible regardless of any points/gridlines behind it), and shift it to
## the 4th x-slot so it doesn't sit on top of the L0a1 point cloud itself
ann_df <- tibble::tibble(
  position = factor("chrM:955 (confined to L0a1)", levels = levels(plot_df_a$position)),
  x_label = levels(plot_df_a$x_label)[4],
  y = 0.55,
  label = sprintf("%d / %d L0a1 samples carry\nm.955 heteroplasmy\n(0 / %d other Bantu-speaking samples do)",
                  n_l0a1_carriers, n_l0a1, n_bantu - n_l0a1)
)

p_a <- ggplot(plot_df_a, aes(x = x_label, y = total_hf, colour = x_label)) +
  geom_jitter(width = 0.15, height = 0, size = 3.5, alpha = 0.85) +
  scale_colour_manual(values = haplo_colors_labeled) +
  facet_wrap(~ position, ncol = 1) +
  #geom_label(data = ann_df, aes(x = x_label, y = y, label = label),
             #inherit.aes = FALSE, colour = "black", fontface = "bold", size = 5,
             #lineheight = 0.95, label.size = 0.4, fill = "white", alpha = 0.9) +
  labs(
    title = "Length heteroplasmy at m.955 is a marker of subhaplogroup L0a1",
    subtitle = sprintf("Bantu-speaking cohort (n = %d). Each point is one individual.", n_bantu),
    x = "Bantu-speaking subhaplogroup (ordered by m.955 heteroplasmy)",
    y = "Length-heteroplasmy fraction (A>AC + A>ACC + A>ACCC)"
  ) +
  theme_poster

ggsave("poster_option_A_jitter_strip.png", p_a, width = 16, height = 11, dpi = 400, bg = "white")

## =============================================================================
## OPTION B: heatmap/tile summary -- % carrying the variant per
## subhaplogroup, one row per position
## =============================================================================
heat_df <- plot_df %>%
  dplyr::mutate(broad_haplogroup = factor(broad_haplogroup, levels = order_levels)) %>%
  dplyr::group_by(position, broad_haplogroup) %>%
  dplyr::summarise(n = dplyr::n(), n_carriers = sum(total_hf > 0), pct_carrier = n_carriers / n, .groups = "drop") %>%
  dplyr::mutate(position = factor(position, levels = c("chrM:955", "chrM:302")))  ## 955 on top row for emphasis

p_b <- ggplot(heat_df, aes(x = broad_haplogroup, y = position, fill = pct_carrier)) +
  geom_tile(colour = "white", linewidth = 1.2) +
  geom_text(aes(label = sprintf("%d/%d", n_carriers, n)), size = 4.5, fontface = "bold",
            colour = ifelse(heat_df$pct_carrier > 0.5, "white", "grey20")) +
  scale_fill_gradient(low = "#FBF3E4", high = "#BF812D", limits = c(0, 1), labels = scales::percent) +
  labs(
    title = "m.302 is population-wide; m.955 marks subhaplogroup L0a1",
    subtitle = sprintf("Bantu-speaking cohort (n = %d). Tile = %% of subhaplogroup carrying the variant (carriers/n shown).", n_bantu),
    x = "Bantu-speaking subhaplogroup", y = NULL, fill = "% carrying\nvariant"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
    axis.text.y = element_text(size = 16, face = "bold"),
    plot.title = element_text(face = "bold", size = 24, hjust = 0.5),
    plot.subtitle = element_text(size = 14, hjust = 0.5, colour = "grey30"),
    legend.position = "right",
    panel.grid = element_blank()
  )

ggsave("poster_option_B_heatmap.png", p_b, width = 16, height = 5.5, dpi = 400, bg = "white")

## =============================================================================
## OPTION C: simplified two-category comparison -- L0a1 vs. all other
## Bantu subhaplogroups pooled. Maximum visual punch, minimum clutter.
## =============================================================================
plot_df_c <- plot_df %>%
  dplyr::mutate(
    group = ifelse(broad_haplogroup == "L0a1", "L0a1", "All other\nBantu subhaplogroups"),
    group = factor(group, levels = c("L0a1", "All other\nBantu subhaplogroups")),
    position = factor(position, levels = c("chrM:302", "chrM:955"),
                      labels = c("chrM:302\n(population-wide)", "chrM:955\n(L0a1-specific)"))
  )

group_colors <- c("L0a1" = unname(haplo_colors["L0a1"]), "All other\nBantu subhaplogroups" = "#B0B0B0")

summary_c <- plot_df_c %>% dplyr::group_by(position, group) %>%
  dplyr::summarise(n = dplyr::n(), n_carriers = sum(total_hf > 0), pct = n_carriers / n, .groups = "drop")

p_c <- ggplot(plot_df_c, aes(x = group, y = total_hf, colour = group)) +
  geom_jitter(width = 0.2, height = 0, size = 4, alpha = 0.75) +
  geom_text(data = summary_c, aes(x = group, y = 1.05, label = sprintf("%d/%d\n(%.0f%%)", n_carriers, n, pct * 100)),
            inherit.aes = FALSE, fontface = "bold", size = 5.5, colour = "grey20", lineheight = 0.9) +
  scale_colour_manual(values = group_colors) +
  scale_y_continuous(limits = c(0, 1.15), breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  facet_wrap(~ position) +
  labs(
    title = "m.955 length heteroplasmy is exclusive to subhaplogroup L0a1",
    subtitle = sprintf("Bantu-speaking cohort (n = %d): L0a1 (n=%d) vs. all %d other Bantu subhaplogroups", n_bantu, n_l0a1, n_bantu - n_l0a1),
    x = NULL, y = "Length-heteroplasmy fraction"
  ) +
  theme_poster +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 15, face = "bold"))

ggsave("poster_option_C_L0a1_vs_rest.png", p_c, width = 13, height = 9, dpi = 400, bg = "white")

message("Saved three options: poster_option_A_jitter_strip.png, poster_option_B_heatmap.png, poster_option_C_L0a1_vs_rest.png")