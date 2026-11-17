# ============================================================
# HETEROPLASMY BURDEN FIGURES
# ResisTB = Bantu
# LTBI    = SAC
# CAAPA   = CAAPA
# ============================================================

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(scales)

# ------------------------------------------------------------
# 1. Check heteroplasmy classes
# ------------------------------------------------------------

table(
  variant_calls_annotated_dedup_for_burden$hf_class,
  useNA = "ifany"
)

# ------------------------------------------------------------
# 2. Create per-sample heteroplasmy burden
# ------------------------------------------------------------

# This assumes hf_class contains labels including the words
# "low" and "intermediate".
# Using grepl() makes the code robust to labels such as
# "Low-level", "low", "Intermediate-level", etc.

het_counts <- variant_calls_annotated_dedup_for_burden %>%
  filter(
    sample_pass_qc == TRUE,
    call_pass_qc == TRUE,
    !is.na(hf_class)
  ) %>%
  mutate(
    hf_group = case_when(
      grepl("low", hf_class, ignore.case = TRUE) ~ "low",
      grepl("intermediate", hf_class, ignore.case = TRUE) ~ "intermediate",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(hf_group)) %>%
  count(
    cohort,
    sample_id,
    hf_group,
    name = "n"
  ) %>%
  pivot_wider(
    names_from = hf_group,
    values_from = n,
    values_fill = 0
  )


# ------------------------------------------------------------
# 3. Start from QC-passed metadata
#    This preserves people with ZERO heteroplasmies
# ------------------------------------------------------------

plot_data <- sample_metadata_QC_passed %>%
  filter(sample_pass_qc == TRUE) %>%
  select(
    sample_id,
    cohort,
    age,
    sex,
    haplogroup,
    broad_haplogroup,
    macrohaplogroup,
    mean_depth,
    mtcn
  ) %>%
  left_join(
    het_counts,
    by = c("sample_id", "cohort")
  ) %>%
  mutate(
    low = replace_na(low, 0),
    intermediate = replace_na(intermediate, 0),
    
    total_heteroplasmy = low + intermediate,
    
    cohort_label = dplyr::recode(
      cohort,
      "Bantu" = "ResisTB",
      "SAC"   = "LTBI",
      "CAAPA" = "CAAPA"
    ),
    
    sex = factor(sex),
    
    macrohaplogroup = factor(macrohaplogroup),
    broad_haplogroup = factor(broad_haplogroup)
  )


# Quick sanity checks
plot_data %>%
  count(cohort_label)

plot_data %>%
  group_by(cohort_label) %>%
  summarise(
    n = n(),
    total_hets = sum(total_heteroplasmy),
    low_hets = sum(low),
    intermediate_hets = sum(intermediate),
    mean_burden = mean(total_heteroplasmy),
    median_burden = median(total_heteroplasmy)
  )

