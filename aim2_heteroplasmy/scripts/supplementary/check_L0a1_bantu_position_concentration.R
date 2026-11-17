## =============================================================================
## check_L0a1_bantu_position_concentration.R
## =============================================================================

library(dplyr); library(readr)

qc_passed_calls <- readRDS("mtDNA_QC_outputs/variant_calls_QC_passed_hgflagged.rds")
sample_meta <- readr::read_csv("mtDNA_QC_outputs/association_analysis/per_sample_heteroplasmic_burden_with_haplogroups.csv", show_col_types = FALSE)

l0a1_ids <- sample_meta %>% dplyr::filter(cohort == "Bantu", broad_haplogroup == "L0a1") %>% dplyr::pull(sample_id)
message(sprintf("N L0a1 Bantu individuals: %d", length(l0a1_ids)))

het_calls <- qc_passed_calls %>%
  dplyr::filter(sample_id %in% l0a1_ids, hf_class %in% c("low_het", "intermediate_het")) %>%
  dplyr::filter(is.na(haplogroup_defining) | !haplogroup_defining)

n_total <- nrow(het_calls)
message(sprintf("Total heteroplasmic calls among L0a1 individuals: %d", n_total))

## -----------------------------------------------------------------------
## Version A: single exact positions only (302, 955) -- UNDERCOUNTS
## -----------------------------------------------------------------------
n_single_pos <- het_calls %>% dplyr::filter(POS %in% c(302, 955)) %>% nrow()
message(sprintf("\n[Version A -- single positions only] At exactly m.302 or m.955: %d (%.0f%%)",
                n_single_pos, 100 * n_single_pos / n_total))

## -----------------------------------------------------------------------
## Version B: full homopolymer tract ranges -- matches manuscript's "95 (94%)"
## -----------------------------------------------------------------------
m302_tract <- c(302, 315)
m955_tract <- c(955, 965)

in_m302_tract <- het_calls %>% dplyr::filter(POS >= m302_tract[1], POS <= m302_tract[2])
in_m955_tract <- het_calls %>% dplyr::filter(POS >= m955_tract[1], POS <= m955_tract[2])
n_tract <- nrow(in_m302_tract) + nrow(in_m955_tract)

message(sprintf("[Version B -- full tract ranges] At m.302 tract (%d-%d) or m.955 tract (%d-%d): %d (%.0f%%)",
                m302_tract[1], m302_tract[2], m955_tract[1], m955_tract[2], n_tract, 100 * n_tract / n_total))
message("  Distinct positions contributing within each tract:")
message(sprintf("    m.302 tract: %s", paste(sort(unique(in_m302_tract$POS)), collapse = ", ")))
message(sprintf("    m.955 tract: %s", paste(sort(unique(in_m955_tract$POS)), collapse = ", ")))

## -----------------------------------------------------------------------
## Per-position penetrance: how many of the 15 L0a1 individuals carry each
## -----------------------------------------------------------------------
n_955_carriers <- het_calls %>% dplyr::filter(POS == 955) %>% dplyr::distinct(sample_id) %>% nrow()
n_302_carriers <- het_calls %>% dplyr::filter(POS == 302) %>% dplyr::distinct(sample_id) %>% nrow()

message(sprintf("\nm.955 present in %d of %d L0a1 individuals", n_955_carriers, length(l0a1_ids)))
message(sprintf("m.302 present in %d of %d L0a1 individuals", n_302_carriers, length(l0a1_ids)))

write.csv(het_calls, "het_calls__Bantu-speakingL0a1_individuals.csv", row.names = FALSE)
