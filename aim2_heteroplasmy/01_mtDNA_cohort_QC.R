## =============================================================================
## 01_mtDNA_cohort_QC.R
##
## Full QC pipeline for mtDNA heteroplasmy calls across three WGS cohorts
## (Bantu, SAC, CAAPA).
##
## ARTIFACT-FILTER POSITIONS -- FINAL DECISION (see project notes for full
## chain of evidence):
##   Four positions (302, 567, 955, 16182) are poly-C-tract heteroplasmic
##   INDEL hotspots explicitly named in Chen et al. 2023, Nature ("Genome-wide
##   detection of heteroplasmic indels in the UK Biobank", chrM:302/567/
##   955/16182 in non-coding mtDNA; fewer indels at protein-coding poly-C
##   tracts, consistent with negative selection -- a biological signature no
##   sequencing artifact would produce). Sibling-pair heteroplasmy
##   concordance for individual length-alleles at m.955 (their Fig. 4c,
##   confirmed via user-provided screenshot) shows these are heritable, real
##   variants -- each specific tract length (A>AC, A>ACC, A>ACCC, ...) is
##   analysed as ITS OWN variant in that paper, not collapsed into one
##   presence/absence indicator, so this pipeline does the same: no ladder
##   collapsing.
##   These 4 positions are EXEMPTED from BOTH artifact filters (position-
##   hotspot AND, where applicable, the per-allele homopolymer-recurrence
##   filter -- m.955 sits inside a homopolymer_regions window and was found
##   to be only PARTIALLY reinstated by an earlier version of this script
##   that patched just one of the two filters; fixed here).
##   Four OTHER positions gnomAD's own pipeline calls "artifact_prone_site"
##   (301, 302, 310, 316, 3107, 16182 per Laricchia et al. 2022) overlap this
##   list at 302/16182 (now reinstated per the stronger, mechanism-specific
##   evidence above) but 301, 310, 316, 3107 remain EXCLUDED -- no paper has
##   named these as genuine recurrent heteroplasmy, only as sites gnomAD's
##   SNV pipeline can't reliably call.
##   16179/16183/16189 (same D-loop tract as 16182) remain EXCLUDED --
##   plausible by the same poly-C mechanism but not directly confirmed by
##   either source; treated as an open question, not extended by inference.
##
## Other consolidated fixes:
##   - vcfR2tidy(single_frame = TRUE) VCF parsing (manual FIX/GT joins
##     silently duplicate rows at multi-record positions)
##   - revised homopolymer region boundaries (302 not 303; added 955-965 and
##     8272-8290; widened 16179-16194)
##   - two-tier artifact filter: per-allele recurrence within homopolymer
##     regions, AND per-position recurrence (any allele) genome-wide
##   - haplogroup-defining variant flagging (Haplogrep2) + Ts/Tv diagnostics
## =============================================================================

library(vcfR); library(dplyr); library(tidyr); library(stringr); library(purrr)
library(readr); library(tibble); library(forcats); library(ggplot2)

## -----------------------------------------------------------------------
## 1. CONFIG
## -----------------------------------------------------------------------
cohorts <- list(
  Bantu = list(vcf = "../Outputs/Bantu_mutect2.03.merge.vcf.gz", meta = "../Outputs/Bantu_metadata.csv"),
  SAC   = list(vcf = "../Outputs/SAC_mutect2.03.merge.vcf.gz",   meta = "../Outputs/SAC_metadata.csv"),
  CAAPA = list(vcf = "../Outputs/CAAPA_mutect2.03.merge.vcf",    meta = "../Outputs/CAAPA_metadata.csv")
)
meta_cols <- list(
  sample_id = "sample_id", age = "age", sex = "sex", haplogroup = "haplogroup",
  mtcn = "copy_number", mean_depth = "mean_depth", contamination = "contamination"
)

haplogrep_files <- list(
  Bantu = "../Outputs/bantu_haplogrep_output.txt",
  SAC   = "../Outputs/sac_haplogrep_output.txt",
  CAAPA = "../Outputs/caapa_haplogrep_output.txt"
)
hg_cols <- list(sample_id = "SampleID", found_polys = "Found_Polys")

## sample-level thresholds
min_mean_depth <- 100
max_contamination <- 0.02
cn_outlier_method <- "iqr"
cn_outlier_k <- 3

## call-level thresholds
min_site_depth <- 100
min_alt_reads <- 5
min_hf <- 0.05
homoplasmy_hf <- 0.95
het_low_hf_max <- 0.10
require_pass_filter <- TRUE

## artifact filters
recurrent_het_freq_threshold <- 0.02    ## per-ALLELE recurrence, within homopolymer regions
position_hotspot_freq_threshold <- 0.05 ## per-POSITION recurrence, any allele, cohort-wide

homopolymer_regions <- list(
  c(66, 71), c(302, 315), c(515, 525),
  c(955, 965),      ## 12S rRNA poly-C tract
  c(3106, 3107),
  c(8272, 8290),    ## COII/tRNA-Lys 9bp deletion region
  c(12418, 12425),
  c(16179, 16194)
)

## Positions confirmed via external evidence to be genuine common
## heteroplasmic INDEL hotspots (poly-C tract length heteroplasmy) rather
## than technical artifact, despite high recurrence -- see header notes for
## the full citation chain (Chen et al. 2023 Nature; sibling-concordance
## evidence at m.955). Exempted from BOTH the position-hotspot filter and
## the per-allele homopolymer-recurrence filter. Do not add positions here
## without equivalent direct, positive, mechanism-specific evidence.
confirmed_common_heteroplasmy_positions <- c(302, 567, 955, 16182)

out_dir <- "mtDNA_QC_outputs"
dir.create(out_dir, showWarnings = FALSE)
dir.create(file.path(out_dir, "plots"), showWarnings = FALSE)

## -----------------------------------------------------------------------
## 2. Helper functions
## -----------------------------------------------------------------------
in_homopolymer <- function(pos, regions = homopolymer_regions) {
  purrr::map_lgl(pos, function(p) any(purrr::map_lgl(regions, ~ p >= .x[1] & p <= .x[2])))
}

split_multiallelic <- function(df) {
  has_multi <- stringr::str_detect(df$ALT, ",")
  if (!any(has_multi, na.rm = TRUE)) return(df)
  message(sprintf("Splitting %d multiallelic records...", sum(has_multi, na.rm = TRUE)))
  single <- df[!has_multi | is.na(has_multi), ]
  multi  <- df[has_multi %in% TRUE, ]
  expand_row <- function(row) {
    alts <- stringr::str_split(row$ALT, ",")[[1]]; n <- length(alts)
    out <- row[rep(1, n), ]; out$ALT <- alts
    if ("gt_AF" %in% names(row) && stringr::str_detect(as.character(row$gt_AF), ",")) {
      afs <- stringr::str_split(as.character(row$gt_AF), ",")[[1]]
      out$gt_AF <- suppressWarnings(as.numeric(afs[seq_len(n)]))
    }
    out
  }
  expanded <- multi %>% dplyr::rowwise() %>% dplyr::group_split() %>% purrr::map_dfr(expand_row)
  dplyr::bind_rows(single, expanded)
}

## single_frame = TRUE lets vcfR resolve the FIX/GT join internally at the
## record level -- a manual join on (ChromKey, POS) silently duplicates
## genotype rows at any position with more than one VCF record, inflating
## every downstream count by ~10%. Do not revert to a manual join.
read_cohort_vcf <- function(vcf_path, cohort_name) {
  message(sprintf("[%s] reading VCF: %s", cohort_name, vcf_path))
  vcf <- vcfR::read.vcfR(vcf_path, verbose = FALSE)
  tidy <- vcfR::vcfR2tidy(vcf, single_frame = TRUE, verbose = FALSE, info_only = FALSE)
  long <- tidy$dat %>%
    dplyr::mutate(cohort = cohort_name, variant_key = paste(CHROM, POS, REF, ALT, sep = "_"))
  split_multiallelic(long)
}

classify_hf <- function(af) {
  dplyr::case_when(
    is.na(af) ~ NA_character_,
    af <  min_hf ~ "below_threshold",
    af >= min_hf & af <  het_low_hf_max ~ "low_het",
    af >= het_low_hf_max & af < homoplasmy_hf ~ "intermediate_het",
    af >= homoplasmy_hf ~ "homoplasmic",
    TRUE ~ NA_character_
  )
}

flag_sample_qc <- function(meta_df) {
  meta_df %>%
    dplyr::mutate(
      flag_low_depth = mean_depth < min_mean_depth,
      flag_contam = contamination > max_contamination,
      flag_missing_hg = is.na(haplogroup) | haplogroup == ""
    ) %>%
    dplyr::group_by(cohort) %>%
    dplyr::mutate(
      cn_outlier = {
        if (cn_outlier_method == "iqr") {
          q <- stats::quantile(mtcn, probs = c(.25, .75), na.rm = TRUE); iqr <- diff(q)
          mtcn < (q[1] - cn_outlier_k * iqr) | mtcn > (q[2] + cn_outlier_k * iqr)
        } else {
          m <- mean(mtcn, na.rm = TRUE); s <- stats::sd(mtcn, na.rm = TRUE)
          abs(mtcn - m) > cn_outlier_k * s
        }
      }
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(sample_pass_qc = !flag_low_depth & !flag_contam & !flag_missing_hg)
}

parse_found_polys <- function(x) {
  if (is.na(x) || x == "") {
    return(tibble(pos = numeric(0), allele = character(0), is_indel = logical(0),
                  back_mutation = logical(0), raw = character(0)))
  }
  toks <- str_split(x, "\\s+")[[1]]; toks <- toks[toks != ""]
  purrr::map_dfr(toks, function(tok) {
    back_mutation <- str_detect(tok, "!$")
    tok_clean <- str_remove(tok, "!$")
    is_indel <- str_detect(tok_clean, "\\.") | str_detect(tok_clean, "[dD]$") | str_detect(tok_clean, "-")
    pos_match <- str_extract(tok_clean, "^[0-9]+(\\.[0-9]+)?")
    allele <- str_remove(tok_clean, "^[0-9]+(\\.[0-9]+)?")
    tibble(pos = suppressWarnings(as.numeric(pos_match)), allele = allele,
           is_indel = is_indel, back_mutation = back_mutation, raw = tok)
  })
}

compute_tstv <- function(df) {
  df <- df %>% filter(nchar(REF) == 1, nchar(ALT) == 1)
  purine <- c("A", "G"); pyrimidine <- c("C", "T")
  df <- df %>% mutate(is_transition = (REF %in% purine & ALT %in% purine) | (REF %in% pyrimidine & ALT %in% pyrimidine))
  tibble(n_variants = nrow(df), n_transitions = sum(df$is_transition, na.rm = TRUE),
         n_transversions = sum(!df$is_transition, na.rm = TRUE),
         ts_tv_ratio = sum(df$is_transition, na.rm = TRUE) / sum(!df$is_transition, na.rm = TRUE))
}

## -----------------------------------------------------------------------
## 3. Sample-level QC
## -----------------------------------------------------------------------
read_one_metadata <- function(path, cohort_name) {
  df <- readr::read_csv(path, show_col_types = FALSE)
  missing <- setdiff(unlist(meta_cols), names(df))
  if (length(missing) > 0) stop(sprintf("[%s] metadata missing column(s): %s", cohort_name, paste(missing, collapse = ", ")))
  df %>%
    dplyr::rename(sample_id = !!meta_cols$sample_id, age = !!meta_cols$age, sex = !!meta_cols$sex,
                  haplogroup = !!meta_cols$haplogroup, mtcn = !!meta_cols$mtcn,
                  mean_depth = !!meta_cols$mean_depth, contamination = !!meta_cols$contamination) %>%
    dplyr::mutate(cohort = cohort_name)
}

metadata_all <- purrr::imap(cohorts, ~ read_one_metadata(.x$meta, .y)) %>% dplyr::bind_rows() %>% flag_sample_qc()
readr::write_csv(metadata_all, file.path(out_dir, "sample_QC_summary.csv"))

message("\n=== Sample-level QC summary ===")
metadata_all %>% dplyr::group_by(cohort) %>%
  dplyr::summarise(n_total = dplyr::n(), n_low_depth = sum(flag_low_depth), n_high_contam = sum(flag_contam),
                   n_missing_hg = sum(flag_missing_hg), n_cn_outlier = sum(cn_outlier),
                   n_pass_qc = sum(sample_pass_qc)) %>% print(n = Inf)

## -----------------------------------------------------------------------
## 4. Load VCFs, join sample flags
## -----------------------------------------------------------------------
vcf_long_all <- purrr::imap(cohorts, ~ read_cohort_vcf(.x$vcf, .y)) %>% dplyr::bind_rows() %>%
  dplyr::rename(sample_id = Indiv) %>%
  dplyr::mutate(
    AF = suppressWarnings(as.numeric(gt_AF)),
    site_DP = if ("DP" %in% names(.)) suppressWarnings(as.numeric(dplyr::coalesce(gt_DP, DP))) else suppressWarnings(as.numeric(gt_DP)),
    alt_reads = round(AF * site_DP),
    pos_class = dplyr::if_else(in_homopolymer(POS), "homopolymer", "standard")
  ) %>%
  dplyr::left_join(
    metadata_all %>% dplyr::select(sample_id, cohort, sample_pass_qc, flag_low_depth, flag_contam, cn_outlier),
    by = c("sample_id", "cohort")
  )

## -----------------------------------------------------------------------
## 5. Call-level QC flags + artifact filters
## -----------------------------------------------------------------------
vcf_long_all <- vcf_long_all %>%
  dplyr::mutate(
    flag_filter_fail = if (require_pass_filter) !(FILTER %in% c("PASS", ".", NA)) else FALSE,
    flag_low_site_depth = site_DP < min_site_depth,
    flag_low_alt_reads = alt_reads < min_alt_reads,
    flag_homopolymer = pos_class == "homopolymer",
    hf_class = classify_hf(AF),
    is_confirmed_common_het = POS %in% confirmed_common_heteroplasmy_positions
  )

het_calls <- vcf_long_all %>% dplyr::filter(sample_pass_qc %in% TRUE, hf_class %in% c("low_het", "intermediate_het"))
n_pass_by_cohort <- metadata_all %>% dplyr::filter(sample_pass_qc) %>% dplyr::count(cohort, name = "n_samples_pass")

recurrent_het_flags <- het_calls %>%
  dplyr::group_by(cohort, POS, REF, ALT) %>% dplyr::summarise(n_het = dplyr::n(), .groups = "drop") %>%
  dplyr::left_join(n_pass_by_cohort, by = "cohort") %>%
  dplyr::mutate(het_freq = n_het / n_samples_pass, flag_recurrent_het_raw = het_freq > recurrent_het_freq_threshold)

position_recurrence <- het_calls %>%
  dplyr::group_by(cohort, POS) %>% dplyr::summarise(n_samples_with_call = dplyr::n_distinct(sample_id), .groups = "drop") %>%
  dplyr::left_join(n_pass_by_cohort, by = "cohort") %>%
  dplyr::mutate(pct_samples = n_samples_with_call / n_samples_pass,
                flag_position_hotspot_raw = pct_samples > position_hotspot_freq_threshold,
                is_confirmed_common_het = POS %in% confirmed_common_heteroplasmy_positions,
                flag_position_hotspot = flag_position_hotspot_raw & !is_confirmed_common_het,
                override_applied = flag_position_hotspot_raw & is_confirmed_common_het)

message("\n=== Position-level hotspot audit ===")
message(sprintf("(%d position(s) exempted from BOTH artifact filters via confirmed-common-heteroplasmy override: %s)",
                length(confirmed_common_heteroplasmy_positions),
                paste(confirmed_common_heteroplasmy_positions, collapse = ", ")))
position_recurrence %>% dplyr::filter(flag_position_hotspot_raw) %>%
  dplyr::arrange(cohort, dplyr::desc(pct_samples)) %>%
  dplyr::select(cohort, POS, n_samples_with_call, n_samples_pass, pct_samples, flag_position_hotspot, override_applied) %>%
  print(n = Inf)

readr::write_csv(position_recurrence %>% dplyr::arrange(dplyr::desc(pct_samples)),
                 file.path(out_dir, "position_hotspot_audit.csv"))

vcf_long_all <- vcf_long_all %>%
  dplyr::left_join(recurrent_het_flags %>% dplyr::select(cohort, POS, REF, ALT, het_freq, flag_recurrent_het_raw),
                   by = c("cohort", "POS", "REF", "ALT")) %>%
  dplyr::left_join(position_recurrence %>% dplyr::select(cohort, POS, pct_samples, flag_position_hotspot),
                   by = c("cohort", "POS")) %>%
  dplyr::mutate(flag_recurrent_het_raw = tidyr::replace_na(flag_recurrent_het_raw, FALSE),
                flag_position_hotspot = tidyr::replace_na(flag_position_hotspot, FALSE),
                ## the per-allele filter is ALSO overridden for confirmed positions -- this is the
                ## fix for the earlier bug where m.955 was only partially reinstated because only
                ## flag_position_hotspot was overridden, leaving the per-allele filter to still
                ## exclude most of its length-allele ladder in Bantu
                flag_recurrent_het = flag_recurrent_het_raw & !is_confirmed_common_het)

vcf_long_all <- vcf_long_all %>%
  dplyr::mutate(
    call_pass_qc = sample_pass_qc %in% TRUE &
      !flag_filter_fail & !flag_low_site_depth & !flag_low_alt_reads &
      !(flag_homopolymer & flag_recurrent_het) &
      !flag_position_hotspot &
      hf_class != "below_threshold"
  )

readr::write_csv(vcf_long_all %>% dplyr::select(-dplyr::any_of("gt_GT_alleles")),
                 file.path(out_dir, "variant_call_QC_flags.csv"))

message("\n=== Variant-level QC summary (post sample-QC, PASS calls only) ===")
vcf_long_all %>% dplyr::filter(call_pass_qc) %>% dplyr::group_by(cohort, hf_class) %>%
  dplyr::summarise(n_calls = dplyr::n(), n_variant_positions = dplyr::n_distinct(POS), .groups = "drop") %>%
  print(n = Inf)

message("\n=== Confirmed common heteroplasmy positions: calls passing QC after full override ===")
vcf_long_all %>% dplyr::filter(call_pass_qc, is_confirmed_common_het) %>%
  dplyr::count(cohort, POS, REF, ALT, hf_class, name = "n_calls") %>%
  dplyr::arrange(POS, cohort, dplyr::desc(n_calls)) %>% print(n = Inf)

## -----------------------------------------------------------------------
## 6. Haplogroup-defining variant flagging (Haplogrep2) + Ts/Tv diagnostic
## -----------------------------------------------------------------------
read_haplogrep <- function(path, cohort_name) {
  df <- readr::read_tsv(path, show_col_types = FALSE)
  missing <- setdiff(unlist(hg_cols), names(df))
  if (length(missing) > 0) stop(sprintf("[%s] Haplogrep output missing column(s): %s", cohort_name, paste(missing, collapse = ", ")))
  df %>% dplyr::rename(sample_id = !!hg_cols$sample_id, found_polys = !!hg_cols$found_polys) %>% dplyr::mutate(cohort = cohort_name)
}

haplogrep_all <- purrr::imap(haplogrep_files, read_haplogrep) %>% dplyr::bind_rows()
hg_markers <- haplogrep_all %>%
  dplyr::mutate(parsed = purrr::map(found_polys, parse_found_polys)) %>%
  dplyr::select(sample_id, cohort, parsed) %>% tidyr::unnest(parsed)
hg_markers_snv <- hg_markers %>% dplyr::filter(!is_indel, !back_mutation) %>%
  dplyr::mutate(match_key = paste(sample_id, pos, allele, sep = "_"))
readr::write_csv(hg_markers, file.path(out_dir, "haplogroup_marker_tokens_all.csv"))

qc_passed_calls <- vcf_long_all %>% dplyr::filter(call_pass_qc) %>%
  dplyr::mutate(match_key = paste(sample_id, POS, ALT, sep = "_"),
                haplogroup_defining = match_key %in% hg_markers_snv$match_key)

message("\n=== Haplogroup-defining variant flagging summary ===")
qc_passed_calls %>% dplyr::filter(hf_class == "homoplasmic") %>% dplyr::group_by(cohort) %>%
  dplyr::summarise(n_homoplasmic_calls = dplyr::n(), n_flagged_haplogroup_marker = sum(haplogroup_defining),
                   n_remaining = n_homoplasmic_calls - n_flagged_haplogroup_marker, .groups = "drop") %>% print(n = Inf)

homoplasmic_calls <- qc_passed_calls %>% dplyr::filter(hf_class == "homoplasmic")
tstv_incl <- homoplasmic_calls %>% dplyr::group_by(cohort) %>% dplyr::group_modify(~ compute_tstv(.x)) %>%
  dplyr::ungroup() %>% dplyr::mutate(marker_set = "including_haplogroup_markers")
tstv_excl <- homoplasmic_calls %>% dplyr::filter(!haplogroup_defining) %>% dplyr::group_by(cohort) %>%
  dplyr::group_modify(~ compute_tstv(.x)) %>% dplyr::ungroup() %>% dplyr::mutate(marker_set = "excluding_haplogroup_markers")
tstv_summary <- dplyr::bind_rows(tstv_incl, tstv_excl) %>% dplyr::arrange(cohort, marker_set)

message("\n=== Ts/Tv ratio, homoplasmic variants (per cohort) ===")
print(tstv_summary, n = Inf)
readr::write_csv(tstv_summary, file.path(out_dir, "homoplasmic_tstv_summary.csv"))

## -----------------------------------------------------------------------
## 7. Diagnostic plots
## -----------------------------------------------------------------------
## cohort colour palette (soft/natural tones), keyed to this script's short
## cohort names -- consistent with heteroplasmy_diagnostic_plots_by_cohort.R
cohort_colors_raw <- c(CAAPA = "#3B6FB6", BantuSpeakingAfrican = "#4F8A3F", SouthAfricanColoured = "#D77A2D")
cohort_name_map <- c(Bantu = "BantuSpeakingAfrican", SAC = "SouthAfricanColoured", CAAPA = "CAAPA")
cohort_colors <- stats::setNames(cohort_colors_raw[cohort_name_map], names(cohort_name_map))

theme_set(theme_bw(base_size = 12))
p_depth <- ggplot(metadata_all, aes(cohort, mean_depth, fill = cohort)) +
  geom_boxplot(outlier.shape = NA, alpha = .7) + geom_jitter(width = .15, size = .8, alpha = .5) +
  geom_hline(yintercept = min_mean_depth, linetype = "dashed", colour = "red") +
  scale_fill_manual(values = cohort_colors) +
  labs(title = "Mean mtDNA read depth per sample", y = "Mean depth", x = NULL) + theme(legend.position = "none")
p_contam <- ggplot(metadata_all, aes(cohort, contamination, fill = cohort)) +
  geom_boxplot(outlier.shape = NA, alpha = .7) + geom_jitter(width = .15, size = .8, alpha = .5) +
  geom_hline(yintercept = max_contamination, linetype = "dashed", colour = "red") +
  scale_fill_manual(values = cohort_colors) +
  labs(title = "Estimated contamination per sample", y = "Contamination", x = NULL) + theme(legend.position = "none")
p_cn <- ggplot(metadata_all, aes(cohort, mtcn, fill = cohort)) +
  geom_boxplot(outlier.shape = NA, alpha = .7) + geom_jitter(aes(colour = cn_outlier), width = .15, size = .8) +
  scale_fill_manual(values = cohort_colors) +
  scale_colour_manual(values = c(`TRUE` = "red", `FALSE` = "grey30")) +
  labs(title = "mtDNA copy number per sample", y = "mtDNA-CN", x = NULL, colour = "CN outlier") + theme(legend.position = "bottom")
variants_per_sample <- qc_passed_calls %>% dplyr::group_by(cohort, sample_id, hf_class) %>% dplyr::summarise(n = dplyr::n(), .groups = "drop")
p_vps <- ggplot(variants_per_sample, aes(hf_class, n, fill = cohort)) +
  geom_boxplot(outlier.shape = NA, alpha = .7, position = position_dodge()) +
  geom_jitter(position = position_jitterdodge(jitter.width = .15), size = .6, alpha = .4) +
  scale_fill_manual(values = cohort_colors) +
  labs(title = "Variants per sample by heteroplasmy class (QC-passed calls)", x = NULL, y = "Number of variants") +
  theme(legend.position = "bottom")
plots <- list(depth = p_depth, contamination = p_contam, cn = p_cn, variants_per_sample = p_vps)
purrr::iwalk(plots, function(p, nm) ggsave(file.path(out_dir, "plots", paste0(nm, ".png")), p, width = 7, height = 5, dpi = 300))


## -----------------------------------------------------------------------
## 8. Save outputs for downstream scripts
## -----------------------------------------------------------------------
qc_passed_samples <- metadata_all %>% dplyr::filter(sample_pass_qc)
saveRDS(vcf_long_all, file.path(out_dir, "variant_calls_all_flagged.rds"))
saveRDS(qc_passed_calls, file.path(out_dir, "variant_calls_QC_passed_hgflagged.rds"))
saveRDS(metadata_all, file.path(out_dir, "sample_metadata_all_flagged.rds"))
saveRDS(qc_passed_samples, file.path(out_dir, "sample_metadata_QC_passed.rds"))

message(sprintf(
  "\nDone. %d / %d samples passed sample-level QC. %d variant calls passed QC across %d samples.\nOutputs in: %s",
  nrow(qc_passed_samples), nrow(metadata_all), nrow(qc_passed_calls),
  dplyr::n_distinct(qc_passed_calls$sample_id), normalizePath(out_dir)
))