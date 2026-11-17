## =============================================================================
## 03_variant_annotation_mitomap.R
##
## Annotates QC-passed variants (gene, region type, functional consequence)
## using the MITOMAP-documented human mitochondrial genetic code
## (https://www.mitomap.org/foswiki/bin/view/MITOMAP/HumanMitoCode), and
## screens against MitoMap pathogenicity data.
##
## Consequence/annotation approach: the codon-translation engine (rCRS-based)
## is the sole annotation source -- MitoMap does not provide a clean
## bulk-download API for its general variant catalogue, so this is used
## directly rather than approximated. Two caveats documented on the
## HumanMitoCode page are surfaced as explicit flag columns (see section 3):
## contested AGA/AGG stop-codon status, and documented alternative
## initiation codons for MT-ND2/MT-ATP6.
##
## NEW: indel calls at the four positions confirmed (in script 01) as
## genuine poly-C length-heteroplasmy hotspots (302, 567, 955, 16182; see
## Chen et al. 2023 Nature) are labelled "poly_C_length_heteroplasmy"
## instead of the generic "indel_manual_review" -- these don't need manual
## review, their mechanism is documented, and each distinct tract length is
## treated as its own variant (not collapsed), matching how the source
## paper itself analyses them (each length-allele plotted as its own point
## in their sibling-concordance figure).
##
## Run after 01_mtDNA_cohort_QC.R and 02_heteroplasmy_diagnostic_plots_by_cohort.R
## =============================================================================

library(dplyr); library(tidyr); library(stringr); library(purrr)
library(readr); library(tibble); library(ggplot2)

qc_dir <- "mtDNA_QC_outputs"
qc_passed_calls <- readRDS(file.path(qc_dir, "variant_calls_QC_passed_hgflagged.rds"))

## -----------------------------------------------------------------------
## 1. CONFIG
## -----------------------------------------------------------------------
pathogenic_path <- "mitomap_confirmed_pathogenic.csv"

## must match script 01's confirmed_common_heteroplasmy_positions
confirmed_common_heteroplasmy_positions <- c(302, 567, 955, 16182)

## -----------------------------------------------------------------------
## 2. rCRS reference sequence (cached)
## -----------------------------------------------------------------------
rcrs_cache <- file.path(qc_dir, "rCRS_NC_012920.1.fasta")
get_rcrs_sequence <- function(cache_path = rcrs_cache) {
  if (!file.exists(cache_path)) {
    message("Downloading rCRS reference sequence (NC_012920.1)...")
    url <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_012920.1&rettype=fasta&retmode=text"
    download.file(url, cache_path, quiet = TRUE)
  }
  lines <- readLines(cache_path); lines <- lines[!str_starts(lines, ">")]
  toupper(paste(lines, collapse = ""))
}
rcrs <- get_rcrs_sequence()
stopifnot(nchar(rcrs) == 16569)

complement_base <- function(b) c(A = "T", T = "A", C = "G", G = "C")[b]
revcomp <- function(seq) paste(rev(complement_base(str_split(seq, "")[[1]])), collapse = "")

## -----------------------------------------------------------------------
## 3. Gene/region coordinate table + codon translation engine
## -----------------------------------------------------------------------
gene_table <- tribble(
  ~gene,   ~start, ~end,   ~strand, ~type,
  "D-LOOP_a", 16024, 16569, "+", "control_region",
  "D-LOOP_b",     1,   576, "+", "control_region",
  "MT-TF",      577,   647, "+", "tRNA",
  "MT-RNR1",    648,  1601, "+", "rRNA",
  "MT-TV",     1602,  1670, "+", "tRNA",
  "MT-RNR2",   1671,  3229, "+", "rRNA",
  "MT-TL1",    3230,  3304, "+", "tRNA",
  "MT-ND1",    3307,  4262, "+", "protein_coding",
  "MT-TI",     4263,  4331, "+", "tRNA",
  "MT-TQ",     4329,  4400, "-", "tRNA",
  "MT-TM",     4402,  4469, "+", "tRNA",
  "MT-ND2",    4470,  5511, "+", "protein_coding",
  "MT-TW",     5512,  5579, "+", "tRNA",
  "MT-TA",     5587,  5655, "-", "tRNA",
  "MT-TN",     5657,  5729, "-", "tRNA",
  "OL",        5721,  5798, "+", "OL",
  "MT-TC",     5761,  5826, "-", "tRNA",
  "MT-TY",     5826,  5891, "-", "tRNA",
  "MT-CO1",    5904,  7445, "+", "protein_coding",
  "MT-TS1",    7446,  7514, "-", "tRNA",
  "MT-TD",     7518,  7585, "+", "tRNA",
  "MT-CO2",    7586,  8269, "+", "protein_coding",
  "MT-TK",     8295,  8364, "+", "tRNA",
  "MT-ATP8",   8366,  8572, "+", "protein_coding",
  "MT-ATP6",   8527,  9207, "+", "protein_coding",
  "MT-CO3",    9207,  9990, "+", "protein_coding",
  "MT-TG",     9991, 10058, "+", "tRNA",
  "MT-ND3",   10059, 10404, "+", "protein_coding",
  "MT-TR",    10405, 10469, "+", "tRNA",
  "MT-ND4L",  10470, 10766, "+", "protein_coding",
  "MT-ND4",   10760, 12137, "+", "protein_coding",
  "MT-TH",    12138, 12206, "+", "tRNA",
  "MT-TS2",   12207, 12265, "+", "tRNA",
  "MT-TL2",   12266, 12336, "+", "tRNA",
  "MT-ND5",   12337, 14148, "+", "protein_coding",
  "MT-ND6",   14149, 14673, "-", "protein_coding",
  "MT-TE",    14674, 14742, "-", "tRNA",
  "MT-CYB",   14747, 15887, "+", "protein_coding",
  "MT-TT",    15888, 15953, "+", "tRNA",
  "MT-TP",    15956, 16023, "-", "tRNA"
)
find_overlapping_genes <- function(pos) gene_table %>% filter(pos >= start, pos <= end)

## Vertebrate mitochondrial genetic code, per MITOMAP's HumanMitoCode page
## (https://www.mitomap.org/foswiki/bin/view/MITOMAP/HumanMitoCode).
## Differences from the universal/nuclear code: UGA=Trp (not stop),
## AUA=Met (not Ile), AGA/AGG=Stop (not Arg -- see caveat below).
mito_codon_table <- c(
  TTT="F",TTC="F",TTA="L",TTG="L",CTT="L",CTC="L",CTA="L",CTG="L",
  ATT="I",ATC="I",ATA="M",ATG="M",GTT="V",GTC="V",GTA="V",GTG="V",
  TCT="S",TCC="S",TCA="S",TCG="S",CCT="P",CCC="P",CCA="P",CCG="P",
  ACT="T",ACC="T",ACA="T",ACG="T",GCT="A",GCC="A",GCA="A",GCG="A",
  TAT="Y",TAC="Y",TAA="*",TAG="*",CAT="H",CAC="H",CAA="Q",CAG="Q",
  AAT="N",AAC="N",AAA="K",AAG="K",GAT="D",GAC="D",GAA="E",GAG="E",
  TGT="C",TGC="C",TGA="W",TGG="W",CGT="R",CGC="R",CGA="R",CGG="R",
  AGT="S",AGC="S",AGA="*",AGG="*",GGT="G",GGC="G",GGA="G",GGG="G"
)
translate_codon <- function(codon) unname(mito_codon_table[codon])

## AGA/AGG stop-codon status is explicitly flagged as under active review on
## MITOMAP's HumanMitoCode page (Temperley et al. 2010 argue human
## mitochondria use only UAA/UAG as true stops; MITOMAP retains AGA/AGG as
## STOP pending resolution). Any consequence call resting on an AGA/AGG stop
## assignment is flagged rather than silently trusted.
uncertain_stop_codons <- c("AGA", "AGG")

## Documented alternative initiation codons (MITOMAP HumanMitoCode): ND2 can
## initiate at ATT (read as Met at that position, not the standard Ile);
## ATP6 has been reported to initiate at GUG (Dubot et al. 2004). These only
## affect the literal first codon of the gene, not translation elsewhere.
alternative_start_genes <- c("MT-ND2", "MT-ATP6")

annotate_in_gene <- function(pos, ref, alt, gene_row) {
  if (nchar(ref) != 1 || nchar(alt) != 1) {
    consequence <- if (pos %in% confirmed_common_heteroplasmy_positions) "poly_C_length_heteroplasmy" else "indel_manual_review"
    return(tibble(gene = gene_row$gene, region_type = gene_row$type, codon_number = NA_integer_,
                  ref_aa = NA_character_, alt_aa = NA_character_, consequence = consequence,
                  ref_base_mismatch = NA, mitomap_stop_uncertainty_flag = FALSE, alternative_start_flag = FALSE))
  }
  if (gene_row$type != "protein_coding") {
    return(tibble(gene = gene_row$gene, region_type = gene_row$type, codon_number = NA_integer_,
                  ref_aa = NA_character_, alt_aa = NA_character_, consequence = paste0(gene_row$type, "_variant"),
                  ref_base_mismatch = NA, mitomap_stop_uncertainty_flag = FALSE, alternative_start_flag = FALSE))
  }
  strand <- gene_row$strand; gstart <- gene_row$start; gend <- gene_row$end
  offset0 <- if (strand == "+") pos - gstart else gend - pos
  codon_number <- (offset0 %/% 3) + 1; pos_in_codon <- (offset0 %% 3) + 1
  if (strand == "+") {
    codon_start <- gstart + (codon_number - 1) * 3; codon_end <- codon_start + 2
    ref_codon <- substr(rcrs, codon_start, codon_end)
  } else {
    codon_end <- gend - (codon_number - 1) * 3; codon_start <- codon_end - 2
    ref_codon <- revcomp(substr(rcrs, codon_start, codon_end))
  }
  if (nchar(ref_codon) != 3) {
    return(tibble(gene = gene_row$gene, region_type = gene_row$type, codon_number = codon_number,
                  ref_aa = NA_character_, alt_aa = NA_character_, consequence = "boundary_codon_manual_review",
                  ref_base_mismatch = NA, mitomap_stop_uncertainty_flag = FALSE, alternative_start_flag = FALSE))
  }
  ref_base_check <- if (strand == "+") substr(rcrs, pos, pos) else complement_base(substr(rcrs, pos, pos))
  ref_mismatch <- !is.na(ref_base_check) && ref_base_check != ref
  alt_codon <- ref_codon
  alt_allele_in_codon <- if (strand == "+") alt else complement_base(alt)
  substr(alt_codon, pos_in_codon, pos_in_codon) <- alt_allele_in_codon
  ref_aa <- translate_codon(ref_codon); alt_aa <- translate_codon(alt_codon)
  consequence <- dplyr::case_when(
    is.na(ref_aa) | is.na(alt_aa) ~ "translation_error_manual_review",
    ref_aa == alt_aa ~ "synonymous", ref_aa == "*" ~ "stop_lost", alt_aa == "*" ~ "stop_gained",
    TRUE ~ "missense"
  )
  
  involves_uncertain_stop <- (ref_codon %in% uncertain_stop_codons && ref_aa == "*") ||
    (alt_codon %in% uncertain_stop_codons && alt_aa == "*")
  mitomap_stop_uncertainty_flag <- involves_uncertain_stop && consequence %in% c("stop_gained", "stop_lost")
  alternative_start_flag <- codon_number == 1 && gene_row$gene %in% alternative_start_genes
  
  tibble(gene = gene_row$gene, region_type = gene_row$type, codon_number = codon_number,
         ref_aa = ref_aa, alt_aa = alt_aa, consequence = consequence, ref_base_mismatch = ref_mismatch,
         mitomap_stop_uncertainty_flag = mitomap_stop_uncertainty_flag,
         alternative_start_flag = alternative_start_flag)
}

annotate_variant_codon <- function(pos, ref, alt) {
  overlaps <- find_overlapping_genes(pos)
  if (nrow(overlaps) == 0) {
    consequence <- if (nchar(ref) != 1 || nchar(alt) != 1) {
      if (pos %in% confirmed_common_heteroplasmy_positions) "poly_C_length_heteroplasmy" else "intergenic_indel_manual_review"
    } else "intergenic_variant"
    return(tibble(gene = "intergenic", region_type = "intergenic", codon_number = NA_integer_,
                  ref_aa = NA_character_, alt_aa = NA_character_, consequence = consequence,
                  ref_base_mismatch = NA, mitomap_stop_uncertainty_flag = FALSE, alternative_start_flag = FALSE))
  }
  purrr::pmap_dfr(overlaps, function(...) annotate_in_gene(pos, ref, alt, tibble(...)))
}

## -----------------------------------------------------------------------
## 4. Annotate distinct QC-passed variants
## -----------------------------------------------------------------------
distinct_variants <- qc_passed_calls %>% distinct(POS, REF, ALT)
message(sprintf("Annotating %d distinct variants via codon-translation engine...", nrow(distinct_variants)))
annotated_codon <- distinct_variants %>%
  mutate(ann = pmap(list(POS, REF, ALT), annotate_variant_codon)) %>% unnest(ann)

n_mismatch <- sum(annotated_codon$ref_base_mismatch, na.rm = TRUE)
if (n_mismatch > 0) warning(sprintf("%d records have REF not matching rCRS -- check reference genome/polarity before trusting their consequence calls.", n_mismatch))

n_stop_uncertain <- sum(annotated_codon$mitomap_stop_uncertainty_flag, na.rm = TRUE)
if (n_stop_uncertain > 0) message(sprintf("%d variant(s) have a stop_gained/stop_lost call resting on the contested AGA/AGG stop-codon assignment (flagged `mitomap_stop_uncertainty_flag`).", n_stop_uncertain))

n_alt_start <- sum(annotated_codon$alternative_start_flag, na.rm = TRUE)
if (n_alt_start > 0) message(sprintf("%d variant(s) affect the first codon of MT-ND2 or MT-ATP6 (flagged `alternative_start_flag`).", n_alt_start))

n_poly_c <- sum(annotated_codon$consequence == "poly_C_length_heteroplasmy", na.rm = TRUE)
message(sprintf("%d distinct length-alleles labelled 'poly_C_length_heteroplasmy' at the 4 confirmed common-heteroplasmy positions.", n_poly_c))

annotated <- annotated_codon

## -----------------------------------------------------------------------
## 5. Join annotation back onto every QC-passed call; build de-duplicated
## version for per-sample burden counting (overlapping-gene regions like
## ATP8/ATP6 and ND4L/ND4 legitimately produce two annotation rows; poly-C
## length-heteroplasmy alleles are NOT collapsed here -- each distinct
## length-allele a sample carries counts as its own call, matching how
## Chen et al. 2023 treat them, so a sample with calls at both m.955:A>AC
## and m.955:A>ACC contributes two rows, not one).
## -----------------------------------------------------------------------
qc_passed_calls_annotated <- qc_passed_calls %>%
  dplyr::left_join(annotated, by = c("POS", "REF", "ALT"), relationship = "many-to-many")

region_priority <- c(protein_coding = 1, tRNA = 2, rRNA = 2, OL = 3, control_region = 4, intergenic = 5)
qc_passed_calls_annotated_dedup <- qc_passed_calls_annotated %>%
  dplyr::group_by(cohort, sample_id, POS, REF, ALT) %>%
  dplyr::arrange(region_priority[region_type], .by_group = TRUE) %>%
  dplyr::summarise(
    dplyr::across(-c(gene, region_type, codon_number, ref_aa, alt_aa, consequence,
                     dplyr::any_of(c("ref_base_mismatch", "mitomap_stop_uncertainty_flag", "alternative_start_flag"))),
                  dplyr::first),
    gene = paste(unique(gene), collapse = ","),
    region_type = dplyr::first(region_type),
    consequence = paste(unique(consequence), collapse = ","),
    mitomap_stop_uncertainty_flag = any(mitomap_stop_uncertainty_flag, na.rm = TRUE),
    alternative_start_flag = any(alternative_start_flag, na.rm = TRUE),
    .groups = "drop"
  )

saveRDS(qc_passed_calls_annotated, file.path(qc_dir, "variant_calls_annotated.rds"))
readr::write_csv(qc_passed_calls_annotated, file.path(qc_dir, "variant_calls_annotated.csv"))
saveRDS(qc_passed_calls_annotated_dedup, file.path(qc_dir, "variant_calls_annotated_dedup_for_burden.rds"))
readr::write_csv(qc_passed_calls_annotated_dedup, file.path(qc_dir, "variant_calls_annotated_dedup_for_burden.csv"))

message("\n=== Variants by region type (distinct positions, per cohort) ===")
qc_passed_calls_annotated %>% distinct(cohort, POS, region_type) %>% count(cohort, region_type) %>% print(n = Inf)

message("\n=== Functional consequence breakdown (protein-coding, distinct positions) ===")
qc_passed_calls_annotated %>% filter(region_type == "protein_coding") %>%
  distinct(cohort, POS, gene, consequence) %>% count(cohort, consequence) %>% print(n = Inf)

message("\n=== Poly-C length-heteroplasmy calls by position and allele (QC-passed) ===")
qc_passed_calls_annotated %>% filter(consequence == "poly_C_length_heteroplasmy") %>%
  count(cohort, POS, REF, ALT, hf_class, name = "n_calls") %>% arrange(POS, cohort, desc(n_calls)) %>% print(n = Inf)

## -----------------------------------------------------------------------
## 6. MitoMap CONFIRMED PATHOGENIC variant screening
## -----------------------------------------------------------------------
patho_raw <- read_csv(pathogenic_path, show_col_types = FALSE, col_names = TRUE)
names(patho_raw) <- c("index", "locus_type", "locus", "diseases", "allele", "position",
                      "aa_or_rna", "status", "last_update")

parse_mitomap_allele <- function(allele_str) {
  m <- str_match(allele_str, "^m\\.(\\d+)([ACGT])>([ACGT])$")
  tibble(pos = suppressWarnings(as.numeric(m[,2])), ref = m[,3], alt = m[,4], is_simple_snv = !is.na(m[,2]))
}
patho_parsed <- patho_raw %>% mutate(parsed = map(allele, parse_mitomap_allele)) %>% unnest(parsed) %>%
  mutate(pos = if_else(is.na(pos), position, pos))
readr::write_csv(patho_parsed, file.path(qc_dir, "mitomap_pathogenic_parsed.csv"))

n_complex <- sum(!patho_parsed$is_simple_snv)
message(sprintf("%d of %d MitoMap pathogenic entries are non-simple-SNV -- matched by POSITION ONLY, manual review needed for those hits.",
                n_complex, nrow(patho_parsed)))

exact_matches <- qc_passed_calls_annotated %>%
  inner_join(patho_parsed %>% filter(is_simple_snv) %>% select(pos, ref, alt, locus, diseases, aa_or_rna, status, last_update),
             by = c("POS" = "pos", "REF" = "ref", "ALT" = "alt")) %>% mutate(match_type = "exact_snv_match")
position_matches <- qc_passed_calls_annotated %>%
  inner_join(patho_parsed %>% filter(!is_simple_snv) %>% select(pos, locus, diseases, aa_or_rna, status, last_update, allele),
             by = c("POS" = "pos")) %>% mutate(match_type = "position_only_manual_review")
pathogenic_hits <- bind_rows(exact_matches, position_matches %>% rename(mitomap_allele = allele))
readr::write_csv(pathogenic_hits, file.path(qc_dir, "mitomap_pathogenic_hits.csv"))

message("\n=== MitoMap confirmed pathogenic variant hits ===")
if (nrow(pathogenic_hits) == 0) {
  message("No individuals carry a variant matching the MitoMap confirmed pathogenic list.")
} else {
  pathogenic_hits %>% select(cohort, sample_id, POS, REF, ALT, AF, hf_class, locus, diseases, status, match_type) %>%
    arrange(match_type, cohort, POS) %>% print(n = Inf)
  message("NOTE: 'position_only_manual_review' hits are matched by genomic position only -- verify the actual")
  message("allele before reporting these as pathogenic. Only 'exact_snv_match' rows are directly confirmed.")
}

message(sprintf("\nDone. Outputs in: %s", normalizePath(qc_dir)))
