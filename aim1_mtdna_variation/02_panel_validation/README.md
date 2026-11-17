# 02_panel_validation: mtDNA Imputation Reference Panel Validation

**Objective:** To evaluate the performance of the mtDNA imputation reference panel across validation cohorts, mitochondrial lineages and genotyping platforms, using sequence-derived mitochondrial genotypes and haplogroups as reference data.

The panel constructed in [`../01_reference_panel/`](../01_reference_panel/) was evaluated by comparing original array genotypes, newly imputed mtSNPs and complete mitochondrial sequence-derived calls. Validation considers **variant recovery**, **genotype concordance** and **haplogroup concordance**. The supplied scripts implement the downstream comparisons and visualisations; the commands used to create the in silico arrays and run IMPUTE2 were not included among the supplied files.

## Pipeline summary

```text
Complete mtDNA sequencing VCFs (truth set) + matched samples
                       │
                       ├───────────────────────────────────────┐
                       ▼                                       │
Original genotyping-array VCFs (or in silico arrays)            │
                       │                                       │
                       ▼                                       │
IMPUTE2 reference panel from ../01_reference_panel/            │
                       │                                       │
                       ▼                                       │
Imputed mtDNA VCFs (after relevant QC / INFO filtering)        │
                       │                                       │
                       ├───────────────────────────────────────┤
                       ▼                                       ▼
[mtDNA_imputation_concordance.R]               Haplogroup calls from sequence,
 1. Retain biallelic mitochondrial SNPs         array and imputed datasets
 2. Exclude all original array positions                      │
 3. Match sequence and imputed variants                       ▼
    using POS + REF + ALT                 [strict_haplogroup_concordance.R]
 4. Match samples; standardise haploid GT      Exact haplogroup agreement
 5. Calculate variant recovery, per-site       PhyloTree 17 hierarchical P/R/F
    TP/TN/FP/FN, MCC and concordance           Paired comparisons + uncertainty
                       │                                       │
                       ▼                                       ▼
 mtdna_imputation_concordance/           Haplogroup concordance CSVs + plots
                       │
                       ▼
[venn_diagrams_updated.R]
 Corrected allele-aware recovery counts supplied as CSV
                       │
                       ▼
 Venn diagrams + recovery summaries (PNG/PDF/CSV)
```

## Directory contents

Suggested repository layout (the input VCFs, haplogroup CSVs and generated results are **not** included):

```text
02_panel_validation/
├── README.md
├── scripts/
│   ├── mtDNA_imputation_concordance.R     # Primary untyped-mtSNP genotype validation
│   ├── strict_haplogroup_concordance.R    # Primary exact and hierarchical haplogroup validation
│   └──  venn_diagrams_updated.R           # Variant recovery visualisation
├── resources/
│   └── phylotree17_parents.csv           # Haplogroup-parent relationships (Build 17)
├── input/                                # User-provided VCFs and haplogroup/count CSVs
└── results/                              # Generated outputs, organised by cohort/platform
```

**Path note:** The R scripts retain their original input/output path settings. Moving them into this suggested structure does **not** automatically update those paths. Set the paths in each script (including `tree_file`) before running, or run from an appropriate working directory.

## Methods

### Validation datasets and imputation

The validation design requires matched complete mitochondrial sequencing data, original array genotypes (including any in silico genotyping-array datasets) and corresponding post-imputation genotypes. The sequence-derived calls are treated as the reference or **truth set**. Performance can be stratified by cohort, genotyping platform and mitochondrial lineage. The supplied scripts include platform labels for Affymetrix 6.0, H3Africa, Infinium GDA, Infinium Omni5-4, Multi-Ethnic AMR/AFR and Multi-Ethnic Global, but do not establish which cohorts/platforms were actually run or their final sample sizes.

The generation of array subsets, IMPUTE2 execution, pre-imputation filtering, post-imputation filtering and haplogroup assignment are upstream prerequisites; their complete commands are not present in the supplied validation scripts. Document the corresponding commands, software versions, thresholds and input manifests separately when available. One example input filename in `mtDNA_imputation_concordance.R` indicates an `INFO0.7`-filtered imputed VCF, but the script itself does not apply that threshold.

### Variant recovery and genotype concordance (`mtDNA_imputation_concordance.R`)

The primary genotype-validation script reads original array, imputed and sequence-derived VCFs, restricts analysis to mitochondrial contigs (`MT`, `chrM`, `chrMT` or `M`) and retains biallelic SNP records. It converts homozygous diploid genotype encodings to haploid allele calls (`0` or `1`); mixed/heterozygous calls and missing genotypes are excluded from genotype comparisons.

**Directly interrogated array positions are excluded by genomic position**, regardless of whether REF/ALT representations differ between datasets. Consequently, only sequence-derived biallelic mtSNPs at positions absent from the original array contribute to the newly imputed variant evaluation. Sequence and imputed variant records are then matched using **position, REF and ALT allele together** (`POS_REF_ALT`), after which shared sample IDs are used for genotype comparison.

Variant recovery is defined as:

```text
Recovery rate = sequence-derived untyped mtSNP records recovered by imputation
                -----------------------------------------------------------
                       all sequence-derived untyped mtSNP records
```

For each recovered site, the script tabulates true positives (TP), true negatives (TN), false positives (FP) and false negatives (FN), calculates genotype concordance and evaluates the Matthews correlation coefficient (MCC):

```text
MCC = (TP × TN − FP × FN) /
      sqrt((TP + FP)(TP + FN)(TN + FP)(TN + FN))
```

MCC is undefined when its denominator is zero. Such sites remain relevant to **variant recovery** but do not contribute to the mean or median MCC across estimable sites. The script also calculates overall genotype concordance from the aggregated counts across compared genotype pairs.

**Outputs** (under the configured `out_dir`, default `mtdna_imputation_concordance/`):

- `imputed_only_site_specific_MCC.csv` — site-specific genotype comparisons and MCC.
- `imputed_only_MCC_summary.csv` — recovery counts/rate, sample counts, MCC summary and overall genotype concordance.
- `recovered_imputed_variant_ids.txt` — allele-aware IDs for recovered variants.

### Variant recovery visualisation (`venn_diagrams_updated.R`)

This script plots corrected two-set Venn diagrams comparing sequence-derived and imputed biallelic mtSNP records at **positions absent from the relevant original array**. The input is a prepared CSV with canonical headers:

```csv
platform,n_imputed,n_seq,n_common
```

The input counts must use the same allele-aware variant matching rule as the genotype-concordance analysis, and **must include recovered variants with undefined MCC**. The plotting script does not read VCFs, identify untyped sites or calculate overlaps itself. It validates the counts, calculates recovery percentages and creates individual and combined figures. Circle areas are intentionally not proportional to counts.

Run from the repository root after preparing the counts:

```bash
Rscript scripts/venn_diagrams_updated.R input/corrected_counts.csv results/venn_diagrams
```

**Outputs:** `venn_<platform>.png/.pdf`, `venn_all_platforms_combined.png/.pdf`, and `variant_recovery_summary.csv` in the chosen output directory.

### Haplogroup concordance (`strict_haplogroup_concordance.R`)

Haplogroup assignment accuracy is evaluated against sequence-derived haplogroups using two complementary metrics:

1. **Exact-match proportion:** the proportion of predicted haplogroup labels matching the sequence-derived label after minimal cleaning of whitespace and uncertainty/parenthetical annotations. Wilson 95% confidence intervals are calculated.
2. **Hierarchical precision, recall and F-measure:** PhyloTree Build 17 parent-child relationships are used to compare the ancestor paths of predicted and sequence-derived haplogroups. Shared ancestry receives partial credit when assignments differ in resolution or placement along the same branch. The root (`mt-MRCA`) is excluded from the similarity calculation.

The strict script **does not progressively truncate haplogroup labels** to find a matching tree entry. Labels absent from the supplied PhyloTree parent table are reported and excluded from hierarchical calculations rather than silently mapped to broader lineages. Exact matching remains a separate calculation.

Within each platform, array-derived and imputed haplogroup calls are compared with the same sequence-derived truth calls using paired samples. The script computes per-sample scores, platform/method summaries and paired imputed-minus-array differences, with bootstrap-based uncertainty and a paired McNemar comparison for exact matches. Its settings specify **2,000 bootstrap replicates** and seed `123`.

Supported input CSV layouts:

```text
Single platform:
Sample_ID,True_Haplo,Array,Imputed

Multiple platforms:
Sample_ID,True_Haplo,Array_Affymetrix,Imputed_Affymetrix,
                     Array_H3Africa,Imputed_H3Africa,...
```

The accompanying `resources/phylotree17_parents.csv` provides `haplogroup,parent` columns. Set `infile`, `tree_file` and `output_prefix` in the script before running.

**Outputs** (in the working directory, prefixed by `output_prefix`):

- `_haplogroup_per_sample_scores.csv`
- `_haplogroup_concordance_summary.csv`
- `_haplogroup_concordance_differences.csv`
- `_unmatched_phylotree_labels.csv`
- `_haplogroup_concordance.png` and `.pdf`

### Earlier analysis scripts (`scripts/legacy/`)

These scripts are retained for traceability rather than used as the primary reported validation pipeline:

- `CalculateMCC_GenotypeConcordance.R` compares array and imputed VCFs separately with sequence-derived genotypes and produces overall/per-variant MCC summaries and a comparison plot. Unlike the primary script, it does **not** explicitly exclude directly typed array positions from its imputed-genotype comparison. **Check its paths before reuse:** the uploaded settings point both array and imputed inputs to `ArraySamples_rename.vcf.gz`.
- `HaplogroupConcordance.R` computes exact and hierarchical haplogroup agreement, but its PhyloTree matching routine progressively truncates unrecognised labels to a valid ancestor. This can inflate hierarchical agreement relative to the strict approach. It defaults to 1,000 bootstrap replicates for summary confidence intervals.
- `GenotypeConcordance_VennDiagrams.R` creates earlier PNG-only Venn diagrams from a count table (`bcftools_stats.csv` by default). It does not enforce the corrected allele-aware, untyped-site count definition used by the updated script.

These are **related versions, not byte-for-byte duplicates**. The preferred scripts should be used consistently for reported results; avoid mixing statistics from different definitions or script versions without explicitly labelling them.

## Running the validation

1. Prepare matched array, imputed and sequence VCFs with harmonised sample IDs and mitochondrial reference coordinates. Set the three VCF paths and `out_dir` in `mtDNA_imputation_concordance.R`.
2. Run `Rscript scripts/mtDNA_imputation_concordance.R` from the directory expected by the configured paths. Repeat for each relevant validation cohort/platform.
3. Assemble `corrected_counts.csv` from the validated **untyped, allele-aware** variant sets for each platform. Run the updated Venn plotting script using the command above.
4. Prepare a haplogroup CSV containing sequence truth, array and imputed calls for matched samples. Set the input, PhyloTree CSV path and output prefix in `strict_haplogroup_concordance.R`, then run `Rscript scripts/strict_haplogroup_concordance.R`.
5. Archive the input manifest, exact scripts, software versions, QC settings and output tables for each cohort/platform.

## Software

- R (`Rscript`; exact R version not specified in the uploaded files).
- Base R — VCF reading, genotype comparisons and MCC calculations in the genotype-validation scripts.
- `data.table`, `ggplot2`, `boot` — haplogroup concordance.
- `VennDiagram`, `gridExtra`, `grid` — variant-recovery plots.
- PhyloTree Build 17 — hierarchical haplogroup relationships, via `phylotree17_parents.csv`.
- IMPUTE2 — upstream imputation (see `../01_reference_panel/`); execution scripts not supplied here.
- HaploGrep — upstream haplogroup assignment; command/version for validation calls not supplied here.

## Interpretation and reproducibility notes

- **Recovery is not genotype accuracy:** a recovered variant can have poor genotype concordance, and MCC may be undefined at sites lacking variation in the evaluated samples.
- **Position exclusion precedes allele matching:** excluding all originally typed positions prevents directly observed array sites from inflating imputation accuracy; matching `POS_REF_ALT` subsequently avoids counting mismatched alternate alleles as recovered.
- **Phylogenetic concordance complements exact concordance:** related but nonidentical haplogroups can have meaningful shared ancestry; strict label matching prevents artificial improvement through forced truncation.

## References

- McInerney, T. et al. (2021). A globally diverse reference alignment and panel for imputation of mitochondrial DNA variants. *BMC Bioinformatics*.
- Howie, B.N. et al. (2009). A flexible and accurate genotype imputation method for the next generation of genome-wide association studies. *PLoS Genetics*.
- van Oven, M. & Kayser, M. (2009). Updated comprehensive phylogenetic tree of global human mitochondrial DNA variation. *Human Mutation*. (PhyloTree; Build 17 used for this analysis.)
- Matthews, B.W. (1975). Comparison of the predicted and observed secondary structure of T4 phage lysozyme. *Biochimica et Biophysica Acta*.
- Wilson, E.B. (1927). Probable inference, the law of succession, and statistical inference. *Journal of the American Statistical Association*.
