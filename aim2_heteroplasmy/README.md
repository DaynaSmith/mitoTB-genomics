# aim2_heteroplasmy: Mitochondrial Heteroplasmy and TB Susceptibility

**Objective:** To investigate the relationship between mitochondrial heteroplasmy and tuberculosis (TB) susceptibility and characterise patterns of heteroplasmic variation in African populations.

This aim analyses mitochondrial variant calls derived from whole-genome sequencing (WGS) across three cohorts: **ResisTB** (internally `Bantu`), **LTBI** (internally `SAC`) and **CAAPA**. It encompasses variant/sample quality control, heteroplasmy classification, mitochondrial genomic annotation, burden analyses, lineage comparisons, recurrent-site analyses and homopolymer length-heteroplasmy characterisation.

## Pipeline summary

```text
WGS-derived mtDNA VCFs + sample metadata + haplogroup assignments
                              │
                              ▼
[01_mtDNA_cohort_QC.R]
  Sample/call QC, heteroplasmy classification, recurrent-site filtering,
  haplogroup-marker flags
             │
             ├────────────────────────────────────────────┐
             ▼                                            ▼
[02_variant_annotation_mitomap.R]          [05_poly_C_length_heteroplasmy_composition.R]
  rCRS-based region and consequence          Per-sample length-allele composition
  annotation; pathogenicity screening       at m.302, m.567, m.955, m.16182
             │                                            │
             ▼                                            └──► focused tract figures
[03_heteroplasmy_haplogroup_association.R]
  Per-sample low/intermediate/total burden (including zero counts)
  Cohort-stratified models: age, sex, mtDNA-CN, depth, haplogroup
             │
             ▼
[04_heteroplasmy_diagnostic_plots_by_cohort.R]
  Six-panel cohort-specific diagnostics
             │
             └──► supplementary/ analyses and figures/ visualisations
```

**Execution order:** Run `01` first, followed by `02` and `03`, then `04`. The poly-C composition stage (`05`) depends on QC outputs and may be run after `01` independently of the association branch. The numbered names represent the recommended order, not a requirement that all branches execute serially.

## Directory contents

```text
aim2_heteroplasmy/
├── README.md
├── scripts/
│   ├── 01_mtDNA_cohort_QC.R
│   ├── 02_variant_annotation_mitomap.R
│   ├── 03_heteroplasmy_haplogroup_association.R
│   ├── 04_heteroplasmy_diagnostic_plots_by_cohort.R
│   ├── 05_poly_C_length_heteroplasmy_composition.R
│   ├── supplementary/
│   │   ├── calculate_percentage_heteroplasmies.R
│   │   ├── check_L0a1_bantu_position_concentration.R
│   │   └── L0a1_m955_concordance.R
│   ├── figures/
│   │   ├── Figure_4_8A.R
│   │   ├── Figure_4_15.R
│   │   ├── Figure_4_17.R
│   │   ├── Heteroplasmy_burden_figures.R
│   │   ├── Six_panel_heteroplasmy_figure.R
│   │   ├── Macrohaplogroup_forest_plots.R
│   │   ├── Subhaplogroups_forest_plots.R
│   │   ├── boxplot_subhaplogroup_burden.R
│   │   ├── region_type_and_consequence_barplots.R
│   │   └── heteroplasmic_region_type_and_consequence_barplots.R
│   └── utils/
│       ├── plotting_themes.R
│       └── Statistical_label_functions.R
└── figures/
```

Some figure scripts expect objects already loaded into an R session, while others contain hardcoded summary counts.

## Methods

### Input data and sample-level quality control (`01_mtDNA_cohort_QC.R`)

The main pipeline reads WGS-derived mitochondrial VCFs (`Bantu_mutect2.03.merge.vcf.gz`, `SAC_mutect2.03.merge.vcf.gz`, `CAAPA_mutect2.03.merge.vcf`), corresponding metadata and cohort-specific HaploGrep output. Metadata includes sample ID, age, sex, haplogroup, mtDNA copy number, mean sequencing depth and contamination estimates. The script uses `vcfR::vcfR2tidy(single_frame = TRUE)` to avoid incorrect joins of fixed and genotype fields at multiallelic positions.

Configured sample-level thresholds are **mean mtDNA depth ≥100×** and **contamination ≤0.02**; mtDNA copy-number outliers are evaluated with an IQR-based procedure (`cn_outlier_k = 3`). Samples without the required haplogroup information are flagged. Inspect the sample-QC flags in the script and `sample_QC_summary.csv` for the precise disposition of each exclusion/diagnostic flag.

### Call-level QC and heteroplasmy classification (`01_mtDNA_cohort_QC.R`)

The configured call-level requirements are **site depth ≥100×**, **≥5 alternate reads**, **heteroplasmy fraction (HF) ≥0.05**, and a passing variant FILTER status. Calls are classified as:

| Class | Heteroplasmy fraction |
|---|---|
| Below threshold | HF < 0.05 |
| Low-level heteroplasmy (`low_het`) | 0.05 ≤ HF < 0.10 |
| Intermediate-level heteroplasmy (`intermediate_het`) | 0.10 ≤ HF < 0.95 |
| Homoplasmic (`homoplasmic`) | HF ≥ 0.95 |

Haplogroup-defining variants are flagged using HaploGrep marker information. Burden calculations exclude haplogroup-defining homoplasmic markers; the analysis distinguishes heteroplasmic variation from fixed maternal-lineage markers.

### Recurrent-site and homopolymer artifact filtering (`01_mtDNA_cohort_QC.R`)

Two recurrence screens are implemented: a **per-allele frequency threshold >2%** within defined homopolymer regions, and a **per-position frequency threshold >5%** across each cohort. The pipeline explicitly exempts **m.302, m.567, m.955 and m.16182** from both screens because they are treated as established recurrent poly-C length-heteroplasmy sites. Distinct alleles at the same position are retained as separate variant records rather than collapsed into a single binary event.

Other artifact-prone sites and homopolymer intervals are specified directly in the QC script. The exemption of the four recurrent sites is a deliberate analysis choice and should be considered when interpreting the observed concentration of heteroplasmy in control-region/poly-C tracts.

### Functional and genomic annotation (`02_variant_annotation_mitomap.R`)

QC-passed variants are annotated by mitochondrial region/gene and, where applicable, coding consequence using the **rCRS (NC_012920.1)** and the human mitochondrial genetic code documented by MITOMAP. The script obtains or caches the rCRS sequence, distinguishes synonymous and nonsynonymous substitutions, flags relevant mitochondrial-code caveats and screens against an external `mitomap_confirmed_pathogenic.csv` file. Length variants at the four specified poly-C hotspots are labelled `poly_C_length_heteroplasmy` rather than generically designated for indel manual review. The script produces annotated and burden-ready variant data, including `variant_calls_annotated_dedup_for_burden.rds`.

### Heteroplasmic burden and demographic/genomic correlates (`03_heteroplasmy_haplogroup_association.R`)

Per-sample burden is defined as the number of QC-passed low-level, intermediate-level or total heteroplasmic variant calls. The sample-level dataset retains **individuals with zero heteroplasmies**, rather than restricting analyses to people with detected heteroplasmy.

Regression is performed **separately within each cohort**, reflecting the confounding between cohort and sequencing depth. For each cohort, models are fitted at both macrohaplogroup and subhaplogroup resolution for three outcomes (`n_low_het`, `n_intermediate_het`, `n_total_het`):

```text
heteroplasmic_burden ~ age + sex + mtcn + mean_depth + haplogroup
```

`MASS::glm.nb` is used for negative binomial regression; the script falls back to a quasi-Poisson GLM when the negative binomial fit fails. Rare haplogroup categories are grouped as `Other`, with the reference category chosen within each cohort. Regression outputs include exponentiated coefficients (incidence rate ratios for negative binomial fits), Wald confidence intervals and p-values. The script documents Bonferroni correction within cohort-specific planned comparison families, and also reports exploratory haplogroup comparisons. **Check the model-family indicator before interpreting fallback quasi-Poisson estimates as negative-binomial IRRs.**

Outputs include `per_sample_heteroplasmic_burden_with_haplogroups.csv`, `cohort_data_prepared.rds`, `negbin_regression_results_by_cohort.csv`, `model_comparison_by_cohort.csv` and model/plot files under `association_analysis/`.

### Diagnostic figures and descriptive tests (`04_heteroplasmy_diagnostic_plots_by_cohort.R`)

Six-panel figures are generated separately for each cohort, showing heteroplasmic burden against **sex, age, mean read depth, mtDNA copy number, macrohaplogroup and subhaplogroup**. Sex comparisons use Wilcoxon rank-sum tests; age, depth and copy number use Spearman rank correlations; haplogroup comparisons use Kruskal–Wallis tests. Continuous-variable panels show observations without fitted smoothers. The script requires `association_analysis/cohort_data_prepared.rds` from the association pipeline.

`Heteroplasmy_burden_figures.R`, `Six_panel_heteroplasmy_figure.R`, `Figure_4_8A.R`, `plotting_themes.R` and `Statistical_label_functions.R` provide supplementary burden summaries, plotting functions, consistent colour palettes and statistical annotations. Some are interactive-session scripts that expect objects such as `plot_data` or `variant_calls_annotated_dedup_for_burden` already to exist.

### Genomic distribution, functional consequences and recurrent sites

`region_type_and_consequence_barplots.R` visualises distinct variant positions across mitochondrial region types and synonymous/missense coding substitutions for **all QC-passed variants**. `heteroplasmic_region_type_and_consequence_barplots.R` presents the analogous summaries **restricted to low- and intermediate-level heteroplasmic variants**. Both scripts contain **hardcoded summary counts** derived from earlier outputs; rerunning QC or annotation does not automatically update those plots.

`calculate_percentage_heteroplasmies.R` counts **distinct QC-passed individuals per mitochondrial position per cohort**, divides by the total number of QC-passed cohort members and plots recurrent sites, including m.302, m.567, m.955 and m.16182. Counting unique individuals avoids inflating site recurrence when one individual carries several length alleles at a position.

### Haplogroup and maternal-lineage comparisons

The association script compares burden at macrohaplogroup and finer subhaplogroup resolutions. `boxplot_subhaplogroup_burden.R` tidies per-cohort subhaplogroup burden plots; `Macrohaplogroup_forest_plots.R` and `Subhaplogroups_forest_plots.R` regenerate labelled forest plots from the saved regression CSV without changing fitted statistics. The subhaplogroup forest plots are explicitly exploratory and should not be presented as Bonferroni-significant results unless supported by the underlying analyses.

`check_L0a1_bantu_position_concentration.R` evaluates how much of the heteroplasmic call burden among ResisTB/Bantu L0a1 individuals falls at **the full m.302 (302–315) and m.955 (955–965) tract ranges**, rather than counting only exact anchor positions. It also counts distinct sample carriers at selected positions.

### Poly-C tract length heteroplasmy (`05_poly_C_length_heteroplasmy_composition.R`)

Per-sample allele-composition profiles are generated for **m.302, m.567, m.955 and m.16182**. For each target position, the script combines passing alternate-allele fractions, an estimate from failing calls and a residual reference fraction. Frequently observed length alleles are tracked individually; other alleles and small composition components (<5%) are assigned to `Other`.

The tracked allele sets are **position-specific**, not assumed identical across poly-C tracts. In particular, the recurrent m.567 alleles are longer insertions, while the m.16182 tracked states include a deletion. The script writes position-specific composition CSVs and plots.

**Depth limitation:** The input VCF does not provide depth for every sample at positions without a variant call. Consequently, the source-paper-style minimum-depth criterion at each target position cannot be verified for all zero-call samples without independent per-position coverage data. Such samples are included by default in this composition workflow.

### Focused m.955 and L0a1 analyses (`L0a1_m955_concordance.R`, `Figure_4_15.R`, `Figure_4_17.R`)

`L0a1_m955_concordance.R` compares the distribution of m.302 and m.955 heteroplasmy across ResisTB/Bantu subhaplogroups using a per-sample jitter plot, a carrier-percentage heatmap and an L0a1-versus-other-lineages comparison.

`Figure_4_15.R` displays multi-allelic m.955 length heteroplasmy in L0a1 individuals using the composition output and QC-passed variant calls. `Figure_4_17.R` visualises **read-level** m.955 poly-C tract lengths in ResisTB L0a1b1a1 individuals from the separately prepared `m955_consensus_vs_FASTQ.tsv` file; states supported by at least 5% of pure-C reads are plotted alongside consensus tract-length calls. The supplied R scripts **plot** this read-level table but do not generate it from FASTQ/BAM files.

## Running the analysis

Run the core stages from the project working directory after checking input paths and installing the required R packages:

```bash
Rscript scripts/01_mtDNA_cohort_QC.R
Rscript scripts/02_variant_annotation_mitomap.R
Rscript 'scripts/03_heteroplasmy_haplogroup_association.R'
Rscript scripts/04_heteroplasmy_diagnostic_plots_by_cohort.R
Rscript scripts/05_poly_C_length_heteroplasmy_composition.R
```

Scripts in `scripts/supplementary/` and `scripts/figures/` can then be run when their respective intermediate data files or in-memory objects are available. Some require manual path adjustments or sourcing of shared plotting functions; they are **not** all independent command-line pipelines.

## Software and resources

- **R**, with packages used across scripts including `vcfR`, `dplyr`, `tidyr`, `stringr`, `purrr`, `readr`, `tibble`, `forcats`, `ggplot2`, `patchwork`, `MASS`, `broom`, `ggrepel` and `scales` (package requirements differ by script).
- **HaploGrep** outputs for haplogroup assignment and haplogroup-defining marker identification.
- **rCRS** mitochondrial reference sequence, accession **NC_012920.1**.
- **MITOMAP** human mitochondrial genetic-code documentation and confirmed-pathogenicity data.
- WGS-derived mitochondrial variant calls produced upstream (the supplied R scripts begin from existing VCFs).

## Scope and reproducibility notes

1. **TB susceptibility:** Although TB susceptibility is a stated aim, the supplied scripts do **not** clearly implement a dedicated case–control model of heteroplasmic burden versus TB status. Do not infer such a model from the demographic/lineage regression; add the relevant script and its phenotype definition if available.
2. **Upstream calling:** Alignment, mitochondrial read extraction, heteroplasmy variant calling and mtDNA copy-number estimation occur before these R scripts and are not reproduced here.
3. **Output paths:** Some scripts use `mtDNA_QC_outputs/`, while the m.955 poster script also defines `analysis_outputs/` and uses explicit QC-output paths. Confirm the working directory before running.
4. **Hardcoded plots:** The region/consequence barplot scripts contain fixed counts and must be regenerated or manually updated after data changes.
5. **Recurrent-site exceptions:** Four biologically motivated poly-C sites bypass recurrence-based artifact exclusions; this choice materially affects burden and recurrence summaries.
6. **Count units:** Distinguish **allele-specific call counts**, **distinct mitochondrial positions**, and **numbers of individuals carrying a variant**; these are not interchangeable.
7. **Statistical interpretation:** Cohort-stratified results should not be read as pooled estimates; exploratory lineage contrasts should be distinguished from corrected planned comparisons.
8. **Read-level confirmation:** `Figure_4_17.R` requires a prepared FASTQ-derived summary table. The read-extraction/counting script was not among the supplied files.

## References and external resources

- Chen et al. (2023), *Nature*: mitochondrial heteroplasmic indels and poly-C length-heteroplasmy composition (referenced in the QC and composition scripts; verify full bibliographic details against the dissertation reference list).
- Laricchia et al. (2022): mitochondrial variant QC and artifact-prone sites (referenced in `01_mtDNA_cohort_QC.R`; verify full citation details).
- [MITOMAP: Human Mitochondrial Genetic Code](https://www.mitomap.org/foswiki/bin/view/MITOMAP/HumanMitoCode).
- NCBI RefSeq: [NC_012920.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_012920.1).

**Citation note:** Full references, versions and upstream WGS processing details should be checked against the final thesis methods and bibliography before public release.
