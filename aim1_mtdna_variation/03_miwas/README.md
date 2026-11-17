# 03_miwas: Inherited mtDNA Variation and TB Susceptibility

**Objective:** To characterise mitochondrial genetic variation and population structure in South African tuberculosis (TB) case–control cohorts, and test associations between mitochondrial haplogroups, individual mtSNPs, and TB susceptibility.

This stage applies the mtDNA imputation reference panel developed in `01_reference_panel/` and evaluated in `02_panel_validation/`. The principal analysis is implemented in `miWAS_cohesive_pipeline.R`; `PCA.R` is an earlier, standalone PCA and visualisation script, not a separate association pipeline.

## Pipeline summary

```
Final QC'd merged imputed mtDNA VCF
  + sample_data.csv (phenotype, demographics, cohort, haplogroups)
  + global_ancestry_proportions.csv (autosomal ancestry)
       │
       ▼
[miWAS_cohesive_pipeline.R]
 1. Characterise cases/controls (age, sex, cohort, haplogroups)
 2. Describe haplogroup distributions and cohort differences
 3. mtDNA PCA: 10 axes, scree plot, haplogroup plots
 4. Assess whether PCs track haplogroup and/or cohort
 5. Recode mtDNA calls to haploid 0/1; uncertain mixed calls → NA
 6. Merge covariates and global autosomal ancestry by SampleID
 7. Test macrohaplogroup association with TB (logistic models/LRT)
 8. Primary mtSNP miWAS: Firth logistic regression
 9. Sensitivity miWAS: add mitochondrial macrohaplogroup
10. Bonferroni and BH-FDR adjustment; QQ plots and lambda
       │
       ▼
miWAS_outputs/ (tables and figures)

[PCA.R] — supplementary standalone mtDNA PCA / cohort plots
```

## Directory contents

Suggested repository layout (input datasets and outputs should only be committed when sharing permissions allow):

```
03_miwas/
├── README.md
├── scripts/
│   ├── miWAS_cohesive_pipeline.R   # Main analysis
│   └── PCA.R                          # Standalone exploratory PCA
├── data/                              # Input files; not included here
│   ├── sample_data.csv
│   ├── global_ancestry_proportions.csv
│   └── <final_qc_imputed_mtDNA>.vcf
└── miWAS_outputs/                     # Created by main script at runtime
```

**Note:** The supplied R scripts currently expect input files in the **working directory**, not `data/`. Either place the inputs in the working directory or edit the paths at the top of the scripts before running. The directory tree above is a proposed organisational structure, not a claim that all these files were supplied.

## Methods

### Input data and cohort characterisation (`miWAS_cohesive_pipeline.R`)

The main script reads `sample_data.csv`, requiring `SampleID`, `Age`, `Sex`, `Phenotype`, `Cohort`, `Haplogroup`, `Macrohaplogroup`, and `BroadHaplogroup`. Phenotype is encoded as `Control` or `Case`; sex as `Male` or `Female`. It reports cohort, phenotype, sex and age summaries; compares ages using a Wilcoxon test (with Shapiro–Wilk normality checks); and assesses sex distributions with a chi-squared test. Haplogroup distributions by phenotype and cohort are compared using Fisher's exact tests with simulated p-values (2,000 replicates). The script also examines phenotype distribution across cohorts and exports descriptive plots.

### mtDNA principal component analysis (`miWAS_cohesive_pipeline.R`; `PCA.R`)

The final QC'd imputed mtDNA VCF is read with `vcfR` and converted to a `genlight` object. `adegenet::glPca()` computes ten principal components. The cohesive pipeline exports a scree plot, sample scores joined to haplogroup/cohort metadata, and PC1–PC2, PC2–PC3 and PC1–PC3 visualisations by mitochondrial lineage. For PCs 1–5, linear models calculate R² for haplogroup alone, cohort alone, and both jointly, to assess the extent to which mtDNA axes reflect lineage or cohort structure.

The separate `PCA.R` script produces `PCA_haplogroups.png` and `PCA_cohorts.png` from the VCF and `sample_data.csv`. Its PCA and plotting functions substantially overlap with the main pipeline, so it is retained for provenance or independent exploratory plots rather than treated as an additional required step.

**Interpretation note:** Both scripts divide the eigenvalues returned by `glPca(nf=10)` by their sum; the resulting percentages describe the **retained axes**, not necessarily the total variance across every possible axis. Although `PCA.R` comments on MAF >1% filtering, it does not itself implement that filter.

### mtDNA genotype recoding

Genotypes are extracted from the VCF and represented as haploid binary calls: `0` for reference, `1` for alternate, and `NA` for missing or uncertain heterozygous/heteroplasmic calls. Homozygous diploid encodings (`0/0`, `1/1`) are mapped to their corresponding haploid calls; homozygous non-reference multiallelic calls are also mapped to `1`. The genotype matrix is transposed to samples × variants. The script prints alternate-allele carrier counts.

### Autosomal ancestry and sample alignment

`global_ancestry_proportions.csv` must contain `SampleID`, `KhoeSan`, `Bantu.speakingAfrican`, `European`, and `EastAsian`. The three latter components enter regression models, with **KhoeSan omitted as the reference component** to avoid including all compositional proportions simultaneously. Ancestry patterns by cohort and mitochondrial haplogroup are summarised. Phenotype, covariates, ancestry and genotype matrices are explicitly aligned using sample identifiers; only overlapping samples are retained for ancestry-adjusted analyses.

### Mitochondrial macrohaplogroup association

The script compares two binomial logistic models:

```
Null:  TB ~ Age + Sex + Cohort + autosomal ancestry
Full:  TB ~ Macrohaplogroup + Age + Sex + Cohort + autosomal ancestry
```

A likelihood-ratio test evaluates the overall macrohaplogroup term. The script exports model AIC/deviance comparisons, haplogroup model coefficients, odds ratios, approximate 95% confidence intervals, p-values, and variance inflation diagnostics (VIF/GVIF). Here the association variable is constructed from `Macrohaplogroup`, distinct from the `BroadHaplogroup` used in several descriptive figures.

### Variant-level mitochondrial-wide association study (miWAS)

Each mtSNP is evaluated using Firth penalised logistic regression (`logistf`) with variant-specific complete-case data. Variants must be polymorphic and have **at least five alternate-allele carriers** among complete cases. Failed or ineligible models are skipped.

**Primary model:**

```
TB ~ mtSNP + Age + Sex + Cohort
     + Bantu.speakingAfrican + European + EastAsian
```

**Haplogroup-adjusted sensitivity model:**

```
TB ~ mtSNP + Age + Sex + Cohort
     + Bantu.speakingAfrican + European + EastAsian
     + Macrohaplogroup
```

The second model tests whether variant-level associations persist after accounting for mitochondrial lineage. Per-variant outputs include effect size (beta), standard error, odds ratio, profile-likelihood confidence interval, p-value, sample count and carrier/non-carrier counts. Because the sensitivity model includes an additional covariate, its tested variant set may differ from the primary analysis.

### Multiple testing and calibration

Within each model's successfully tested variants, p-values are adjusted using Bonferroni correction and the Benjamini–Hochberg false discovery rate (BH-FDR). The script calculates a Bonferroni threshold of `0.05 / number_of_tested_variants`. It also estimates genomic inflation factor lambda from the median one-degree-of-freedom chi-squared statistic derived from valid p-values and produces QQ plots. These are descriptive calibration diagnostics and should be interpreted with care for linked mitochondrial variants.

## Key outputs

All files below are written by the **main script** to `miWAS_outputs/` (subject to successful execution):

| Output | Description |
|---|---|
| `Age_distribution_cases_controls.png`, `Sex_distribution_cases_controls.png`, `Age_and_sex_distribution.png` | Demographic summaries |
| `Macrohaplogroup_counts.png`, `Macrohaplogroup_distribution_by_phenotype.png`, `Macrohaplogroup_distribution_by_cohort.png` | Lineage composition |
| `Phenotype_distribution_by_cohort.png` | Cohort case/control distribution |
| `mtDNA_PCA_variance_explained.csv`, `mtDNA_PCA_scree_plot.png`, `mtDNA_PCA_by_macrohaplogroup.png` | mtDNA PCA |
| `mtDNA_PC_structure_by_haplogroup_and_cohort.csv` | R² of PC structure models |
| `Global_ancestry_by_cohort.csv`, `Global_ancestry_by_macrohaplogroup.png` | Autosomal ancestry summaries |
| `Haplogroup_model_comparison.csv`, `Haplogroup_association_results.csv`, `Haplogroup_model_VIF.csv` | Haplogroup association and diagnostics; VIF CSV is conditional on returned format |
| `miWAS_primary_results.csv` | Primary mtSNP association results |
| `miWAS_haplogroup_adjusted_sensitivity_results.csv` | Sensitivity results, if models succeed |
| `QQ_primary_miWAS.png`, `QQ_haplogroup_adjusted_sensitivity.png` | QQ plots; sensitivity plot conditional |
| `miWAS_lambda_results.csv` | Genomic inflation factors |

## Running the analysis

Edit `covariate_file`, `ancestry_file`, `vcf_file` and `outdir` at the beginning of `miWAS_cohesive_pipeline.R`, then run from the directory containing the input files:

```bash
Rscript 'scripts/miWAS_cohesive_pipeline.R'
```

For the standalone PCA, update its hard-coded VCF and metadata paths and run:

```bash
Rscript scripts/PCA.R
```

The main script's default VCF filename is `Merged_mtDNA_imputed_geno0.05_mind0.1_INFO0.7_SetHeteroplasmiesToMissing_ExcludeMissingData.recode.vcf`. This name documents intended upstream filtering (including INFO 0.7 and setting heteroplasmies to missing), but **the uploaded scripts do not implement the upstream VCF QC, imputation or haplogroup calling steps**. Preserve the scripts and parameters for those steps separately to reproduce the complete analysis.

## Software

- R (version not specified in the supplied scripts)
- `tidyverse`, `ggplot2`, `ggpubr`, `patchwork`
- `adegenet`, `vcfR`
- `car`, `logistf`, `qqman`
- `RColorBrewer`, `devtools` (loaded by standalone `PCA.R`)

For reproducibility, record package versions using `sessionInfo()` or `renv` and retain the exact input VCF and metadata versions.

## Reproducibility notes

- Sample IDs must agree across VCF, phenotype/covariate and ancestry files. The main pipeline checks duplicates and explicitly matches sample order.
- The script assumes phenotype labels `Case`/`Control` and sex labels `Male`/`Female` exactly as specified.
- Results depend on the final VCF QC, the definition of macrohaplogroups and the available ancestry estimates; those upstream processes are not documented by the two uploaded R scripts.
- The cohort-level Fisher tests with simulated p-values and any stochastic components should use a recorded random seed if exact reruns are required; the main script does not set a global seed.
- Firth regression uses `pl = TRUE`, enabling profile-likelihood confidence intervals. Failed fits are skipped, so inspect run logs and the number of successfully tested variants.
- The main script does not perform functional annotation of associated variants; that analysis, if used, requires additional documentation/scripts.

## References

- Jombart, T. (2008). `adegenet`: a R package for the multivariate analysis of genetic markers. *Bioinformatics*.
- Knaus, B.J. & Grünwald, N.J. (2017). `vcfR`: a package to manipulate and visualize variant call format data in R. *Molecular Ecology Resources*.
- Heinze, G. & Schemper, M. (2002). A solution to the problem of separation in logistic regression. *Statistics in Medicine*.
- Benjamini, Y. & Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. *Journal of the Royal Statistical Society, Series B*.
