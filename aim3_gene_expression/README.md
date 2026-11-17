# Aim 3: Mitochondrial Gene Expression and TB Status

## Overview

This directory contains the computational scripts and reference resources used to investigate mitochondrial gene expression in relation to tuberculosis (TB) status, mitochondrial haplogroup and genetic variation.

The analyses integrate RNA-sequencing (RNA-seq) data processing, mitochondrial genome reconstruction, personalised mitochondrial read mapping, mitochondrial and nuclear gene-expression quantification, differential-expression analysis, expression quantitative trait locus (eQTL) mapping, nuclear–mitochondrial co-expression, functional enrichment and sensitivity analyses.

The workflow includes paired-end (PE; GC6-74) and single-end (SE; TANDEM) RNA-seq datasets. These datasets are generally analysed separately, allowing cross-cohort comparison and replication while accounting for differences in sequencing design and study composition.

Only analytical scripts, documentation and approved reference resources or figures are intended for public distribution. Individual-level sequencing data, clinical metadata, genotype matrices, expression matrices and restricted analytical outputs are excluded in accordance with applicable ethical approvals and data-access agreements.

## Research objectives

1. Identify mtDNA-encoded and nuclear-encoded mitochondrial genes in RNA-sequencing datasets from TB cases and controls.
2. Conduct differential gene-expression analyses to identify mitochondrial genes associated with TB status.
3. Characterise biological pathways and functional processes represented by differentially expressed mitochondrial genes.
4. Investigate whether mitochondrial gene expression differs according to mitochondrial haplogroup.
5. Evaluate whether previously reported differences in mitochondrial gene expression between African L and non-L mitochondrial haplogroups are replicated in southern African populations.

Additional analyses investigate associations between mitochondrial gene expression and mitochondrial or nuclear genetic variation, nuclear–mitochondrial co-expression, and the sensitivity of results to technical and statistical choices.

---

## Repository structure

Scripts are grouped by analytical function. The directory numbers are organisational labels and do not imply that every script must be executed sequentially.

```text
aim3_gene_expression/
├── README.md
│
├── scripts/
│   ├── 01_rnaseq_processing/
│   │   ├── align_PE.pbs
│   │   ├── align_SE.pbs
│   │   ├── mtdna_coverage.pbs
│   │   └── aggregate_qc.py
│   │
│   ├── 02_mitogenome_reconstruction/
│   │   ├── mitgard_run.pbs
│   │   ├── check_mitgard_output.py
│   │   ├── check_haplogroup_quality.py
│   │   ├── personalized_mtdna_remap.pbs
│   │   └── extract_personalised_MT-RNR2.py
│   │
│   ├── 03_variant_calling/
│   │   ├── gatk_rnaseq_varcall.pbs
│   │   ├── gatk_gvcf_call.pbs
│   │   ├── gatk_joint_genotype.pbs
│   │   ├── gatk_genotype_parallel.pbs
│   │   ├── check_gatk_output.py
│   │   ├── build_ok_sample_list.py
│   │   ├── call_mtdna_snps.pbs
│   │   ├── merge_mtdna_snps.pbs
│   │   ├── dedup_ndna_vcf.pbs
│   │   ├── remove_chrm_from_ndna.pbs
│   │   ├── filtering_rna_editing_sites.pbs
│   │   ├── merge_and_filter.pbs
│   │   ├── extract_mtdna_genotypes.pbs
│   │   ├── extract_ndna_genotypes.pbs
│   │   ├── make_mito_target_bed.py
│   │   └── split_bed_for_genotyping.py
│   │
│   ├── 04_expression_quantification/
│   │   ├── aggregate_mtdna_counts.py
│   │   └── aggregate_nuclear_counts.py
│   │
│   ├── 05_metadata_and_qc/
│   │   ├── build_pe_metadata.py
│   │   ├── build_se_metadata.py
│   │   ├── batch_effect_PCA.R
│   │   └── check_TB_covariate_balance.R
│   │
│   ├── 06_differential_expression/
│   │   ├── deseq2_TB_case_control.R
│   │   ├── deseq2_PE_L_vs_nonL.R
│   │   ├── deseq2_SE_L_vs_nonL.R
│   │   ├── deseq2_L_subgroups.R
│   │   └── deseq2_L_subgroups_covariate_adjusted.R
│   │
│   ├── 07_eqtl_analysis/
│   │   ├── eQTL_pipeline.R
│   │   ├── eQTL_pipeline_withSexCovariates.R
│   │   └── annotate_eQTL_hits.R
│   │
│   ├── 08_coexpression_and_enrichment/
│   │   ├── nuclear_coexpression.pbs
│   │   ├── nuclear_coexpression.R
│   │   ├── coexpression_clean.R
│   │   ├── apply_effect_size_filter.R
│   │   └── go_enrichment.R
│   │
│   └── 09_sensitivity_and_validation/
│       └── mtDNA_genomewide_normalisation_sensitivity.R
│
├── resources/
│   ├── chrM_genes.gtf
│   └── mito_nuclear_genes.bed
│
└── figures/
```

---

## Workflow summary

```text
RNA-seq FASTQ files
        |
        v
STAR alignment (PE / SE)
        |
        +-----------------------------+
        |                             |
        v                             v
Alignment QC                   MITGARD reconstruction
        |                             |
        |                      Mitogenome QC
        |                             |
        |                      Haplogroup assignment
        |                             |
        |                      Personalised mtDNA remapping
        |                             |
        v                             v
Nuclear gene counts            Mitochondrial gene counts
        |                             |
        +--------------+--------------+
                       |
                       v
             Unified sample metadata
                       |
          +------------+-------------+
          |                          |
          v                          v
Differential expression        Variant calling
          |                          |
          |                    mtDNA / nDNA genotypes
          |                          |
          |                          v
          |                     eQTL analyses
          |
          +------------+-------------+
                       |
                       v
              Co-expression analysis
                       |
                       v
                GO enrichment

Additional branches:
- Batch-effect PCA
- TB covariate-balance checks
- Normalisation sensitivity
- Targeted MT-RNR2 validation
- Haplogroup subgroup analyses
```

This diagram summarises the analytical relationships. It is not an executable workflow specification.

---

## 1. RNA-seq processing and quality control

### `align_PE.pbs` and `align_SE.pbs`

PBS scripts for RNA-seq alignment of paired-end and single-end datasets.

Separate scripts accommodate differences in sequencing-library structure and alignment configuration.

Alignment outputs support subsequent gene-expression quantification, mitochondrial coverage assessment, mitochondrial genome reconstruction and RNA-derived variant calling.

The original STAR configuration, reference genome, gene annotation and library strandedness should be preserved when reproducing the analyses.

### `mtdna_coverage.pbs`

Generates mitochondrial coverage information from aligned sequencing data.

Coverage information supports evaluation of mitochondrial sequence representation and identification of samples with insufficient coverage for downstream analyses.

### `aggregate_qc.py`

Aggregates STAR alignment statistics and mitochondrial coverage information.

Reported metrics include read counts, uniquely mapped reads, multimapped reads, mitochondrial read counts, mean mitochondrial coverage and coverage breadth at specified depth thresholds.

The script supports identification of samples meeting the required mitochondrial coverage criteria.

---

## 2. Mitochondrial genome reconstruction

### `mitgard_run.pbs`

Runs MITGARD to reconstruct mitochondrial genomes from RNA-seq data.

MITGARD is an external pipeline for mitochondrial genome assembly from transcriptomic sequencing data.

Original repository:

https://github.com/pedronachtigall/MITGARD

The reconstructed mitochondrial sequences provide inputs for mitochondrial haplogroup assignment and personalised mitochondrial mapping.

### `check_mitgard_output.py`

Evaluates mitochondrial assembly outputs, including missing or empty FASTA files, unexpected numbers of sequences and assemblies shorter than a specified proportion of the human mitochondrial reference genome.

The reference mitochondrial genome length is 16,569 bp.

### `check_haplogroup_quality.py`

Summarises mitochondrial haplogroup assignments and quality scores generated using HaploGrep3.

The default quality threshold is 0.90.

The resulting QC information supports sample inclusion and downstream mitochondrial haplogroup comparisons.

### `personalized_mtdna_remap.pbs`

Implements personalised mitochondrial read remapping using sample-specific mitochondrial sequences.

Personalised remapping accounts for differences between individual mitochondrial sequences and the reference genome and provides inputs for mitochondrial expression quantification.

---

## 3. Variant calling and genotype preparation

Mitochondrial and nuclear variants were processed from RNA-seq data for downstream genotype–expression analyses.

### GATK variant-calling scripts

The following scripts implement components of the GATK workflow:

- `gatk_rnaseq_varcall.pbs`
- `gatk_gvcf_call.pbs`
- `gatk_joint_genotype.pbs`
- `gatk_genotype_parallel.pbs`

These represent related variant-calling and genotyping procedures. Their precise invocation order and whether particular scripts constitute alternative implementations should be determined from the PBS commands and original execution records.

### `check_gatk_output.py`

Performs quality control of GATK variant-calling outputs.

Checks include the presence and readability of VCF files, variant counts and PASS-filtered variant counts.

### `build_ok_sample_list.py`

Extracts samples classified as successfully processed from the GATK QC summary.

This sample list supports downstream restriction to successfully genotyped samples.

### Mitochondrial variant-processing scripts

| Script | Function |
|---|---|
| `call_mtdna_snps.pbs` | Mitochondrial variant calling |
| `merge_mtdna_snps.pbs` | Consolidation of mitochondrial variant calls |
| `filtering_rna_editing_sites.pbs` | Filtering of RNA-editing-associated sites |
| `extract_mtdna_genotypes.pbs` | Extraction of mitochondrial genotype data |

### Nuclear variant-processing scripts

| Script | Function |
|---|---|
| `dedup_ndna_vcf.pbs` | Nuclear VCF deduplication |
| `remove_chrm_from_ndna.pbs` | Exclusion of mitochondrial chromosome records |
| `merge_and_filter.pbs` | Variant consolidation and filtering |
| `extract_ndna_genotypes.pbs` | Extraction of nuclear genotype data |
| `make_mito_target_bed.py` | Construction of nuclear-encoded mitochondrial gene intervals |
| `split_bed_for_genotyping.py` | Splitting of target intervals for parallel genotyping |

### Nuclear-encoded mitochondrial gene selection

`make_mito_target_bed.py` uses mitochondrial gene symbols and GENCODE genomic annotations to construct a BED file of nuclear-encoded mitochondrial genes.

The supplied workflow uses the MitoCarta3.0 human mitochondrial gene inventory.

`split_bed_for_genotyping.py` divides the resulting intervals into smaller regions to support parallel genotyping.

### Interpretation of RNA-derived genotypes

Variant calls derived from RNA-seq data may be influenced by expression-dependent coverage, allele-specific expression, RNA editing, reference-mapping bias and related technical effects.

Consequently, RNA-derived genotypes should not be assumed to be equivalent to genotypes derived from matched whole-genome sequencing.

---

## 4. Gene-expression quantification

### `aggregate_mtdna_counts.py`

Aggregates mitochondrial gene counts from personalised mitochondrial remapping outputs.

The supplied workflow combines counts from two alignment passes:

- A main alignment pass allowing up to eight mismatches, used for mitochondrial protein-coding and ribosomal RNA genes.
- A stringent tRNA alignment pass allowing zero mismatches, used for mitochondrial transfer RNA genes.

Counts are extracted from the reverse-stranded column of STAR's `ReadsPerGene.out.tab` files.

The resulting mitochondrial expression matrix contains gene identifiers, gene names, gene types and per-sample counts.

Expected output:

```text
mtdna_expression_matrix.tsv
```

The appropriateness of reverse-stranded counting should be verified against the original library preparation and alignment settings before applying the script to new datasets.

### `aggregate_nuclear_counts.py`

Aggregates genome-wide nuclear gene counts from the original STAR alignment outputs.

This produces the nuclear expression matrix used in downstream co-expression analyses.

Expected output:

```text
nuclear_expression_matrix.tsv
```

Both expression matrices contain sample-level data and are excluded from public distribution.

---

## 5. Metadata integration and exploratory quality control

### `build_pe_metadata.py`

Constructs a unified metadata table for the paired-end dataset by integrating sample identifiers, clinical information, mitochondrial haplogroup assignments, quality scores and exclusion information.

Expected output:

```text
pe_metadata_unified.tsv
```

### `build_se_metadata.py`

Constructs the corresponding metadata table for the single-end dataset.

The script maps eligible TB-only samples to cases and healthy controls to controls for the TB case–control comparison.

Other clinical categories are not assigned to the primary TB case–control comparison, although they may be retained for other analyses where appropriate.

Expected output:

```text
se_metadata_unified.tsv
```

### `batch_effect_PCA.R`

Evaluates whether mitochondrial gene-expression profiles cluster according to sequencing cohort.

The script:

1. Combines mitochondrial expression data from PE and SE samples.
2. Applies a DESeq2 variance-stabilising transformation.
3. Performs principal component analysis.
4. Visualises cohort separation.
5. Evaluates the relationship between cohort and the leading principal components.

This analysis supports assessment of cross-cohort technical differences and the decision to perform separate discovery and replication analyses rather than pooling the datasets without adjustment.

### `check_TB_covariate_balance.R`

Evaluates covariate balance between TB cases and controls.

This analysis informs the interpretation of the primary TB differential-expression model and the potential inclusion of mitochondrial L/non-L status as an additional covariate.

The supplied TB differential-expression script defaults to excluding L status from the model. That setting should be interpreted as a configuration choice, not evidence that an adjusted model was necessarily executed.

---

## 6. Differential gene-expression analyses

Differential-expression analyses were performed using DESeq2.

The analyses investigate mitochondrial gene expression in relation to TB status and mitochondrial haplogroup, with separate cohort analyses and replication where appropriate.

### 6.1 TB case–control differential expression

**Script:** `deseq2_TB_case_control.R`

This script evaluates mitochondrial gene-expression differences between TB cases and controls.

The primary model is:

```r
design = ~ sex + age_centered + condition
```

Here, `condition` distinguishes TB cases from controls, with controls used as the reference group.

The script also permits adjustment for mitochondrial L/non-L status:

```r
design = ~ sex + age_centered + L_status + condition
```

This alternative is controlled by:

```r
INCLUDE_L_STATUS_COVARIATE <- FALSE
```

The supplied configuration therefore uses the sex- and age-adjusted model without L status.

The script runs the analysis separately in the PE and SE cohorts and compares the resulting effect estimates and adjusted p-values.

Statistical significance is evaluated using Benjamini–Hochberg-adjusted p-values, with a threshold of 0.05.

Outputs include cohort-specific differential-expression tables, normalised expression counts and a PE–SE comparison table.

### 6.2 African L versus non-L differential expression: PE cohort

**Script:** `deseq2_PE_L_vs_nonL.R`

This script evaluates mitochondrial gene-expression differences between individuals assigned to African L mitochondrial haplogroups and individuals assigned to non-L haplogroups.

L status is derived from the mitochondrial haplogroup assignment.

Two models are implemented.

**Model A: Full available sample**

```r
design = ~ L_status
```

This unadjusted model maximises the number of eligible samples with mitochondrial expression data.

**Model B: Clinically restricted sample**

```r
design = ~ sex + age_centered + L_status
```

This model uses samples with the required clinical and demographic information.

The two analyses allow assessment of whether L/non-L expression differences are sensitive to sample restriction and covariate adjustment.

The non-L group is the reference category.

### 6.3 African L versus non-L differential expression: SE replication

**Script:** `deseq2_SE_L_vs_nonL.R`

Applies corresponding L/non-L differential-expression analyses to the SE cohort.

Both unadjusted and sex- and age-adjusted models are implemented.

The script also compares the SE results with the PE discovery results, including effect direction and adjusted p-values.

This analysis directly addresses replication of previously reported mitochondrial expression differences between African L and non-L lineages.

### 6.4 Mitochondrial haplogroup subgroup analysis

**Script:** `deseq2_L_subgroups.R`

Extends the binary L/non-L comparison to individual African L subgroups.

In the PE dataset, the script evaluates:

- L0 versus non-L
- L1 versus non-L
- L2 versus non-L
- L3 versus non-L
- L4 versus non-L

L4 is treated as exploratory because of its smaller sample size.

L5 and L6 are excluded from formal testing because the available subgroup sample sizes are insufficient.

In the SE dataset, formal subgroup testing is restricted to L0 versus non-L because other L subgroups have very small sample sizes.

The script uses:

```r
design = ~ L_subgroup
```

It implements an omnibus likelihood-ratio test and individual subgroup contrasts against the non-L reference group.

### 6.5 Covariate-adjusted haplogroup subgroup analysis

**Script:** `deseq2_L_subgroups_covariate_adjusted.R`

Provides an additional subgroup differential-expression analysis incorporating covariate adjustment.

This analysis is intended to assess whether subgroup-associated expression differences persist after accounting for available demographic or clinical covariates.

The exact model formula and eligibility criteria should be retained from the script rather than inferred from the unadjusted subgroup analysis.

### Interpretation

Haplogroup-associated expression differences should be interpreted in the context of mitochondrial phylogeny, population structure, uneven subgroup sample sizes and possible technical differences between cohorts.

The L/non-L grouping represents broad maternal-lineage structure and should not be interpreted as a direct measure of ancestry independent of other genetic and demographic characteristics.

---

## 7. Expression quantitative trait locus analyses

The eQTL analyses investigate associations between genetic variation and mitochondrial gene expression.

Two distinct genotype classes are considered:

1. Mitochondrial variants.
2. Nuclear variants.

The primary combined workflow is implemented in `analysis1_eqtl_clean.R`.

### 7.1 Primary combined eQTL workflow

**Script:** `eQTL_pipeline.R`

This script implements mitochondrial and nuclear eQTL analyses using Matrix eQTL, with separate PE and SE analyses.

#### Expression preprocessing

Mitochondrial expression counts are normalised using DESeq2 size factors.

Normalised counts are transformed as:

```text
log2(normalised count + 1)
```

Genes must be detected in at least 90% of samples within the relevant cohort.

#### Genotype quality control

The supplied script applies:

| Parameter | Threshold |
|---|---|
| Minimum minor allele count | 10 |
| Maximum major allele frequency for mtDNA | 95% |
| Minimum expression presence | 90% |
| Number of nuclear genotype PCs | 5 |
| Multiple-testing correction | Bonferroni |
| Significance threshold | 0.05 after correction |

#### Mitochondrial eQTL model

Mitochondrial genotype–expression associations are evaluated using the categorical Matrix eQTL model:

```r
useModel = modelANOVA
```

This treats genotype categories as groups rather than imposing a linear dosage relationship.

The approach is appropriate for evaluating categorical allele assignments in the supplied mitochondrial genotype representation.

The script does not apply nuclear genotype principal components to the mitochondrial eQTL model.

#### Nuclear eQTL model

Nuclear genotype–expression associations are evaluated using the additive linear model:

```r
useModel = modelLINEAR
```

Nuclear genotypes are represented using diploid allele dosages.

Sex and nuclear genotype principal components are included as available covariates.

The principal components are derived from filtered nuclear genotype data and provide a proxy for genetic population structure.

#### Multiple testing and replication

The script calculates the Bonferroni significance threshold from the number of tested genotype–gene pairs.

Analyses are conducted separately in PE and SE samples.

Discovery associations identified in PE are compared against corresponding associations in SE, with replication assessed using the cohort-specific significance indicators.

The script also exports genotype principal components for downstream robustness analyses.

#### Deviations from the comparison study

The script documents several differences between the original comparison study and the available datasets:

- Matched DNA-derived mitochondrial copy-number measurements were unavailable.
- Laboratory-of-origin covariates were unavailable.
- Nuclear genotype principal components were used as an ancestry-structure proxy in the nuclear model.
- The nuclear mitochondrial-gene selection was based on MitoCarta3.0 rather than the narrower candidate-gene list used in the comparison study.
- RNA-derived nuclear variants were not restricted to variants catalogued in the 1000 Genomes Project.

These deviations are important when comparing results across studies.

**Implementation caution:** The script's minor allele count and major allele frequency calculations use numeric genotype coding and a ploidy-based denominator. Multiallelic mitochondrial sites require particular care, because these calculations may not fully capture allele frequencies across all alternate alleles. The genotype representation and filtering logic should be verified before interpreting multiallelic associations.

### 7.2 Nuclear eQTL refinement

**Script:** `eQTL_pipeline_withSexCovariates.R`

Implements a separate nuclear eQTL workflow incorporating stringent nuclear genotype QC, NUMT-related exclusions, sex covariates, genetic principal components, discovery–replication comparisons and mitochondrial haplogroup sensitivity analyses.

This script does not rerun the mitochondrial eQTL analysis.

It should therefore be documented as a complementary nuclear eQTL analysis rather than a replacement for the combined mitochondrial and nuclear workflow.

### 7.3 Annotation of replicating eQTL associations

**Script:** `annotate_eQTL_hits.R`

Annotates genotype–expression associations identified in the combined eQTL workflow.

The script distinguishes associations that were significant in PE and merely present in SE from those meeting the significance criterion in both cohorts.

It maps variant coordinates to mitochondrial gene annotations or nuclear mitochondrial-gene intervals.

The annotation workflow uses:

```text
chrM_genes.gtf
mito_nuclear_genes.bed
mtdna_expression_matrix.tsv
```

It also reconciles mitochondrial reference-sequence naming conventions when interpreting mitochondrial variant identifiers.

The supplied script expects the combined eQTL outputs in the `analysis1_eqtl_v3` output directory, which is also the default output directory used by `analysis1_eqtl_clean.R`.

---

## 8. Nuclear–mitochondrial co-expression

Co-expression analyses investigate relationships between mitochondrial gene expression and nuclear gene-expression profiles.

### `nuclear_coexpression.pbs`

PBS script associated with the nuclear–mitochondrial co-expression workflow.

The precise commands and software dependencies are defined in the PBS script.

### `nuclear_coexpression.R`

Performs genome-wide nuclear–mitochondrial co-expression analysis.

The script:

1. Loads mitochondrial and nuclear expression matrices.
2. Normalises expression data.
3. Calculates Pearson correlations between mitochondrial and nuclear gene-expression profiles.
4. Performs analyses separately in PE and SE samples.
5. Evaluates replication between cohorts.

False-discovery-rate correction is applied to identify statistically significant associations.

### `apply_effect_size_filter.R`

Applies an additional effect-size requirement to co-expression replication results.

The supplied script requires an absolute Pearson correlation coefficient of at least 0.5 in both PE and SE for the principal filtered analysis.

It also evaluates alternative correlation thresholds from 0.1 to 0.5.

This script filters existing results rather than recalculating correlations.

### `coexpression_clean.R`

Implements a separate integrated analysis linking mitochondrial–nuclear co-expression to mitochondrial haplogroup-associated nuclear expression.

The analysis comprises three stages.

**Stage 1: Co-expression screening**

Pearson correlations are calculated between 15 mitochondrial mRNA/rRNA genes and nuclear genes.

Bonferroni correction is applied based on the number of tested associations.

Significant co-expression is required in both PE and SE datasets.

**Stage 2: L/non-L nuclear expression analysis**

Nuclear genes passing the co-expression stage are tested for expression differences between African L and non-L mitochondrial haplogroup backgrounds.

The analysis uses Matrix eQTL with sex as a covariate.

**Stage 3: Functional enrichment**

Genes meeting the relevant co-expression and differential-expression criteria are assessed for functional enrichment.

This script should be distinguished from `nuclear_coexpression.R` because the analyses differ in their gene-selection strategy, significance criteria and downstream objectives.

---

## 9. Functional enrichment

### `go_enrichment.R`

Performs Gene Ontology enrichment analysis using the g:Profiler API.

The script accepts a filtered list of nuclear genes identified through the co-expression analysis.

Three Gene Ontology domains are queried:

- Biological Process (`GO:BP`)
- Molecular Function (`GO:MF`)
- Cellular Component (`GO:CC`)

The script uses `wget` to communicate with the g:Profiler API and the R package `jsonlite` to process the returned results.

An internet connection is required.

The functional enrichment implemented within `coexpression_clean.R` represents a related but distinct analytical branch.

Enrichment results should be interpreted with reference to the precise input gene list and statistical background used in the relevant analysis.

---

## 10. Sensitivity and validation analyses

### `mtDNA_genomewide_normalisation_sensitivity.R`

Evaluates the sensitivity of mitochondrial gene-expression analyses to the normalisation strategy.

The script compares:

1. Mitochondrial-only expression normalisation.
2. Normalisation using size factors derived from genome-wide STAR gene counts.

This comparison assesses whether the choice of normalisation approach influences mitochondrial expression results.

The script uses personalised mitochondrial expression counts and original STAR nuclear gene-count outputs.

The supplied implementation contains absolute HPC paths that must be adapted before reuse.

Its SE processing configuration should also be checked against the original alignment and library strandedness settings.

### Additional validation and diagnostic analyses

The following scripts provide complementary checks:

| Script | Purpose |
|---|---|
| `batch_effect_PCA.R` | Evaluates cohort-associated expression structure |
| `check_TB_covariate_balance.R` | Evaluates covariate balance in TB case–control comparisons |
| `extract_personalised_MT-RNR2.py` | Validates personalised mitochondrial sequence extraction |
| `check_mitgard_output.py` | Evaluates mitochondrial assembly completeness |
| `check_haplogroup_quality.py` | Evaluates haplogroup assignment quality |
| `check_gatk_output.py` | Evaluates variant-calling outputs |

---

## 11. Reference resources

### `chrM_genes.gtf`

Mitochondrial gene annotation used for mitochondrial expression quantification and variant annotation.

The annotation must be compatible with the reference mitochondrial genome and chromosome naming conventions used in the analysis.

### `mito_nuclear_genes.bed`

Genomic intervals corresponding to nuclear-encoded mitochondrial genes.

The supplied workflow generates these intervals from MitoCarta3.0 gene symbols and GENCODE gene annotations.

### External references

- MITGARD: https://github.com/pedronachtigall/MITGARD
- Human MitoCarta3.0: https://www.broadinstitute.org/mitocarta/mitocarta30-inventory-mammalian-mitochondrial-proteins-and-pathways
- GENCODE: https://www.gencodegenes.org/
- HaploGrep: https://haplogrep.i-med.ac.at/
- g:Profiler: https://biit.cs.ut.ee/gprofiler/

---

## 12. Software and computational environment

The analyses were developed for execution in a high-performance computing environment using PBS job submission.

The supplied scripts use or reference the following software and packages.

| Software | Application |
|---|---|
| STAR | RNA-seq alignment and gene quantification |
| MITGARD | Mitochondrial genome reconstruction |
| HaploGrep3 | Mitochondrial haplogroup assignment |
| GATK | RNA-seq variant calling and genotyping |
| SAMtools / BCFtools | Alignment and variant-file processing |
| BEDTools | Genomic interval processing |
| Python | Workflow automation and data processing |
| Biopython | Sequence alignment and extraction |
| R | Statistical analysis |
| DESeq2 | Normalisation and differential expression |
| MatrixEQTL | Genotype–expression association testing |
| ggplot2 | Visualisation |
| g:Profiler | Functional enrichment |
| jsonlite | JSON processing for enrichment results |

Exact versions should be obtained from the original analysis environment or recorded module configurations.

### PBS execution

PBS scripts are submitted using the cluster scheduler, for example:

```bash
qsub align_PE.pbs
```

The submission commands, resource requests and software modules may require adaptation to another HPC environment.

### Python execution

Python scripts are generally executed after their required upstream files have been generated:

```bash
python aggregate_qc.py --help
```

Where supported, command-line arguments can be inspected using `--help`.

### R execution

Statistical analyses can be run using:

```bash
Rscript deseq2_TB_case_control.R
Rscript deseq2_PE_L_vs_nonL.R
Rscript deseq2_SE_L_vs_nonL.R
Rscript analysis1_eqtl_clean.R
```

These commands assume that the configured input paths are valid and all required dependencies have been installed.

The SE L/non-L replication script additionally requires the PE results produced by `deseq2_PE_L_vs_nonL.R`.

---

## 13. Execution dependencies

The workflow contains several branches that share common upstream inputs.

| Analysis | Required upstream products |
|---|---|
| Mitochondrial assembly QC | MITGARD assembly outputs |
| Haplogroup QC | Haplogroup assignment outputs |
| Personalised mtDNA mapping | Reconstructed mitochondrial sequences and RNA-seq alignments |
| Mitochondrial count aggregation | Personalised remapping outputs and mitochondrial annotation |
| Nuclear count aggregation | Original STAR gene-count outputs |
| Metadata integration | Clinical metadata, manifests, haplogroup assignments and exclusion information |
| Differential expression | Mitochondrial expression matrix and cohort metadata |
| Combined eQTL analysis | Mitochondrial expression matrix, mitochondrial and nuclear genotype matrices, cohort metadata |
| Nuclear eQTL refinement | Strict-QC nuclear genotype matrix, mitochondrial expression matrix and metadata |
| Co-expression | Mitochondrial and nuclear expression matrices |
| Effect-size filtering | Existing co-expression replication results |
| Functional enrichment | Filtered gene lists from upstream analyses |
| Normalisation sensitivity | Personalised mitochondrial counts and original STAR gene counts |

The repository is not configured as a single executable workflow manager.

Where scripts implement alternative or supplementary approaches, users should consult the relevant analytical section rather than execute every script indiscriminately.

---

## 14. Reproducibility and methodological considerations

### Sample identifiers and restricted intermediate files

Several scripts expect files such as:

```text
bam_manifest.txt
mtdna_expression_matrix.tsv
nuclear_expression_matrix.tsv
mtdna_genotype_matrix.tsv
ndna_genotype_matrix.tsv
pe_metadata_unified.tsv
se_metadata_unified.tsv
gatk_qc_summary.tsv
```

These are examples of intermediate files required by the workflow.

They are not publicly distributed because they may contain individual-level or restricted information.

### File paths

Some scripts contain absolute HPC paths, while others expect input files in the current working directory.

Paths must be adapted before reuse.

After moving scripts into the proposed subdirectories, any relative references should be reviewed.

### Multiple analytical implementations

The repository includes related but non-identical implementations of several analyses.

In particular:

- `analysis1_eqtl_clean.R` is the selected primary combined mitochondrial and nuclear eQTL script.
- `Final_eQTL_pipeline.R` is a separate nuclear eQTL refinement workflow.
- `nuclear_coexpression.R` and `coexpression_clean.R` implement different co-expression strategies.
- `apply_effect_size_filter.R` performs downstream filtering rather than recalculating correlations.
- The GATK PBS scripts include multiple variant-calling and genotyping procedures.

These scripts should not be assumed to produce identical results.

### Cohort differences

The PE and SE datasets differ in sequencing-library structure and may differ in sample composition and clinical characteristics.

Cross-cohort replication is therefore interpreted with consideration of technical and biological heterogeneity.

### Multiple testing

Different analytical branches use different multiple-testing approaches.

These include:

- Benjamini–Hochberg false-discovery-rate correction for differential-expression analyses.
- Bonferroni correction for the primary combined eQTL workflow.
- False-discovery-rate correction for genome-wide co-expression screening.
- Bonferroni correction in the integrated co-expression workflow.

Significance thresholds and correction methods should be reported separately for each analysis.

### Mitochondrial population structure

Mitochondrial variants are phylogenetically linked and mitochondrial haplogroups represent shared maternal ancestry.

Associations involving mitochondrial genetic variation may therefore reflect broader haplotypic or population structure rather than independent variant-specific effects.

Haplogroup sensitivity analyses and cross-cohort comparisons provide important context for interpretation.

---

## 15. Data availability and ethical considerations

The analyses described in this repository were conducted using human RNA-sequencing and associated clinical and genetic data.

Individual-level sequencing files, clinical metadata, sample identifiers, genotype matrices, expression matrices and restricted intermediate results are not publicly distributed.

The repository is intended to contain only analytical scripts, documentation and reference resources or figures approved for public release, in accordance with applicable Human Research Ethics Committee (HREC) approvals and data-access agreements.

Access to underlying datasets is subject to the relevant study-specific data-access procedures and ethical requirements.

The absence of restricted input data means that the repository does not independently support full reproduction of the analyses without authorised access to the original datasets.

---

## 16. Acknowledgements and citation

The analyses were conducted as part of doctoral research in mitochondrial genetics and tuberculosis susceptibility at Stellenbosch University, South Africa.

The computational workflows use established open-source bioinformatics software and externally developed reference resources.

Users should cite the relevant software and resource publications when reproducing or adapting the analyses.

### MITGARD

Nachtigall PG, Grazziotin FG, Junqueira-de-Azevedo ILM. 2021. MITGARD: an automated pipeline for mitochondrial genome assembly in eukaryotic species using RNA-seq data. *Briefings in Bioinformatics*. 22(5):bbaa429.

https://doi.org/10.1093/bib/bbaa429

Original software repository:

https://github.com/pedronachtigall/MITGARD

Additional relevant tools include STAR, GATK, HaploGrep, MitoCarta3.0, DESeq2, Matrix eQTL and g:Profiler.

---

## 17. Scope and limitations

This repository documents the computational methods and statistical analyses associated with Aim 3 of the dissertation.

It is intended to improve transparency, traceability and reproducibility while respecting restrictions on human genetic and clinical data sharing.

The scripts were developed for a specific HPC environment and may require changes to paths, software installations and reference resources before execution elsewhere.

The public repository does not include restricted datasets or guarantee complete reproduction without access to the original computational environment and approved data.

All reference files, figures and scripts should be reviewed for data-sharing compliance before public release.
