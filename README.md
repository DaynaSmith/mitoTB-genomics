# Mitochondrial DNA Variation, Heteroplasmy, and Gene Expression in TB Susceptibility

This repository contains the code, analysis pipelines, and intermediate files associated with my PhD research investigating the role of mitochondrial DNA (mtDNA) variation, heteroplasmy, and gene expression in tuberculosis (TB) susceptibility in South African and broader African populations.

## Overview

Tuberculosis remains a major public health burden in southern Africa, and host genetic factors — including mitochondrial genetic variation — are thought to contribute to differences in individual susceptibility. This project combines mtDNA imputation, haplogroup inference, heteroplasmy characterisation, and mitochondrial gene expression analysis to build a comprehensive picture of how mitochondrial genomic variation relates to TB risk in African populations.

## Research Aims

### Aim 1: Inherited mtDNA variation and TB susceptibility
To investigate the association between inherited mitochondrial DNA variation and TB susceptibility in South African populations.

1. Develop an mtDNA imputation reference panel using publicly available mitochondrial genome sequences.
2. Evaluate the performance of the mtDNA imputation reference panel in southern African populations representing diverse mitochondrial lineages and across different genotyping platforms.
3. Use mtDNA imputation to recover untyped mitochondrial variants from existing South African TB case-control genotype datasets.
4. Infer mitochondrial haplogroups and characterise patterns of mtDNA variation within the TB case-control datasets.
5. Conduct a mitochondrial-wide association study (miWAS) to investigate associations between mtDNA variants and TB susceptibility.
6. Investigate associations between mitochondrial haplogroups and TB susceptibility.
7. Evaluate the predicted functional significance of mtDNA variants identified through the association analyses.

### Aim 2: Mitochondrial heteroplasmy and TB susceptibility
To investigate the relationship between mitochondrial heteroplasmy and TB susceptibility and characterise patterns of heteroplasmic variation in African populations.

1. Identify and characterise heteroplasmic mtDNA variants from whole-genome sequencing data.
2. Investigate whether mitochondrial heteroplasmic burden is associated with TB susceptibility.
3. Investigate the relationship between heteroplasmic burden and demographic and mitochondrial genomic characteristics, including age, mitochondrial DNA copy number, sequencing depth, and mitochondrial haplogroup.
4. Characterise the genomic distribution of heteroplasmic variation and identify recurrent heteroplasmic sites across study populations.
5. Investigate whether patterns of heteroplasmy differ among mitochondrial haplogroups and maternal lineages.
6. Investigate the contribution of local mtDNA sequence context, including homopolymeric regions, to recurrent and lineage-associated patterns of heteroplasmy.

### Aim 3: Mitochondrial gene expression and TB status
To investigate mitochondrial gene expression in relation to TB status and mitochondrial genetic variation.

1. Identify mtDNA-encoded and nuclear-encoded mitochondrial genes in RNA-sequencing datasets from TB cases and controls.
2. Conduct differential gene expression analyses to identify mitochondrial genes associated with TB status.
3. Characterise biological pathways and functional processes represented by differentially expressed mitochondrial genes.
4. Investigate whether mitochondrial gene expression differs according to mitochondrial haplogroup.
5. Evaluate whether previously reported differences in mitochondrial gene expression between African L and non-L mitochondrial haplogroups are replicated in southern African populations.

## Repository Structure

```
.
├── aim1_mtdna_variation/          # Imputation panel development, miWAS, haplogroup association
├── aim2_heteroplasmy/             # Heteroplasmy detection, burden analysis, sequence context
├── aim3_gene_expression/          # RNA-seq processing, differential expression, pathway analysis
├── data/                          # Raw data pointers / metadata (no controlled-access data committed)
├── docs/                          # Supplementary documentation, notes, protocols
├── envs/                          # Conda/renv/environment specification files
├── scripts/                       # Shared utility scripts
└── results/                       # Summary outputs, figures, tables (non-sensitive only)
```

*Note: Controlled-access genotype, sequencing, and clinical data are not stored in this repository; only code, intermediate/derived non-identifiable files, and documentation are included.*

## Requirements

Software and package dependencies for each analysis module are listed in the `envs/` directory (to be added). Broadly, this project makes use of tools for:

- mtDNA sequence retrieval and alignment
- Genotype imputation (mtDNA-specific reference panel construction and evaluation)
- Mitochondrial haplogroup classification
- Association testing (miWAS)
- Heteroplasmy detection from whole-genome sequencing data
- RNA-seq processing and differential expression analysis
- Population genetic and phylogenetic analysis

## Usage

Each aim-specific directory contains its own README with details on pipeline steps, input/output file formats, and how to reproduce individual analyses. Start with the top-level workflow diagram (to be added) for an overview of how the modules connect.

## Data Availability

This project uses publicly available mitochondrial genome sequences alongside South African TB case-control genotype, whole-genome sequencing, and RNA-sequencing datasets. Access to controlled datasets is subject to the relevant data access agreements and ethics approvals; details are provided in `docs/`.

## Contact

Dayna Croock, Stellenbosch University: dayna@sun.ac.za
