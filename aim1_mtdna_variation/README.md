# Aim 1: Inherited mtDNA Variation and TB Susceptibility

**Aim:** To investigate the association between inherited mitochondrial DNA variation and TB susceptibility in South African populations.

This aim is organised into three stages, each with its own directory:

| Stage | Directory | Description |
|-------|-----------|-------------|
| 1. Reference panel construction | [`01_reference_panel/`](01_reference_panel) | Building an mtDNA imputation reference panel from publicly available mitochondrial genome sequences |
| 2. Reference panel validation | [`02_panel_validation/`](02_panel_validation) | Evaluating imputation performance across southern African populations, mitochondrial lineages, and genotyping platforms |
| 3. Application: miWAS | [`03_miwas/`](03_miwas) | Applying the validated panel to impute untyped variants in South African TB case-control datasets and testing for association with TB susceptibility |

Haplogroup inference and characterisation of mtDNA variation within the case-control datasets, haplogroup-level association testing, and functional annotation of significant variants are carried out as part of the `03_miwas/` workflow (see that directory's README for details).

## Workflow overview

```
Reference panel construction  →  Reference panel validation  →  Imputation + miWAS
     (01_reference_panel)          (02_panel_validation)           (03_miwas)
```

Each stage's directory README documents the specific methods, scripts, inputs, and outputs for that stage.
