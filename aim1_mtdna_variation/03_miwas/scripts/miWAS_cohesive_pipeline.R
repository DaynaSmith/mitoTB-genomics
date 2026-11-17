# =============================================================================
# Cohort characterisation, mtDNA population structure, haplogroup association,
# and mitochondrial-wide association analysis (miWAS) of TB susceptibility
#
# Primary miWAS:
#   TB ~ mtSNP + Age + Sex + Cohort + global autosomal ancestry
#
# Sensitivity miWAS:
#   TB ~ mtSNP + Age + Sex + Cohort + global autosomal ancestry
#        + mitochondrial macrohaplogroup
#
# Notes:
#   - mtDNA genotypes are treated as haploid for association testing:
#       0 = reference allele
#       1 = alternate allele
#       NA = missing/uncertain heterozygous call
#   - KhoeSan is the omitted/reference autosomal ancestry component.
#   - Sample alignment between covariates and genotype matrices is always
#     performed explicitly by sample ID.
# =============================================================================


# =============================================================================
# 0. PACKAGES AND FILE PATHS
# =============================================================================

library(tidyverse)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(adegenet)
library(vcfR)
library(car)
library(logistf)
library(qqman)

covariate_file <- "sample_data.csv"
ancestry_file  <- "global_ancestry_proportions.csv"
vcf_file       <- "Merged_mtDNA_imputed_geno0.05_mind0.1_INFO0.7_SetHeteroplasmiesToMissing_ExcludeMissingData.recode.vcf"

# Create output directory if required
outdir <- "miWAS_outputs"
dir.create(outdir, showWarnings = FALSE)


# =============================================================================
# 1. COLOUR PALETTES
# =============================================================================

haplo_colors <- c(
  # L0a branch
  "L0a'g" = "#F6D89A",
  "L0a"    = "#DFC27D",
  "L0a2"   = "#BF812D",

  # L0d branch
  "L0d1" = "#E9A03B",
  "L0d2" = "#F5C46B",
  "L0d3" = "#FEE090",

  # Other L0
  "L0f" = "#BFA26A",
  "L0k" = "#AEBF7A",

  # L1 branch
  "L1b" = "#4D9221",
  "L1c" = "#2B7A0B",

  # L2 branch
  "L2b" = "#66C2A4",
  "L2c" = "#31A354",

  # L3 branch
  "L3"  = "#35978F",
  "L3b" = "#80CDC1",
  "L3d" = "#1C9099",
  "L3e" = "#7BC8C3",
  "L3f" = "#9EDBD7",
  "L3h" = "#4EB3A5",

  # L4 / L5
  "L4b" = "#2B8CBE",
  "L5a" = "#756BB1",

  # Eurasian lineages
  "M" = "#C475E6",
  "N" = "#E67597",
  "R" = "#C51B7D"
)

ancestry_colors <- c(
  "KhoeSan" = "#D92300",
  "Bantu.speakingAfrican" = "#12A802",
  "European" = "#875DC9",
  "EastAsian" = "#2282BF"
)


# =============================================================================
# 2. LOAD AND FORMAT SAMPLE-LEVEL DATA
# =============================================================================

covariate_data <- read.csv(
  covariate_file,
  header = TRUE,
  stringsAsFactors = FALSE
)

required_covariates <- c(
  "SampleID", "Age", "Sex", "Phenotype", "Cohort",
  "Haplogroup", "Macrohaplogroup", "BroadHaplogroup"
)

stopifnot(all(required_covariates %in% colnames(covariate_data)))
stopifnot(!anyDuplicated(covariate_data$SampleID))

# Keep human-readable factor versions for descriptive analyses
covariate_data$Phenotype <- factor(
  covariate_data$Phenotype,
  levels = c("Control", "Case")
)

covariate_data$Sex <- factor(
  covariate_data$Sex,
  levels = c("Male", "Female")
)

covariate_data$Cohort <- factor(covariate_data$Cohort)
covariate_data$BroadHaplogroup <- factor(covariate_data$BroadHaplogroup)
covariate_data$Macrohaplogroup <- factor(covariate_data$Macrohaplogroup)


# =============================================================================
# 3. DESCRIPTIVE CHARACTERISATION OF THE ANALYTICAL COHORT
# =============================================================================

cat("\n=== Sample counts by phenotype ===\n")
print(table(covariate_data$Phenotype, useNA = "ifany"))

cat("\n=== Sample counts by sex ===\n")
print(table(covariate_data$Sex, useNA = "ifany"))

cat("\n=== Sample counts by cohort ===\n")
print(table(covariate_data$Cohort, useNA = "ifany"))

cat("\n=== Age summary: full cohort ===\n")
print(summary(covariate_data$Age))

cat("\n=== Age summary: cases ===\n")
print(summary(covariate_data$Age[covariate_data$Phenotype == "Case"]))

cat("\n=== Age summary: controls ===\n")
print(summary(covariate_data$Age[covariate_data$Phenotype == "Control"]))


# -----------------------------------------------------------------------------
# 3.1 Age differences between cases and controls
# -----------------------------------------------------------------------------

shapiro_case <- shapiro.test(
  covariate_data$Age[covariate_data$Phenotype == "Case"]
)

shapiro_control <- shapiro.test(
  covariate_data$Age[covariate_data$Phenotype == "Control"]
)

cat("\nShapiro-Wilk p-value, cases:", shapiro_case$p.value, "\n")
cat("Shapiro-Wilk p-value, controls:", shapiro_control$p.value, "\n")

age_test <- wilcox.test(
  Age ~ Phenotype,
  data = covariate_data,
  exact = FALSE
)

print(age_test)

age_plot <- ggplot(
  covariate_data,
  aes(x = Phenotype, y = Age)
) +
  geom_boxplot(fill = "grey") +
  theme_classic() +
  labs(
    x = "TB status",
    y = "Age (years)"
  )

ggsave(
  file.path(outdir, "Age_distribution_cases_controls.png"),
  plot = age_plot,
  width = 5,
  height = 5,
  dpi = 300
)


# -----------------------------------------------------------------------------
# 3.2 Sex distribution between cases and controls
# -----------------------------------------------------------------------------

sex_table <- table(
  covariate_data$Sex,
  covariate_data$Phenotype
)

sex_test <- chisq.test(sex_table)

cat("\n=== Sex by phenotype ===\n")
print(sex_table)
print(sex_test)

sex_counts <- covariate_data %>%
  count(Phenotype, Sex, name = "Count")

sex_plot <- ggplot(
  covariate_data,
  aes(x = Phenotype, fill = Sex)
) +
  geom_bar(position = "dodge") +
  geom_text(
    data = sex_counts,
    aes(y = Count, label = Count, group = Sex),
    position = position_dodge(width = 0.9),
    vjust = -0.4
  ) +
  theme_classic() +
  labs(
    x = "TB status",
    y = "Number of individuals",
    fill = "Sex"
  ) +
  theme(legend.position = "top")

ggsave(
  file.path(outdir, "Sex_distribution_cases_controls.png"),
  plot = sex_plot,
  width = 7,
  height = 5,
  dpi = 300
)

ggsave(
  file.path(outdir, "Age_and_sex_distribution.png"),
  plot = sex_plot | age_plot,
  width = 10,
  height = 5,
  dpi = 300
)


# -----------------------------------------------------------------------------
# 3.3 Mitochondrial macrohaplogroup distributions
# -----------------------------------------------------------------------------

haplogroup_counts <- covariate_data %>%
  count(BroadHaplogroup, name = "Count")

haplogroup_plot <- ggplot(
  haplogroup_counts,
  aes(x = BroadHaplogroup, y = Count, fill = BroadHaplogroup)
) +
  geom_col() +
  geom_text(aes(label = Count), vjust = -0.4) +
  scale_fill_manual(values = haplo_colors) +
  theme_classic() +
  theme(legend.position = "none") +
  labs(
    x = "Mitochondrial macrohaplogroup",
    y = "Number of individuals"
  )

ggsave(
  file.path(outdir, "Macrohaplogroup_counts.png"),
  plot = haplogroup_plot,
  width = 10,
  height = 6,
  dpi = 300
)

# Haplogroup distribution by phenotype
haplo_phenotype_table <- table(
  covariate_data$Macrohaplogroup,
  covariate_data$Phenotype
)

haplo_phenotype_test <- fisher.test(
  haplo_phenotype_table,
  simulate.p.value = TRUE,
  B = 2000
)

cat("\n=== Macrohaplogroup distribution by phenotype ===\n")
print(haplo_phenotype_table)
print(haplo_phenotype_test)

haplo_by_phenotype <- ggplot(
  covariate_data,
  aes(x = Phenotype, fill = BroadHaplogroup)
) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = haplo_colors) +
  theme_classic() +
  labs(
    x = "TB status",
    y = "Proportion",
    fill = "Mitochondrial macrohaplogroup"
  )

ggsave(
  file.path(outdir, "Macrohaplogroup_distribution_by_phenotype.png"),
  plot = haplo_by_phenotype,
  width = 10,
  height = 7,
  dpi = 300
)

# Haplogroup distribution by cohort
haplo_cohort_table <- table(
  covariate_data$Macrohaplogroup,
  covariate_data$Cohort
)

haplo_cohort_test <- fisher.test(
  haplo_cohort_table,
  simulate.p.value = TRUE,
  B = 2000
)

cat("\n=== Macrohaplogroup distribution by cohort ===\n")
print(haplo_cohort_table)
print(haplo_cohort_test)

haplo_by_cohort <- ggplot(
  covariate_data,
  aes(x = Cohort, fill = BroadHaplogroup)
) +
  geom_bar(position = "fill") +
  scale_fill_manual(values = haplo_colors) +
  theme_classic() +
  labs(
    x = "Study cohort",
    y = "Proportion",
    fill = "Mitochondrial macrohaplogroup"
  )

ggsave(
  file.path(outdir, "Macrohaplogroup_distribution_by_cohort.png"),
  plot = haplo_by_cohort,
  width = 10,
  height = 7,
  dpi = 300
)


# -----------------------------------------------------------------------------
# 3.4 Case/control distribution across cohorts
# -----------------------------------------------------------------------------

case_control_table <- table(
  covariate_data$Phenotype,
  covariate_data$Cohort
)

case_control_test <- fisher.test(case_control_table)

cat("\n=== Phenotype distribution by cohort ===\n")
print(case_control_table)
print(case_control_test)

phenotype_by_cohort <- ggplot(
  covariate_data,
  aes(x = Cohort, fill = Phenotype)
) +
  geom_bar(position = "dodge") +
  theme_classic() +
  labs(
    x = "Study cohort",
    y = "Number of individuals",
    fill = "TB status"
  ) +
  theme(legend.position = "top")

ggsave(
  file.path(outdir, "Phenotype_distribution_by_cohort.png"),
  plot = phenotype_by_cohort,
  width = 7,
  height = 5,
  dpi = 300
)


# =============================================================================
# 4. LOAD FINAL IMPUTED mtDNA VCF
# =============================================================================

vcf <- read.vcfR(vcf_file, verbose = FALSE)

variant_info <- data.frame(
  CHROM = getCHROM(vcf),
  POS   = getPOS(vcf),
  REF   = getREF(vcf),
  ALT   = getALT(vcf),
  stringsAsFactors = FALSE
)

cat("\nVariants in final mtDNA VCF:", nrow(variant_info), "\n")
cat("Samples in final mtDNA VCF:", ncol(vcf@gt) - 1, "\n")


# =============================================================================
# 5. mtDNA PRINCIPAL COMPONENT ANALYSIS
# =============================================================================

genlight_obj <- vcfR2genlight(vcf)
pca_res <- glPca(genlight_obj, nf = 10)

# Important:
# glPca() returns eigenvalues for retained axes. Percentage variance is
# calculated relative to the eigenvalues available in pca_res.
variance_explained <- 100 * pca_res$eig / sum(pca_res$eig)

pca_variance <- data.frame(
  PC = paste0("PC", seq_along(variance_explained)),
  Variance = variance_explained
)

write.csv(
  pca_variance,
  file.path(outdir, "mtDNA_PCA_variance_explained.csv"),
  row.names = FALSE
)

scree_plot <- ggplot(
  pca_variance,
  aes(x = factor(PC, levels = PC), y = Variance, group = 1)
) +
  geom_line() +
  geom_point(size = 2) +
  theme_classic() +
  labs(
    x = "Principal component",
    y = "Percentage of retained PCA variance"
  )

ggsave(
  file.path(outdir, "mtDNA_PCA_scree_plot.png"),
  plot = scree_plot,
  width = 8,
  height = 5,
  dpi = 300
)

pca_scores <- as.data.frame(pca_res$scores)
pca_scores$Sample <- rownames(pca_res$scores)

# Rename first 10 score columns explicitly
score_cols <- setdiff(colnames(pca_scores), "Sample")
colnames(pca_scores)[match(score_cols, colnames(pca_scores))] <-
  paste0("PC", seq_along(score_cols))

# Merge PCA scores with sample metadata using sample IDs
pca_metadata <- covariate_data %>%
  transmute(
    Sample = SampleID,
    Haplogroup = Haplogroup,
    Macrohaplogroup = BroadHaplogroup,
    Cohort = Cohort
  )

pca_df <- merge(
  pca_scores,
  pca_metadata,
  by = "Sample",
  all.x = TRUE,
  sort = FALSE
)

stopifnot(!any(is.na(pca_df$Macrohaplogroup)))
stopifnot(!any(is.na(pca_df$Cohort)))

# PC plots
pc1_var <- round(variance_explained[1], 2)
pc2_var <- round(variance_explained[2], 2)
pc3_var <- round(variance_explained[3], 2)

PC1_2_haps <- ggplot(
  pca_df,
  aes(x = PC1, y = PC2, colour = Macrohaplogroup)
) +
  geom_point(size = 2.5, alpha = 0.8) +
  scale_colour_manual(values = haplo_colors) +
  theme_classic() +
  labs(
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC2 (", pc2_var, "%)"),
    colour = "Macrohaplogroup"
  ) + theme(legend.position = "none")

PC2_3_haps <- ggplot(
  pca_df,
  aes(x = PC2, y = PC3, colour = Macrohaplogroup)
) +
  geom_point(size = 2.5, alpha = 0.8) +
  scale_colour_manual(values = haplo_colors) +
  theme_classic() +
  labs(
    x = paste0("PC2 (", pc2_var, "%)"),
    y = paste0("PC3 (", pc3_var, "%)"),
    colour = "Macrohaplogroup"
  )+ theme(legend.position = "none")

PC1_3_haps <- ggplot(
  pca_df,
  aes(x = PC1, y = PC3, colour = Macrohaplogroup)
) +
  geom_point(size = 2.5, alpha = 0.8) +
  scale_colour_manual(values = haplo_colors) +
  theme_classic() +
  labs(
    x = paste0("PC1 (", pc1_var, "%)"),
    y = paste0("PC3 (", pc3_var, "%)"),
    colour = "Macrohaplogroup"
  )

combined_pca_plot <- (PC1_2_haps | PC2_3_haps) / PC1_3_haps

ggsave(
  file.path(outdir, "mtDNA_PCA_by_macrohaplogroup.png"),
  plot = combined_pca_plot,
  width = 13,
  height = 9,
  dpi = 300
)


# =============================================================================
# 6. ASSESS WHETHER mtDNA PCs TRACK HAPLOGROUP OR COHORT
# =============================================================================
#
# IMPORTANT FIX:
# Models are fitted using pca_df itself so the phenotype/covariate order
# is guaranteed to match the PC scores. Do not refer back to covariate_data
# by row position here.

pc_structure_results <- lapply(1:5, function(i) {

  pc_name <- paste0("PC", i)

  fit_haplo <- lm(
    reformulate("Macrohaplogroup", response = pc_name),
    data = pca_df
  )

  fit_cohort <- lm(
    reformulate("Cohort", response = pc_name),
    data = pca_df
  )

  fit_both <- lm(
    reformulate(c("Macrohaplogroup", "Cohort"), response = pc_name),
    data = pca_df
  )

  data.frame(
    PC = pc_name,
    R2_haplogroup = summary(fit_haplo)$r.squared,
    R2_cohort = summary(fit_cohort)$r.squared,
    R2_both = summary(fit_both)$r.squared,
    Residual_variance = 1 - summary(fit_both)$r.squared
  )
})

pc_structure_results <- bind_rows(pc_structure_results)

print(pc_structure_results)

write.csv(
  pc_structure_results,
  file.path(outdir, "mtDNA_PC_structure_by_haplogroup_and_cohort.csv"),
  row.names = FALSE
)


# =============================================================================
# 7. EXTRACT AND RECODE mtDNA GENOTYPES
# =============================================================================

gt_character <- extract.gt(
  vcf,
  element = "GT"
)

cat("\nOriginal genotype values in VCF:\n")
print(sort(unique(as.vector(gt_character))))


recode_mt_call <- function(x) {
  
  # Missing genotype
  if (
    is.na(x) ||
    x %in% c(".", "./.", ".|.")
  ) {
    return(NA_real_)
  }
  
  # Convert phased separator to /
  x <- gsub("\\|", "/", x)
  
  # ---------------------------------------------------------
  # Haploid representation
  # ---------------------------------------------------------
  if (x == "0") {
    return(0)
  }
  
  if (x == "1") {
    return(1)
  }
  
  # Other alternate alleles at multiallelic sites
  if (grepl("^[2-9][0-9]*$", x)) {
    return(1)
  }
  
  # ---------------------------------------------------------
  # Diploid representation used for imputation compatibility
  # ---------------------------------------------------------
  if (x == "0/0") {
    return(0)
  }
  
  if (x == "1/1") {
    return(1)
  }
  
  # Homozygous non-reference multiallelic call, e.g. 2/2
  alleles <- strsplit(x, "/", fixed = TRUE)[[1]]
  
  if (
    length(alleles) == 2 &&
    alleles[1] == alleles[2] &&
    alleles[1] != "0"
  ) {
    return(1)
  }
  
  # Any heterozygous mitochondrial call is uncertain
  # and is therefore treated as missing
  if (length(alleles) == 2 && alleles[1] != alleles[2]) {
    return(NA_real_)
  }
  
  # Anything unexpected
  return(NA_real_)
}


gt_binary <- apply(
  gt_character,
  c(1, 2),
  recode_mt_call
)

# Samples as rows, variants as columns
gt_t <- t(gt_binary)

rownames(gt_t) <- colnames(gt_character)
colnames(gt_t) <- rownames(gt_character)


cat("\nRecoded genotype values:\n")
print(sort(unique(as.vector(gt_t))))

cat("\nCarrier count summary:\n")

carrier_counts <- colSums(
  gt_t == 1,
  na.rm = TRUE
)

print(summary(carrier_counts))

cat(
  "\nVariants with >=5 alternate-allele carriers:",
  sum(carrier_counts >= 5),
  "\n"
)


# =============================================================================
# 8. PREPARE ASSOCIATION COVARIATES
# =============================================================================

association_covariates <- covariate_data %>%
  transmute(
    Sample = SampleID,
    Age = Age,
    Sex = factor(Sex),
    Phenotype = ifelse(Phenotype == "Case", 1, 0),
    Haplogroup = factor(Macrohaplogroup),
    Cohort = factor(Cohort)
  )

stopifnot(!anyDuplicated(association_covariates$Sample))

# =============================================================================
# 9. LOAD AND MERGE GLOBAL AUTOSOMAL ANCESTRY
# =============================================================================

ancestry_data <- read.csv(
  ancestry_file,
  header = TRUE,
  stringsAsFactors = FALSE
)

stopifnot("SampleID" %in% colnames(ancestry_data))
stopifnot(!anyDuplicated(ancestry_data$SampleID))

ancestry_cols <- c(
  "KhoeSan",
  "Bantu.speakingAfrican",
  "European",
  "EastAsian"
)

stopifnot(all(ancestry_cols %in% colnames(ancestry_data)))

# Explicitly use KhoeSan as the omitted/reference ancestry component.
reference_ancestry <- "KhoeSan"
ancestry_covariates <- setdiff(ancestry_cols, reference_ancestry)

missing_ancestry <- setdiff(
  association_covariates$Sample,
  ancestry_data$SampleID
)

if (length(missing_ancestry) > 0) {
  cat(
    "\nSamples without global ancestry data:",
    length(missing_ancestry),
    "\n"
  )
  print(missing_ancestry)
}

analysis_data <- merge(
  association_covariates,
  ancestry_data,
  by.x = "Sample",
  by.y = "SampleID",
  all = FALSE,
  sort = FALSE
)

cat(
  "\nIndividuals retained after ancestry merge:",
  nrow(analysis_data),
  "\n"
)


# -----------------------------------------------------------------------------
# 9.1 Descriptive ancestry analyses
# -----------------------------------------------------------------------------

ancestry_by_cohort_results <- lapply(
  ancestry_cols,
  function(component) {

    fit <- lm(
      reformulate("Cohort", response = component),
      data = analysis_data
    )

    anova_result <- anova(fit)

    data.frame(
      Component = component,
      R2 = summary(fit)$r.squared,
      p_value = anova_result$`Pr(>F)`[1]
    )
  }
) %>%
  bind_rows()

print(ancestry_by_cohort_results)

write.csv(
  ancestry_by_cohort_results,
  file.path(outdir, "Global_ancestry_by_cohort.csv"),
  row.names = FALSE
)

ancestry_long <- analysis_data %>%
  select(Sample, Haplogroup, all_of(ancestry_cols)) %>%
  pivot_longer(
    cols = all_of(ancestry_cols),
    names_to = "Component",
    values_to = "Proportion"
  )

ancestry_by_haplogroup <- ggplot(
  ancestry_long,
  aes(x = Haplogroup, y = Proportion, fill = Component)
) +
  geom_boxplot(outlier.size = 0.5) +
  scale_fill_manual(values = ancestry_colors) +
  theme_classic() +
  labs(
    x = "Mitochondrial macrohaplogroup",
    y = "Global ancestry proportion",
    fill = "Autosomal ancestry"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

ggsave(
  file.path(outdir, "Global_ancestry_by_macrohaplogroup.png"),
  plot = ancestry_by_haplogroup,
  width = 12,
  height = 6,
  dpi = 300
)


# =============================================================================
# 10. ALIGN mtDNA GENOTYPES AND COVARIATES BY SAMPLE ID
# =============================================================================

common_samples <- intersect(
  rownames(gt_t),
  analysis_data$Sample
)

gt_aligned <- gt_t[
  common_samples,
  ,
  drop = FALSE
]

data_aligned <- analysis_data[
  match(common_samples, analysis_data$Sample),
  ,
  drop = FALSE
]

stopifnot(
  identical(
    rownames(gt_aligned),
    data_aligned$Sample
  )
)

cat("\nSamples in genotype matrix:", nrow(gt_t), "\n")
cat("Samples with phenotype/covariate/ancestry data:", nrow(analysis_data), "\n")
cat("Samples used for ancestry-adjusted analyses:", length(common_samples), "\n")


# =============================================================================
# 11. MITOCHONDRIAL MACROHAPLOGROUP ASSOCIATION WITH TB
# =============================================================================

ancestry_terms <- paste(
  ancestry_covariates,
  collapse = " + "
)

haplo_null_formula <- as.formula(
  paste(
    "Phenotype ~ Age + Sex + Cohort +",
    ancestry_terms
  )
)

haplo_full_formula <- as.formula(
  paste(
    "Phenotype ~ Haplogroup + Age + Sex + Cohort +",
    ancestry_terms
  )
)

fit_haplo_null <- glm(
  haplo_null_formula,
  data = data_aligned,
  family = binomial(link = "logit")
)

fit_haplo_full <- glm(
  haplo_full_formula,
  data = data_aligned,
  family = binomial(link = "logit")
)

cat("\n=== Overall mitochondrial macrohaplogroup association ===\n")

haplogroup_lrt <- anova(
  fit_haplo_null,
  fit_haplo_full,
  test = "Chisq"
)

print(haplogroup_lrt)

haplogroup_model_comparison <- data.frame(
  Model = c(
    "Age + Sex + Cohort + Ancestry",
    "Haplogroup + Age + Sex + Cohort + Ancestry"
  ),
  AIC = c(
    AIC(fit_haplo_null),
    AIC(fit_haplo_full)
  ),
  Residual_deviance = c(
    deviance(fit_haplo_null),
    deviance(fit_haplo_full)
  )
)

write.csv(
  haplogroup_model_comparison,
  file.path(outdir, "Haplogroup_model_comparison.csv"),
  row.names = FALSE
)

# Haplogroup-specific estimates
haplo_coef <- summary(fit_haplo_full)$coefficients

haplogroup_results <- data.frame(
  Term = rownames(haplo_coef),
  Beta = haplo_coef[, "Estimate"],
  SE = haplo_coef[, "Std. Error"],
  OR = exp(haplo_coef[, "Estimate"]),
  CI_lower = exp(
    haplo_coef[, "Estimate"] -
      1.96 * haplo_coef[, "Std. Error"]
  ),
  CI_upper = exp(
    haplo_coef[, "Estimate"] +
      1.96 * haplo_coef[, "Std. Error"]
  ),
  p_value = haplo_coef[, "Pr(>|z|)"]
)

write.csv(
  haplogroup_results,
  file.path(outdir, "Haplogroup_association_results.csv"),
  row.names = FALSE
)

# Multicollinearity diagnostic
vif_results <- car::vif(fit_haplo_full)

cat("\n=== VIF/GVIF ===\n")
print(vif_results)

# Convert GVIF output to a tidy table where applicable
if (is.matrix(vif_results)) {
  vif_table <- as.data.frame(vif_results)
  vif_table$Term <- rownames(vif_table)

  if (all(c("GVIF", "Df") %in% colnames(vif_table))) {
    vif_table$GVIF_scaled <- vif_table$GVIF^(1 / (2 * vif_table$Df))
  }

  write.csv(
    vif_table,
    file.path(outdir, "Haplogroup_model_VIF.csv"),
    row.names = FALSE
  )
}


# =============================================================================
# 12. VARIANT-LEVEL FIRTH LOGISTIC REGRESSION FUNCTION
# =============================================================================

run_mt_variant <- function(
    variant_idx,
    genotypes,
    data,
    ancestry_covariates,
    adjust_haplogroup = FALSE,
    min_carriers = 5) {

  df <- data.frame(
    phenotype = data$Phenotype,
    genotype = genotypes[, variant_idx],
    age = data$Age,
    sex = factor(data$Sex),
    cohort = factor(data$Cohort),
    haplogroup = factor(data$Haplogroup)
  )

  df <- cbind(
    df,
    data[, ancestry_covariates, drop = FALSE]
  )

  required_variables <- c(
    "phenotype",
    "genotype",
    "age",
    "sex",
    "cohort",
    ancestry_covariates
  )

  if (adjust_haplogroup) {
    required_variables <- c(
      required_variables,
      "haplogroup"
    )
  }

  # Variant-specific complete-case analysis
  df <- df[
    complete.cases(df[, required_variables, drop = FALSE]),
    ,
    drop = FALSE
  ]

  # Variant must be polymorphic
  if (length(unique(df$genotype)) < 2) {
    return(NULL)
  }

  n_carriers <- sum(df$genotype == 1)
  n_noncarriers <- sum(df$genotype == 0)

  if (n_carriers < min_carriers) {
    return(NULL)
  }

  model_terms <- c(
    "genotype",
    "age",
    "sex",
    "cohort",
    ancestry_covariates
  )

  if (adjust_haplogroup) {
    model_terms <- c(
      model_terms,
      "haplogroup"
    )
  }

  model_formula <- as.formula(
    paste(
      "phenotype ~",
      paste(model_terms, collapse = " + ")
    )
  )

  tryCatch({

    fit <- logistf(
      model_formula,
      data = df,
      pl = TRUE
    )

    coef_idx <- which(
      names(fit$coefficients) == "genotype"
    )

    if (length(coef_idx) != 1) {
      return(NULL)
    }

    beta <- unname(
      fit$coefficients[coef_idx]
    )

    se <- unname(
      sqrt(diag(vcov(fit)))[coef_idx]
    )

    p_value <- unname(
      fit$prob[coef_idx]
    )

    data.frame(
      beta = beta,
      se = se,
      OR = exp(beta),
      CI_lower = exp(unname(fit$ci.lower[coef_idx])),
      CI_upper = exp(unname(fit$ci.upper[coef_idx])),
      p_value = p_value,
      n_samples = nrow(df),
      n_carriers = n_carriers,
      n_noncarriers = n_noncarriers
    )

  }, error = function(e) {

    message(
      "Variant ",
      variant_idx,
      " failed: ",
      conditionMessage(e)
    )

    return(NULL)
  })
}


# =============================================================================
# 13. PRIMARY miWAS
#
# TB ~ mtSNP + Age + Sex + Cohort + global autosomal ancestry
# =============================================================================

cat("\n=== Running primary miWAS ===\n")

primary_list <- lapply(
  seq_len(ncol(gt_aligned)),
  function(i) {

    if (i %% 100 == 0) {
      cat("Variant", i, "of", ncol(gt_aligned), "\n")
    }

    run_mt_variant(
      variant_idx = i,
      genotypes = gt_aligned,
      data = data_aligned,
      ancestry_covariates = ancestry_covariates,
      adjust_haplogroup = FALSE,
      min_carriers = 5
    )
  }
)

valid_primary <- !vapply(
  primary_list,
  is.null,
  logical(1)
)

if (!any(valid_primary)) {
  stop("No variants passed the primary miWAS analysis.")
}

primary_results <- cbind(
  variant_info[valid_primary, , drop = FALSE],
  bind_rows(primary_list[valid_primary])
)

primary_results <- primary_results %>%
  arrange(p_value)


# =============================================================================
# 14. HAPLOGROUP-ADJUSTED SENSITIVITY miWAS
#
# TB ~ mtSNP + Age + Sex + Cohort + global ancestry + macrohaplogroup
# =============================================================================

cat("\n=== Running haplogroup-adjusted sensitivity miWAS ===\n")

sensitivity_list <- lapply(
  seq_len(ncol(gt_aligned)),
  function(i) {

    if (i %% 100 == 0) {
      cat("Variant", i, "of", ncol(gt_aligned), "\n")
    }

    run_mt_variant(
      variant_idx = i,
      genotypes = gt_aligned,
      data = data_aligned,
      ancestry_covariates = ancestry_covariates,
      adjust_haplogroup = TRUE,
      min_carriers = 5
    )
  }
)

valid_sensitivity <- !vapply(
  sensitivity_list,
  is.null,
  logical(1)
)

if (!any(valid_sensitivity)) {
  warning("No variants passed the haplogroup-adjusted sensitivity analysis.")
  sensitivity_results <- NULL
} else {

  sensitivity_results <- cbind(
    variant_info[valid_sensitivity, , drop = FALSE],
    bind_rows(sensitivity_list[valid_sensitivity])
  ) %>%
    arrange(p_value)
}


# =============================================================================
# 15. MULTIPLE-TESTING CORRECTION
# =============================================================================

primary_results <- primary_results %>%
  mutate(
    p_Bonferroni = p.adjust(p_value, method = "bonferroni"),
    p_FDR = p.adjust(p_value, method = "BH")
  )

primary_bonferroni_threshold <- 0.05 / nrow(primary_results)

cat("\nPrimary variants tested:", nrow(primary_results), "\n")
cat(
  "Primary Bonferroni threshold:",
  format(primary_bonferroni_threshold, scientific = TRUE),
  "\n"
)

if (!is.null(sensitivity_results)) {

  sensitivity_results <- sensitivity_results %>%
    mutate(
      p_Bonferroni = p.adjust(p_value, method = "bonferroni"),
      p_FDR = p.adjust(p_value, method = "BH")
    )

  sensitivity_bonferroni_threshold <- 0.05 / nrow(sensitivity_results)

  cat(
    "Sensitivity variants tested:",
    nrow(sensitivity_results),
    "\n"
  )

  cat(
    "Sensitivity Bonferroni threshold:",
    format(sensitivity_bonferroni_threshold, scientific = TRUE),
    "\n"
  )
}


# =============================================================================
# 16. ASSOCIATION TEST CALIBRATION
# =============================================================================

calculate_lambda <- function(p_values) {

  p_values <- p_values[
    !is.na(p_values) &
      p_values > 0 &
      p_values < 1
  ]

  chi_sq <- qchisq(
    1 - p_values,
    df = 1
  )

  median(chi_sq) /
    qchisq(0.5, df = 1)
}

lambda_primary <- calculate_lambda(
  primary_results$p_value
)

cat(
  "\nPrimary miWAS lambda:",
  round(lambda_primary, 3),
  "\n"
)

png(
  file.path(outdir, "QQ_primary_miWAS.png"),
  width = 1800,
  height = 1800,
  res = 300
)

qq(
  primary_results$p_value,
  main = paste0(
    "Primary miWAS (lambda = ",
    round(lambda_primary, 3),
    ")"
  )
)

dev.off()


if (!is.null(sensitivity_results)) {

  lambda_sensitivity <- calculate_lambda(
    sensitivity_results$p_value
  )

  cat(
    "Sensitivity miWAS lambda:",
    round(lambda_sensitivity, 3),
    "\n"
  )

  png(
    file.path(outdir, "QQ_haplogroup_adjusted_sensitivity.png"),
    width = 1800,
    height = 1800,
    res = 300
  )

  qq(
    sensitivity_results$p_value,
    main = paste0(
      "Haplogroup-adjusted sensitivity (lambda = ",
      round(lambda_sensitivity, 3),
      ")"
    )
  )

  dev.off()

  lambda_table <- data.frame(
    Model = c(
      "Primary",
      "Haplogroup-adjusted sensitivity"
    ),
    Lambda = c(
      lambda_primary,
      lambda_sensitivity
    )
  )

} else {

  lambda_table <- data.frame(
    Model = "Primary",
    Lambda = lambda_primary
  )
}

write.csv(
  lambda_table,
  file.path(outdir, "miWAS_lambda_results.csv"),
  row.names = FALSE
)


# =============================================================================
# 17. SAVE FINAL RESULTS
# =============================================================================

write.csv(
  primary_results,
  file.path(outdir, "miWAS_primary_results.csv"),
  row.names = FALSE
)

if (!is.null(sensitivity_results)) {

  write.csv(
    sensitivity_results,
    file.path(outdir, "miWAS_haplogroup_adjusted_sensitivity_results.csv"),
    row.names = FALSE
  )
}

cat("\n============================================================\n")
cat("Pipeline complete.\n")
cat("Outputs written to:", outdir, "\n")
cat("============================================================\n")
