## PCA 
#install.packages("adegenet")
library(adegenet)
#install.packages("vcfR")
library(vcfR)
library(devtools)
library(RColorBrewer)
library(patchwork)
library(ggplot2)
## Make PCA plots with QC'd VCF (filter for MAF > 1%)
vcf <- read.vcfR("Merged_mtDNA_imputed_geno0.05_mind0.1_INFO0.7_SetHeteroplasmiesToMissing_ExcludeMissingData.recode.vcf")
gen <- vcfR2genlight(vcf)
pca_res <- glPca(gen, nf=10)
print(pca_res)

# Compute percentage variance explained
variance_explained <- 100 * pca_res$eig / sum(pca_res$eig)
print(variance_explained)

df_var <- data.frame(
  PC = paste0("PC", seq_along(variance_explained)),
  Variance = variance_explained
)

df_var <- df_var[1:10, ]  # show only first 10 PCs
df_var$PC <- factor(df_var$PC, levels = paste0("PC", seq_along(variance_explained)))

scree_plot <- ggplot(df_var, aes(x = PC, y = Variance, group = 1)) +
  geom_line(color="darkgrey") +
  geom_point(shape=21, color="darkgrey", fill="darkgrey", size=2) +
  theme_classic() +
  labs(x = "Principal Component",
       y = "Percentage of Variance Explained"
  )
scree_plot

# Extract PC1 - PC10 variance
pc1_var <- round(variance_explained[1], 2)
pc2_var <- round(variance_explained[2], 2)
pc3_var <- round(variance_explained[3], 2)
pc4_var <- round(variance_explained[4], 2)
pc5_var <- round(variance_explained[5], 2)
pc6_var <- round(variance_explained[6], 2)
pc7_var <- round(variance_explained[7], 2)
pc8_var <- round(variance_explained[8], 2)
pc9_var <- round(variance_explained[9], 2)
pc10_var <- round(variance_explained[10], 2)

# Create a dataframe for plotting
pca_df <- data.frame(
  Sample = rownames(pca_res$scores),
  PC1 = pca_res$scores[, 1],
  PC2 = pca_res$scores[, 2],
  PC3 = pca_res$scores[, 3],
  PC4 = pca_res$scores[, 4],
  PC5 = pca_res$scores[, 5],
  PC6 = pca_res$scores[, 6],
  PC7 = pca_res$scores[, 7],
  PC8 = pca_res$scores[, 8],
  PC9 = pca_res$scores[, 9],
  PC10 = pca_res$scores[, 10]
)
head(pca_df)

# Make a dataframe for PCs and Haplogroup
meta <- read.csv("sample_data.csv", header=TRUE)
haplogroups <- meta[,c("SampleID", "BroadHaplogroup", "Macrohaplogroup")]
colnames(haplogroups) <- c("Sample", "Haplogroup", "Macrohaplogroup")
head(haplogroups)

pca_haplogroups_df <- merge(pca_df, haplogroups, by="Sample")
head(pca_haplogroups_df)

# Make a dataframe for PCs and Batch
batches <- meta[,c("SampleID", "Cohort")]
colnames(batches) <- c("Sample", "Cohort")
head(batches)

pca_batches_df <- merge(pca_df, batches, by="Sample")
head(pca_batches_df)

haplo_colors <- c(
  
  # L0a branch
  "L0a'g"     = "#F6D89A",
  "L0a"           = "#DFC27D",
  "L0a2"          = "#BF812D",
  
  # L0d branch
  "L0d1"          = "#E9A03B",
  "L0d2"          = "#F5C46B",
  "L0d3"          = "#FEE090",
  
  # Other L0
  "L0f"           = "#BFA26A",
  "L0k"           = "#AEBF7A",
  
  # L1 branch
  "L1b"          = "#4D9221",
  "L1c"          = "#2B7A0B",
  
  # L2 branch
  "L2b"           = "#66C2A4",
  "L2c"           = "#31A354",
  
  # L3 branch
  "L3"            = "#35978F",
  "L3b"           = "#80CDC1",
  "L3d"           = "#1C9099",
  "L3e"           = "#7BC8C3",
  "L3f"           = "#9EDBD7",
  "L3h"           = "#4EB3A5",
  
  # L4
  "L4b"           = "#2B8CBE",
  
  # L5
  "L5a"           = "#756BB1",
  
  # Eurasian lineages
  "M"             = "#C475E6",
  "N"             = "#E67597",
  "R"             = "#C51B7D"
)

# Look at PCA for haplogroups
PC1_2_haps <- ggplot(pca_haplogroups_df, aes(x = PC1, y = PC2, color = Haplogroup)) +
  geom_point(size = 3, alpha = 0.8) +  # Adjust size and transparency
  labs(
    title = "PC1 versus PC2",
    x = paste0("PC1 (", pc1_var, "% variance)"),
    y = paste0("PC2 (", pc2_var, "% variance)"),
    color = "Haplogroup"
  ) +
  scale_color_manual(values = haplo_colors) +
  theme_classic() +
  theme(legend.position = "left")

PC2_3_haps <- ggplot(pca_haplogroups_df, aes(x = PC2, y = PC3, color = Haplogroup)) +
  geom_point(size = 3, alpha = 0.8) +  # Adjust size and transparency
  labs(
    title = "PC2 versus PC3",
    x = paste0("PC2 (", pc2_var, "% variance)"),
    y = paste0("PC3 (", pc3_var, "% variance)"),
    color = "Haplogroup"
  ) +
  scale_color_manual(values = haplo_colors) +
  theme_classic() +
  theme(legend.position = "none")

PC1_3_haps <- ggplot(pca_haplogroups_df, aes(x = PC1, y = PC3, color = Haplogroup)) +
  geom_point(size = 3, alpha = 0.8) +  # Adjust size and transparency
  labs(
    title = "PC1 versus PC3",
    x = paste0("PC1 (", pc1_var, "% variance)"),
    y = paste0("PC3 (", pc3_var, "% variance)"),
    color = "Haplogroup"
  ) +
  scale_color_manual(values = haplo_colors) +
  theme_classic() +
  theme(legend.position = "none")

combined_plot_haps <- PC1_2_haps + PC2_3_haps / PC1_3_haps
combined_plot_haps
ggsave("PCA_haplogroups.png", plot = combined_plot_haps, height = 8, width = 15, dpi=300)

cohort_colors <- c("SAC_MEGA" = "#1F77B4", 
                   "SAC_TANDEM" = "#FF7F0E")

# Look at PCA for cohorts
PC1_2_cohort <- ggplot(pca_batches_df, aes(x = PC1, y = PC2, color = Cohort)) +
  geom_point(size = 3, alpha = 0.8) +  # Adjust size and transparency
  labs(
    title = "PC1 versus PC2",
    x = paste0("PC1 (", pc1_var, "% variance)"),
    y = paste0("PC2 (", pc2_var, "% variance)"),
    color = "Cohort"
  ) +
  scale_color_manual(values=cohort_colors) +
  theme_classic() +
  theme(legend.position = "left")

PC2_3_cohort <- ggplot(pca_batches_df, aes(x = PC2, y = PC3, color = Cohort)) +
  geom_point(size = 3, alpha = 0.8) +  # Adjust size and transparency
  labs(
    title = "PC2 versus PC3",
    x = paste0("PC2 (", pc2_var, "% variance)"),
    y = paste0("PC3 (", pc3_var, "% variance)"),
    color = "Cohort"
  ) +
  theme_classic() +
  scale_color_manual(values=cohort_colors) +
  theme(legend.position = "none")

PC1_3_cohort <- ggplot(pca_batches_df, aes(x = PC1, y = PC3, color = Cohort)) +
  geom_point(size = 3, alpha = 0.8) +  # Adjust size and transparency
  labs(
    title = "PC1 versus PC3",
    x = paste0("PC1 (", pc1_var, "% variance)"),
    y = paste0("PC3 (", pc3_var, "% variance)"),
    color = "Cohort"
  ) +
  theme_classic() +
  scale_color_manual(values=cohort_colors) +
  theme(legend.position = "none")

combined_plot_cohort <- PC1_2_cohort + PC2_3_cohort / PC1_3_cohort
combined_plot_cohort
ggsave("PCA_cohorts.png", plot = combined_plot_cohort, height = 8, width = 15, dpi=300)
