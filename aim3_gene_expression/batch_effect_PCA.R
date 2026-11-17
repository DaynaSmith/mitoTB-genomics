#!/usr/bin/env Rscript
# =============================================================================
# Batch-effect PCA check: does mtDNA gene expression cluster by PE vs SE?
#
# For documentation purposes, since PE and SE are already being analyzed as
# separate discovery/replication cohorts rather than pooled -- this confirms
# empirically that separation was the right call, the same check that led
# the original paper's authors to abandon cross-study pooling.
#
# Usage:
#   Rscript batch_pca_check.R
# =============================================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
})

# ---- EDIT THESE PATHS ----
EXPR_MATRIX <- "mtdna_expression_matrix.tsv"
PE_METADATA <- "pe_metadata_unified.tsv"
SE_METADATA <- "se_metadata_unified.tsv"
OUT_PLOT    <- "batch_pca_check.pdf"
# ---------------------------

expr <- read.delim(EXPR_MATRIX, stringsAsFactors = FALSE, check.names = FALSE)
pe_meta <- read.delim(PE_METADATA, stringsAsFactors = FALSE)
se_meta <- read.delim(SE_METADATA, stringsAsFactors = FALSE)

meta <- rbind(
  data.frame(sample_id = pe_meta$sample_id, cohort = "PE"),
  data.frame(sample_id = se_meta$sample_id, cohort = "SE")
)

count_cols <- setdiff(colnames(expr), c("gene_id", "gene_name", "gene_type"))
mat <- as.matrix(expr[, count_cols])
rownames(mat) <- expr$gene_id
mode(mat) <- "numeric"

# Keep only samples present in both the expression matrix and metadata, with complete counts
common_samples <- intersect(meta$sample_id, colnames(mat))
mat <- mat[, common_samples]
complete <- colSums(is.na(mat)) == 0
mat <- mat[, complete]
meta <- meta[match(colnames(mat), meta$sample_id), ]

cat("PCA on", ncol(mat), "samples (", sum(meta$cohort == "PE"), "GC6-74,",
    sum(meta$cohort == "SE"), "TANDEM )\n")

storage.mode(mat) <- "integer"
col_data <- data.frame(cohort = factor(meta$cohort), row.names = colnames(mat))

dds <- DESeqDataSetFromMatrix(countData = mat, colData = col_data, design = ~1)
vsd <- varianceStabilizingTransformation(dds, blind = TRUE)

pca_data <- plotPCA(vsd, intgroup = "cohort", returnData = TRUE)
percent_var <- round(100 * attr(pca_data, "percentVar"))

library(dplyr)
library(stringr)

# Syntax
pca_data <- pca_data %>% 
  mutate(group = str_replace_all(group, "PE", "GC6-74"))
pca_data <- pca_data %>% 
  mutate(group = str_replace_all(group, "SE", "TANDEM"))


p <- ggplot(pca_data, aes(PC1, PC2, color = group)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = c("#F0B27A", "#7FB3D5")) +
  labs(x = paste0("PC1: ", percent_var[1], "% variance"),
       y = paste0("PC2: ", percent_var[2], "% variance"),
       title = "mtDNA expression PCA: GC6-74 vs TANDEM") +
  theme_classic(base_size = 12)
p
ggsave(OUT_PLOT, p, width = 7, height = 6)
cat("Wrote", OUT_PLOT, "\n")

# Quick statistical check: does cohort explain a significant share of PC1/PC2 variance?
cat("\nPC1 ~ cohort:\n")
print(summary(aov(PC1 ~ cohort, data = pca_data)))
cat("\nPC2 ~ cohort:\n")
print(summary(aov(PC2 ~ cohort, data = pca_data)))

