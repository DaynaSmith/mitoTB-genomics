#!/usr/bin/env Rscript
# =============================================================================
# Check covariate balance (age, sex, L-status) between TB Case and Control
# before building Analysis 3's DE model -- specifically checking the two
# concerns flagged in the project checklist: (1) does age need a nonlinear
# term or stratification rather than a simple linear covariate, and (2) is
# L-haplogroup status balanced between Case/Control (if not, Analysis 2's
# findings mean it could confound a naive Case/Control comparison).
#
# Usage:
#   Rscript check_TB_covariate_balance.R
# =============================================================================

PE_METADATA <- "pe_metadata_unified.tsv"
SE_METADATA <- "se_metadata_unified.tsv"

pe_meta <- read.delim(PE_METADATA, stringsAsFactors = FALSE)
se_meta <- read.delim(SE_METADATA, stringsAsFactors = FALSE)

pe_meta$L_status <- ifelse(startsWith(pe_meta$haplogroup, "L"), "L", "non_L")
se_meta$L_status <- ifelse(startsWith(se_meta$haplogroup, "L"), "L", "non_L")

check_cohort <- function(meta, label) {
  cat("\n=====", label, "=====\n")
  meta_cc <- meta[meta$condition %in% c("Case", "Control"), ]
  meta_cc$age <- suppressWarnings(as.numeric(meta_cc$age))
  cat("N:", nrow(meta_cc), " (Case:", sum(meta_cc$condition == "Case"),
      ", Control:", sum(meta_cc$condition == "Control"), ")\n")
  
  cat("\nAge by condition:\n")
  print(aggregate(age ~ condition, data = meta_cc,
                  FUN = function(x) c(mean = mean(x, na.rm=TRUE), sd = sd(x, na.rm=TRUE),
                                      min = min(x, na.rm=TRUE), max = max(x, na.rm=TRUE))))
  age_test <- t.test(age ~ condition, data = meta_cc)
  cat("Age t-test p-value (Case vs Control):", age_test$p.value, "\n")
  
  cat("\nAge distribution shape (quartiles) by condition:\n")
  print(tapply(meta_cc$age, meta_cc$condition, quantile, na.rm = TRUE))
  
  cat("\nSex by condition:\n")
  sex_tab <- table(meta_cc$condition, meta_cc$sex)
  print(sex_tab)
  if (all(dim(sex_tab) == c(2,2))) {
    cat("Sex chi-square p-value:", chisq.test(sex_tab)$p.value, "\n")
  }
  
  cat("\nL-status by condition:\n")
  L_tab <- table(meta_cc$condition, meta_cc$L_status)
  print(L_tab)
  if (all(dim(L_tab) == c(2,2))) {
    cat("L-status chi-square p-value:", chisq.test(L_tab)$p.value, "\n")
  }
  
  invisible(meta_cc)
}

check_cohort(pe_meta, "PE")
check_cohort(se_meta, "SE")
