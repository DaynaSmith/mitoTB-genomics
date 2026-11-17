# ============================================================
# STATISTICAL LABEL FUNCTIONS
# ============================================================

format_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) {
    format(p, scientific = TRUE, digits = 3)
  } else {
    format(round(p, 4), nsmall = 4)
  }
}


# Sex: Wilcoxon rank-sum test
get_wilcox_label <- function(data, outcome = "total_heteroplasmy") {
  
  d <- data %>%
    filter(
      !is.na(.data[[outcome]]),
      !is.na(sex)
    )
  
  if (n_distinct(d$sex) != 2) {
    return("Wilcoxon rank-sum p = NA")
  }
  
  test <- wilcox.test(
    d[[outcome]] ~ d$sex,
    exact = FALSE
  )
  
  paste0(
    "Wilcoxon rank-sum p = ",
    format_p(test$p.value)
  )
}


# Continuous variables: Spearman correlation
get_spearman_label <- function(data, xvar,
                               outcome = "total_heteroplasmy") {
  
  d <- data %>%
    filter(
      !is.na(.data[[xvar]]),
      !is.na(.data[[outcome]])
    )
  
  test <- cor.test(
    d[[xvar]],
    d[[outcome]],
    method = "spearman",
    exact = FALSE
  )
  
  rho <- unname(test$estimate)
  
  paste0(
    "Spearman rho = ",
    round(rho, 4),
    ", p = ",
    format_p(test$p.value)
  )
}


# Haplogroup comparisons: Kruskal-Wallis
get_kw_label <- function(data, groupvar,
                         outcome = "total_heteroplasmy") {
  
  d <- data %>%
    filter(
      !is.na(.data[[groupvar]]),
      !is.na(.data[[outcome]])
    )
  
  if (n_distinct(d[[groupvar]]) < 2) {
    return("Kruskal-Wallis p = NA")
  }
  
  formula_test <- reformulate(
    groupvar,
    response = outcome
  )
  
  test <- kruskal.test(
    formula_test,
    data = d
  )
  
  paste0(
    "Kruskal-Wallis p = ",
    format_p(test$p.value)
  )
}

