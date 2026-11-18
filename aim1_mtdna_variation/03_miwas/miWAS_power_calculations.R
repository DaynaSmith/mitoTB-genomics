# miWAS simulation-based power / minimum detectable odds ratio
# Run AFTER sections 0-10 of miWAS_cohesive_pipeline(2).R
# Uses real aligned covariates and real variant carrier patterns.
# No original association tests are rerun; simulated TB outcomes only.

stopifnot(exists('data_aligned'), exists('gt_aligned'), exists('ancestry_covariates'))
stopifnot(requireNamespace('logistf', quietly = TRUE))
set.seed(20261118)
alpha <- 2.15e-4
n_sim <- 300L  # exploratory; increase to >=1000 for final estimates
n_cases_expected <- 412L
stopifnot(nrow(data_aligned) == 815L)
stopifnot(sum(data_aligned$Phenotype == 1) == n_cases_expected)

covars <- c('Age', 'Sex', 'Cohort', ancestry_covariates)
stopifnot(!anyNA(data_aligned[, c('Phenotype', covars)]))
base_df <- data.frame(y = as.integer(data_aligned$Phenotype),
                      age = data_aligned$Age,
                      sex = factor(data_aligned$Sex),
                      cohort = factor(data_aligned$Cohort),
                      data_aligned[, ancestry_covariates, drop = FALSE],
                      check.names = FALSE)
base_formula <- reformulate(c('age', 'sex', 'cohort', ancestry_covariates), response='y')
base_fit <- glm(base_formula, family=binomial(), data=base_df)
base_lp <- as.numeric(predict(base_fit, type='link'))

# Choose representative REAL variants by carrier count, not arbitrary independent
# Bernoulli genotypes (which would erase mtDNA-ancestry correlation).
carrier_n <- colSums(gt_aligned == 1, na.rm=TRUE)
valid <- which(carrier_n >= 5 & carrier_n < nrow(gt_aligned))
if (!length(valid)) stop('No variants with >=5 carriers')
target_carriers <- c(5, 10, 20, 40, 80, 160)
representative <- unique(vapply(target_carriers, function(k) {
  valid[which.min(abs(carrier_n[valid] - k))]
}, integer(1)))

# For a specified OR, recalibrate intercept to keep expected total TB cases = 412.
# This is an unconditional logistic outcome simulation with fixed EXPECTED case
# count, not exact case count; it preserves observed covariate relationships.
make_prob <- function(genotype, log_or, baseline_lp) {
  lp <- baseline_lp + log_or * genotype
  intercept_shift <- uniroot(function(a) sum(plogis(lp + a)) - n_cases_expected,
                             c(-30, 30))$root
  plogis(lp + intercept_shift)
}

one_power <- function(j, OR, nsim=n_sim) {
  keep <- !is.na(gt_aligned[, j])
  df <- base_df[keep, , drop=FALSE]
  df$genotype <- as.numeric(gt_aligned[keep, j])
  if (sum(df$genotype == 1) < 5 || length(unique(df$genotype)) < 2) return(NA_real_)
  # Calibrate expected case fraction for variant-complete subset.
  lp <- base_lp[keep] + log(OR) * df$genotype
  target_cases_subset <- sum(base_df$y[keep])
  shift <- uniroot(function(a) sum(plogis(lp + a)) - target_cases_subset,
                   c(-30, 30))$root
  p <- plogis(lp + shift)
  full_formula <- reformulate(c('genotype', 'age', 'sex', 'cohort', ancestry_covariates), response='y')
  detected <- replicate(nsim, {
    df$y <- rbinom(nrow(df), size=1L, prob=p)
    fit <- tryCatch(suppressWarnings(logistf::logistf(full_formula, data=df,
                                                      pl=TRUE)), error=function(e) NULL)
    if (is.null(fit)) {
      NA_real_
    } else {
      pv <- unname(fit$prob['genotype'])
      if (!length(pv) || !is.finite(pv)) NA_real_ else as.numeric(pv < alpha)
    }
  })
  if (mean(is.na(detected)) > 0.05) warning('>5% model failures for variant ', j)
  if (all(is.na(detected))) return(NA_real_)
  mean(detected, na.rm=TRUE)
}

# Exploratory power curve; profile-likelihood fitting is computationally expensive.
OR_grid <- c(1, 1.5, 2, 3, 4, 6, 8, 12)
res <- do.call(rbind, lapply(representative, function(j) {
  data.frame(variant=colnames(gt_aligned)[j], carriers=carrier_n[j],
             OR=OR_grid, power=vapply(OR_grid, function(or) one_power(j, or), numeric(1)))
}))
write.csv(res, 'miWAS_Firth_power_curves.csv', row.names=FALSE)
print(res)

# Minimum detectable OR (MDE) via monotone interpolation on exploratory grid.
# If 80% is not reached, report >max(OR_grid), not a fabricated value.
MDE <- do.call(rbind, lapply(split(res, res$variant), function(d) {
  d <- d[order(d$OR), ]
  smoothed <- cummax(d$power)
  ix <- which(smoothed >= 0.80)[1]
  mde <- if (is.na(ix)) NA_real_ else if (ix == 1) d$OR[1] else
    exp(approx(smoothed[c(ix-1,ix)], log(d$OR[c(ix-1,ix)]), xout=0.8)$y)
  data.frame(variant=d$variant[1], carriers=d$carriers[1],
             MDE_OR=mde, status=if(is.na(mde)) paste0('> ', max(OR_grid)) else 'within grid')
}))
write.csv(MDE, 'miWAS_Firth_MDE_exploratory.csv', row.names=FALSE)
print(MDE)

# IMPORTANT:
# * This is conditional on observed carrier patterns and a fitted baseline model.
# * Simulated outcomes, not observed TB outcomes, are tested.
# * Monte Carlo uncertainty can be substantial with 300 replicates.
# * MDE values are exploratory until refined with >=1000 simulations and a finer OR grid.
# * The observed variants may be linked; this is per-variant power at fixed Bonferroni alpha.
# * For an exact fixed 412/404 case-control design, use weighted sampling with
#   probabilities calibrated under a retrospective ascertainment model.

