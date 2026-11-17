#!/usr/bin/env Rscript

# ==============================================================================
# STRICT mtDNA HAPLOGROUP CONCORDANCE
#
# Metrics:
#   1. Exact haplogroup match
#   2. Hierarchical Precision, Recall and F-measure using PhyloTree Build 17
#
# Design:
#   - Array and Imputed calls are evaluated against sequence-derived truth.
#   - Comparisons are made on the SAME paired samples within each platform.
#   - Haplogroup labels are NOT truncated to force a PhyloTree match.
#   - Unmatched PhyloTree labels are reported and excluded only from
#     hierarchical calculations.
#   - Exact matches remain genuine exact matches after minimal label cleaning.
#
# Input formats supported:
#
# Single platform:
#   Sample_ID | True_Haplo | Array | Imputed
#
# Multiple platforms:
#   Sample_ID | True_Haplo |
#   Array_Affymetrix | Imputed_Affymetrix |
#   Array_H3Africa   | Imputed_H3Africa   | ...
#
# ==============================================================================


# ==============================================================================
# 0. SETUP
# ==============================================================================

required_pkgs <- c(
  "data.table",
  "ggplot2",
  "boot"
)

new_pkgs <- required_pkgs[
  !(required_pkgs %in% installed.packages()[, "Package"])
]

if (length(new_pkgs) > 0) {
  install.packages(new_pkgs)
}

library(data.table)
library(ggplot2)
library(boot)


# ==============================================================================
# 1. USER SETTINGS
# ==============================================================================

# Change for each validation cohort
infile <- "Real_cohort_haplogroups.csv"

# PhyloTree Build 17 parent table
tree_file <- "../../Scripts/Analyses/phylotree17_parents.csv"

# Output prefix
output_prefix <- "Experimental_MT-TB"

# Number of bootstrap replicates
n_boot <- 2000

set.seed(123)


# ==============================================================================
# 2. READ INPUT DATA
# ==============================================================================

dt <- fread(
  infile,
  na.strings = c(
    "",
    "NA",
    "N/A",
    "."
  )
)

required_cols <- c(
  "Sample_ID",
  "True_Haplo"
)

missing_cols <- setdiff(
  required_cols,
  names(dt)
)

if (length(missing_cols) > 0) {
  stop(
    "Missing required columns: ",
    paste(missing_cols, collapse = ", ")
  )
}

# Check duplicate sample IDs
if (anyDuplicated(dt$Sample_ID)) {
  warning(
    "Duplicate Sample_ID values detected in input file."
  )
}


# ==============================================================================
# 3. DETECT PLATFORM STRUCTURE
# ==============================================================================

detect_platform_map <- function(dt) {
  
  cn <- names(dt)
  
  # Single-platform structure
  if (all(c("Array", "Imputed") %in% cn)) {
    
    return(
      data.table(
        Platform = "Overall",
        Array_col = "Array",
        Imputed_col = "Imputed"
      )
    )
  }
  
  # Multi-platform structure
  array_cols <- grep(
    "^Array_",
    cn,
    value = TRUE
  )
  
  imputed_cols <- grep(
    "^Imputed_",
    cn,
    value = TRUE
  )
  
  platforms_array <- sub(
    "^Array_",
    "",
    array_cols
  )
  
  platforms_imputed <- sub(
    "^Imputed_",
    "",
    imputed_cols
  )
  
  common_platforms <- intersect(
    platforms_array,
    platforms_imputed
  )
  
  if (length(common_platforms) == 0) {
    
    stop(
      "Could not identify matching Array_<Platform> and ",
      "Imputed_<Platform> columns."
    )
  }
  
  data.table(
    Platform = common_platforms,
    Array_col = paste0(
      "Array_",
      common_platforms
    ),
    Imputed_col = paste0(
      "Imputed_",
      common_platforms
    )
  )
}


platform_map <- detect_platform_map(dt)

is_single_platform <- (
  nrow(platform_map) == 1 &&
    platform_map$Platform[1] == "Overall"
)

cat("\nDetected platform structure:\n")
print(platform_map)


# ==============================================================================
# 4. LOAD AND VALIDATE PHYLOTREE BUILD 17
# ==============================================================================

tree_tbl <- fread(
  tree_file,
  na.strings = c("", "NA")
)

if (!all(
  c("haplogroup", "parent") %in%
  names(tree_tbl)
)) {
  
  stop(
    "PhyloTree table must contain columns ",
    "'haplogroup' and 'parent'."
  )
}

# Remove duplicate haplogroup definitions
if (anyDuplicated(tree_tbl$haplogroup)) {
  
  stop(
    "Duplicate haplogroup names detected ",
    "in PhyloTree parent table."
  )
}

parent_map <- setNames(
  tree_tbl$parent,
  tree_tbl$haplogroup
)

valid_names <- tree_tbl$haplogroup

ROOT <- "mt-MRCA"


# ==============================================================================
# 5. STRICT HAPLOGROUP LABEL CLEANING
# ==============================================================================

# Minimal cleaning only.
#
# We remove:
#   - surrounding whitespace
#   - HaploGrep uncertainty annotations following ! or ?
#   - parenthetical annotations
#
# IMPORTANT:
# We DO NOT progressively truncate a haplogroup until it matches PhyloTree.

clean_haplogroup <- function(x) {
  
  if (
    is.na(x) ||
    trimws(x) == ""
  ) {
    return(NA_character_)
  }
  
  x <- trimws(x)
  
  # Remove uncertainty / annotation suffixes
  x <- sub(
    "[!?].*$",
    "",
    x
  )
  
  # Remove parenthetical annotation
  x <- sub(
    "\\s*\\(.*\\)\\s*$",
    "",
    x
  )
  
  trimws(x)
}


# ==============================================================================
# 6. STRICT PHYLOTREE MATCHING
# ==============================================================================

match_tree_haplogroup <- function(x) {
  
  x <- clean_haplogroup(x)
  
  if (is.na(x)) {
    return(NA_character_)
  }
  
  if (x %in% valid_names) {
    return(x)
  }
  
  # No forced truncation
  NA_character_
}


# ==============================================================================
# 7. ANCESTOR PATH FUNCTION
# ==============================================================================

ancestor_cache <- new.env(
  parent = emptyenv()
)


get_ancestors <- function(hg_raw) {
  
  cleaned <- clean_haplogroup(hg_raw)
  
  if (is.na(cleaned)) {
    return(NA_character_)
  }
  
  if (
    exists(
      cleaned,
      envir = ancestor_cache,
      inherits = FALSE
    )
  ) {
    
    return(
      get(
        cleaned,
        envir = ancestor_cache,
        inherits = FALSE
      )
    )
  }
  
  hg <- match_tree_haplogroup(cleaned)
  
  if (is.na(hg)) {
    
    assign(
      cleaned,
      NA_character_,
      envir = ancestor_cache
    )
    
    return(NA_character_)
  }
  
  chain <- character(0)
  
  cur <- hg
  
  visited <- character(0)
  
  repeat {
    
    # Detect malformed loops in parent table
    if (cur %in% visited) {
      
      stop(
        "Cycle detected in PhyloTree parent table at: ",
        cur
      )
    }
    
    visited <- c(
      visited,
      cur
    )
    
    chain <- c(
      chain,
      cur
    )
    
    nxt <- parent_map[[cur]]
    
    if (
      is.null(nxt) ||
      is.na(nxt) ||
      nxt == ""
    ) {
      break
    }
    
    cur <- nxt
  }
  
  # Root itself should not contribute to similarity
  chain <- chain[
    chain != ROOT
  ]
  
  assign(
    cleaned,
    chain,
    envir = ancestor_cache
  )
  
  chain
}


# ==============================================================================
# 8. HIERARCHICAL PRECISION / RECALL / F
# ==============================================================================

hier_prf <- function(
    true_hg,
    pred_hg
) {
  
  anc_true <- get_ancestors(
    true_hg
  )
  
  anc_pred <- get_ancestors(
    pred_hg
  )
  
  # Invalid tree label
  if (
    length(anc_true) == 1 &&
    is.na(anc_true)
  ) {
    
    return(
      c(
        hP = NA_real_,
        hR = NA_real_,
        hF = NA_real_
      )
    )
  }
  
  if (
    length(anc_pred) == 1 &&
    is.na(anc_pred)
  ) {
    
    return(
      c(
        hP = NA_real_,
        hR = NA_real_,
        hF = NA_real_
      )
    )
  }
  
  if (
    length(anc_true) == 0 ||
    length(anc_pred) == 0
  ) {
    
    return(
      c(
        hP = NA_real_,
        hR = NA_real_,
        hF = NA_real_
      )
    )
  }
  
  shared <- length(
    intersect(
      anc_true,
      anc_pred
    )
  )
  
  hP <- shared /
    length(anc_pred)
  
  hR <- shared /
    length(anc_true)
  
  hF <- if (
    hP + hR == 0
  ) {
    
    0
    
  } else {
    
    2 * hP * hR /
      (hP + hR)
  }
  
  c(
    hP = hP,
    hR = hR,
    hF = hF
  )
}


# ==============================================================================
# 9. WILSON 95% CI
# ==============================================================================

wilson_ci <- function(
    x,
    n,
    conf = 0.95
) {
  
  if (n == 0) {
    
    return(
      c(
        low = NA_real_,
        high = NA_real_
      )
    )
  }
  
  z <- qnorm(
    1 - (1 - conf) / 2
  )
  
  phat <- x / n
  
  denom <- 1 +
    z^2 / n
  
  center <- (
    phat +
      z^2 / (2 * n)
  ) / denom
  
  half <- (
    z *
      sqrt(
        phat *
          (1 - phat) /
          n +
          z^2 /
          (4 * n^2)
      )
  ) / denom
  
  c(
    low = max(
      0,
      center - half
    ),
    high = min(
      1,
      center + half
    )
  )
}


# ==============================================================================
# 10. BOOTSTRAP FUNCTIONS
# ==============================================================================

bootstrap_mean_ci <- function(
    x,
    R = n_boot
) {
  
  x <- x[
    !is.na(x)
  ]
  
  if (length(x) < 2) {
    
    return(
      c(
        low = NA_real_,
        high = NA_real_
      )
    )
  }
  
  b <- boot(
    data = x,
    statistic = function(d, i) {
      mean(d[i])
    },
    R = R
  )
  
  ci <- tryCatch(
    
    boot.ci(
      b,
      type = "perc"
    )$percent[4:5],
    
    error = function(e) {
      c(
        NA_real_,
        NA_real_
      )
    }
  )
  
  c(
    low = ci[1],
    high = ci[2]
  )
}


bootstrap_paired_difference <- function(
    array,
    imputed,
    R = n_boot
) {
  
  keep <- (
    !is.na(array) &
      !is.na(imputed)
  )
  
  dat <- data.frame(
    Array = array[keep],
    Imputed = imputed[keep]
  )
  
  if (nrow(dat) < 2) {
    
    return(
      c(
        diff = NA_real_,
        low = NA_real_,
        high = NA_real_,
        p = NA_real_
      )
    )
  }
  
  statistic <- function(
    d,
    i
  ) {
    
    mean(
      d$Imputed[i]
    ) -
      mean(
        d$Array[i]
      )
  }
  
  observed_diff <- statistic(
    dat,
    seq_len(nrow(dat))
  )
  
  b <- boot(
    data = dat,
    statistic = statistic,
    R = R
  )
  
  ci <- tryCatch(
    
    boot.ci(
      b,
      type = "perc"
    )$percent[4:5],
    
    error = function(e) {
      c(
        NA_real_,
        NA_real_
      )
    }
  )
  
  # Add-one correction prevents p = 0
  n_le_zero <- sum(
    b$t <= 0,
    na.rm = TRUE
  )
  
  n_ge_zero <- sum(
    b$t >= 0,
    na.rm = TRUE
  )
  
  p <- min(
    1,
    2 * min(
      (n_le_zero + 1) /
        (length(b$t) + 1),
      (n_ge_zero + 1) /
        (length(b$t) + 1)
    )
  )
  
  c(
    diff = observed_diff,
    low = ci[1],
    high = ci[2],
    p = p
  )
}


# ==============================================================================
# 11. BUILD PER-SAMPLE RESULTS
# ==============================================================================

per_sample <- data.table()


for (
  i in seq_len(
    nrow(platform_map)
  )
) {
  
  p <- platform_map$Platform[i]
  
  array_col <- platform_map$Array_col[i]
  
  imputed_col <- platform_map$Imputed_col[i]
  
  
  # ---------------------------------------------------------------------------
  # CRITICAL:
  # Restrict to samples with truth + BOTH Array and Imputed calls
  # ---------------------------------------------------------------------------
  
  sub <- dt[
    !is.na(True_Haplo) &
      !is.na(get(array_col)) &
      !is.na(get(imputed_col)),
    .(
      Sample_ID,
      True_Haplo,
      Array = get(array_col),
      Imputed = get(imputed_col)
    )
  ]
  
  
  if (nrow(sub) == 0) {
    
    warning(
      "No paired samples for platform: ",
      p
    )
    
    next
  }
  
  
  for (
    method in c(
      "Array",
      "Imputed"
    )
  ) {
    
    temp <- sub[
      ,
      .(
        Sample_ID,
        True_Haplo,
        Call = get(method)
      )
    ]
    
    
    # -------------------------------------------------------------------------
    # Exact match
    # -------------------------------------------------------------------------
    
    temp[
      ,
      True_clean := vapply(
        True_Haplo,
        clean_haplogroup,
        character(1)
      )
    ]
    
    temp[
      ,
      Call_clean := vapply(
        Call,
        clean_haplogroup,
        character(1)
      )
    ]
    
    temp[
      ,
      ExactMatch :=
        !is.na(True_clean) &
        !is.na(Call_clean) &
        True_clean == Call_clean
    ]
    
    
    # -------------------------------------------------------------------------
    # Hierarchical metrics
    # -------------------------------------------------------------------------
    
    scores <- t(
      mapply(
        hier_prf,
        temp$True_Haplo,
        temp$Call
      )
    )
    
    temp[
      ,
      `:=`(
        hP = scores[, "hP"],
        hR = scores[, "hR"],
        hF = scores[, "hF"],
        Platform = p,
        Method = method
      )
    ]
    
    per_sample <- rbind(
      per_sample,
      temp,
      fill = TRUE
    )
  }
}


# ==============================================================================
# 12. IDENTIFY UNMATCHED PHYLOTREE LABELS
# ==============================================================================

all_labels <- unique(
  c(
    per_sample$True_clean,
    per_sample$Call_clean
  )
)

all_labels <- all_labels[
  !is.na(all_labels)
]

unmatched_labels <- all_labels[
  !(all_labels %in% valid_names)
]

unmatched_dt <- data.table(
  haplogroup = sort(
    unique(unmatched_labels)
  )
)

fwrite(
  unmatched_dt,
  paste0(
    output_prefix,
    "_unmatched_phylotree_labels.csv"
  )
)

cat(
  "\nUnique labels not found in PhyloTree Build 17: ",
  nrow(unmatched_dt),
  "\n"
)

if (nrow(unmatched_dt) > 0) {
  print(unmatched_dt)
}


# ==============================================================================
# 13. SAVE PER-SAMPLE SCORES
# ==============================================================================

fwrite(
  per_sample,
  paste0(
    output_prefix,
    "_haplogroup_per_sample_scores.csv"
  )
)


# ==============================================================================
# 14. SUMMARY STATISTICS
# ==============================================================================

summary_dt <- per_sample[
  ,
  {
    
    exact_valid <- ExactMatch[
      !is.na(ExactMatch)
    ]
    
    hF_valid <- hF[
      !is.na(hF)
    ]
    
    
    n_exact <- length(
      exact_valid
    )
    
    n_correct <- sum(
      exact_valid
    )
    
    exact_ci <- wilson_ci(
      n_correct,
      n_exact
    )
    
    
    hF_ci <- bootstrap_mean_ci(
      hF_valid
    )
    
    
    list(
      
      N_exact = n_exact,
      
      ExactMatch_n = n_correct,
      
      ExactMatch_Prop =
        n_correct /
        n_exact,
      
      Exact_CI_low =
        exact_ci["low"],
      
      Exact_CI_high =
        exact_ci["high"],
      
      N_hF =
        length(hF_valid),
      
      hP =
        mean(
          hP,
          na.rm = TRUE
        ),
      
      hR =
        mean(
          hR,
          na.rm = TRUE
        ),
      
      hF =
        mean(
          hF_valid,
          na.rm = TRUE
        ),
      
      hF_CI_low =
        hF_ci["low"],
      
      hF_CI_high =
        hF_ci["high"]
    )
  },
  
  by = .(
    Platform,
    Method
  )
]


cat(
  "\n=== Haplogroup concordance summary ===\n"
)

print(summary_dt)


fwrite(
  summary_dt,
  paste0(
    output_prefix,
    "_haplogroup_concordance_summary.csv"
  )
)


# ==============================================================================
# 15. PAIRED ARRAY vs IMPUTED COMPARISONS
# ==============================================================================

diff_dt <- data.table()


for (
  p in platform_map$Platform
) {
  
  sub_p <- per_sample[
    Platform == p
  ]
  
  
  # ---------------------------------------------------------------------------
  # Exact-match paired comparison
  # ---------------------------------------------------------------------------
  
  wide_exact <- dcast(
    sub_p,
    Sample_ID + True_Haplo ~ Method,
    value.var = "ExactMatch"
  )
  
  wide_exact <- wide_exact[
    !is.na(Array) &
      !is.na(Imputed)
  ]
  
  
  exact_diff <- bootstrap_paired_difference(
    as.numeric(wide_exact$Array),
    as.numeric(wide_exact$Imputed)
  )
  
  
  # Exact McNemar test using discordant pairs
  b <- sum(
    wide_exact$Array == TRUE &
      wide_exact$Imputed == FALSE
  )
  
  c <- sum(
    wide_exact$Array == FALSE &
      wide_exact$Imputed == TRUE
  )
  
  
  if (b + c > 0) {
    
    mcnemar_exact_p <- binom.test(
      min(b, c),
      b + c,
      p = 0.5,
      alternative = "two.sided"
    )$p.value
    
  } else {
    
    mcnemar_exact_p <- 1
  }
  
  
  # ---------------------------------------------------------------------------
  # Hierarchical F paired comparison
  # ---------------------------------------------------------------------------
  
  wide_hf <- dcast(
    sub_p,
    Sample_ID + True_Haplo ~ Method,
    value.var = "hF"
  )
  
  wide_hf <- wide_hf[
    !is.na(Array) &
      !is.na(Imputed)
  ]
  
  
  hF_diff <- bootstrap_paired_difference(
    wide_hf$Array,
    wide_hf$Imputed
  )
  
  
  diff_dt <- rbind(
    diff_dt,
    
    data.table(
      
      Platform = p,
      
      N_exact_paired =
        nrow(wide_exact),
      
      Exact_Diff =
        exact_diff["diff"],
      
      Exact_Diff_CI_low =
        exact_diff["low"],
      
      Exact_Diff_CI_high =
        exact_diff["high"],
      
      Exact_boot_p =
        exact_diff["p"],
      
      Exact_McNemar_p =
        mcnemar_exact_p,
      
      N_hF_paired =
        nrow(wide_hf),
      
      hF_Diff =
        hF_diff["diff"],
      
      hF_Diff_CI_low =
        hF_diff["low"],
      
      hF_Diff_CI_high =
        hF_diff["high"],
      
      hF_boot_p =
        hF_diff["p"]
    ),
    
    fill = TRUE
  )
}


cat(
  "\n=== Imputed minus Array paired differences ===\n"
)

print(diff_dt)


fwrite(
  diff_dt,
  paste0(
    output_prefix,
    "_haplogroup_concordance_differences.csv"
  )
)


# ==============================================================================
# 16. SANITY CHECKS
# ==============================================================================

cat(
  "\n=== SANITY CHECKS ===\n"
)


# hF should fall between 0 and 1
bad_hf <- per_sample[
  !is.na(hF) &
    (
      hF < 0 |
        hF > 1
    )
]

cat(
  "hF values outside [0,1]: ",
  nrow(bad_hf),
  "\n"
)


# Exact matches with valid tree labels must have hF = 1
exact_not_one <- per_sample[
  ExactMatch == TRUE &
    !is.na(hF) &
    abs(
      hF - 1
    ) > 1e-10
]

cat(
  "Exact matches with hF != 1: ",
  nrow(exact_not_one),
  "\n"
)


# Compare exact-match and hierarchical F
audit_dt <- summary_dt[
  ,
  .(
    Platform,
    Method,
    N_exact,
    N_hF,
    ExactMatch_Prop,
    hF,
    hF_minus_exact =
      hF -
      ExactMatch_Prop
  )
]

cat(
  "\nExact match vs hierarchical F:\n"
)

print(audit_dt)


# Sample counts by platform
cat(
  "\nPaired sample counts:\n"
)

print(
  diff_dt[
    ,
    .(
      Platform,
      N_exact_paired,
      N_hF_paired
    )
  ]
)


# ==============================================================================
# 17. PLOT
# ==============================================================================

combined_long <- rbind(
  
  summary_dt[
    ,
    .(
      Platform,
      Method,
      Value =
        ExactMatch_Prop,
      CI_low =
        Exact_CI_low,
      CI_high =
        Exact_CI_high,
      Metric =
        "Exact match"
    )
  ],
  
  summary_dt[
    ,
    .(
      Platform,
      Method,
      Value =
        hF,
      CI_low =
        hF_CI_low,
      CI_high =
        hF_CI_high,
      Metric =
        "Hierarchical F-measure"
    )
  ]
)


combined_long[
  ,
  Metric := factor(
    Metric,
    levels = c(
      "Exact match",
      "Hierarchical F-measure"
    )
  )
]


method_colors <- c(
  "Array" = "#BDBDBD",
  "Imputed" = "#35978F"
)


p_bars <- ggplot(
  combined_long,
  aes(
    x = Platform,
    y = Value,
    fill = Method
  )
) +
  
  geom_col(
    position =
      position_dodge(
        width = 0.72
      ),
    width = 0.62,
    colour = "grey30",
    linewidth = 0.2
  ) +
  
  geom_errorbar(
    aes(
      ymin = CI_low,
      ymax = CI_high
    ),
    position =
      position_dodge(
        width = 0.72
      ),
    width = 0.14,
    linewidth = 0.4
  ) +
  
  facet_wrap(
    ~ Metric,
    ncol = 1
  ) +
  
  scale_fill_manual(
    values = method_colors
  ) +
  
  scale_y_continuous(
    limits = c(
      0,
      1
    ),
    breaks = seq(
      0,
      1,
      0.25
    ),
    expand = expansion(
      mult = c(
        0,
        0.03
      )
    )
  ) +
  
  labs(
    x = if (
      is_single_platform
    ) {
      NULL
    } else {
      "Genotyping platform"
    },
    y = "Concordance",
    fill = "Method"
  ) +
  
  theme_classic(
    base_size = 11
  ) +
  
  theme(
    
    strip.background =
      element_rect(
        fill = "white",
        colour = "grey30",
        linewidth = 0.4
      ),
    
    strip.text =
      element_text(
        size = 10
      ),
    
    legend.position =
      "right",
    
    axis.text.x =
      element_text(
        angle =
          if (is_single_platform) 0 else 35,
        hjust =
          if (is_single_platform) 0.5 else 1
      )
  )


print(p_bars)

ggsave(
  paste0(
    output_prefix,
    "_haplogroup_concordance.png"
  ),
  p_bars,
  width =
    if (is_single_platform) 6 else 8,
  height = 6,
  dpi = 600
)


ggsave(
  paste0(
    output_prefix,
    "_haplogroup_concordance.pdf"
  ),
  p_bars,
  width =
    if (is_single_platform) 6 else 8,
  height = 6
)


# ==============================================================================
# 18. FINISHED
# ==============================================================================

cat(
  "\nDone.\n",
  "\nOutputs:\n",
  " - ",
  output_prefix,
  "_haplogroup_per_sample_scores.csv\n",
  " - ",
  output_prefix,
  "_haplogroup_concordance_summary.csv\n",
  " - ",
  output_prefix,
  "_haplogroup_concordance_differences.csv\n",
  " - ",
  output_prefix,
  "_unmatched_phylotree_labels.csv\n",
  " - ",
  output_prefix,
  "_haplogroup_concordance.png\n",
  " - ",
  output_prefix,
  "_haplogroup_concordance.pdf\n",
  sep = ""
)

