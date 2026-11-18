
library(dplyr)

resistb_results <- data.frame(
  comparison = c(
    "Active TB vs LTBI",
    "LTBI vs resisters",
    "Active TB vs resisters"
  ),
  IRR = c(0.82, 1.08, 0.84),
  CI_lower = c(0.53, 0.75, 0.57),
  CI_upper = c(1.26, 1.55, 1.24)
)

alpha <- 0.05
target_power <- 0.80

z_alpha <- qnorm(1 - alpha / 2)
z_power <- qnorm(target_power)

resistb_power <- resistb_results %>%
  mutate(
    SE_logIRR = (
      log(CI_upper) - log(CI_lower)
    ) / (2 * qnorm(0.975)),
    
    MDE_logIRR = (z_alpha + z_power) * SE_logIRR,
    
    MDE_increase = exp(MDE_logIRR),
    MDE_decrease = exp(-MDE_logIRR)
  )

print(resistb_power)

write.csv(
  resistb_power,
  "ResisTB_approximate_MDE.csv",
  row.names = FALSE
)
