library(dplyr)
library(ggplot2)

if (!exists("results")) {
  load("revised_simulation_workspace.RData")
}

results <- results %>%
  mutate(
    abs_bias = abs(bias),
    total_prediction_error = E_m + E_h
  )

method_summary <- results %>%
  group_by(method) %>%
  summarise(
    n = n(),
    mean_E_m = mean(E_m),
    mean_E_h = mean(E_h),
    mean_D_score = mean(D_score),
    bias = mean(bias),
    abs_bias = mean(abs_bias),
    RMSE = sqrt(mean(squared_error)),
    empirical_SE = sd(tau_hat),
    mean_model_SE = mean(se),
    coverage = mean(coverage),
    .groups = "drop"
  )

print(
  as.data.frame(method_summary)
)

write.csv(
  method_summary,
  "baseline_method_summary.csv",
  row.names = FALSE
)


correlations <- results %>%
  filter(method != "Oracle") %>%
  group_by(method) %>%
  summarise(
    cor_Em_absbias =
      cor(E_m, abs_bias),
    cor_Eh_absbias =
      cor(E_h, abs_bias),
    cor_D_absbias =
      cor(D_score, abs_bias),
    .groups = "drop"
  )

print(
  as.data.frame(correlations)
)

write.csv(
  correlations,
  "baseline_correlations.csv",
  row.names = FALSE
)


regression_results <- list()

for (m in unique(results$method)) {
  
  dat <- results %>%
    filter(method == m)
  
  if (m == "Oracle") {
    next
  }
  
  fit_E <- lm(
    abs_bias ~ E_m + E_h,
    data = dat
  )
  
  fit_D <- lm(
    abs_bias ~ D_score,
    data = dat
  )
  
  regression_results[[m]] <- data.frame(
    method = m,
    
    R2_E =
      summary(fit_E)$r.squared,
    
    p_E_m =
      coef(summary(fit_E))["E_m", "Pr(>|t|)"],
    
    p_E_h =
      coef(summary(fit_E))["E_h", "Pr(>|t|)"],
    
    R2_D =
      summary(fit_D)$r.squared,
    
    beta_D =
      coef(fit_D)["D_score"],
    
    p_D =
      coef(summary(fit_D))["D_score", "Pr(>|t|)"]
  )
}

regression_summary <- bind_rows(
  regression_results
)

print(
  as.data.frame(regression_summary)
)

write.csv(
  regression_summary,
  "baseline_regression_summary.csv",
  row.names = FALSE
)


p_bias <- ggplot(
  method_summary,
  aes(
    x = reorder(method, abs_bias),
    y = abs_bias
  )
) +
  geom_col() +
  labs(
    x = NULL,
    y = "Mean absolute causal bias"
  ) +
  theme_bw()

ggsave(
  "baseline_bias_by_method.png",
  p_bias,
  width = 8,
  height = 6,
  dpi = 300
)


p_error <- ggplot(
  results %>%
    filter(method != "Oracle"),
  aes(
    x = total_prediction_error,
    y = abs_bias
  )
) +
  geom_point(
    alpha = 0.3
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE
  ) +
  facet_wrap(
    ~ method,
    scales = "free"
  ) +
  labs(
    x = "Overall nuisance prediction error",
    y = "Absolute causal bias"
  ) +
  theme_bw()

ggsave(
  "baseline_error_vs_bias.png",
  p_error,
  width = 10,
  height = 7,
  dpi = 300
)


p_D <- ggplot(
  results %>%
    filter(method != "Oracle"),
  aes(
    x = D_score,
    y = abs_bias
  )
) +
  geom_point(
    alpha = 0.3
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE
  ) +
  facet_wrap(
    ~ method,
    scales = "free"
  ) +
  labs(
    x = "|E(delta_m delta_h)|",
    y = "Absolute causal bias"
  ) +
  theme_bw()

ggsave(
  "baseline_Dscore_vs_bias.png",
  p_D,
  width = 10,
  height = 7,
  dpi = 300
)


cat(
  "\nBaseline analysis finished.\n"
)

cat(
  "Files saved in:\n",
  getwd(),
  "\n"
)