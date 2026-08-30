set.seed(20260827)

N_ORACLE <- 1000

oracle_results <- vector(
  "list",
  N_ORACLE
)

for (s in seq_len(N_ORACLE)) {
  
  dat <- generate_data(
    n = n_eval,
    p = p,
    tau = tau_true,
    rho = rho,
    sigma_x = sigma_x,
    sigma_y = sigma_y,
    clustered = FALSE,
    n_cluster = n_cluster
  )
  
  fit <- estimate_tau(
    X = dat$X,
    Y = dat$Y,
    m_hat = dat$m0,
    h_hat = dat$h0
  )
  
  oracle_results[[s]] <- data.frame(
    tau_hat = fit$tau,
    se = fit$se,
    bias = fit$tau - tau_true,
    squared_error =
      (fit$tau - tau_true)^2,
    coverage =
      tau_true >= fit$lower &&
      tau_true <= fit$upper
  )
}

oracle_results <- bind_rows(
  oracle_results
)

oracle_summary <- oracle_results %>%
  summarise(
    bias = mean(bias),
    RMSE = sqrt(mean(squared_error)),
    empirical_SE = sd(tau_hat),
    mean_model_SE = mean(se),
    coverage = mean(coverage)
  )

print(
  oracle_summary,
  width = Inf
)

oracle_summary <- oracle_results %>%
  summarise(
    bias = mean(bias),
    RMSE = sqrt(mean(squared_error)),
    empirical_SE = sd(tau_hat),
    mean_model_SE = mean(se),
    coverage = mean(coverage)
  )

print(as.data.frame(oracle_summary))

N_SIM <- 500

results_list <- vector(
  "list",
  N_SIM
)

cat(
  "Running scenario 2 with",
  N_SIM,
  "replicates\n"
)

for (s in seq_len(N_SIM)) {
  
  cat(
    "Replicate",
    s,
    "of",
    N_SIM,
    "\n"
  )
  
  results_list[[s]] <- run_one(
    seed = 20000 + s
  )
  
  if (s %% 50 == 0) {
    
    save(
      results_list,
      file = "scenario2_checkpoint.RData"
    )
    
    cat(
      "Checkpoint saved at replicate",
      s,
      "\n"
    )
  }
}


results <- bind_rows(
  results_list,
  .id = "replicate"
)


results$replicate <- as.integer(
  results$replicate
)


summary_results <- results %>%
  group_by(method) %>%
  summarise(
    n = n(),
    mean_E_m = mean(E_m),
    mean_E_h = mean(E_h),
    mean_D_score = mean(D_score),
    bias = mean(bias),
    abs_bias = mean(abs(bias)),
    RMSE = sqrt(mean(squared_error)),
    empirical_SE = sd(tau_hat),
    mean_model_SE = mean(se),
    coverage = mean(coverage),
    .groups = "drop"
  )


print(
  as.data.frame(summary_results)
)


cat(
  "\nMissing values:\n"
)

print(
  colSums(
    is.na(results)
  )
)


cat(
  "\nNon-finite values:\n"
)

print(
  sapply(
    results,
    function(x) {
      if (is.numeric(x)) {
        sum(!is.finite(x))
      } else {
        0
      }
    }
  )
)


write.csv(
  results,
  "scenario2_simulation_results.csv",
  row.names = FALSE
)


write.csv(
  summary_results,
  "scenario2_simulation_summary.csv",
  row.names = FALSE
)


save(
  results,
  results_list,
  summary_results,
  file = "scenario2_simulation_workspace.RData"
)


cat(
  "\nScenario 2 finished.\n"
)

cat(
  "Results saved in:\n",
  getwd(),
  "\n"
)

library(dplyr)
library(ggplot2)



if (!exists("results")) {
  load("scenario2_simulation_workspace.RData")
}

results <- results %>%
  mutate(
    method = factor(
      method,
      levels = c(
        "Oracle",
        "OLS",
        "GAM",
        "XGBoost",
        "DML-XGBoost"
      )
    )
  )

cat("Rows:\n")
print(nrow(results))

cat("\nReplicates per method:\n")
print(
  results %>%
    count(method)
)

cat("\nMissing values:\n")
print(
  colSums(is.na(results))
)

cat("\nNon-finite values:\n")
print(
  sapply(
    results,
    function(x) {
      if (is.numeric(x)) {
        sum(!is.finite(x))
      } else {
        0
      }
    }
  )
)



performance_table <- results %>%
  group_by(method) %>%
  summarise(
    N = n(),
    Bias = mean(bias),
    Abs_Bias = mean(abs(bias)),
    RMSE = sqrt(mean(squared_error)),
    Empirical_SE = sd(tau_hat),
    Mean_Model_SE = mean(se),
    Coverage = mean(coverage),
    .groups = "drop"
  )

print(
  as.data.frame(performance_table),
  row.names = FALSE
)

write.csv(
  performance_table,
  "scenario2_final_performance_table.csv",
  row.names = FALSE
)



prediction_table <- results %>%
  group_by(method) %>%
  summarise(
    N = n(),
    Mean_E_m = mean(E_m),
    Mean_E_h = mean(E_h),
    Mean_D_score = mean(D_score),
    .groups = "drop"
  )

print(
  as.data.frame(prediction_table),
  row.names = FALSE
)

write.csv(
  prediction_table,
  "scenario2_final_prediction_table.csv",
  row.names = FALSE
)



coverage_check <- performance_table %>%
  mutate(
    Coverage_Error = Coverage - 0.95
  )

print(
  as.data.frame(coverage_check),
  row.names = FALSE
)

write.csv(
  coverage_check,
  "scenario2_final_coverage_check.csv",
  row.names = FALSE
)


plot_data <- results %>%
  filter(method != "Oracle") %>%
  mutate(
    abs_bias = abs(bias),
    total_prediction_error = E_m + E_h
  )

p1 <- ggplot(
  plot_data,
  aes(
    x = total_prediction_error,
    y = abs_bias
  )
) +
  geom_point(alpha = 0.35) +
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
  "scenario2_prediction_error_vs_causal_bias.png",
  p1,
  width = 10,
  height = 7,
  dpi = 300
)

p2 <- ggplot(
  plot_data,
  aes(
    x = D_score,
    y = abs_bias
  )
) +
  geom_point(alpha = 0.35) +
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
  "scenario2_joint_error_vs_causal_bias.png",
  p2,
  width = 10,
  height = 7,
  dpi = 300
)



correlations <- results %>%
  filter(method != "Oracle") %>%
  group_by(method) %>%
  summarise(
    cor_Em_absbias = cor(
      E_m,
      abs(bias)
    ),
    cor_Eh_absbias = cor(
      E_h,
      abs(bias)
    ),
    cor_D_absbias = cor(
      D_score,
      abs(bias)
    ),
    .groups = "drop"
  )

print(
  as.data.frame(correlations),
  row.names = FALSE
)

write.csv(
  correlations,
  "scenario2_final_correlations.csv",
  row.names = FALSE
)


oracle_results <- results %>%
  filter(method == "Oracle")

oracle_summary <- oracle_results %>%
  summarise(
    N = n(),
    Bias = mean(bias),
    Abs_Bias = mean(abs(bias)),
    RMSE = sqrt(mean(squared_error)),
    Empirical_SE = sd(tau_hat),
    Mean_Model_SE = mean(se),
    Coverage = mean(coverage)
  )

print(
  as.data.frame(oracle_summary),
  row.names = FALSE
)

write.csv(
  oracle_summary,
  "scenario2_oracle_sanity_check.csv",
  row.names = FALSE
)



files_needed <- c(
  "scenario2_simulation_results.csv",
  "scenario2_simulation_summary.csv",
  "scenario2_simulation_workspace.RData",
  "scenario2_final_performance_table.csv",
  "scenario2_final_prediction_table.csv",
  "scenario2_final_coverage_check.csv",
  "scenario2_final_correlations.csv",
  "scenario2_oracle_sanity_check.csv",
  "scenario2_prediction_error_vs_causal_bias.png",
  "scenario2_joint_error_vs_causal_bias.png"
)

print(
  data.frame(
    file = files_needed,
    exists = file.exists(files_needed)
  ),
  row.names = FALSE
)

cat("\nFinal QC finished.\n")
cat("Working directory:\n")
cat(getwd(), "\n")