library(dplyr)
library(openxlsx)



scenario1 <- read.csv(
  "revised_simulation_summary.csv"
)

scenario2 <- read.csv(
  "scenario2_simulation_summary.csv"
)

scenario1 <- scenario1 %>%
  mutate(
    Scenario = "Scenario 1"
  )

scenario2 <- scenario2 %>%
  mutate(
    Scenario = "Scenario 2"
  )

combined_summary <- bind_rows(
  scenario1,
  scenario2
) %>%
  select(
    Scenario,
    method,
    everything()
  )

write.csv(
  combined_summary,
  "final_simulation_summary_all_scenarios.csv",
  row.names = FALSE
)

scenario1_perf <- read.csv(
  "revised_simulation_summary.csv"
)

scenario2_perf <- read.csv(
  "scenario2_final_performance_table.csv"
)

scenario1_perf <- scenario1_perf %>%
  mutate(
    Scenario = "Scenario 1"
  )

scenario2_perf <- scenario2_perf %>%
  mutate(
    Scenario = "Scenario 2"
  )

performance_all <- bind_rows(
  scenario1_perf,
  scenario2_perf
) %>%
  select(
    Scenario,
    everything()
  )

write.csv(
  performance_all,
  "final_performance_all_scenarios.csv",
  row.names = FALSE
)

scenario2_prediction <- read.csv(
  "scenario2_final_prediction_table.csv"
)

scenario2_coverage <- read.csv(
  "scenario2_final_coverage_check.csv"
)

scenario2_correlations <- read.csv(
  "scenario2_final_correlations.csv"
)

oracle_check <- read.csv(
  "scenario2_oracle_sanity_check.csv"
)

wb <- createWorkbook()

addWorksheet(
  wb,
  "All Results"
)

writeData(
  wb,
  "All Results",
  combined_summary
)

addWorksheet(
  wb,
  "Performance"
)

writeData(
  wb,
  "Performance",
  performance_all
)

addWorksheet(
  wb,
  "Scenario 2 Prediction"
)

writeData(
  wb,
  "Scenario 2 Prediction",
  scenario2_prediction
)

addWorksheet(
  wb,
  "Scenario 2 Coverage"
)

writeData(
  wb,
  "Scenario 2 Coverage",
  scenario2_coverage
)

addWorksheet(
  wb,
  "Scenario 2 Correlations"
)

writeData(
  wb,
  "Scenario 2 Correlations",
  scenario2_correlations
)

addWorksheet(
  wb,
  "Oracle Check"
)

writeData(
  wb,
  "Oracle Check",
  oracle_check
)

saveWorkbook(
  wb,
  "final_simulation_results.xlsx",
  overwrite = TRUE
)

cat(
  "\nFiles created:\n"
)

print(
  c(
    "final_simulation_summary_all_scenarios.csv",
    "final_performance_all_scenarios.csv",
    "final_simulation_results.xlsx"
  )
)

cat(
  "\nFinal simulation tables are ready.\n"
)

rm(list = setdiff(ls(), c(
  "results",
  "results_list",
  "summary_results"
)))

library(dplyr)
library(ggplot2)
library(openxlsx)

options(stringsAsFactors = FALSE)

 

scenario1_summary <- read.csv(
  "revised_simulation_summary.csv"
)

scenario1_results <- read.csv(
  "revised_simulation_results.csv"
)

scenario2_summary <- read.csv(
  "scenario2_simulation_summary.csv"
)

scenario2_results <- read.csv(
  "scenario2_simulation_results.csv"
)

scenario1_summary <- scenario1_summary %>%
  mutate(
    Scenario = "Scenario 1"
  )

scenario2_summary <- scenario2_summary %>%
  mutate(
    Scenario = "Scenario 2"
  )

scenario1_results <- scenario1_results %>%
  mutate(
    Scenario = "Scenario 1"
  )

scenario2_results <- scenario2_results %>%
  mutate(
    Scenario = "Scenario 2"
  )

all_results <- bind_rows(
  scenario1_results,
  scenario2_results
)

all_summary <- bind_rows(
  scenario1_summary,
  scenario2_summary
)

 

print(
  all_results %>%
    count(
      Scenario,
      method
    )
)

cat("\nTotal rows:\n")
print(nrow(all_results))

cat("\nMissing values:\n")
print(
  colSums(
    is.na(all_results)
  )
)

cat("\nNon-finite numeric values:\n")
print(
  sapply(
    all_results,
    function(x) {
      if (is.numeric(x)) {
        sum(!is.finite(x))
      } else {
        0
      }
    }
  )
)

 

performance_all <- all_results %>%
  group_by(
    Scenario,
    method
  ) %>%
  summarise(
    N = n(),
    Bias = mean(bias),
    Abs_Bias = mean(abs(bias)),
    RMSE = sqrt(
      mean(squared_error)
    ),
    Empirical_SE = sd(tau_hat),
    Mean_Model_SE = mean(se),
    Coverage = mean(coverage),
    .groups = "drop"
  )

print(
  as.data.frame(performance_all),
  row.names = FALSE
)

write.csv(
  performance_all,
  "FINAL_performance_all_scenarios.csv",
  row.names = FALSE
)
 

prediction_all <- all_results %>%
  group_by(
    Scenario,
    method
  ) %>%
  summarise(
    N = n(),
    Mean_E_m = mean(E_m),
    Mean_E_h = mean(E_h),
    Mean_D_score = mean(D_score),
    .groups = "drop"
  )

print(
  as.data.frame(prediction_all),
  row.names = FALSE
)

write.csv(
  prediction_all,
  "FINAL_prediction_error_all_scenarios.csv",
  row.names = FALSE
)

 

inference_all <- all_results %>%
  group_by(
    Scenario,
    method
  ) %>%
  summarise(
    N = n(),
    Bias = mean(bias),
    RMSE = sqrt(
      mean(squared_error)
    ),
    Empirical_SE = sd(tau_hat),
    Mean_Model_SE = mean(se),
    Coverage = mean(coverage),
    .groups = "drop"
  ) %>%
  mutate(
    SE_Ratio =
      Mean_Model_SE / Empirical_SE,
    Coverage_Error =
      Coverage - 0.95
  )

print(
  as.data.frame(inference_all),
  row.names = FALSE
)

write.csv(
  inference_all,
  "FINAL_inference_all_scenarios.csv",
  row.names = FALSE
)

rm(list = ls())

library(dplyr)
library(ggplot2)
library(readr)

setwd("C:/Users/congc/Downloads/Projects/Nuisance error")



scenario1_file <- "revised_simulation_results.csv"
scenario2_file <- "scenario2_simulation_results.csv"

if (!file.exists(scenario1_file)) {
  stop("Missing Scenario 1 results file.")
}

if (!file.exists(scenario2_file)) {
  stop("Missing Scenario 2 results file.")
}

scenario1 <- read.csv(
  scenario1_file,
  stringsAsFactors = FALSE
)

scenario2 <- read.csv(
  scenario2_file,
  stringsAsFactors = FALSE
)

scenario1$Scenario <- "Scenario 1"
scenario2$Scenario <- "Scenario 2"

results_all <- bind_rows(
  scenario1,
  scenario2
)

results_all <- results_all %>%
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

cat("Total rows:\n")
print(nrow(results_all))

cat("\nRows by scenario and method:\n")

print(
  results_all %>%
    count(Scenario, method)
)



required_columns <- c(
  "Scenario",
  "method",
  "E_m",
  "E_h",
  "D_score",
  "tau_hat",
  "se",
  "bias",
  "squared_error",
  "coverage"
)

missing_columns <- setdiff(
  required_columns,
  names(results_all)
)

if (length(missing_columns) > 0) {
  stop(
    paste(
      "Missing columns:",
      paste(missing_columns, collapse = ", ")
    )
  )
}

print(
  colSums(
    is.na(results_all)
  )
)

cat("\nNon-finite values:\n")

print(
  sapply(
    results_all,
    function(x) {
      if (is.numeric(x)) {
        sum(!is.finite(x))
      } else {
        0
      }
    }
  )
)



performance_table <- results_all %>%
  group_by(Scenario, method) %>%
  summarise(
    n = n(),
    Bias = mean(bias),
    Abs_Bias = mean(abs(bias)),
    RMSE = sqrt(mean(squared_error)),
    Empirical_SE = sd(tau_hat),
    Model_SE = mean(se),
    Coverage = mean(coverage),
    SE_Ratio = Model_SE / Empirical_SE,
    .groups = "drop"
  )

performance_table <- as.data.frame(
  performance_table
)

print(
  performance_table,
  row.names = FALSE
)

write.csv(
  performance_table,
  "final_performance_table.csv",
  row.names = FALSE
)


prediction_table <- results_all %>%
  group_by(Scenario, method) %>%
  summarise(
    mean_E_m = mean(E_m),
    mean_E_h = mean(E_h),
    mean_D_score = mean(D_score),
    .groups = "drop"
  )

prediction_table <- as.data.frame(
  prediction_table
)

print(
  prediction_table,
  row.names = FALSE
)

write.csv(
  prediction_table,
  "final_prediction_table.csv",
  row.names = FALSE
)





coverage_table <- performance_table %>%
  select(
    Scenario,
    method,
    Empirical_SE,
    Model_SE,
    SE_Ratio,
    Coverage
  )

print(
  coverage_table,
  row.names = FALSE
)

write.csv(
  coverage_table,
  "final_coverage_table.csv",
  row.names = FALSE
)



oracle_check <- performance_table %>%
  filter(method == "Oracle") %>%
  mutate(
    Bias_OK = abs(Bias) < 0.02,
    Coverage_OK = Coverage >= 0.90 & Coverage <= 0.99,
    SE_OK = SE_Ratio >= 0.90 & SE_Ratio <= 1.10
  )

print(
  oracle_check,
  row.names = FALSE
)

write.csv(
  oracle_check,
  "oracle_sanity_check.csv",
  row.names = FALSE
)



plot_data <- results_all %>%
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
  geom_point(
    alpha = 0.30
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE
  ) +
  facet_grid(
    Scenario ~ method,
    scales = "free"
  ) +
  labs(
    x = "Overall nuisance prediction error",
    y = "Absolute causal bias"
  ) +
  theme_bw()

ggsave(
  "final_prediction_error_vs_causal_bias.png",
  p1,
  width = 13,
  height = 8,
  dpi = 300
)

cat("Saved: final_prediction_error_vs_causal_bias.png\n")

p2 <- ggplot(
  plot_data,
  aes(
    x = D_score,
    y = abs_bias
  )
) +
  geom_point(
    alpha = 0.30
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE
  ) +
  facet_grid(
    Scenario ~ method,
    scales = "free"
  ) +
  labs(
    x = "|E(delta_m delta_h)|",
    y = "Absolute causal bias"
  ) +
  theme_bw()

ggsave(
  "final_joint_error_vs_causal_bias.png",
  p2,
  width = 13,
  height = 8,
  dpi = 300
)


correlations <- results_all %>%
  group_by(Scenario, method) %>%
  summarise(
    cor_Em_absbias = ifelse(
      sd(E_m) == 0 | sd(abs(bias)) == 0,
      NA_real_,
      cor(E_m, abs(bias))
    ),
    cor_Eh_absbias = ifelse(
      sd(E_h) == 0 | sd(abs(bias)) == 0,
      NA_real_,
      cor(E_h, abs(bias))
    ),
    cor_D_absbias = ifelse(
      sd(D_score) == 0 | sd(abs(bias)) == 0,
      NA_real_,
      cor(D_score, abs(bias))
    ),
    .groups = "drop"
  )

correlations <- as.data.frame(
  correlations
)

print(
  correlations,
  row.names = FALSE
)

write.csv(
  correlations,
  "final_correlations.csv",
  row.names = FALSE
)



regression_results <- list()

for (sc in unique(results_all$Scenario)) {
  
  dat_sc <- results_all %>%
    filter(Scenario == sc)
  
  for (m in unique(dat_sc$method)) {
    
    dat_m <- dat_sc %>%
      filter(method == m)
    
    if (
      nrow(dat_m) >= 10 &&
      sd(dat_m$E_m) > 0 &&
      sd(dat_m$E_h) > 0
    ) {
      
      fit_E <- lm(
        abs(bias) ~ E_m + E_h,
        data = dat_m
      )
      
      fit_D <- NULL
      
      if (sd(dat_m$D_score) > 0) {
        fit_D <- lm(
          abs(bias) ~ D_score,
          data = dat_m
        )
      }
      
      regression_results[[length(regression_results) + 1]] <-
        data.frame(
          Scenario = sc,
          method = as.character(m),
          R2_prediction =
            summary(fit_E)$r.squared,
          p_E_m =
            coef(summary(fit_E))["E_m", "Pr(>|t|)"],
          p_E_h =
            coef(summary(fit_E))["E_h", "Pr(>|t|)"],
          R2_D =
            ifelse(
              is.null(fit_D),
              NA,
              summary(fit_D)$r.squared
            ),
          p_D =
            ifelse(
              is.null(fit_D),
              NA,
              coef(summary(fit_D))[
                "D_score",
                "Pr(>|t|)"
              ]
            )
        )
    }
  }
}

regression_table <- bind_rows(
  regression_results
)

print(
  regression_table,
  row.names = FALSE
)

write.csv(
  regression_table,
  "final_regression_checks.csv",
  row.names = FALSE
)



ranking <- performance_table %>%
  group_by(Scenario) %>%
  mutate(
    Rank_RMSE =
      rank(RMSE, ties.method = "min"),
    Rank_Abs_Bias =
      rank(Abs_Bias, ties.method = "min")
  ) %>%
  ungroup()

print(
  ranking,
  row.names = FALSE
)

write.csv(
  ranking,
  "final_method_ranking.csv",
  row.names = FALSE
)



manuscript_table <- performance_table %>%
  mutate(
    Bias = round(Bias, 4),
    Abs_Bias = round(Abs_Bias, 4),
    RMSE = round(RMSE, 4),
    Empirical_SE = round(Empirical_SE, 4),
    Model_SE = round(Model_SE, 4),
    SE_Ratio = round(SE_Ratio, 3),
    Coverage = round(Coverage, 3)
  ) %>%
  select(
    Scenario,
    method,
    Bias,
    Abs_Bias,
    RMSE,
    Empirical_SE,
    Model_SE,
    SE_Ratio,
    Coverage
  )

print(
  manuscript_table,
  row.names = FALSE
)

write.csv(
  manuscript_table,
  "manuscript_simulation_table.csv",
  row.names = FALSE
)



save(
  results_all,
  performance_table,
  prediction_table,
  coverage_table,
  oracle_check,
  correlations,
  regression_table,
  ranking,
  manuscript_table,
  file = "final_simulation_analysis.RData"
)

cat(
  "Saved: final_simulation_analysis.RData\n"
)



files_needed <- c(
  "revised_simulation_results.csv",
  "scenario2_simulation_results.csv",
  "final_performance_table.csv",
  "final_prediction_table.csv",
  "final_coverage_table.csv",
  "oracle_sanity_check.csv",
  "final_correlations.csv",
  "final_regression_checks.csv",
  "final_method_ranking.csv",
  "manuscript_simulation_table.csv",
  "final_prediction_error_vs_causal_bias.png",
  "final_joint_error_vs_causal_bias.png",
  "final_simulation_analysis.RData"
)

file_check <- data.frame(
  file = files_needed,
  exists = file.exists(files_needed)
)

print(
  file_check,
  row.names = FALSE
)

if (all(file_check$exists)) {
  
  cat(
    "\nALL FINAL FILES ARE PRESENT.\n"
  )
  
} else {
  
  cat(
    "\nWARNING: SOME FINAL FILES ARE MISSING.\n"
  )
  
}

cat(
  "\nWorking directory:\n",
  getwd(),
  "\n"
)

cat(
  "\nFINAL ANALYSIS FINISHED.\n"
)