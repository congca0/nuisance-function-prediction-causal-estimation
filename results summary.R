    
# FINAL TABLE FIX
 
    

library(dplyr)

cat("\n   \n")
cat("REBUILDING FINAL TABLES FROM EXISTING RESULTS\n")
cat("   \n\n")

   
# 1. Load existing Scenario 1
   

scenario1_candidates <- c(
  "revised_simulation_results.csv",
  "baseline_simulation_results.csv",
  "scenario1_simulation_results.csv"
)

scenario1_file <- scenario1_candidates[
  file.exists(scenario1_candidates)
][1]

if (is.na(scenario1_file)) {
  stop("No Scenario 1 results file found.")
}

scenario1 <- read.csv(
  scenario1_file,
  stringsAsFactors = FALSE
)

scenario1$Scenario <- "Scenario 1"

scenario1$replicate <- as.integer(
  scenario1$replicate
)


   
# 2. Load existing Scenario 2
   

if (!file.exists("scenario2_simulation_results.csv")) {
  stop("scenario2_simulation_results.csv not found.")
}

scenario2 <- read.csv(
  "scenario2_simulation_results.csv",
  stringsAsFactors = FALSE
)

scenario2$Scenario <- "Scenario 2"

scenario2$replicate <- as.integer(
  scenario2$replicate
)


   
# 3. Check Scenario 1
   

cat("Scenario 1 file:\n")
cat(scenario1_file, "\n\n")

cat("Scenario 1 rows by method:\n")

print(
  scenario1 %>%
    count(method)
)


   
# 4. Check Scenario 2
   

cat("\nScenario 2 rows by method:\n")

print(
  scenario2 %>%
    count(method)
)


   
# 5. HARD CHECK
   

s1_count <- scenario1 %>%
  count(method)

s2_count <- scenario2 %>%
  count(method)


if (!all(s1_count$n == 500)) {
  stop(
    "STOP: Selected Scenario 1 file is NOT the 500-replicate Scenario 1 result."
  )
}

if (!all(s2_count$n == 500)) {
  stop(
    "STOP: Scenario 2 does not contain 500 replicates for every method."
  )
}

if (!"Oracle" %in% unique(scenario1$method)) {
  stop(
    "STOP: Scenario 1 is missing Oracle."
  )
}

if (!"Oracle" %in% unique(scenario2$method)) {
  stop(
    "STOP: Scenario 2 is missing Oracle."
  )
}


cat(
  "\nScenario 1 and Scenario 2 replicate checks PASSED.\n"
)


   
# 6. Standardize columns
   

required_columns <- c(
  "replicate",
  "method",
  "E_m",
  "E_h",
  "D_score",
  "tau_hat",
  "se",
  "bias",
  "squared_error",
  "coverage",
  "Scenario"
)

missing_s1 <- setdiff(
  required_columns,
  names(scenario1)
)

missing_s2 <- setdiff(
  required_columns,
  names(scenario2)
)

if (length(missing_s1) > 0) {
  stop(
    paste(
      "Scenario 1 missing:",
      paste(missing_s1, collapse = ", ")
    )
  )
}

if (length(missing_s2) > 0) {
  stop(
    paste(
      "Scenario 2 missing:",
      paste(missing_s2, collapse = ", ")
    )
  )
}


scenario1 <- scenario1 %>%
  select(all_of(required_columns))

scenario2 <- scenario2 %>%
  select(all_of(required_columns))


   
# 7. Combine
   

all_results <- bind_rows(
  scenario1,
  scenario2
)


cat("\nCombined rows:\n")
cat(nrow(all_results), "\n")


if (nrow(all_results) != 5000) {
  stop(
    paste(
      "STOP: Expected 5000 total rows, found",
      nrow(all_results)
    )
  )
}


   
# 8. Final Scenario × method count
   

final_count <- all_results %>%
  count(
    Scenario,
    method
  )

cat("\nFinal Scenario × method counts:\n")

print(final_count)


if (!all(final_count$n == 500)) {
  stop(
    "STOP: Every Scenario × method combination must contain 500 observations."
  )
}


 # 9. Missing-value check
 
na_check <- sapply(
  all_results,
  function(x) sum(is.na(x))
)

cat("\nMissing values:\n")

print(na_check)

if (any(na_check > 0)) {
  stop(
    "STOP: Missing values detected."
  )
}


 
# 10. Non-finite check
 

nonfinite_check <- sapply(
  all_results,
  function(x) {
    if (is.numeric(x)) {
      sum(!is.finite(x))
    } else {
      0
    }
  }
)

cat("\nNon-finite values:\n")

print(nonfinite_check)

if (any(nonfinite_check > 0)) {
  stop(
    "STOP: Non-finite values detected."
  )
}


#    
# 11. FINAL MANUSCRIPT PERFORMANCE TABLE
#    

manuscript_table <- all_results %>%
  group_by(
    Scenario,
    method
  ) %>%
  summarise(
    
    n = n(),
    
    Bias = mean(bias),
    
    Abs_Bias = mean(abs(bias)),
    
    RMSE = sqrt(
      mean(squared_error)
    ),
    
    Empirical_SE = sd(tau_hat),
    
    Model_SE = mean(se),
    
    SE_Ratio =
      Model_SE / Empirical_SE,
    
    Coverage = mean(coverage),
    
    .groups = "drop"
  )


cat(
    
)

cat(
  "FINAL MANUSCRIPT PERFORMANCE TABLE\n"
)

cat(
  "   \n\n"
)

print(
  as.data.frame(manuscript_table)
)


 # 12. Oracle sanity check
 
oracle_check <- manuscript_table %>%
  filter(
    method == "Oracle"
  ) %>%
  mutate(
    
    Bias_OK =
      abs(Bias) < 0.02,
    
    Coverage_OK =
      Coverage >= 0.90 &
      Coverage <= 0.99,
    
    SE_OK =
      SE_Ratio >= 0.90 &
      SE_Ratio <= 1.10
    
  )


cat(
    
)

cat(
  "ORACLE SANITY CHECK\n"
)

cat(
  "   \n\n"
)

print(
  as.data.frame(oracle_check)
)


if (
  !all(
    oracle_check$Bias_OK &
    oracle_check$Coverage_OK &
    oracle_check$SE_OK
  )
) {
  
  stop(
    "STOP: Oracle sanity check failed."
  )
}


cat(
  "\nORACLE SANITY CHECK PASSED.\n"
)


#    
# 13. FINAL PREDICTION ERROR TABLE
#    

prediction_table <- all_results %>%
  group_by(
    Scenario,
    method
  ) %>%
  summarise(
    
    n = n(),
    
    Mean_E_m = mean(E_m),
    
    Mean_E_h = mean(E_h),
    
    Mean_D_score = mean(D_score),
    
    .groups = "drop"
  )


cat(
    
)

cat(
  "FINAL PREDICTION ERROR TABLE\n"
)

cat(
  "   \n\n"
)

print(
  as.data.frame(prediction_table)
)


#    
# 14. FINAL COVERAGE / SE TABLE
#    

coverage_table <- manuscript_table %>%
  select(
    Scenario,
    method,
    n,
    Empirical_SE,
    Model_SE,
    SE_Ratio,
    Coverage
  )


cat(
    
)

cat(
  "FINAL COVERAGE / SE TABLE\n"
)

cat(
  "   \n\n"
)

print(
  as.data.frame(coverage_table)
)


#    
# 15. CORRELATION TABLE
#    

correlations <- all_results %>%
  group_by(
    Scenario,
    method
  ) %>%
  summarise(
    
    cor_Em_absbias =
      if (
        sd(E_m) == 0 ||
        sd(abs(bias)) == 0
      ) {
        NA_real_
      } else {
        cor(
          E_m,
          abs(bias)
        )
      },
    
    cor_Eh_absbias =
      if (
        sd(E_h) == 0 ||
        sd(abs(bias)) == 0
      ) {
        NA_real_
      } else {
        cor(
          E_h,
          abs(bias)
        )
      },
    
    cor_D_absbias =
      if (
        sd(D_score) == 0 ||
        sd(abs(bias)) == 0
      ) {
        NA_real_
      } else {
        cor(
          D_score,
          abs(bias)
        )
      },
    
    .groups = "drop"
  )


cat(
    
)

cat(
  "FINAL CORRELATIONS\n"
)

cat(
  "   \n\n"
)

print(
  as.data.frame(correlations)
)


#    
# 16. SAVE ALL FINAL TABLES
#    

write.csv(
  all_results,
  "FINAL_all_simulation_results.csv",
  row.names = FALSE
)

write.csv(
  manuscript_table,
  "FINAL_manuscript_performance_table.csv",
  row.names = FALSE
)

write.csv(
  prediction_table,
  "FINAL_prediction_error_table.csv",
  row.names = FALSE
)

write.csv(
  coverage_table,
  "FINAL_coverage_SE_table.csv",
  row.names = FALSE
)

write.csv(
  correlations,
  "FINAL_correlations.csv",
  row.names = FALSE
)

write.csv(
  oracle_check,
  "FINAL_oracle_sanity_check.csv",
  row.names = FALSE
)


#    
# 17. FINAL FILE CHECK
#    

files_needed <- c(
  "FINAL_all_simulation_results.csv",
  "FINAL_manuscript_performance_table.csv",
  "FINAL_prediction_error_table.csv",
  "FINAL_coverage_SE_table.csv",
  "FINAL_correlations.csv",
  "FINAL_oracle_sanity_check.csv"
)


file_check <- data.frame(
  file = files_needed,
  exists = file.exists(files_needed)
)


 

cat(
  "FINAL FILE CHECK\n"
)

cat(
  "   \n\n"
)

print(
  file_check,
  row.names = FALSE
)


if (!all(file_check$exists)) {
  stop(
    "STOP: One or more final files were not created."
  )
}


#    
# 18. FINAL SUMMARY
 
 

cat(
  "FINAL ANALYSIS COMPLETE\n"
)

cat(
  "   \n\n"
)

cat(
  "Scenario 1: 500 replicates\n"
)

cat(
  "Scenario 2: 500 replicates\n"
)

cat(
  "Methods per scenario: 5\n"
)

cat(
  "Total rows: ",
  nrow(all_results),
  "\n",
  sep = ""
)

cat(
  "\nScenario 1 was NOT rerun.\n"
)

cat(
  "Scenario 2 was NOT rerun.\n"
)

cat(
  "Only existing simulation results were analyzed.\n"
)

cat(
  "\nFINAL OUTPUTS:\n"
)

print(
  files_needed
)

cat(
  "\nWorking directory:\n"
)

cat(
  getwd(),
  "\n"
)

 