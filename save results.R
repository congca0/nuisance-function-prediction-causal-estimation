 results <- dplyr::bind_rows(
  results_list,
  .id = "replicate"
)
 
summary_results <- results %>%
  dplyr::group_by(method) %>%
  dplyr::summarise(
    mean_E_m = mean(E_m),
    mean_E_h = mean(E_h),
    mean_D_score = mean(D_score),
    bias = mean(bias),
    RMSE = sqrt(mean(squared_error)),
    empirical_SE = sd(tau_hat),
    mean_model_SE = mean(se),
    coverage = mean(coverage),
    .groups = "drop"
  )
 
write.csv(
  results,
  "baseline_simulation_results.csv",
  row.names = FALSE
)

write.csv(
  summary_results,
  "baseline_simulation_summary.csv",
  row.names = FALSE
)

 
cat("Results saved to:\n")
cat(getwd(), "\n\n")

cat("Files:\n")
print(list.files(
  pattern = "baseline_simulation",
  ignore.case = TRUE
))

 print(summary_results)