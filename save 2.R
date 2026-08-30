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

write.csv(
  plot_data,
  "baseline_simulation_plot_data.csv",
  row.names = FALSE
)

ggsave(
  "prediction_error_vs_causal_bias.png",
  p1,
  width = 10,
  height = 7,
  dpi = 300
)

ggsave(
  "joint_error_vs_causal_bias.png",
  p2,
  width = 10,
  height = 7,
  dpi = 300
)

cat("\nSaved to:\n")
cat(getwd(), "\n\n")

print(list.files(
  pattern = "baseline_simulation|prediction_error|joint_error",
  ignore.case = TRUE
))

cat("\nSummary:\n")
print(summary_results)