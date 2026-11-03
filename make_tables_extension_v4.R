# make_tables_extension_v4.R
# Builds Tables 4-6 (Scenarios 3-6) from v4_output/S*_results.csv,
# written by extension_scenarios_v4.R. Run from the folder containing v4_output/.

suppressMessages(library(dplyr))
dir <- "v4_output"

summ_one <- function(sc) {
  f <- file.path(dir, paste0(sc, "_results.csv"))
  if (!file.exists(f)) return(NULL)
  read.csv(f) %>%
    group_by(method) %>%
    summarise(
      r_E    = if (all(E_m == 0)) NA else cor(E_m + E_h, abs(bias)),
      r_pred = if (all(E_m == 0)) NA else cor(pred_bias, bias),
      r_ref  = if (all(is.na(ref_pred_bias))) NA else cor(ref_pred_bias, bias),
      r_gap  = if (all(is.na(gap_vs_ref))) NA else cor(gap_vs_ref, bias),
      R_rep = n(), Bias = mean(bias), MCSE_bias = sd(tau_hat) / sqrt(n()),
      RMSE = sqrt(mean(squared_error)), SEratio = mean(se) / sd(tau_hat),
      Coverage = mean(coverage), E_m = mean(E_m), E_h = mean(E_h),
      pred_bias = mean(pred_bias), cs_bound = mean(cs_bound),
      ref_pred_bias = mean(ref_pred_bias), gap = mean(gap_vs_ref),
      cv_m = mean(cv_mse_m), cv_h = mean(cv_mse_h),
      depth_m = paste(sort(unique(na.omit(depth_m))), collapse = "/"),
      .groups = "drop") %>%
    mutate(Scenario = sub("S", "", sc))
}

all <- bind_rows(lapply(c("S3", "S4", "S5", "S6"), summ_one))
ord <- c("Oracle", "OLS", "GAM", "XGBoost", "DML-XGBoost", "XGBoost-tuned", "DML-XGBoost-tuned")
all <- all %>% mutate(method = factor(method, levels = ord)) %>% arrange(Scenario, method)
write.csv(all, file.path(dir, "extension_numbers_for_text.csv"), row.names = FALSE)
print(as.data.frame(all), digits = 3, width = 250)

f2 <- function(x) ifelse(is.na(x), "--", sprintf("%.2f", x))
f3 <- function(x) ifelse(is.na(x), "--", sprintf("%.3f", x))
f4 <- function(x) ifelse(is.na(x), "--", sprintf("%.4f", x))
mid <- function(rows, sc) {
  br <- which(diff(as.integer(sc)) != 0)
  out <- character(0)
  for (i in seq_along(rows)) { out <- c(out, rows[i]); if (i %in% br) out <- c(out, "\\midrule") }
  out
}

# Table 4: S3-S5 performance
t4 <- all %>% filter(Scenario %in% c("3", "4", "5"))
rows4 <- with(t4, sprintf("%s & %-11s & %s & %s & %s & %s & %s & %s & %s & %s & %s \\\\",
  Scenario, as.character(method),
  ifelse(method == "Oracle", "--", f3(E_m)), ifelse(method == "Oracle", "--", f3(E_h)),
  ifelse(method == "Oracle", "--", f2(cs_bound)), ifelse(method == "Oracle", "--", f4(pred_bias)),
  f4(Bias), f4(RMSE), f3(SEratio), f3(Coverage), f2(r_E)))
writeLines(mid(rows4, t4$Scenario), file.path(dir, "table4_rows.tex"))

# Table 5: observable diagnostics, S3-S6 (learners compared with the DML-XGBoost reference)
refb <- all %>% filter(method == "DML-XGBoost") %>% select(Scenario, Bias_ref = Bias)
t5 <- all %>% filter(!is.na(ref_pred_bias)) %>% left_join(refb, by = "Scenario")
rows5 <- with(t5, sprintf("%s & %-11s & %s & %s & %s & %s & %s & %s \\\\",
  Scenario, as.character(method), f4(Bias), f4(Bias_ref), f4(gap), f2(r_gap),
  f4(ref_pred_bias), f2(r_ref)))
writeLines(mid(rows5, t5$Scenario), file.path(dir, "table5_rows.tex"))

# Table 6: S6 tuning
t6 <- all %>% filter(Scenario == "6", method != "Oracle")
rows6 <- with(t6, sprintf("%-17s & %s & %s & %s & %s & %s & %s & %s & %s \\\\",
  as.character(method), f3(cv_m), f3(cv_h), f3(E_m), f3(E_h), f4(Bias), f4(RMSE),
  f3(SEratio), f3(Coverage)))
writeLines(rows6, file.path(dir, "table6_rows.tex"))

mc <- all %>% group_by(Scenario) %>% summarise(R = max(R_rep),
  mcse_bias_lo = min(MCSE_bias), mcse_bias_hi = max(MCSE_bias),
  mcse_cov = max(sqrt(Coverage * (1 - Coverage) / R_rep)[method != "OLS"]))
print(as.data.frame(mc), digits = 3)
