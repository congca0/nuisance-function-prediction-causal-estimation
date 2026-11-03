# v4: build Tables 1-3, the numbers quoted in the text, and Figures 1-2
# directly from the saved v3 results (no hand copying).
# v4 changes:
#  - correlations and the remainder-predicted bias are computed from the
#    replicate-level values; in v3, E_m had already been replaced by its mean
#    inside summarise(), so pred_bias was mean(R)/(Var + mean(E_m)) and r_pred
#    was effectively cor(R, bias). The differences are small but the table
#    should be regenerated.
#  - new column r_E: within-method correlation between E_m + E_h and |bias|
#  - Table 3: method ranks by prediction error and by causal performance
# Run from the folder that contains v3_output/ (or from inside v3_output/).

library(dplyr)
library(ggplot2)

dir <- if (dir.exists("v3_output")) "v3_output" else "."

load_scenario <- function(file, label, s2x) {
  e <- new.env()
  load(file.path(dir, file), envir = e)
  e$results %>%
    mutate(Scenario = label, s2x = s2x, replicate = as.integer(replicate))
}

res <- bind_rows(
  load_scenario("revised_simulation_workspace.RData", "1", 1),          # Var(X|C) = 1
  load_scenario("scenario2_final_workspace.RData",    "2", 1 + 0.8^2)   # Var(X|C) = 1.64
)

method_order <- c("Oracle", "OLS", "GAM", "XGBoost", "DML-XGBoost")

summ <- res %>%
  group_by(Scenario, method) %>%
  summarise(
    cor_pred_bias = if (method[1] == "Oracle") NA else cor(R / (s2x + E_m), bias),
    cor_D_absbias = if (method[1] == "Oracle") NA else cor(D_score, abs(bias)),
    cor_E_absbias = if (method[1] == "Oracle") NA else cor(E_m + E_h, abs(bias)),
    pred_bias = mean(R / (s2x + E_m)),
    R_rep = n(),
    Bias = mean(bias),
    RMSE = sqrt(mean(squared_error)),
    EmpSE = sd(tau_hat),
    ModelSE = mean(se),
    SEratio = mean(se) / sd(tau_hat),
    Coverage = mean(coverage),
    MCSE_cov = sqrt(mean(coverage) * (1 - mean(coverage)) / n()),
    E_m = mean(E_m),
    E_h = mean(E_h),
    D_joint = mean(D_score),
    .groups = "drop"
  ) %>%
  mutate(method = factor(method, levels = method_order)) %>%
  arrange(Scenario, method)

write.csv(summ, file.path(dir, "v3_numbers_for_text.csv"), row.names = FALSE)
print(as.data.frame(summ), digits = 4)

#  Table 1 (LaTeX, booktabs)  
f4 <- function(x) sprintf("%.4f", x)
f3 <- function(x) sprintf("%.3f", x)
rows <- with(summ, sprintf("%s & %-11s & %7s & %s & %s & %s & %s & %s \\\\",
                           Scenario, as.character(method), f4(Bias), f4(RMSE),
                           f4(EmpSE), f4(ModelSE), f3(SEratio), f3(Coverage)))
mid <- which(summ$Scenario == "2")[1]
rows <- c(rows[seq_len(mid - 1)], "\\midrule", rows[mid:length(rows)])
mcse <- sprintf("%.3f", max(summ$MCSE_cov[summ$method != "OLS"], na.rm = TRUE))

tex <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Causal estimation performance across simulation scenarios.}",
  "\\label{tab:causal_performance}",
  "\\small",
  "\\begin{tabular}{llrrrrrr}",
  "\\toprule",
  "Scenario & Method & Bias & RMSE & Empirical SE & Mean model SE & SE ratio & Coverage \\\\",
  "\\midrule",
  rows,
  "\\bottomrule",
  "\\end{tabular}",
  "",
  "\\vspace{0.15cm}",
  "",
  "\\begin{minipage}{0.95\\textwidth}",
  "\\small",
  "\\textit{Note.}",
  sprintf(paste0(
    "All estimators used the training sample ($n=1{,}000$); %d replicates per scenario. ",
    "Bias is the Monte Carlo mean of $\\widehat{\\tau}-\\tau$; RMSE is the root mean squared error. ",
    "Empirical SE is the standard deviation of $\\widehat{\\tau}$ across replicates; mean model SE is the ",
    "average estimated standard error; SE ratio is mean model SE divided by empirical SE. Coverage is the ",
    "proportion of nominal 95\\%% confidence intervals containing $\\tau$ (Monte Carlo SE at most %s). ",
    "Scenario~2 uses cluster-level cross-fitting and cluster-robust variance estimation."),
    max(summ$R_rep), mcse),
  "\\end{minipage}",
  "\\end{table}"
)
writeLines(tex, file.path(dir, "table1.tex"))

#   Table 2 (LaTeX): prediction error, remainder-predicted bias, observed bias  
f2 <- function(x) sprintf("%.2f", x)
t2 <- summ %>% filter(method != "Oracle")
rows2 <- with(t2, sprintf("%s & %-11s & %s & %s & %s & %s & %s & %s & %s & %s \\\\",
                          Scenario, as.character(method), f3(E_m), f3(E_h), f3(D_joint),
                          f4(pred_bias), f4(Bias), f2(cor_E_absbias), f2(cor_pred_bias),
                          f2(cor_D_absbias)))
mid2 <- which(t2$Scenario == "2")[1]
rows2 <- c(rows2[seq_len(mid2 - 1)], "\\midrule", rows2[mid2:length(rows2)])
mcse_bias <- t2$EmpSE / sqrt(t2$R_rep)

tex2 <- c(
  "\\begin{table}[htbp]",
  "\\centering",
  "\\caption{Nuisance-function prediction error, remainder-predicted bias, and observed bias.}",
  "\\label{tab:prediction}",
  "\\small",
  "\\begin{tabular}{llrrrrrrrr}",
  "\\toprule",
  "Scenario & Method & $\\mathcal{E}_m$ & $\\mathcal{E}_h$ & $\\mathcal{D}_{\\mathrm{joint}}$ & Predicted bias & Observed bias & $r_{\\mathcal{E}}$ & $r_{\\mathrm{pred}}$ & $r_{\\mathcal{D}}$ \\\\",
  "\\midrule",
  rows2,
  "\\bottomrule",
  "\\end{tabular}",
  "",
  "\\vspace{0.15cm}",
  "",
  "\\begin{minipage}{0.95\\textwidth}",
  "\\small",
  "\\textit{Note.}",
  sprintf(paste0(
    "Means over %d replicates. $\\mathcal{E}_m$, $\\mathcal{E}_h$, and $\\mathcal{D}_{\\mathrm{joint}}$ are computed on the evaluation sample. ",
    "Predicted bias is the mean of $R/\\{\\mathrm{Var}(X\\mid C)+\\mathcal{E}_m\\}$ from~\\eqref{eq:remainder}; ",
    "observed bias is the mean of $\\widehat\\tau-\\tau$ (Monte Carlo SE %.4f to %.4f). ",
    "$r_{\\mathrm{pred}}$ is the correlation across replicates between predicted bias and $\\widehat\\tau-\\tau$; ",
    "$r_{\\mathcal{D}}$ is the correlation between $\\mathcal{D}_{\\mathrm{joint}}$ and $|\\widehat\\tau-\\tau|$."),
    max(t2$R_rep), min(mcse_bias), max(mcse_bias)),
  "\\end{minipage}",
  "\\end{table}"
)
writeLines(tex2, file.path(dir, "table2.tex"))

#   v4: Table 3, ranks (1 = best) among the non-oracle methods  
t3 <- summ %>% filter(method != "Oracle") %>% group_by(Scenario) %>%
  mutate(rk_Em = rank(E_m), rk_Eh = rank(E_h), rk_tot = rank(E_m + E_h),
         rk_bias = rank(abs(Bias)), rk_rmse = rank(RMSE),
         rk_cov = rank(abs(Coverage - 0.95))) %>% ungroup()
print(as.data.frame(t3[, c("Scenario", "method", "rk_Em", "rk_Eh", "rk_tot",
                           "rk_bias", "rk_rmse", "rk_cov")]))

#   Figures (same file names as in the manuscript)  
plot_dat <- res %>%
  filter(method != "Oracle") %>%
  mutate(method = factor(method, levels = method_order[-1]),
         Scenario = paste("Scenario", Scenario),
         abs_err = abs(tau_hat - 1),
         total_pred = E_m + E_h)

p1 <- ggplot(plot_dat, aes(total_pred, abs_err)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.6) +
  facet_grid(Scenario ~ method, scales = "free_x") +
  labs(x = expression("Total nuisance prediction error " * (E[m] + E[h])),
       y = expression("Absolute estimation error " * group("|", hat(tau) - tau, "|"))) +
  theme_bw()
ggsave(file.path(dir, "final_prediction_error_vs_causal_bias.png"), p1,
       width = 11, height = 5.5, dpi = 300)

p2 <- ggplot(plot_dat, aes(D_score, abs_err)) +
  geom_point(alpha = 0.3, size = 0.8) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.6) +
  facet_grid(Scenario ~ method, scales = "free_x") +
  labs(x = expression(D[joint] == group("|", E(delta[m] * delta[h]), "|")),
       y = expression("Absolute estimation error " * group("|", hat(tau) - tau, "|"))) +
  theme_bw()
ggsave(file.path(dir, "final_joint_error_vs_causal_bias.png"), p2,
       width = 11, height = 5.5, dpi = 300)

cat("\nWrote table1.tex, table2.tex, v3_numbers_for_text.csv and both figures to", normalizePath(dir), "\n")
