# v3 changes (Oct 2026)

# FINAL SCENARIO 2
# TRUE CLUSTERED DGP
# 500 REPLICATES
 


rm(list = ls())

set.seed(20260828)

library(mgcv)
library(xgboost)
library(dplyr)
library(ggplot2)


 
# 1. SETTINGS
 

N_SIM <- 500

n_train <- 1000
n_eval  <- 10000

p <- 10

cluster_size_train <- 10
cluster_size_eval  <- 10

tau_true <- 1

rho <- 0.3

sigma_x <- 1
sigma_y <- 1

sigma_cluster_x <- 0.8
sigma_cluster_y <- 0.8

nfold <- 5

# v3: Var(X | C), used for the remainder-predicted bias R / (Var(X | C) + E_m)
sigma2_x_given_C <- sigma_x^2 + sigma_cluster_x^2

if (basename(getwd()) != "v3_output") {
  dir.create("v3_output", showWarnings = FALSE)
  setwd("v3_output")
}


 
# 2. COVARIATES
 

make_covariates <- function(
    n,
    p,
    rho = 0.3
) {
  
  Sigma <- matrix(
    rho,
    p,
    p
  )
  
  diag(Sigma) <- 1
  
  Z <- matrix(
    rnorm(n * p),
    nrow = n,
    ncol = p
  )
  
  C <- Z %*% chol(Sigma)
  
  colnames(C) <- paste0(
    "C",
    1:p
  )
  
  C
}


 


m0_fun <- function(C) {
  
  0.8 * C[, 1] +
    0.5 * C[, 2] +
    0.4 * sin(C[, 3]) +
    0.3 * C[, 4] * C[, 5]
}


g0_fun <- function(C) {
  
  0.7 * C[, 1] -
    0.5 * C[, 2] +
    0.4 * C[, 3]^2 +
    0.3 * C[, 4] * C[, 5]
}





generate_clustered_data <- function(
    n,
    p,
    cluster_size,
    tau = 1,
    rho = 0.3,
    sigma_x = 1,
    sigma_y = 1,
    sigma_cluster_x = 0.8,
    sigma_cluster_y = 0.8
) {

  n_cluster <- ceiling(n / cluster_size)
  cluster_id <- rep(seq_len(n_cluster), each = cluster_size)
  cluster_id <- cluster_id[seq_len(n)]

  C <- make_covariates(n = n, p = p, rho = rho)

  m0 <- m0_fun(C)
  g0 <- g0_fun(C)
  h0 <- tau * m0 + g0          # v3 (C2): E(Y | C); v2 used g0 here

  # Cluster-level random effects
  u_x <- rnorm(n_cluster, mean = 0, sd = sigma_cluster_x)
  u_y <- rnorm(n_cluster, mean = 0, sd = sigma_cluster_y)

  # Individual-level errors
  epsilon_x <- rnorm(n, mean = 0, sd = sigma_x)
  epsilon_y <- rnorm(n, mean = 0, sd = sigma_y)

  # Treatment / exposure
  X <- m0 + u_x[cluster_id] + epsilon_x

  # Outcome (unchanged DGP: Y = tau X + g0 + u_y + e_y)
  Y <- tau * X + g0 + u_y[cluster_id] + epsilon_y

  list(
    C = C,
    X = X,
    Y = Y,
    m0 = m0,
    g0 = g0,
    h0 = h0,
    cluster_id = cluster_id,
    u_x = u_x,
    u_y = u_y
  )
}


 


prediction_metrics <- function(
    pred_m,
    pred_h,
    true_m,
    true_h
) {

  delta_m <- pred_m - true_m
  delta_h <- pred_h - true_h

  list(
    E_m = mean(delta_m^2),
    E_h = mean(delta_h^2),
    D_score = abs(mean(delta_m * delta_h)),              # v2 measure, kept
    cross = mean(delta_m * delta_h),                     # v3: signed
    R = mean(delta_m * (delta_h - tau_true * delta_m)),  # v3: remainder
    delta_m = delta_m,
    delta_h = delta_h
  )
}


 


fit_ols <- function(
    train,
    eval
) {
  
  train_m <- data.frame(
    X = train$X,
    train$C
  )
  
  train_h <- data.frame(
    Y = train$Y,
    train$C
  )
  
  eval_df <- data.frame(
    eval$C
  )
  
  
  fit_m <- lm(
    X ~ .,
    data = train_m
  )
  
  fit_h <- lm(
    Y ~ .,
    data = train_h
  )
  
  
  pred_m <- predict(
    fit_m,
    newdata = eval_df
  )
  
  pred_h <- predict(
    fit_h,
    newdata = eval_df
  )
  
  
  list(
    m = as.numeric(pred_m),
    h = as.numeric(pred_h)
  )
}


 


fit_gam <- function(
    train,
    eval
) {
  
  train_m <- data.frame(
    X = train$X,
    train$C
  )
  
  train_h <- data.frame(
    Y = train$Y,
    train$C
  )
  
  eval_df <- data.frame(
    eval$C
  )
  
  
  form_m <- as.formula(
    paste(
      "X ~",
      paste0(
        "s(C",
        1:p,
        ")",
        collapse = " + "
      )
    )
  )
  
  
  form_h <- as.formula(
    paste(
      "Y ~",
      paste0(
        "s(C",
        1:p,
        ")",
        collapse = " + "
      )
    )
  )
  
  
  fit_m <- mgcv::gam(
    form_m,
    data = train_m,
    method = "REML"
  )
  
  fit_h <- mgcv::gam(
    form_h,
    data = train_h,
    method = "REML"
  )
  
  
  pred_m <- predict(
    fit_m,
    newdata = eval_df
  )
  
  pred_h <- predict(
    fit_h,
    newdata = eval_df
  )
  
  
  list(
    m = as.numeric(pred_m),
    h = as.numeric(pred_h)
  )
}


 

 

fit_xgb <- function(
    train,
    eval
) {
  
  dtrain_m <- xgboost::xgb.DMatrix(
    data = train$C,
    label = train$X
  )
  
  dtrain_h <- xgboost::xgb.DMatrix(
    data = train$C,
    label = train$Y
  )
  
  dtest <- xgboost::xgb.DMatrix(
    data = eval$C
  )
  
  
  params <- list(
    objective = "reg:squarederror",
    eta = 0.05,
    max_depth = 4,
    subsample = 0.8,
    colsample_bytree = 0.8
  )
  
  
  fit_m <- xgboost::xgb.train(
    params = params,
    data = dtrain_m,
    nrounds = 300,
    verbose = 0
  )
  
  
  fit_h <- xgboost::xgb.train(
    params = params,
    data = dtrain_h,
    nrounds = 300,
    verbose = 0
  )
  
  
  pred_m <- predict(
    fit_m,
    dtest
  )
  
  pred_h <- predict(
    fit_h,
    dtest
  )
  
  
  list(
    m = as.numeric(pred_m),
    h = as.numeric(pred_h)
  )
}


 

 

estimate_tau_cluster <- function(
    X,
    Y,
    pred_m,
    pred_h,
    cluster_id
) {
  
  X_res <- X -
    pred_m
  
  Y_res <- Y -
    pred_h
  
  
  tau_hat <- sum(
    X_res * Y_res
  ) /
    sum(
      X_res^2
    )
  
  
  psi <- X_res *
    (
      Y_res -
        tau_hat * X_res
    )
  
  
  denom <- mean(
    X_res^2
  )
  
  
  cluster_sums <- tapply(
    psi,
    cluster_id,
    sum
  )
  
  
  n <- length(X)
  
  G <- length(
    cluster_sums
  )
  
  
  variance <- sum(
    cluster_sums^2
  ) /
    (
      n^2 *
        denom^2
    )
  
  
  finite_sample_correction <- G /
    (G - 1)
  
  
  variance <- variance *
    finite_sample_correction
  
  
  se <- sqrt(
    variance
  )
  
  
  ci_low <- tau_hat -
    1.96 * se
  
  ci_high <- tau_hat +
    1.96 * se
  
  
  list(
    tau = tau_hat,
    se = se,
    ci_low = ci_low,
    ci_high = ci_high
  )
}


 
 

dml_xgb_cluster <- function(
    dat,
    nfold = 5,
    eval = NULL
) {

  cluster_ids <- unique(dat$cluster_id)
  cluster_folds <- sample(rep(seq_len(nfold), length.out = length(cluster_ids)))
  names(cluster_folds) <- cluster_ids
  folds <- cluster_folds[as.character(dat$cluster_id)]

  m_hat <- numeric(length(dat$X))
  h_hat <- numeric(length(dat$X))
  fold_metrics <- vector("list", nfold)

  for (k in seq_len(nfold)) {

    train_idx <- which(folds != k)
    test_idx <- which(folds == k)

    train <- list(
      C = dat$C[train_idx, , drop = FALSE],
      X = dat$X[train_idx],
      Y = dat$Y[train_idx]
    )

    # v3 (C3): one fit per fold, predicting the held-out fold AND the
    # evaluation sample, so prediction error is measured out of sample
    newC <- dat$C[test_idx, , drop = FALSE]
    if (!is.null(eval)) newC <- rbind(newC, eval$C)

    pred <- fit_xgb(train, list(C = newC))

    nt <- length(test_idx)
    m_hat[test_idx] <- pred$m[seq_len(nt)]
    h_hat[test_idx] <- pred$h[seq_len(nt)]

    if (!is.null(eval)) {
      pm <- prediction_metrics(pred$m[-seq_len(nt)], pred$h[-seq_len(nt)],
                               eval$m0, eval$h0)
      fold_metrics[[k]] <- unlist(pm[c("E_m", "E_h", "D_score", "cross", "R")])
    }
  }

  estimate <- estimate_tau_cluster(X = dat$X, Y = dat$Y, pred_m = m_hat,
                                   pred_h = h_hat, cluster_id = dat$cluster_id)

  eval_metrics <- if (is.null(eval)) NULL else colMeans(do.call(rbind, fold_metrics))

  list(
    tau = estimate$tau,
    se = estimate$se,
    ci_low = estimate$ci_low,
    ci_high = estimate$ci_high,
    m_hat = m_hat,
    h_hat = h_hat,
    folds = folds,
    eval_metrics = eval_metrics
  )
}


 

 


make_row <- function(method, metrics, est) {
  data.frame(
    method = method,
    E_m = metrics[["E_m"]],
    E_h = metrics[["E_h"]],
    D_score = metrics[["D_score"]],
    cross = metrics[["cross"]],
    R = metrics[["R"]],
    tau_hat = est$tau,
    se = est$se,
    bias = est$tau - tau_true,
    squared_error = (est$tau - tau_true)^2,
    coverage = tau_true >= est$ci_low && tau_true <= est$ci_high
  )
}

run_one <- function(
    seed
) {

  set.seed(seed)

  train <- generate_clustered_data(
    n = n_train, p = p, cluster_size = cluster_size_train, tau = tau_true,
    rho = rho, sigma_x = sigma_x, sigma_y = sigma_y,
    sigma_cluster_x = sigma_cluster_x, sigma_cluster_y = sigma_cluster_y
  )

  eval <- generate_clustered_data(
    n = n_eval, p = p, cluster_size = cluster_size_eval, tau = tau_true,
    rho = rho, sigma_x = sigma_x, sigma_y = sigma_y,
    sigma_cluster_x = sigma_cluster_x, sigma_cluster_y = sigma_cluster_y
  )


  oracle_metrics <- prediction_metrics(eval$m0, eval$h0, eval$m0, eval$h0)
  oracle_tau <- estimate_tau_cluster(
    X = train$X, Y = train$Y, pred_m = train$m0, pred_h = train$h0,
    cluster_id = train$cluster_id
  )
  rows <- list(make_row("Oracle", oracle_metrics, oracle_tau))


  both <- list(C = rbind(train$C, eval$C))
  tr <- seq_len(nrow(train$C))

  for (method in c("OLS", "GAM", "XGBoost")) {

    pred <- switch(method,
      OLS = fit_ols(train, both),
      GAM = fit_gam(train, both),
      XGBoost = fit_xgb(train, both)
    )

    metrics <- prediction_metrics(pred$m[-tr], pred$h[-tr], eval$m0, eval$h0)

    est <- estimate_tau_cluster(
      X = train$X, Y = train$Y, pred_m = pred$m[tr], pred_h = pred$h[tr],
      cluster_id = train$cluster_id
    )

    rows[[length(rows) + 1]] <- make_row(method, metrics, est)
  }

  # v3 (C3): DML prediction error on the evaluation sample, averaged over folds
  dml <- dml_xgb_cluster(train, nfold = nfold, eval = eval)
  rows[[length(rows) + 1]] <- make_row("DML-XGBoost", dml$eval_metrics, dml)

  bind_rows(rows)
}



 


required_functions <- c(
  "make_covariates",
  "m0_fun",
  "g0_fun",
  "generate_clustered_data",
  "prediction_metrics",
  "fit_ols",
  "fit_gam",
  "fit_xgb",
  "estimate_tau_cluster",
  "dml_xgb_cluster",
  "run_one"
)


function_check <- sapply(
  required_functions,
  exists,
  mode = "function"
)


print(
  function_check
)


if (!all(function_check)) {
  
  stop(
    "Required function missing. Simulation NOT started."
  )
}


cat(
  "\nAll required functions are available.\n"
)


 
# 13. ONE-REPLICATE TEST
 

 


test_result <- run_one(
  seed = 20001
)


print(
  as.data.frame(test_result)
)


if (
  nrow(test_result) != 5 ||
  any(
    !is.finite(
      test_result$tau_hat
    )
  ) ||
  any(
    !is.finite(
      test_result$se
    )
  )
) {
  
  stop(
    "One-replicate test failed. 500-replicate simulation NOT started."
  )
}


cat(
  "\nOne-replicate test PASSED.\n"
)


 
# 14. FINAL 500-REPLICATE SIMULATION
 

results_list <- vector(
  "list",
  N_SIM
)


 

cat(
  "STARTING FINAL SCENARIO 2\n"
)

cat(
  "TRUE CLUSTERED DGP\n"
)

cat(
  "N_SIM = ",
  N_SIM,
  "\n",
  sep = ""
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
  
  
  if (
    s %% 50 == 0 ||
    s == N_SIM
  ) {
    
    save(
      results_list,
      file =
        "scenario2_final_checkpoint.RData"
    )
    
    cat(
      "Checkpoint saved at replicate",
      s,
      "\n"
    )
  }
}


 
# 15. COMBINE RESULTS
 

results <- bind_rows(
  results_list,
  .id = "replicate"
)


results$replicate <- as.integer(
  results$replicate
)


 
# 16. SUMMARY
 

summary_results <- results %>%
  group_by(method) %>%
  summarise(
    n = n(),
    mean_E_m = mean(E_m),
    mean_E_h = mean(E_h),
    mean_D_score = mean(D_score),
    mean_cross = mean(cross),
    mean_R = mean(R),
    pred_bias = mean(R / (sigma2_x_given_C + E_m)),  # v3: remainder-predicted bias
    abs_bias = mean(abs(bias)),                      # v3: before bias is summarised
    mcse_bias = sd(tau_hat) / sqrt(n()),
    bias = mean(bias),
    RMSE = sqrt(mean(squared_error)),
    empirical_SE = sd(tau_hat),
    mean_model_SE = mean(se),
    mcse_coverage = sqrt(mean(coverage) * (1 - mean(coverage)) / n()),
    coverage = mean(coverage),
    .groups = "drop"
  )

 
cat(
  "FINAL SCENARIO 2 SUMMARY\n"
)

 

print(
  as.data.frame(summary_results)
)


 
# 17. DATA QUALITY CHECK
 
 


cat(
  "Missing values:\n"
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
        
        sum(
          !is.finite(x)
        )
        
      } else {
        
        0
      }
    }
  )
)


 
# 18. SAVE RAW RESULTS
 

write.csv(
  results,
  "scenario2_final_results.csv",
  row.names = FALSE
)


write.csv(
  summary_results,
  "scenario2_final_summary.csv",
  row.names = FALSE
)


save(
  results,
  results_list,
  summary_results,
  file =
    "scenario2_final_workspace.RData"
)


 
# 19. PLOT DATA
 

plot_data <- results %>%
  
  mutate(
    
    abs_bias =
      abs(bias),
    
    total_prediction_error =
      E_m + E_h
  )


 
# 20. PLOT 1
 

p1 <- ggplot(
  
  plot_data,
  
  aes(
    x =
      total_prediction_error,
    
    y =
      abs_bias
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
    
    x =
      "Overall nuisance prediction error",
    
    y =
      "Absolute causal bias"
  ) +
  
  theme_bw()


ggsave(
  "scenario2_final_prediction_error_vs_bias.png",
  p1,
  width = 10,
  height = 7,
  dpi = 300
)


 
# 21. PLOT 2
 

p2 <- ggplot(
  
  plot_data,
  
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
    
    x =
      "|E(delta_m delta_h)|",
    
    y =
      "Absolute causal bias"
  ) +
  
  theme_bw()


ggsave(
  "scenario2_final_Dscore_vs_bias.png",
  p2,
  width = 10,
  height = 7,
  dpi = 300
)


 

 

bias_summary <- plot_data %>%
  
  group_by(method) %>%
  
  summarise(
    
    mean_abs_bias =
      mean(abs_bias),
    
    .groups = "drop"
  )


p3 <- ggplot(
  
  bias_summary,
  
  aes(
    x =
      reorder(
        method,
        mean_abs_bias
      ),
    
    y =
      mean_abs_bias
  )
) +
  
  geom_col() +
  
  labs(
    
    x = NULL,
    
    y =
      "Mean absolute causal bias"
  ) +
  
  theme_bw()


ggsave(
  "scenario2_final_bias_by_method.png",
  p3,
  width = 9,
  height = 6,
  dpi = 300
)


 
# 23. CORRELATIONS
 

correlations <- results %>%
  
  group_by(method) %>%
  
  summarise(
    
    cor_Em_absbias =
      ifelse(
        
        sd(E_m) == 0 ||
          sd(abs(bias)) == 0,
        
        NA_real_,
        
        cor(
          E_m,
          abs(bias)
        )
      ),
    
    
    cor_Eh_absbias =
      ifelse(
        
        sd(E_h) == 0 ||
          sd(abs(bias)) == 0,
        
        NA_real_,
        
        cor(
          E_h,
          abs(bias)
        )
      ),
    
    
    cor_D_absbias =
      ifelse(
        
        sd(D_score) == 0 ||
          sd(abs(bias)) == 0,
        
        NA_real_,
        
        cor(
          D_score,
          abs(bias)
        )
      ),
    
    .groups = "drop"
  )


print(
  as.data.frame(correlations)
)


write.csv(
  correlations,
  "scenario2_final_correlations.csv",
  row.names = FALSE
)


 
# 24. REGRESSION: E_m + E_h
 

lm_E_results <- results %>%
  
  group_by(method) %>%
  
  group_modify(
    ~ {
      
      if (
        nrow(.x) < 10 ||
        sd(.x$E_m) == 0 ||
        sd(.x$E_h) == 0
      ) {
        
        return(
          data.frame(
            intercept = NA_real_,
            E_m = NA_real_,
            E_h = NA_real_,
            R2 = NA_real_,
            p_value = NA_real_
          )
        )
      }
      
      
      fit <- lm(
        abs(bias) ~ E_m + E_h,
        data = .x
      )
      
      
      sm <- summary(fit)
      
      
      data.frame(
        
        intercept =
          coef(fit)[1],
        
        E_m =
          coef(fit)[2],
        
        E_h =
          coef(fit)[3],
        
        R2 =
          sm$r.squared,
        
        p_value =
          pf(
            sm$fstatistic[1],
            sm$fstatistic[2],
            sm$fstatistic[3],
            lower.tail = FALSE
          )
      )
    }
  )


print(
  as.data.frame(lm_E_results)
)


write.csv(
  lm_E_results,
  "scenario2_final_lm_prediction_error.csv",
  row.names = FALSE
)


 
# 25. REGRESSION: D SCORE
 

lm_D_results <- results %>%
  
  group_by(method) %>%
  
  group_modify(
    ~ {
      
      if (
        nrow(.x) < 10 ||
        sd(.x$D_score) == 0
      ) {
        
        return(
          data.frame(
            intercept = NA_real_,
            D_score = NA_real_,
            R2 = NA_real_,
            p_value = NA_real_
          )
        )
      }
      
      
      fit <- lm(
        abs(bias) ~ D_score,
        data = .x
      )
      
      
      sm <- summary(fit)
      
      
      data.frame(
        
        intercept =
          coef(fit)[1],
        
        D_score =
          coef(fit)[2],
        
        R2 =
          sm$r.squared,
        
        p_value =
          pf(
            sm$fstatistic[1],
            sm$fstatistic[2],
            sm$fstatistic[3],
            lower.tail = FALSE
          )
      )
    }
  )


print(
  as.data.frame(lm_D_results)
)


write.csv(
  lm_D_results,
  "scenario2_final_lm_Dscore.csv",
  row.names = FALSE
)


 

 
performance_table <- summary_results %>%
  
  mutate(
    
    Bias =
      round(
        bias,
        4
      ),
    
    Abs_Bias =
      round(
        abs_bias,
        4
      ),
    
    RMSE =
      round(
        RMSE,
        4
      ),
    
    Empirical_SE =
      round(
        empirical_SE,
        4
      ),
    
    Model_SE =
      round(
        mean_model_SE,
        4
      ),
    
    SE_Ratio =
      round(
        mean_model_SE /
          empirical_SE,
        3
      ),
    
    Coverage =
      round(
        coverage,
        3
      )
  ) %>%
  
  select(
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
  as.data.frame(
    performance_table
  )
)


write.csv(
  performance_table,
  "scenario2_final_performance_table.csv",
  row.names = FALSE
)


 

prediction_table <- summary_results %>%
  
  mutate(
    
    E_m =
      round(
        mean_E_m,
        4
      ),
    
    E_h =
      round(
        mean_E_h,
        4
      ),
    
    D_score =
      round(
        mean_D_score,
        4
      )
  ) %>%
  
  select(
    method,
    E_m,
    E_h,
    D_score
  )


print(
  as.data.frame(
    prediction_table
  )
)


write.csv(
  prediction_table,
  "scenario2_final_prediction_table.csv",
  row.names = FALSE
)



 
coverage_check <- summary_results %>%
  
  mutate(
    
    SE_Ratio =
      mean_model_SE /
      empirical_SE,
    
    Coverage_OK =
      coverage >= 0.90 &
      coverage <= 0.99,
    
    SE_OK =
      SE_Ratio >= 0.80 &
      SE_Ratio <= 1.20
  ) %>%
  
  select(
    method,
    n,
    bias,
    abs_bias,
    RMSE,
    empirical_SE,
    mean_model_SE,
    SE_Ratio,
    coverage,
    Coverage_OK,
    SE_OK
  )


print(
  as.data.frame(
    coverage_check
  )
)


write.csv(
  coverage_check,
  "scenario2_final_coverage_check.csv",
  row.names = FALSE
)



 

oracle_check <- summary_results %>%
  
  filter(
    method == "Oracle"
  ) %>%
  
  mutate(
    
    Bias_OK =
      abs(bias) < 0.02,
    
    Coverage_OK =
      coverage >= 0.90 &
      coverage <= 0.99,
    
    SE_Ratio =
      mean_model_SE /
      empirical_SE,
    
    SE_OK =
      SE_Ratio >= 0.80 &
      SE_Ratio <= 1.20
  )


print(
  as.data.frame(
    oracle_check
  )
)


write.csv(
  oracle_check,
  "scenario2_final_oracle_sanity_check.csv",
  row.names = FALSE
)


 
files_needed <- c(
  
  "scenario2_final_results.csv",
  
  "scenario2_final_summary.csv",
  
  "scenario2_final_workspace.RData",
  
  "scenario2_final_performance_table.csv",
  
  "scenario2_final_prediction_table.csv",
  
  "scenario2_final_coverage_check.csv",
  
  "scenario2_final_correlations.csv",
  
  "scenario2_final_lm_prediction_error.csv",
  
  "scenario2_final_lm_Dscore.csv",
  
  "scenario2_final_oracle_sanity_check.csv",
  
  "scenario2_final_prediction_error_vs_bias.png",
  
  "scenario2_final_Dscore_vs_bias.png",
  
  "scenario2_final_bias_by_method.png"
)


file_check <- data.frame(
  
  file =
    files_needed,
  
  exists =
    file.exists(
      files_needed
    )
)


 
cat(
  "FINAL FILE CHECK\n"
)
 


print(
  file_check,
  row.names = FALSE
)


 # 31. FINAL QC
 
if (
  any(
    is.na(
      results
    )
  )
) {
  
  warning(
    "There are missing values in results."
  )
}


numeric_nonfinite <- sapply(
  
  results,
  
  function(x) {
    
    if (is.numeric(x)) {
      
      sum(
        !is.finite(x)
      )
      
    } else {
      
      0
    }
  }
)


if (
  any(
    numeric_nonfinite > 0
  )
) {
  
  warning(
    "Non-finite values detected."
  )
}


if (
  !all(
    file_check$exists
  )
) {
  
  warning(
    "One or more expected output files are missing."
  )
}


 

cat(
  "FINAL SCENARIO 2 ANALYSIS FINISHED.\n"
)

 


cat(
  "Working directory:\n"
)

cat(
  getwd(),
  "\n\n"
)


cat(
  "Number of replicates:\n"
)

print(
  unique(
    summary_results$n
  )
)


cat(
  "\nMethods:\n"
)

print(
  summary_results$method
)


cat(
  "\nAll final outputs were written to the working directory.\n"
)

cat(
  "Scenario 1 was NOT rerun.\n"
)

cat(
  "Scenario 2 used the TRUE CLUSTERED DGP.\n"
)
 
