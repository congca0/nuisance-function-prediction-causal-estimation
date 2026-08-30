rm(list = ls())
set.seed(20260827)

library(mgcv)
library(xgboost)
library(dplyr)
library(ggplot2)


TEST_MODE <- TRUE

N_SIM <- if (TEST_MODE) 5 else 500

n_train <- if (TEST_MODE) 1000 else 2000
n_eval <- 10000

p <- 10
tau_true <- 1

rho <- 0.3
sigma_x <- 1
sigma_y <- 1

nfold <- 5

clustered <- TRUE
n_cluster <- 100


make_covariates <- function(n, p, rho = 0.3) {
  
  Sigma <- matrix(rho, p, p)
  diag(Sigma) <- 1
  
  Z <- matrix(
    rnorm(n * p),
    nrow = n,
    ncol = p
  )
  
  C <- Z %*% chol(Sigma)
  
  colnames(C) <- paste0("C", 1:p)
  
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


generate_data <- function(
    n,
    p,
    tau = 1,
    rho = 0.3,
    sigma_x = 1,
    sigma_y = 1,
    clustered = FALSE,
    n_cluster = 100) {
  
  C <- make_covariates(
    n = n,
    p = p,
    rho = rho
  )
  
  m0 <- m0_fun(C)
  g0 <- g0_fun(C)
  
  if (clustered) {
    
    cluster_id <- sample(
      seq_len(n_cluster),
      size = n,
      replace = TRUE
    )
    
    cluster_effect_x <- rnorm(
      n_cluster,
      0,
      0.4
    )
    
    cluster_effect_y <- rnorm(
      n_cluster,
      0,
      0.4
    )
    
    X <- m0 +
      cluster_effect_x[cluster_id] +
      rnorm(n, 0, sigma_x)
    
    Y <- tau * X +
      g0 +
      cluster_effect_y[cluster_id] +
      rnorm(n, 0, sigma_y)
    
  } else {
    
    cluster_id <- rep(
      NA_integer_,
      n
    )
    
    X <- m0 +
      rnorm(n, 0, sigma_x)
    
    Y <- tau * X +
      g0 +
      rnorm(n, 0, sigma_y)
  }
  
  h0 <- tau * m0 + g0
  
  list(
    C = C,
    X = X,
    Y = Y,
    m0 = m0,
    h0 = h0,
    g0 = g0,
    cluster_id = cluster_id
  )
}


prediction_metrics <- function(
    pred_m,
    pred_h,
    true_m,
    true_h,
    C) {
  
  delta_m <- pred_m - true_m
  delta_h <- pred_h - true_h
  
  region <- ifelse(
    C[, 1] >= median(C[, 1]),
    "High_C1",
    "Low_C1"
  )
  
  regional <- data.frame(
    region = region,
    delta_m = delta_m,
    delta_h = delta_h
  ) %>%
    group_by(region) %>%
    summarise(
      E_m = mean(delta_m^2),
      E_h = mean(delta_h^2),
      D_score = abs(mean(delta_m * delta_h)),
      .groups = "drop"
    )
  
  list(
    E_m = mean(delta_m^2),
    E_h = mean(delta_h^2),
    D_score = abs(mean(delta_m * delta_h)),
    regional = regional
  )
}


fit_ols <- function(train, newdata) {
  
  train_m <- data.frame(
    X = train$X,
    train$C
  )
  
  train_h <- data.frame(
    Y = train$Y,
    train$C
  )
  
  newdata <- data.frame(newdata)
  
  fit_m <- lm(
    X ~ .,
    data = train_m
  )
  
  fit_h <- lm(
    Y ~ .,
    data = train_h
  )
  
  list(
    m = as.numeric(
      predict(fit_m, newdata)
    ),
    h = as.numeric(
      predict(fit_h, newdata)
    )
  )
}


fit_gam <- function(train, newdata) {
  
  train_m <- data.frame(
    X = train$X,
    train$C
  )
  
  train_h <- data.frame(
    Y = train$Y,
    train$C
  )
  
  newdata <- data.frame(newdata)
  
  form_m <- as.formula(
    paste(
      "X ~",
      paste0(
        "s(C", 1:p, ")",
        collapse = " + "
      )
    )
  )
  
  form_h <- as.formula(
    paste(
      "Y ~",
      paste0(
        "s(C", 1:p, ")",
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
  
  list(
    m = as.numeric(
      predict(fit_m, newdata)
    ),
    h = as.numeric(
      predict(fit_h, newdata)
    )
  )
}


fit_xgb <- function(train, newdata) {
  
  dtrain_m <- xgboost::xgb.DMatrix(
    data = train$C,
    label = train$X
  )
  
  dtrain_h <- xgboost::xgb.DMatrix(
    data = train$C,
    label = train$Y
  )
  
  dtest <- xgboost::xgb.DMatrix(
    data = newdata
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
  
  list(
    m = as.numeric(
      predict(
        fit_m,
        xgboost::xgb.DMatrix(newdata)
      )
    ),
    h = as.numeric(
      predict(
        fit_h,
        xgboost::xgb.DMatrix(newdata)
      )
    )
  )
}


estimate_tau <- function(
    X,
    Y,
    m_hat,
    h_hat) {
  
  X_res <- X - m_hat
  Y_res <- Y - h_hat
  
  tau_hat <- sum(
    X_res * Y_res
  ) / sum(
    X_res^2
  )
  
  psi <- X_res *
    (Y_res - tau_hat * X_res)
  
  n <- length(X)
  
  J <- mean(X_res^2)
  
  variance <- mean(psi^2) /
    J^2
  
  se <- sqrt(
    variance / n
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


crossfit_xgb <- function(
    dat,
    nfold = 5) {
  
  n <- length(dat$X)
  
  fold_id <- sample(
    rep(
      seq_len(nfold),
      length.out = n
    )
  )
  
  m_hat <- numeric(n)
  h_hat <- numeric(n)
  
  for (k in seq_len(nfold)) {
    
    train_id <- which(
      fold_id != k
    )
    
    test_id <- which(
      fold_id == k
    )
    
    train <- list(
      C = dat$C[
        train_id,
        ,
        drop = FALSE
      ],
      X = dat$X[train_id],
      Y = dat$Y[train_id]
    )
    
    test_C <- dat$C[
      test_id,
      ,
      drop = FALSE
    ]
    
    pred <- fit_xgb(
      train = train,
      newdata = test_C
    )
    
    m_hat[test_id] <- pred$m
    h_hat[test_id] <- pred$h
  }
  
  tau <- estimate_tau(
    X = dat$X,
    Y = dat$Y,
    m_hat = m_hat,
    h_hat = h_hat
  )
  
  list(
    m_hat = m_hat,
    h_hat = h_hat,
    tau = tau
  )
}


evaluate_method <- function(
    method,
    train,
    eval) {
  
  if (method == "OLS") {
    
    pred <- fit_ols(
      train,
      eval$C
    )
    
  } else if (method == "GAM") {
    
    pred <- fit_gam(
      train,
      eval$C
    )
    
  } else if (method == "XGBoost") {
    
    pred <- fit_xgb(
      train,
      eval$C
    )
    
  } else {
    
    stop(
      "Unknown method"
    )
  }
  
  metrics <- prediction_metrics(
    pred_m = pred$m,
    pred_h = pred$h,
    true_m = eval$m0,
    true_h = eval$h0,
    C = eval$C
  )
  
  tau <- estimate_tau(
    X = eval$X,
    Y = eval$Y,
    m_hat = pred$m,
    h_hat = pred$h
  )
  
  list(
    metrics = metrics,
    tau = tau
  )
}


run_one <- function(seed) {
  
  set.seed(seed)
  
  train <- generate_data(
    n = n_train,
    p = p,
    tau = tau_true,
    rho = rho,
    sigma_x = sigma_x,
    sigma_y = sigma_y,
    clustered = clustered,
    n_cluster = n_cluster
  )
  
  eval <- generate_data(
    n = n_eval,
    p = p,
    tau = tau_true,
    rho = rho,
    sigma_x = sigma_x,
    sigma_y = sigma_y,
    clustered = clustered,
    n_cluster = n_cluster
  )
  
  
  out <- list()
  
  
  for (method in c(
    "OLS",
    "GAM",
    "XGBoost"
  )) {
    
    fit <- evaluate_method(
      method,
      train,
      eval
    )
    
    out[[method]] <- data.frame(
      method = method,
      E_m = fit$metrics$E_m,
      E_h = fit$metrics$E_h,
      D_score = fit$metrics$D_score,
      tau_hat = fit$tau$tau,
      se = fit$tau$se,
      bias = fit$tau$tau -
        tau_true,
      squared_error =
        (fit$tau$tau -
           tau_true)^2,
      coverage =
        tau_true >=
        fit$tau$ci_low &&
        tau_true <=
        fit$tau$ci_high
    )
  }
  
  
  dml <- crossfit_xgb(
    dat = train,
    nfold = nfold
  )
  
  dml_eval_metrics <- prediction_metrics(
    pred_m = dml$m_hat,
    pred_h = dml$h_hat,
    true_m = train$m0,
    true_h = train$h0,
    C = train$C
  )
  
  out[["DML-XGBoost"]] <- data.frame(
    method = "DML-XGBoost",
    E_m = dml_eval_metrics$E_m,
    E_h = dml_eval_metrics$E_h,
    D_score =
      dml_eval_metrics$D_score,
    tau_hat = dml$tau$tau,
    se = dml$tau$se,
    bias =
      dml$tau$tau -
      tau_true,
    squared_error =
      (dml$tau$tau -
         tau_true)^2,
    coverage =
      tau_true >=
      dml$tau$ci_low &&
      tau_true <=
      dml$tau$ci_high
  )
  
  
  bind_rows(out)
}


if (file.exists(
  "simulation_checkpoint.RData"
)) {
  
  load(
    "simulation_checkpoint.RData"
  )
  
} else {
  
  results_list <- list()
  start_sim <- 1
}


cat(
  "Running",
  N_SIM,
  "replicates\n"
)


for (s in start_sim:N_SIM) {
  
  cat(
    "Replicate",
    s,
    "of",
    N_SIM,
    "\n"
  )
  
  results_list[[s]] <- run_one(
    seed = 10000 + s
  )
  
  if (
    s %% 5 == 0 ||
    s == N_SIM
  ) {
    
    save(
      results_list,
      s,
      N_SIM,
      file =
        "simulation_checkpoint.RData"
    )
  }
}


results <- bind_rows(
  results_list,
  .id = "replicate"
)


summary_results <- results %>%
  group_by(method) %>%
  summarise(
    mean_E_m = mean(E_m),
    mean_E_h = mean(E_h),
    mean_D_score =
      mean(D_score),
    bias = mean(bias),
    abs_bias =
      mean(abs(bias)),
    RMSE =
      sqrt(mean(squared_error)),
    empirical_SE =
      sd(tau_hat),
    mean_model_SE =
      mean(se),
    coverage =
      mean(coverage),
    .groups = "drop"
  )


regional_results <- bind_rows(
  lapply(
    seq_along(results_list),
    function(i) {
      
      rep_result <- run_one(
        seed = 10000 + i
      )
      
      rep_result$replicate <- i
      
      rep_result
    }
  )
)


plot_data <- results %>%
  mutate(
    abs_bias = abs(bias),
    total_prediction_error =
      E_m + E_h
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

write.csv(
  plot_data,
  "baseline_simulation_plot_data.csv",
  row.names = FALSE
)


save(
  results,
  results_list,
  summary_results,
  plot_data,
  file =
    "baseline_simulation_workspace.RData"
)


p1 <- ggplot(
  plot_data,
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
    x = "|E(delta_m delta_h)|",
    y = "Absolute causal bias"
  ) +
  theme_bw()


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


print(summary_results)

cat(
  "\nSimulation finished.\n"
)

cat(
  "Results saved in:\n",
  getwd(),
  "\n"
)