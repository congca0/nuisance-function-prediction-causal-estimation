rm(list = ls())
set.seed(20260827)

library(mgcv)
library(xgboost)
library(dplyr)
library(ggplot2)

N_SIM <- 5

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


h0_fun <- function(C) {
  
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
  
  cluster_id <- rep(
    seq_len(n_cluster),
    each = cluster_size
  )
  
  cluster_id <- cluster_id[seq_len(n)]
  
  C <- make_covariates(
    n = n,
    p = p,
    rho = rho
  )
  
  m0 <- m0_fun(C)
  h0 <- h0_fun(C)
  
  u_x <- rnorm(
    n_cluster,
    mean = 0,
    sd = sigma_cluster_x
  )
  
  u_y <- rnorm(
    n_cluster,
    mean = 0,
    sd = sigma_cluster_y
  )
  
  epsilon_x <- rnorm(
    n,
    mean = 0,
    sd = sigma_x
  )
  
  epsilon_y <- rnorm(
    n,
    mean = 0,
    sd = sigma_y
  )
  
  X <- m0 +
    u_x[cluster_id] +
    epsilon_x
  
  Y <- tau * X +
    h0 +
    u_y[cluster_id] +
    epsilon_y
  
  list(
    C = C,
    X = X,
    Y = Y,
    m0 = m0,
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
    D_score = abs(mean(delta_m * delta_h)),
    delta_m = delta_m,
    delta_h = delta_h
  )
}


fit_ols <- function(train, eval) {
  
  train_m <- data.frame(
    X = train$X,
    train$C
  )
  
  train_h <- data.frame(
    Y = train$Y,
    train$C
  )
  
  eval_df <- data.frame(eval$C)
  
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


fit_gam <- function(train, eval) {
  
  train_m <- data.frame(
    X = train$X,
    train$C
  )
  
  train_h <- data.frame(
    Y = train$Y,
    train$C
  )
  
  eval_df <- data.frame(eval$C)
  
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


fit_xgb <- function(train, eval) {
  
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
  
  X_res <- X - pred_m
  Y_res <- Y - pred_h
  
  tau_hat <- sum(
    X_res * Y_res
  ) / sum(
    X_res^2
  )
  
  psi <- X_res *
    (Y_res - tau_hat * X_res)
  
  denom <- mean(
    X_res^2
  )
  
  cluster_sums <- tapply(
    psi,
    cluster_id,
    sum
  )
  
  n <- length(X)
  G <- length(cluster_sums)
  
  variance <- sum(
    cluster_sums^2
  ) / (n^2 * denom^2)
  
  finite_sample_correction <- G / (G - 1)
  
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
    nfold = 5
) {
  
  cluster_ids <- unique(
    dat$cluster_id
  )
  
  cluster_folds <- sample(
    rep(
      seq_len(nfold),
      length.out = length(cluster_ids)
    )
  )
  
  names(cluster_folds) <- cluster_ids
  
  folds <- cluster_folds[
    as.character(dat$cluster_id)
  ]
  
  m_hat <- numeric(
    length(dat$X)
  )
  
  h_hat <- numeric(
    length(dat$X)
  )
  
  for (k in seq_len(nfold)) {
    
    train_idx <- which(
      folds != k
    )
    
    test_idx <- which(
      folds == k
    )
    
    train <- list(
      C = dat$C[
        train_idx,
        ,
        drop = FALSE
      ],
      X = dat$X[
        train_idx
      ],
      Y = dat$Y[
        train_idx
      ]
    )
    
    test <- list(
      C = dat$C[
        test_idx,
        ,
        drop = FALSE
      ],
      X = dat$X[
        test_idx
      ],
      Y = dat$Y[
        test_idx
      ]
    )
    
    pred <- fit_xgb(
      train,
      test
    )
    
    m_hat[
      test_idx
    ] <- pred$m
    
    h_hat[
      test_idx
    ] <- pred$h
  }
  
  estimate <- estimate_tau_cluster(
    X = dat$X,
    Y = dat$Y,
    pred_m = m_hat,
    pred_h = h_hat,
    cluster_id = dat$cluster_id
  )
  
  list(
    tau = estimate$tau,
    se = estimate$se,
    ci_low = estimate$ci_low,
    ci_high = estimate$ci_high,
    m_hat = m_hat,
    h_hat = h_hat,
    folds = folds
  )
}


run_one <- function(seed) {
  
  set.seed(seed)
  
  train <- generate_clustered_data(
    n = n_train,
    p = p,
    cluster_size = cluster_size_train,
    tau = tau_true,
    rho = rho,
    sigma_x = sigma_x,
    sigma_y = sigma_y,
    sigma_cluster_x = sigma_cluster_x,
    sigma_cluster_y = sigma_cluster_y
  )
  
  eval <- generate_clustered_data(
    n = n_eval,
    p = p,
    cluster_size = cluster_size_eval,
    tau = tau_true,
    rho = rho,
    sigma_x = sigma_x,
    sigma_y = sigma_y,
    sigma_cluster_x = sigma_cluster_x,
    sigma_cluster_y = sigma_cluster_y
  )
  
  
  oracle_pred <- list(
    m = eval$m0,
    h = eval$h0
  )
  
  oracle_metrics <- prediction_metrics(
    oracle_pred$m,
    oracle_pred$h,
    eval$m0,
    eval$h0
  )
  
  oracle_tau <- estimate_tau_cluster(
    X = eval$X,
    Y = eval$Y,
    pred_m = oracle_pred$m,
    pred_h = oracle_pred$h,
    cluster_id = eval$cluster_id
  )
  
  
  ols_pred <- fit_ols(
    train,
    eval
  )
  
  ols_metrics <- prediction_metrics(
    ols_pred$m,
    ols_pred$h,
    eval$m0,
    eval$h0
  )
  
  ols_tau <- estimate_tau_cluster(
    X = eval$X,
    Y = eval$Y,
    pred_m = ols_pred$m,
    pred_h = ols_pred$h,
    cluster_id = eval$cluster_id
  )
  
  
  gam_pred <- fit_gam(
    train,
    eval
  )
  
  gam_metrics <- prediction_metrics(
    gam_pred$m,
    gam_pred$h,
    eval$m0,
    eval$h0
  )
  
  gam_tau <- estimate_tau_cluster(
    X = eval$X,
    Y = eval$Y,
    pred_m = gam_pred$m,
    pred_h = gam_pred$h,
    cluster_id = eval$cluster_id
  )
  
  
  xgb_pred <- fit_xgb(
    train,
    eval
  )
  
  xgb_metrics <- prediction_metrics(
    xgb_pred$m,
    xgb_pred$h,
    eval$m0,
    eval$h0
  )
  
  xgb_tau <- estimate_tau_cluster(
    X = eval$X,
    Y = eval$Y,
    pred_m = xgb_pred$m,
    pred_h = xgb_pred$h,
    cluster_id = eval$cluster_id
  )
  
  
  dml <- dml_xgb_cluster(
    train,
    nfold = nfold
  )
  
  dml_metrics <- prediction_metrics(
    dml$m_hat,
    dml$h_hat,
    train$m0,
    train$h0
  )
  
  
  bind_rows(
    
    data.frame(
      method = "Oracle",
      E_m = oracle_metrics$E_m,
      E_h = oracle_metrics$E_h,
      D_score = oracle_metrics$D_score,
      tau_hat = oracle_tau$tau,
      se = oracle_tau$se,
      bias = oracle_tau$tau - tau_true,
      squared_error =
        (oracle_tau$tau - tau_true)^2,
      coverage =
        tau_true >= oracle_tau$ci_low &&
        tau_true <= oracle_tau$ci_high
    ),
    
    data.frame(
      method = "OLS",
      E_m = ols_metrics$E_m,
      E_h = ols_metrics$E_h,
      D_score = ols_metrics$D_score,
      tau_hat = ols_tau$tau,
      se = ols_tau$se,
      bias = ols_tau$tau - tau_true,
      squared_error =
        (ols_tau$tau - tau_true)^2,
      coverage =
        tau_true >= ols_tau$ci_low &&
        tau_true <= ols_tau$ci_high
    ),
    
    data.frame(
      method = "GAM",
      E_m = gam_metrics$E_m,
      E_h = gam_metrics$E_h,
      D_score = gam_metrics$D_score,
      tau_hat = gam_tau$tau,
      se = gam_tau$se,
      bias = gam_tau$tau - tau_true,
      squared_error =
        (gam_tau$tau - tau_true)^2,
      coverage =
        tau_true >= gam_tau$ci_low &&
        tau_true <= gam_tau$ci_high
    ),
    
    data.frame(
      method = "XGBoost",
      E_m = xgb_metrics$E_m,
      E_h = xgb_metrics$E_h,
      D_score = xgb_metrics$D_score,
      tau_hat = xgb_tau$tau,
      se = xgb_tau$se,
      bias = xgb_tau$tau - tau_true,
      squared_error =
        (xgb_tau$tau - tau_true)^2,
      coverage =
        tau_true >= xgb_tau$ci_low &&
        tau_true <= xgb_tau$ci_high
    ),
    
    data.frame(
      method = "DML-XGBoost",
      E_m = dml_metrics$E_m,
      E_h = dml_metrics$E_h,
      D_score = dml_metrics$D_score,
      tau_hat = dml$tau,
      se = dml$se,
      bias = dml$tau - tau_true,
      squared_error =
        (dml$tau - tau_true)^2,
      coverage =
        tau_true >= dml$ci_low &&
        tau_true <= dml$ci_high
    )
  )
}


results_list <- vector(
  "list",
  N_SIM
)


cat(
  "Running clustered simulation with",
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
  
  if (
    s %% 10 == 0 ||
    s == N_SIM
  ) {
    
    save(
      results_list,
      file = "clustered_simulation_checkpoint.RData"
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


summary_results <- results %>%
  group_by(method) %>%
  summarise(
    n = n(),
    mean_E_m = mean(E_m),
    mean_E_h = mean(E_h),
    mean_D_score = mean(D_score),
    bias = mean(bias),
    abs_bias = mean(abs(bias)),
    RMSE = sqrt(
      mean(squared_error)
    ),
    empirical_SE = sd(tau_hat),
    mean_model_SE = mean(se),
    coverage = mean(coverage),
    .groups = "drop"
  )


print(
  as.data.frame(summary_results)
)


write.csv(
  results,
  "clustered_simulation_results.csv",
  row.names = FALSE
)


write.csv(
  summary_results,
  "clustered_simulation_summary.csv",
  row.names = FALSE
)


plot_data <- results %>%
  mutate(
    abs_bias = abs(bias),
    total_prediction_error =
      E_m + E_h
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


ggsave(
  "clustered_error_vs_bias.png",
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
  "clustered_Dscore_vs_bias.png",
  p2,
  width = 10,
  height = 7,
  dpi = 300
)


bias_summary <- plot_data %>%
  group_by(method) %>%
  summarise(
    mean_abs_bias = mean(abs_bias),
    .groups = "drop"
  )


p3 <- ggplot(
  bias_summary,
  aes(
    x = reorder(method, mean_abs_bias),
    y = mean_abs_bias
  )
) +
  geom_col() +
  labs(
    x = NULL,
    y = "Mean absolute causal bias"
  ) +
  theme_bw()


ggsave(
  "clustered_bias_by_method.png",
  p3,
  width = 9,
  height = 6,
  dpi = 300
)


correlations <- results %>%
  group_by(method) %>%
  summarise(
    cor_Em_absbias =
      ifelse(
        sd(E_m) == 0 ||
          sd(abs(bias)) == 0,
        NA_real_,
        cor(E_m, abs(bias))
      ),
    cor_Eh_absbias =
      ifelse(
        sd(E_h) == 0 ||
          sd(abs(bias)) == 0,
        NA_real_,
        cor(E_h, abs(bias))
      ),
    cor_D_absbias =
      ifelse(
        sd(D_score) == 0 ||
          sd(abs(bias)) == 0,
        NA_real_,
        cor(D_score, abs(bias))
      ),
    .groups = "drop"
  )


print(
  as.data.frame(correlations)
)


write.csv(
  correlations,
  "clustered_correlations.csv",
  row.names = FALSE
)


cat(
  "\nSimulation finished.\n"
)


cat(
  "Results saved in:\n",
  getwd(),
  "\n"
)


print(
  list.files(
    pattern = "clustered",
    ignore.case = TRUE
  )
)