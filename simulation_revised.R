rm(list = ls())
set.seed(20260827)

library(mgcv)
library(xgboost)
library(dplyr)
library(ggplot2)


N_SIM <- 500

n_train <- 1000
n_eval <- 10000

p <- 10
tau_true <- 1

rho <- 0.3
sigma_x <- 1
sigma_y <- 1

nfold <- 5

n_cluster <- 100

use_clustered <- FALSE


make_covariates <- function(n, p, rho) {
  
  Sigma <- matrix(rho, p, p)
  diag(Sigma) <- 1
  
  Z <- matrix(
    rnorm(n * p),
    nrow = n,
    ncol = p
  )
  
  C <- Z %*% chol(Sigma)
  
  colnames(C) <- paste0("C", seq_len(p))
  
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
    tau,
    rho,
    sigma_x,
    sigma_y,
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
      n,
      replace = TRUE
    )
    
    ux <- rnorm(
      n_cluster,
      0,
      0.5
    )
    
    uy <- rnorm(
      n_cluster,
      0,
      0.5
    )
    
    X <- m0 +
      ux[cluster_id] +
      rnorm(n, 0, sigma_x)
    
    Y <- tau * X +
      g0 +
      uy[cluster_id] +
      rnorm(n, 0, sigma_y)
    
  } else {
    
    cluster_id <- seq_len(n)
    
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
    g0 = g0,
    h0 = h0,
    cluster_id = cluster_id
  )
}


estimate_tau <- function(
    X,
    Y,
    m_hat,
    h_hat) {
  
  x_res <- X - m_hat
  y_res <- Y - h_hat
  
  tau_hat <- sum(
    x_res * y_res
  ) / sum(
    x_res^2
  )
  
  psi <- x_res *
    (y_res - tau_hat * x_res)
  
  J <- mean(
    x_res^2
  )
  
  variance <- mean(
    psi^2
  ) / J^2
  
  se <- sqrt(
    variance / length(X)
  )
  
  list(
    tau = tau_hat,
    se = se,
    lower = tau_hat - 1.96 * se,
    upper = tau_hat + 1.96 * se
  )
}


estimate_tau_cluster <- function(
    X,
    Y,
    m_hat,
    h_hat,
    cluster_id) {
  
  x_res <- X - m_hat
  y_res <- Y - h_hat
  
  tau_hat <- sum(
    x_res * y_res
  ) / sum(
    x_res^2
  )
  
  psi <- x_res *
    (y_res - tau_hat * x_res)
  
  J <- mean(
    x_res^2
  )
  
  cluster_score <- tapply(
    psi,
    cluster_id,
    sum
  )
  
  G <- length(
    cluster_score
  )
  
  variance <- sum(
    cluster_score^2
  ) / (G^2 * J^2)
  
  se <- sqrt(
    variance
  )
  
  list(
    tau = tau_hat,
    se = se,
    lower = tau_hat - 1.96 * se,
    upper = tau_hat + 1.96 * se
  )
}


fit_ols <- function(
    train,
    newdata) {
  
  train_m <- data.frame(
    X = train$X,
    train$C
  )
  
  train_h <- data.frame(
    Y = train$Y,
    train$C
  )
  
  newdata <- data.frame(
    newdata
  )
  
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
      predict(
        fit_m,
        newdata = newdata
      )
    ),
    h = as.numeric(
      predict(
        fit_h,
        newdata = newdata
      )
    )
  )
}


fit_gam <- function(
    train,
    newdata) {
  
  train_m <- data.frame(
    X = train$X,
    train$C
  )
  
  train_h <- data.frame(
    Y = train$Y,
    train$C
  )
  
  newdata <- data.frame(
    newdata
  )
  
  form_m <- as.formula(
    paste(
      "X ~",
      paste0(
        "s(C", seq_len(p), ")",
        collapse = " + "
      )
    )
  )
  
  form_h <- as.formula(
    paste(
      "Y ~",
      paste0(
        "s(C", seq_len(p), ")",
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
      predict(
        fit_m,
        newdata = newdata
      )
    ),
    h = as.numeric(
      predict(
        fit_h,
        newdata = newdata
      )
    )
  )
}


fit_xgb <- function(
    train,
    newdata) {
  
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
        dtest
      )
    ),
    h = as.numeric(
      predict(
        fit_h,
        dtest
      )
    )
  )
}


get_folds <- function(
    cluster_id,
    nfold) {
  
  clusters <- unique(
    cluster_id
  )
  
  fold_cluster <- sample(
    rep(
      seq_len(nfold),
      length.out = length(clusters)
    )
  )
  
  fold_id <- fold_cluster[
    match(
      cluster_id,
      clusters
    )
  ]
  
  fold_id
}


crossfit_xgb <- function(
    dat,
    nfold = 5,
    clustered = FALSE) {
  
  if (clustered) {
    
    fold_id <- get_folds(
      dat$cluster_id,
      nfold
    )
    
  } else {
    
    fold_id <- sample(
      rep(
        seq_len(nfold),
        length.out = nrow(dat$C)
      )
    )
  }
  
  m_hat <- numeric(
    nrow(dat$C)
  )
  
  h_hat <- numeric(
    nrow(dat$C)
  )
  
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
      X = dat$X[
        train_id
      ],
      Y = dat$Y[
        train_id
      ]
    )
    
    pred <- fit_xgb(
      train = train,
      newdata = dat$C[
        test_id,
        ,
        drop = FALSE
      ]
    )
    
    m_hat[
      test_id
    ] <- pred$m
    
    h_hat[
      test_id
    ] <- pred$h
  }
  
  if (clustered) {
    
    tau <- estimate_tau_cluster(
      X = dat$X,
      Y = dat$Y,
      m_hat = m_hat,
      h_hat = h_hat,
      cluster_id = dat$cluster_id
    )
    
  } else {
    
    tau <- estimate_tau(
      X = dat$X,
      Y = dat$Y,
      m_hat = m_hat,
      h_hat = h_hat
    )
  }
  
  list(
    m_hat = m_hat,
    h_hat = h_hat,
    tau = tau,
    fold_id = fold_id
  )
}


prediction_metrics <- function(
    pred_m,
    pred_h,
    true_m,
    true_h) {
  
  delta_m <- pred_m - true_m
  delta_h <- pred_h - true_h
  
  c(
    E_m = mean(delta_m^2),
    E_h = mean(delta_h^2),
    D_score = abs(
      mean(delta_m * delta_h)
    )
  )
}


run_oracle <- function(dat) {
  
  if (use_clustered) {
    
    fit <- estimate_tau_cluster(
      X = dat$X,
      Y = dat$Y,
      m_hat = dat$m0,
      h_hat = dat$h0,
      cluster_id = dat$cluster_id
    )
    
  } else {
    
    fit <- estimate_tau(
      X = dat$X,
      Y = dat$Y,
      m_hat = dat$m0,
      h_hat = dat$h0
    )
  }
  
  data.frame(
    method = "Oracle",
    E_m = 0,
    E_h = 0,
    D_score = 0,
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


run_one <- function(seed) {
  
  set.seed(seed)
  
  train <- generate_data(
    n = n_train,
    p = p,
    tau = tau_true,
    rho = rho,
    sigma_x = sigma_x,
    sigma_y = sigma_y,
    clustered = use_clustered,
    n_cluster = n_cluster
  )
  
  eval <- generate_data(
    n = n_eval,
    p = p,
    tau = tau_true,
    rho = rho,
    sigma_x = sigma_x,
    sigma_y = sigma_y,
    clustered = use_clustered,
    n_cluster = n_cluster
  )
  
  
  oracle <- run_oracle(
    eval
  )
  
  
  method_results <- list(
    oracle
  )
  
  
  for (method in c(
    "OLS",
    "GAM",
    "XGBoost"
  )) {
    
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
      
    } else {
      
      pred <- fit_xgb(
        train,
        eval$C
      )
    }
    
    
    metrics <- prediction_metrics(
      pred_m = pred$m,
      pred_h = pred$h,
      true_m = eval$m0,
      true_h = eval$h0
    )
    
    
    if (use_clustered) {
      
      tau <- estimate_tau_cluster(
        X = eval$X,
        Y = eval$Y,
        m_hat = pred$m,
        h_hat = pred$h,
        cluster_id = eval$cluster_id
      )
      
    } else {
      
      tau <- estimate_tau(
        X = eval$X,
        Y = eval$Y,
        m_hat = pred$m,
        h_hat = pred$h
      )
    }
    
    
    method_results[[length(method_results) + 1]] <-
      data.frame(
        method = method,
        E_m = metrics["E_m"],
        E_h = metrics["E_h"],
        D_score = metrics["D_score"],
        tau_hat = tau$tau,
        se = tau$se,
        bias = tau$tau - tau_true,
        squared_error =
          (tau$tau - tau_true)^2,
        coverage =
          tau_true >= tau$lower &&
          tau_true <= tau$upper
      )
  }
  
  
  dml <- crossfit_xgb(
    dat = train,
    nfold = nfold,
    clustered = use_clustered
  )
  
  
  dml_metrics <- prediction_metrics(
    pred_m = dml$m_hat,
    pred_h = dml$h_hat,
    true_m = train$m0,
    true_h = train$h0
  )
  
  
  method_results[[length(method_results) + 1]] <-
    data.frame(
      method = "DML-XGBoost",
      E_m = dml_metrics["E_m"],
      E_h = dml_metrics["E_h"],
      D_score = dml_metrics["D_score"],
      tau_hat = dml$tau$tau,
      se = dml$tau$se,
      bias =
        dml$tau$tau - tau_true,
      squared_error =
        (dml$tau$tau - tau_true)^2,
      coverage =
        tau_true >= dml$tau$lower &&
        tau_true <= dml$tau$upper
    )
  
  
  bind_rows(
    method_results
  )
}

if (file.exists("revised_baseline_checkpoint.RData")) {
  
  load("revised_baseline_checkpoint.RData")
  
  start_sim <- length(results_list) + 1
  
} else {
  
  results_list <- vector(
    "list",
    N_SIM
  )
  
  start_sim <- 1
}


cat(
  "Running replicates",
  start_sim,
  "to",
  N_SIM,
  "\n"
)


if (start_sim <= N_SIM) {
  
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
    
    if (s %% 5 == 0 || s == N_SIM) {
      
      save(
        results_list,
        file =
          "revised_baseline_checkpoint.RData"
      )
      
      cat(
        "Checkpoint saved at replicate",
        s,
        "\n"
      )
    }
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


write.csv(
  results,
  "revised_simulation_results.csv",
  row.names = FALSE
)

write.csv(
  summary_results,
  "revised_simulation_summary.csv",
  row.names = FALSE
)

save(
  results,
  results_list,
  summary_results,
  file = "revised_simulation_workspace.RData"
)


print(
  summary_results,
  width = Inf
)


cat(
  "\nSimulation finished.\n"
)

cat(
  "Results saved in:\n",
  getwd(),
  "\n"
)