 
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


 
# 3. TRUE NUISANCE FUNCTIONS
 

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


 
# 4. TRUE CLUSTERED DGP
 

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
  
  n_cluster <- ceiling(
    n / cluster_size
  )
  
  cluster_id <- rep(
    seq_len(n_cluster),
    each = cluster_size
  )
  
  cluster_id <- cluster_id[
    seq_len(n)
  ]
  
  C <- make_covariates(
    n = n,
    p = p,
    rho = rho
  )
  
  m0 <- m0_fun(C)
  
  h0 <- h0_fun(C)
  
  
  # Cluster-level random effects
  
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
  
  
  # Individual-level errors
  
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
  
  
  # Treatment / exposure
  
  X <- m0 +
    u_x[cluster_id] +
    epsilon_x
  
  
  # Outcome
  
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


 
# 5. PREDICTION METRICS
 

prediction_metrics <- function(
    pred_m,
    pred_h,
    true_m,
    true_h
) {
  
  delta_m <- pred_m -
    true_m
  
  delta_h <- pred_h -
    true_h
  
  list(
    
    E_m = mean(
      delta_m^2
    ),
    
    E_h = mean(
      delta_h^2
    ),
    
    D_score = abs(
      mean(
        delta_m * delta_h
      )
    ),
    
    delta_m = delta_m,
    delta_h = delta_h
  )
}


 
# 6. OLS
 

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


 
# 7. GAM
 

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


 
# 8. XGBOOST
 

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


 
# 9. CLUSTER-ROBUST TAU ESTIMATION
 

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


 
# 10. CLUSTERED DML-XGBOOST
 

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
      length.out =
        length(cluster_ids)
    )
  )
  
  
  names(cluster_folds) <-
    cluster_ids
  
  
  folds <- cluster_folds[
    as.character(
      dat$cluster_id
    )
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
    
    cluster_id =
      dat$cluster_id
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


 
# 11. FINAL run_one()
 

run_one <- function(
    seed
) {
  
  set.seed(
    seed
  )
  
  
  ##########################################################
  # TRAIN DATA
  ##########################################################
  
  train <- generate_clustered_data(
    
    n = n_train,
    
    p = p,
    
    cluster_size =
      cluster_size_train,
    
    tau = tau_true,
    
    rho = rho,
    
    sigma_x = sigma_x,
    
    sigma_y = sigma_y,
    
    sigma_cluster_x =
      sigma_cluster_x,
    
    sigma_cluster_y =
      sigma_cluster_y
  )
  
   
  # EVALUATION DATA
 
  
  eval <- generate_clustered_data(
    
    n = n_eval,
    
    p = p,
    
    cluster_size =
      cluster_size_eval,
    
    tau = tau_true,
    
    rho = rho,
    
    sigma_x = sigma_x,
    
    sigma_y = sigma_y,
    
    sigma_cluster_x =
      sigma_cluster_x,
    
    sigma_cluster_y =
      sigma_cluster_y
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
    
    pred_m =
      oracle_pred$m,
    
    pred_h =
      oracle_pred$h,
    
    cluster_id =
      eval$cluster_id
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
    
    pred_m =
      ols_pred$m,
    
    pred_h =
      ols_pred$h,
    
    cluster_id =
      eval$cluster_id
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
    
    pred_m =
      gam_pred$m,
    
    pred_h =
      gam_pred$h,
    
    cluster_id =
      eval$cluster_id
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
    
    pred_m =
      xgb_pred$m,
    
    pred_h =
      xgb_pred$h,
    
    cluster_id =
      eval$cluster_id
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
      
      E_m =
        oracle_metrics$E_m,
      
      E_h =
        oracle_metrics$E_h,
      
      D_score =
        oracle_metrics$D_score,
      
      tau_hat =
        oracle_tau$tau,
      
      se =
        oracle_tau$se,
      
      bias =
        oracle_tau$tau -
        tau_true,
      
      squared_error =
        (
          oracle_tau$tau -
            tau_true
        )^2,
      
      coverage =
        tau_true >=
        oracle_tau$ci_low &&
        tau_true <=
        oracle_tau$ci_high
    ),
    
    
    data.frame(
      
      method = "OLS",
      
      E_m =
        ols_metrics$E_m,
      
      E_h =
        ols_metrics$E_h,
      
      D_score =
        ols_metrics$D_score,
      
      tau_hat =
        ols_tau$tau,
      
      se =
        ols_tau$se,
      
      bias =
        ols_tau$tau -
        tau_true,
      
      squared_error =
        (
          ols_tau$tau -
            tau_true
        )^2,
      
      coverage =
        tau_true >=
        ols_tau$ci_low &&
        tau_true <=
        ols_tau$ci_high
    ),
    
    
    data.frame(
      
      method = "GAM",
      
      E_m =
        gam_metrics$E_m,
      
      E_h =
        gam_metrics$E_h,
      
      D_score =
        gam_metrics$D_score,
      
      tau_hat =
        gam_tau$tau,
      
      se =
        gam_tau$se,
      
      bias =
        gam_tau$tau -
        tau_true,
      
      squared_error =
        (
          gam_tau$tau -
            tau_true
        )^2,
      
      coverage =
        tau_true >=
        gam_tau$ci_low &&
        tau_true <=
        gam_tau$ci_high
    ),
    
    
    data.frame(
      
      method = "XGBoost",
      
      E_m =
        xgb_metrics$E_m,
      
      E_h =
        xgb_metrics$E_h,
      
      D_score =
        xgb_metrics$D_score,
      
      tau_hat =
        xgb_tau$tau,
      
      se =
        xgb_tau$se,
      
      bias =
        xgb_tau$tau -
        tau_true,
      
      squared_error =
        (
          xgb_tau$tau -
            tau_true
        )^2,
      
      coverage =
        tau_true >=
        xgb_tau$ci_low &&
        tau_true <=
        xgb_tau$ci_high
    ),
    
    
    data.frame(
      
      method = "DML-XGBoost",
      
      E_m =
        dml_metrics$E_m,
      
      E_h =
        dml_metrics$E_h,
      
      D_score =
        dml_metrics$D_score,
      
      tau_hat =
        dml$tau,
      
      se =
        dml$se,
      
      bias =
        dml$tau -
        tau_true,
      
      squared_error =
        (
          dml$tau -
            tau_true
        )^2,
      
      coverage =
        tau_true >=
        dml$ci_low &&
        tau_true <=
        dml$ci_high
    )
  )
}


 
# 12. FUNCTION CHECK
 

 


required_functions <- c(
  "make_covariates",
  "m0_fun",
  "h0_fun",
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
    
    mean_E_m =
      mean(E_m),
    
    mean_E_h =
      mean(E_h),
    
    mean_D_score =
      mean(D_score),
    
    bias =
      mean(bias),
    
    abs_bias =
      mean(abs(bias)),
    
    RMSE =
      sqrt(
        mean(squared_error)
      ),
    
    empirical_SE =
      sd(tau_hat),
    
    mean_model_SE =
      mean(se),
    
    coverage =
      mean(coverage),
    
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


 
# 22. PLOT 3: BIAS BY METHOD
 

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


 
# 26. FINAL PERFORMANCE TABLE
 
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


 # 27. PREDICTION TABLE
 

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


 # 28. COVERAGE CHECK
 
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


 # 29. ORACLE SANITY CHECK
 

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


 # 30. FINAL FILE CHECK
 
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
 