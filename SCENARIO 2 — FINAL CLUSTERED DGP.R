  
# SCENARIO 2 — FINAL CLUSTERED DGP
# ORACLE SANITY CHECK
 

rm(list = ls())

set.seed(20260828)

library(mgcv)
library(xgboost)
library(dplyr)
library(ggplot2)

 
# 1. FINAL PARAMETERS
 

N_SIM <- 1

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


 
# 5. PREDICTION METRICS
 

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
    D_score = abs(
      mean(delta_m * delta_h)
    )
  )
}


 
# 6. CLUSTER-ROBUST TAU ESTIMATOR
 

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
  
  variance <- variance *
    G / (G - 1)
  
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


 
# 7. GENERATE ONE CLUSTERED DATASET
 

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


 
# 8. ORACLE PREDICTIONS
 

oracle_m <- eval$m0
oracle_h <- eval$h0


 
# 9. ORACLE PREDICTION SANITY
 

oracle_E_m <- mean(
  (oracle_m - eval$m0)^2
)

oracle_E_h <- mean(
  (oracle_h - eval$h0)^2
)

oracle_D <- abs(
  mean(
    (oracle_m - eval$m0) *
      (oracle_h - eval$h0)
  )
)


 
# 10. ORACLE CAUSAL ESTIMATION
 

oracle_est <- estimate_tau_cluster(
  X = eval$X,
  Y = eval$Y,
  pred_m = oracle_m,
  pred_h = oracle_h,
  cluster_id = eval$cluster_id
)


oracle_bias <-
  oracle_est$tau -
  tau_true


 
# 11. STRUCTURAL CHECKS
 

cat("\n")
cat(" \n")
cat("SCENARIO 2 — ORACLE SANITY CHECK\n")
cat(" \n\n")

cat(
  "Number of observations:",
  n_eval,
  "\n"
)

cat(
  "Cluster size:",
  cluster_size_eval,
  "\n"
)

cat(
  "Number of clusters:",
  length(unique(eval$cluster_id)),
  "\n\n"
)

cat(
  "Oracle E_m:",
  format(oracle_E_m, digits = 8),
  "\n"
)

cat(
  "Oracle E_h:",
  format(oracle_E_h, digits = 8),
  "\n"
)

cat(
  "Oracle D_score:",
  format(oracle_D, digits = 8),
  "\n\n"
)

cat(
  "Oracle tau_hat:",
  format(oracle_est$tau, digits = 8),
  "\n"
)

cat(
  "True tau:",
  tau_true,
  "\n"
)

cat(
  "Oracle bias:",
  format(oracle_bias, digits = 8),
  "\n"
)

cat(
  "Oracle cluster-robust SE:",
  format(oracle_est$se, digits = 8),
  "\n"
)

cat(
  "Oracle CI:",
  "[",
  format(oracle_est$ci_low, digits = 8),
  ", ",
  format(oracle_est$ci_high, digits = 8),
  "]\n\n"
)


 
# 12. PASS / FAIL CRITERIA
 

check_E_m <- is.finite(oracle_E_m) &&
  oracle_E_m < 1e-12

check_E_h <- is.finite(oracle_E_h) &&
  oracle_E_h < 1e-12

check_D <- is.finite(oracle_D) &&
  oracle_D < 1e-12

check_tau <- is.finite(oracle_est$tau) &&
  abs(oracle_bias) < 0.15

check_se <- is.finite(oracle_est$se) &&
  oracle_est$se > 0

check_ci <- is.finite(
  oracle_est$ci_low
) &&
  is.finite(
    oracle_est$ci_high
  ) &&
  oracle_est$ci_low < oracle_est$ci_high

oracle_pass <-
  check_E_m &&
  check_E_h &&
  check_D &&
  check_tau &&
  check_se &&
  check_ci


oracle_check <- data.frame(
  
  check = c(
    "Oracle E_m == 0",
    "Oracle E_h == 0",
    "Oracle D_score == 0",
    "|Oracle bias| < 0.15",
    "Oracle SE > 0",
    "Oracle CI valid"
  ),
  
  result = c(
    check_E_m,
    check_E_h,
    check_D,
    check_tau,
    check_se,
    check_ci
  )
)


print(
  oracle_check,
  row.names = FALSE
)


 
# 13. HARD STOP IF STRUCTURE IS BROKEN
 

if (!oracle_pass) {
  
  stop(
    paste0(
      "\nOracle sanity check FAILED.\n",
      "DO NOT run the 500-replicate simulation.\n",
      "Inspect the clustered DGP or estimator first.\n"
    )
  )
  
} else {
  
  cat(
    " "
  )
  
  cat(
    "ORACLE SANITY CHECK PASSED.\n"
  )
  
  cat(
    "The clustered DGP is structurally valid.\n"
  )
  
  cat(
    "You may now run the FINAL 500-replicate Scenario 2.\n"
  )
  
 
}


 
# 14. SAVE SANITY CHECK
 

write.csv(
  oracle_check,
  "scenario2_clustered_oracle_sanity_check.csv",
  row.names = FALSE
)

save(
  oracle_check,
  oracle_est,
  oracle_E_m,
  oracle_E_h,
  oracle_D,
  oracle_bias,
  file =
    "scenario2_clustered_oracle_sanity_check.RData"
)