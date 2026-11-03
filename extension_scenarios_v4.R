# extension_scenarios_v4.R
# Additional scenarios added in revision (v4). Same learners, estimator and
# evaluation as simulation_revised_v3.R, with four changes to the design:
#
#   S3  heterogeneous effect     tau(C) = 1 + 0.5 C1
#   S4  limited overlap          sd(eps_X) = 0.5, so Var(X | C) = 0.25
#   S5  unmeasured confounding   U ~ N(0,1) enters X and Y with loading 0.5
#   S6  prediction-tuned XGBoost Scenario 1 DGP; max_depth and rounds chosen
#                                by 5-fold CV prediction error
#
# Each replicate also computes the reference-learner remainder (an observable
# diagnostic that uses cross-fitted DML-XGBoost as the reference learner).
#
# Usage:   Rscript extension_scenarios_v4.R S3 300
# Backend: XGB_BACKEND=R (default; R package xgboost) or XGB_BACKEND=python
#          (reticulate + Python xgboost; same C++ library).  The results in the
#          revised manuscript were produced with xgboost 3.0.5 via the Python
#          backend, R 4.3.3, mgcv 1.9-1.

args  <- commandArgs(trailingOnly = TRUE)
scen  <- if (length(args) >= 1) args[1] else "S3"
N_SIM <- if (length(args) >= 2) as.integer(args[2]) else 300

suppressMessages({
  library(mgcv)
  library(dplyr)
})

n_train <- 1000
n_eval  <- 10000
p       <- 10
rho     <- 0.3
nfold   <- 5
theta   <- 1   # target parameter: tau (S4-S6) or E{tau(C)} (S3)

cfg <- switch(
  scen,
  S3 = list(sigma_x = 1.0, het = 0.5, u_load = 0.0, tune = FALSE, seed = 30000),
  S4 = list(sigma_x = 0.5, het = 0.0, u_load = 0.0, tune = FALSE, seed = 40000),
  S5 = list(sigma_x = 1.0, het = 0.0, u_load = 0.5, tune = FALSE, seed = 50000),
  S6 = list(sigma_x = 1.0, het = 0.0, u_load = 0.0, tune = TRUE,  seed = 60000),
  stop("unknown scenario")
)

s2x <- cfg$sigma_x^2 + cfg$u_load^2   # Var(X | C), used in the remainder

out_dir <- "v4_output"
dir.create(out_dir, showWarnings = FALSE)

# ---------------------------------------------------------------- XGBoost
XGB_BACKEND <- Sys.getenv("XGB_BACKEND", "R")

if (XGB_BACKEND == "python") {
  library(reticulate)
  if (nzchar(Sys.getenv("XGB_PYTHON"))) use_python(Sys.getenv("XGB_PYTHON"), required = TRUE)
  pyxgb <- import("xgboost", convert = TRUE)
}

xgb_params <- function(max_depth, seed) {
  list(objective = "reg:squarederror", eta = 0.05, max_depth = as.integer(max_depth),
       subsample = 0.8, colsample_bytree = 0.8, nthread = 1L, seed = as.integer(seed))
}

# Fit on (Ctr, y); return predictions for Cnew
xgb_fit_predict <- function(Ctr, y, Cnew, max_depth = 4, nrounds = 300) {
  seed <- sample.int(1e6, 1)
  if (XGB_BACKEND == "python") {
    d   <- pyxgb$DMatrix(r_to_py(Ctr), label = r_to_py(y))
    fit <- pyxgb$train(do.call(dict, xgb_params(max_depth, seed)), d,
                       num_boost_round = as.integer(nrounds), verbose_eval = FALSE)
    as.numeric(fit$predict(pyxgb$DMatrix(r_to_py(Cnew))))
  } else {
    d   <- xgboost::xgb.DMatrix(Ctr, label = y)
    fit <- xgboost::xgb.train(params = xgb_params(max_depth, seed), data = d,
                              nrounds = nrounds, verbose = 0)
    as.numeric(predict(fit, xgboost::xgb.DMatrix(Cnew)))
  }
}

# 5-fold CV with early stopping: returns c(rounds, cv_mse)
# (fixed_rounds: evaluate a fixed number of rounds instead, no early stopping)
xgb_cv <- function(Ctr, y, max_depth, fixed_rounds = NULL) {
  seed <- sample.int(1e6, 1)
  es   <- is.null(fixed_rounds)
  if (XGB_BACKEND == "python") {
    d  <- pyxgb$DMatrix(r_to_py(Ctr), label = r_to_py(y))
    cv <- pyxgb$cv(do.call(dict, xgb_params(max_depth, seed)), d,
                   num_boost_round = if (es) 1000L else as.integer(fixed_rounds),
                   nfold = 5L, early_stopping_rounds = if (es) 25L else NULL,
                   seed = as.integer(seed), verbose_eval = FALSE)
    rmse <- unlist(cv[["test-rmse-mean"]])       # pandas frame arrives as a list
    c(rounds = length(rmse), cv_mse = tail(rmse, 1)^2)
  } else {
    d  <- xgboost::xgb.DMatrix(Ctr, label = y)
    cv <- xgboost::xgb.cv(params = xgb_params(max_depth, seed), data = d,
                          nrounds = if (es) 1000 else fixed_rounds, nfold = 5,
                          early_stopping_rounds = if (es) 25 else NULL, verbose = 0)
    best <- if (!es) fixed_rounds else if (!is.null(cv$best_iteration)) cv$best_iteration else
      cv$early_stop$best_iteration
    c(rounds = best, cv_mse = cv$evaluation_log$test_rmse_mean[best]^2)
  }
}

tune_xgb <- function(Ctr, y, depths = c(2, 4, 6)) {
  res  <- sapply(depths, function(dd) xgb_cv(Ctr, y, dd))
  best <- which.min(res["cv_mse", ])
  list(max_depth = depths[best], nrounds = res["rounds", best], cv_mse = res["cv_mse", best])
}

# ---------------------------------------------------------------- DGP
make_covariates <- function(n) {
  Sigma <- matrix(rho, p, p); diag(Sigma) <- 1
  C <- matrix(rnorm(n * p), n, p) %*% chol(Sigma)
  colnames(C) <- paste0("C", seq_len(p))
  C
}
m0_fun <- function(C) 0.8 * C[, 1] + 0.5 * C[, 2] + 0.4 * sin(C[, 3]) + 0.3 * C[, 4] * C[, 5]
g0_fun <- function(C) 0.7 * C[, 1] - 0.5 * C[, 2] + 0.4 * C[, 3]^2 + 0.3 * C[, 4] * C[, 5]

generate_data <- function(n) {
  C    <- make_covariates(n)
  m0   <- m0_fun(C)
  g0   <- g0_fun(C)
  tauC <- theta + cfg$het * C[, 1]
  U    <- rnorm(n)
  X    <- m0 + cfg$u_load * U + rnorm(n, 0, cfg$sigma_x)
  Y    <- tauC * X + g0 + cfg$u_load * U + rnorm(n, 0, 1)
  list(C = C, X = X, Y = Y, m0 = m0, h0 = tauC * m0 + g0)   # h0 = E(Y | C)
}

# ---------------------------------------------------------------- estimator
estimate_tau <- function(X, Y, m_hat, h_hat) {
  xr <- X - m_hat; yr <- Y - h_hat
  tau <- sum(xr * yr) / sum(xr^2)
  psi <- xr * (yr - tau * xr)
  se  <- sqrt(mean(psi^2) / mean(xr^2)^2 / length(X))
  list(tau = tau, se = se, lower = tau - 1.96 * se, upper = tau + 1.96 * se,
       den = mean(xr^2))
}

prediction_metrics <- function(pred_m, pred_h, true_m, true_h) {
  dm <- pred_m - true_m
  dh <- pred_h - true_h
  dg <- dh - theta * dm                    # error in the implied g = h - theta m
  c(E_m = mean(dm^2), E_h = mean(dh^2), E_g = mean(dg^2),
    cross = mean(dm * dh), D_joint = abs(mean(dm * dh)), R = mean(dm * dg))
}

# ---------------------------------------------------------------- learners
fit_ols <- function(tr, newC) {
  dfm <- data.frame(X = tr$X, tr$C); dfh <- data.frame(Y = tr$Y, tr$C); nd <- data.frame(newC)
  list(m = as.numeric(predict(lm(X ~ ., dfm), nd)), h = as.numeric(predict(lm(Y ~ ., dfh), nd)))
}
fit_gam <- function(tr, newC) {
  rhs <- paste0("s(C", seq_len(p), ")", collapse = " + ")
  dfm <- data.frame(X = tr$X, tr$C); dfh <- data.frame(Y = tr$Y, tr$C); nd <- data.frame(newC)
  fm <- gam(as.formula(paste("X ~", rhs)), data = dfm, method = "REML")
  fh <- gam(as.formula(paste("Y ~", rhs)), data = dfh, method = "REML")
  list(m = as.numeric(predict(fm, nd)), h = as.numeric(predict(fh, nd)))
}
fit_xgb <- function(tr, newC, hp_m = list(max_depth = 4, nrounds = 300),
                    hp_h = hp_m) {
  list(m = xgb_fit_predict(tr$C, tr$X, newC, hp_m$max_depth, hp_m$nrounds),
       h = xgb_fit_predict(tr$C, tr$Y, newC, hp_h$max_depth, hp_h$nrounds))
}

crossfit_xgb <- function(dat, eval, hp_m, hp_h) {
  fold <- sample(rep(seq_len(nfold), length.out = nrow(dat$C)))
  m_hat <- h_hat <- numeric(nrow(dat$C))
  met <- vector("list", nfold)
  for (k in seq_len(nfold)) {
    te <- which(fold == k)
    tr <- list(C = dat$C[-te, , drop = FALSE], X = dat$X[-te], Y = dat$Y[-te])
    pr <- fit_xgb(tr, rbind(dat$C[te, , drop = FALSE], eval$C), hp_m, hp_h)
    nt <- length(te)
    m_hat[te] <- pr$m[seq_len(nt)]
    h_hat[te] <- pr$h[seq_len(nt)]
    met[[k]] <- prediction_metrics(pr$m[-seq_len(nt)], pr$h[-seq_len(nt)], eval$m0, eval$h0)
  }
  list(m_hat = m_hat, h_hat = h_hat,
       tau = estimate_tau(dat$X, dat$Y, m_hat, h_hat),
       metrics = colMeans(do.call(rbind, met)))
}

# Reference-learner remainder: observable proxy of R for learner A, using the
# cross-fitted reference learner B in place of the unknown true functions.
ref_remainder <- function(mA, hA, mB, hB, tauB, den) {
  gA <- hA - tauB * mA
  gB <- hB - tauB * mB
  mean((mA - mB) * (gA - gB)) / den
}

make_row <- function(method, met, fit, extra = list()) {
  data.frame(method = method, as.list(met), tau_hat = fit$tau, se = fit$se,
             bias = fit$tau - theta, squared_error = (fit$tau - theta)^2,
             coverage = theta >= fit$lower && theta <= fit$upper,
             pred_bias = met[["R"]] / (s2x + met[["E_m"]]),
             cs_bound  = sqrt(met[["E_m"]] * met[["E_g"]]) / (s2x + met[["E_m"]]),
             ref_pred_bias = if (is.null(extra$ref)) NA_real_ else extra$ref,
             gap_vs_ref    = if (is.null(extra$gap)) NA_real_ else extra$gap,
             cv_mse_m = if (is.null(extra$cv_m)) NA_real_ else extra$cv_m,
             cv_mse_h = if (is.null(extra$cv_h)) NA_real_ else extra$cv_h,
             depth_m  = if (is.null(extra$d_m)) NA_real_ else extra$d_m,
             depth_h  = if (is.null(extra$d_h)) NA_real_ else extra$d_h,
             check.names = FALSE)
}

zero_met <- c(E_m = 0, E_h = 0, E_g = 0, cross = 0, D_joint = 0, R = 0)

# ---------------------------------------------------------------- one replicate
run_one <- function(seed) {
  set.seed(seed)
  train <- generate_data(n_train)
  eval  <- generate_data(n_eval)
  newC  <- rbind(train$C, eval$C)
  tr    <- seq_len(n_train)

  rows <- list(make_row("Oracle", zero_met, estimate_tau(train$X, train$Y, train$m0, train$h0)))

  default_hp <- list(max_depth = 4, nrounds = 300)

  # Reference learner: cross-fitted XGBoost with default hyperparameters
  dml <- crossfit_xgb(train, eval, default_hp, default_hp)

  learners <- if (cfg$tune) c("XGBoost") else c("OLS", "GAM", "XGBoost")
  for (L in learners) {
    pr <- switch(L, OLS = fit_ols(train, newC), GAM = fit_gam(train, newC),
                 XGBoost = fit_xgb(train, newC, default_hp, default_hp))
    met <- prediction_metrics(pr$m[-tr], pr$h[-tr], eval$m0, eval$h0)
    fit <- estimate_tau(train$X, train$Y, pr$m[tr], pr$h[tr])
    ref <- ref_remainder(pr$m[tr], pr$h[tr], dml$m_hat, dml$h_hat, dml$tau$tau, fit$den)
    rows[[length(rows) + 1]] <- make_row(L, met, fit,
                                         list(ref = ref, gap = fit$tau - dml$tau$tau))
  }
  rows[[length(rows) + 1]] <- make_row("DML-XGBoost", dml$metrics, dml$tau)

  if (cfg$tune) {
    # hyperparameters chosen by 5-fold CV prediction error on the training sample
    tm <- tune_xgb(train$C, train$X)
    th <- tune_xgb(train$C, train$Y)
    cvd_m <- xgb_cv(train$C, train$X, 4, fixed_rounds = 300)   # CV error of the
    cvd_h <- xgb_cv(train$C, train$Y, 4, fixed_rounds = 300)   # default settings
    hp_m <- list(max_depth = tm$max_depth, nrounds = tm$nrounds)
    hp_h <- list(max_depth = th$max_depth, nrounds = th$nrounds)

    pr  <- fit_xgb(train, newC, hp_m, hp_h)
    met <- prediction_metrics(pr$m[-tr], pr$h[-tr], eval$m0, eval$h0)
    fit <- estimate_tau(train$X, train$Y, pr$m[tr], pr$h[tr])
    rows[[length(rows) + 1]] <- make_row("XGBoost-tuned", met, fit,
                                         list(cv_m = tm$cv_mse, cv_h = th$cv_mse,
                                              d_m = tm$max_depth, d_h = th$max_depth))
    dmlt <- crossfit_xgb(train, eval, hp_m, hp_h)
    rows[[length(rows) + 1]] <- make_row("DML-XGBoost-tuned", dmlt$metrics, dmlt$tau,
                                         list(cv_m = tm$cv_mse, cv_h = th$cv_mse,
                                              d_m = tm$max_depth, d_h = th$max_depth))
    rows[[2]]$cv_mse_m <- cvd_m[["cv_mse"]]
    rows[[2]]$cv_mse_h <- cvd_h[["cv_mse"]]
  }
  bind_rows(rows)
}

# ---------------------------------------------------------------- run
ckpt <- file.path(out_dir, paste0(scen, "_checkpoint.rds"))
results_list <- if (file.exists(ckpt)) readRDS(ckpt) else list()
start <- length(results_list) + 1
if (start <= N_SIM) for (s in start:N_SIM) {
  results_list[[s]] <- run_one(cfg$seed + s)
  if (s %% 10 == 0 || s == N_SIM) {
    saveRDS(results_list, ckpt)
    cat(scen, "replicate", s, "of", N_SIM, format(Sys.time()), "\n")
  }
}

results <- bind_rows(results_list, .id = "replicate")
write.csv(results, file.path(out_dir, paste0(scen, "_results.csv")), row.names = FALSE)

summ <- results %>%
  group_by(method) %>%
  summarise(r_E    = if (all(E_m == 0)) NA else cor(E_m + E_h, abs(bias)),
            r_pred = if (all(E_m == 0)) NA else cor(pred_bias, bias),
            r_D    = if (all(E_m == 0)) NA else cor(D_joint, abs(bias)),
            r_ref  = if (all(is.na(ref_pred_bias))) NA else cor(ref_pred_bias, bias),
            r_gap  = if (all(is.na(gap_vs_ref))) NA else cor(gap_vs_ref, bias),
            R_rep = n(), MCSE_bias = sd(tau_hat) / sqrt(n()),
            RMSE = sqrt(mean(squared_error)), EmpSE = sd(tau_hat), ModelSE = mean(se),
            SEratio = mean(se) / sd(tau_hat), Coverage = mean(coverage),
            E_m = mean(E_m), E_h = mean(E_h), E_g = mean(E_g), D_joint = mean(D_joint),
            pred_bias = mean(pred_bias), cs_bound = mean(cs_bound),
            ref_pred_bias = mean(ref_pred_bias), gap_vs_ref = mean(gap_vs_ref),
            cv_mse_m = mean(cv_mse_m), cv_mse_h = mean(cv_mse_h),
            depth_m = paste(names(table(depth_m)), table(depth_m), sep = ":", collapse = " "),
            depth_h = paste(names(table(depth_h)), table(depth_h), sep = ":", collapse = " "),
            Bias = mean(bias),
            .groups = "drop")
write.csv(summ, file.path(out_dir, paste0(scen, "_summary.csv")), row.names = FALSE)
print(as.data.frame(summ), digits = 3)
