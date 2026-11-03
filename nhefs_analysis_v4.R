# nhefs_analysis_v4.R
# Real-data illustration added in revision (v4).
#
# Data: NHEFS (Hernan and Robins, Causal Inference: What If, 2020), the same
# subset used for the continuous-treatment example in their Chapter 12:
# participants with weight measured in 1982 and smokeintensity <= 25 (n = 1,162).
#   X = smkintensity82_71  change in cigarettes per day, 1971-1982
#   Y = wt82_71            weight change in kg, 1971-1982
#   C = sex, race, age, education, smokeintensity, smokeyrs, exercise, active, wt71
# Download nhefs.csv from the book website (or use causaldata::nhefs_complete)
# and place it in the working directory.
#
# Backend for XGBoost as in extension_scenarios_v4.R (XGB_BACKEND=R or python).

suppressMessages({ library(mgcv); library(dplyr) })
set.seed(20261102)

XGB_BACKEND <- Sys.getenv("XGB_BACKEND", "R")
if (XGB_BACKEND == "python") {
  library(reticulate)
  if (nzchar(Sys.getenv("XGB_PYTHON"))) use_python(Sys.getenv("XGB_PYTHON"), required = TRUE)
  pyxgb <- import("xgboost", convert = TRUE)
}
xgb_fit_predict <- function(Ctr, y, Cnew, max_depth = 4, nrounds = 300) {
  seed <- sample.int(1e6, 1)
  par <- list(objective = "reg:squarederror", eta = 0.05, max_depth = as.integer(max_depth),
              subsample = 0.8, colsample_bytree = 0.8, nthread = 1L, seed = as.integer(seed))
  if (XGB_BACKEND == "python") {
    fit <- pyxgb$train(do.call(dict, par), pyxgb$DMatrix(r_to_py(Ctr), label = r_to_py(y)),
                       num_boost_round = as.integer(nrounds), verbose_eval = FALSE)
    as.numeric(fit$predict(pyxgb$DMatrix(r_to_py(Cnew))))
  } else {
    fit <- xgboost::xgb.train(params = par, data = xgboost::xgb.DMatrix(Ctr, label = y),
                              nrounds = nrounds, verbose = 0)
    as.numeric(predict(fit, xgboost::xgb.DMatrix(Cnew)))
  }
}

# 5-fold CV with early stopping (rounds <= 1000): returns c(rounds, cv_mse)
xgb_cv <- function(Ctr, y, max_depth) {
  seed <- sample.int(1e6, 1)
  par <- list(objective = "reg:squarederror", eta = 0.05, max_depth = as.integer(max_depth),
              subsample = 0.8, colsample_bytree = 0.8, nthread = 1L, seed = as.integer(seed))
  if (XGB_BACKEND == "python") {
    cv <- pyxgb$cv(do.call(dict, par), pyxgb$DMatrix(r_to_py(Ctr), label = r_to_py(y)),
                   num_boost_round = 1000L, nfold = 5L, early_stopping_rounds = 25L,
                   seed = as.integer(seed), verbose_eval = FALSE)
    rmse <- unlist(cv[["test-rmse-mean"]])
    c(rounds = length(rmse), cv_mse = tail(rmse, 1)^2)
  } else {
    cv <- xgboost::xgb.cv(params = par, data = xgboost::xgb.DMatrix(Ctr, label = y),
                          nrounds = 1000, nfold = 5, early_stopping_rounds = 25, verbose = 0)
    best <- if (!is.null(cv$best_iteration)) cv$best_iteration else cv$early_stop$best_iteration
    c(rounds = best, cv_mse = cv$evaluation_log$test_rmse_mean[best]^2)
  }
}
tune_xgb <- function(Ctr, y, depths = c(2, 4, 6)) {
  res  <- sapply(depths, function(dd) xgb_cv(Ctr, y, dd))
  best <- which.min(res["cv_mse", ])
  list(max_depth = depths[best], nrounds = unname(res["rounds", best]))
}

# ------------------------------------------------------------------ data
raw <- read.csv(if (file.exists("nhefs.csv")) "nhefs.csv" else "nhefs_complete.csv")
d <- raw %>%
  filter(!is.na(wt82), smokeintensity <= 25) %>%
  transmute(X = smkintensity82_71, Y = wt82_71,
            sex = factor(sex), race = factor(race), education = factor(education),
            exercise = factor(exercise), active = factor(active),
            age, smokeintensity, smokeyrs, wt71)
n <- nrow(d)
cat("n =", n, "\n")

Cmat <- model.matrix(~ sex + race + education + exercise + active + age +
                       smokeintensity + smokeyrs + wt71, data = d)[, -1]
# hyperparameters for XGBoost-tuned: chosen once by 5-fold CV prediction error
hp_m <- tune_xgb(Cmat, d$X)
hp_h <- tune_xgb(Cmat, d$Y)
cat("tuned m: depth", hp_m$max_depth, "rounds", hp_m$nrounds,
    "| tuned h: depth", hp_h$max_depth, "rounds", hp_h$nrounds, "\n")

rhs_lin <- "sex + race + education + exercise + active + age + smokeintensity + smokeyrs + wt71"
rhs_gam <- "sex + race + education + exercise + active + s(age) + s(smokeintensity) + s(smokeyrs) + s(wt71)"

# learner: fit on rows `tr`, predict rows `te`; returns list(m, h)
learn <- function(L, tr, te) {
  if (L == "OLS") {
    fm <- lm(as.formula(paste("X ~", rhs_lin)), data = d[tr, ])
    fh <- lm(as.formula(paste("Y ~", rhs_lin)), data = d[tr, ])
    list(m = predict(fm, d[te, ]), h = predict(fh, d[te, ]))
  } else if (L == "GAM") {
    fm <- gam(as.formula(paste("X ~", rhs_gam)), data = d[tr, ], method = "REML")
    fh <- gam(as.formula(paste("Y ~", rhs_gam)), data = d[tr, ], method = "REML")
    list(m = as.numeric(predict(fm, d[te, ])), h = as.numeric(predict(fh, d[te, ])))
  } else if (L == "XGBoost-tuned") {
    list(m = xgb_fit_predict(Cmat[tr, , drop = FALSE], d$X[tr], Cmat[te, , drop = FALSE],
                             hp_m$max_depth, hp_m$nrounds),
         h = xgb_fit_predict(Cmat[tr, , drop = FALSE], d$Y[tr], Cmat[te, , drop = FALSE],
                             hp_h$max_depth, hp_h$nrounds))
  } else {
    list(m = xgb_fit_predict(Cmat[tr, , drop = FALSE], d$X[tr], Cmat[te, , drop = FALSE]),
         h = xgb_fit_predict(Cmat[tr, , drop = FALSE], d$Y[tr], Cmat[te, , drop = FALSE]))
  }
}

est <- function(m, h) {
  xr <- d$X - m; yr <- d$Y - h
  tau <- sum(xr * yr) / sum(xr^2)
  psi <- xr * (yr - tau * xr)
  c(tau = tau, se = sqrt(mean(psi^2) / mean(xr^2)^2 / n), den = mean(xr^2))
}

learners <- c("OLS", "GAM", "XGBoost", "XGBoost-tuned")
S <- 20     # number of random 5-fold splits for cross-fitting

# ------------------------------------------------------------------ full-sample fits
full <- lapply(learners, function(L) learn(L, seq_len(n), seq_len(n)))
names(full) <- learners
full_est <- t(sapply(full, function(f) est(f$m, f$h)))

# ------------------------------------------------------------------ cross-fitted fits
cf <- lapply(learners, function(L) {
  lapply(seq_len(S), function(s) {
    fold <- sample(rep(1:5, length.out = n))
    m <- h <- numeric(n)
    for (k in 1:5) {
      te <- which(fold == k)
      pr <- learn(L, which(fold != k), te)
      m[te] <- pr$m; h[te] <- pr$h
    }
    list(m = m, h = h, est = est(m, h),
         cv_m = mean((d$X - m)^2), cv_h = mean((d$Y - h)^2))
  })
})
names(cf) <- learners

agg <- function(runs) {   # median aggregation over splits (Chernozhukov et al. 2018)
  tau <- sapply(runs, function(r) r$est[["tau"]])
  se  <- sapply(runs, function(r) r$est[["se"]])
  tm  <- median(tau)
  c(tau = tm, se = median(sqrt(se^2 + (tau - tm)^2)), split_sd = sd(tau),
    cv_m = mean(sapply(runs, `[[`, "cv_m")), cv_h = mean(sapply(runs, `[[`, "cv_h")))
}
cf_est <- t(sapply(cf, agg))

# ------------------------------------------------------------------ diagnostics
# reference learner B = cross-fitted XGBoost; average its predictions over splits
mB <- rowMeans(sapply(cf$XGBoost, `[[`, "m"))
hB <- rowMeans(sapply(cf$XGBoost, `[[`, "h"))
tauB <- cf_est["XGBoost", "tau"]
ref <- sapply(learners, function(L) {
  gA <- full[[L]]$h - tauB * full[[L]]$m
  gB <- hB - tauB * mB
  mean((full[[L]]$m - mB) * (gA - gB)) / full_est[L, "den"]
})

tab <- data.frame(
  learner = learners,
  cv_mse_m = cf_est[, "cv_m"], cv_mse_h = cf_est[, "cv_h"],
  tau_full = full_est[, "tau"], se_full = full_est[, "se"],
  tau_cf = cf_est[, "tau"], se_cf = cf_est[, "se"], split_sd = cf_est[, "split_sd"],
  gap_full_minus_cf = full_est[, "tau"] - cf_est[, "tau"],
  ref_remainder = ref,
  row.names = NULL)
# for XGBoost, ref_remainder compares its full-sample fit with its own cross-fitted fit
print(tab, digits = 3)
dir.create("v4_output", showWarnings = FALSE)
write.csv(tab, "v4_output/nhefs_results.csv", row.names = FALSE)
