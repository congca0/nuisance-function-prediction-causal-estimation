 # nuisance-function-prediction-causal-estimation

R code for my paper "When Prediction Error Is Not Enough: Evaluating Nuisance-Function
Prediction for Causal Estimation" (arXiv:2609.00071).

There are three scripts:

    simulation_revised_v3.R      scenario 1, independent data
    simulation_revised_2_v3.R    scenario 2, clustered data
    make_tables_figures_v3.R     tables 1-2 and figures 1-2

You need mgcv, xgboost, dplyr and ggplot2. Run the two simulation scripts first, then
make_tables_figures_v3.R, all from the same folder:

    Rscript simulation_revised_v3.R
    Rscript simulation_revised_2_v3.R
    Rscript make_tables_figures_v3.R

Each scenario is 500 replicates, so the simulations take a while. Results go into
v3_output/. Scenario 1 saves a checkpoint every 5 replicates and continues from it if it
stops. Scenario 2 saves one every 50 replicates but starts over if you rerun it.

What changed from v2:

- the oracle now uses the training sample, like the other methods
- scenario 2: h0 is now E(Y|C) = tau*m0 + g0 (v2 used g0, which was wrong)
- DML prediction error is computed on the evaluation sample and averaged over folds
- the cluster-robust variance now has the G/(G-1) correction
- the second-order remainder R is saved for each replicate (Table 2)

After these fixes DML-XGBoost has the lowest RMSE of the non-oracle methods, which is not
what v2 found. The v2 code is still in the commit history.

All data are simulated.

Cong Cao, cong.cao@yale.edu
