# When Prediction Error Is Not Enough

Simulation code for the study "When Prediction Error Is Not Enough: Evaluating Nuisance-Function Prediction for Causal Estimation."

## Overview

This study examines the relationship between nuisance-function prediction error and causal estimation performance in a partially linear model.

Methods include:

- OLS
- GAM
- XGBoost
- DML-XGBoost
- Oracle estimation

Performance measures include prediction error, bias, RMSE, standard error, confidence interval coverage, and joint nuisance-function error.

Two simulation settings are considered:

1. Independent observations
2. Clustered observations

Each setting uses 500 Monte Carlo replicates. The clustered analysis uses cluster-level cross-fitting and variance estimation.

## Files

The repository contains R scripts for simulation, analysis, and figure generation, together with saved RData workspaces and checkpoints.

`figures.R` generates the figures reported in the analysis.

Some scripts use previously saved RData objects from the simulation workflow.

## Data

No human participant data are used. All data are simulated from the data-generating mechanisms described in the manuscript.

## Citation

Cao, C. When Prediction Error Is Not Enough: Evaluating Nuisance-Function Prediction for Causal Estimation.

## Author

Cong Cao, PhD  
Department of Biostatistics, Yale University  
New Haven, CT, USA

Contact: cong.cao@yale.edu
