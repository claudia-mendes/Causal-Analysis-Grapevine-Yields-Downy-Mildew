# Causal-Analysis-Grapevine-Yields-Downy-Mildew

This repository contains R scripts developed for a study investigating the causal effects of high precipitation (HP) events on grapevine yields and downy mildew (DM) outbreaks in France. The analyses focus on the growing season months (April to September) from 2000 to 2023 and apply causal inference and machine learning methods to quantify direct and indirect yield impacts.

Three main types of analyses and scripts are included:
- Estimation of HP effects on grapevine yields - Script "Causal_estimate_SP_gbm_dep_vines.R"
- Estimation of HP effects on downy mildew (DM) incidence and severity - Script "dr_pests_binomial_glmm.R"
- Yield modeling using DM and HP as predictors - Script "tab_coef_figure_rmse_mff.R"

# Scripts

# "Causal_estimate_SP_gbm_dep_vines.R"

Purpose: Estimates the causal effects of high precipitation events on grapevine yields using double robust (DR) methods.

Methods: Combines outcome models (GLMM or GBM) with treatment models (logistic regression or GBM), adjusting for covariates and regional effects (French Departments).

Inputs: vines_combined_2000_2023_treat_SP.rda

Outputs: Monthly text files with point estimates and bootstrapped confidence intervals.

Note: Run this script separately for each month from April to September by updating the treatment variable suffix (e.g., "apr", "may") and the corresponding climate variable index (e.g., _4, _5).


# "dr_pests_binomial_glmm.R"

Purpose: Estimates the causal effect of HP on downy mildew incidence and severity using DR models.

Methods: Uses binomial GLMM and GBM models; includes bootstrap resampling (200 iterations) for confidence interval estimation.

Inputs: output_files/vines_combined_2000_2023_treat_SP_pests_selected_departements.rda

Outputs: Monthly results saved in text files under results_DR_pests/.

Note: Note: Run this script separately for each month from April to September by updating the treatment variable suffix (e.g., "apr", "may") and the corresponding climate variable index (e.g., _4, _5).

MFF = downy mildew incidence (DI)

MFI = downy mildew severity (DS)

To switch between MFF and MFI, replace the relevant variable references in the script.


# "tab_coef_figure_rmse_mff.R"

Purpose: Models the effect of downy mildew incidence (MFF = DI) and high precipitation (HP) in May on grapevine yield using both statistical (LMER) and machine learning (GBM) approaches.

Methods: The analysis focuses specifically on May HP events (treat2_may) and May climate variables (e.g., MM_SP_5, MM_ST_5). The yield variable is log-transformed before modeling. Predictions and RMSE values are back-transformed for interpretation. Predicts MFF using both GLMER and GBM models. Applies leave-one-year-out cross-validation to assess model performance.

Fits three yield models:

- LM1: yield ~ HP + climate covariates
- LMd1: yield ~ predicted MFF + climate covariates
- LMd2: yield ~ predicted MFF + HP + climate covariates

Inputs: output_files/vines_combined_2000_2023_treat_SP_pests_selected_departements.rda

Outputs: CSV files: model_plots/proba_glmer_mff.csv, yield_models_mff.csv, rmse_all_models_mff.csv, PNG plots: Predicted vs. observed yield (model_plots/*.png) for each model

Notes:

MFF = downy mildew incidence (DI)

MFI = downy mildew severity (DS)

To analyze MFI instead of MFF, replace variable references accordingly. This script analyzes May only. For other months, use or adapt different scripts.


# Notes on January–March Analyses

Separate scripts, not included in this repository, were used for the months January to March. These scripts excluded the NDT30 variable (monthly number of days with maximum daily temperatures above 30°C), which is only relevant during the active growing season (April to September).
