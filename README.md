# Causal-Analysis-Grapevine-Yields-Downy-Mildew

This repository contains R scripts developed for a study investigating the causal effects of high precipitation (HP) events on grapevine yields and downy mildew (DM) outbreaks in France. The analyses focus on the growing season months (April to September) from 2000 to 2023 and apply causal inference and machine learning methods to quantify direct and indirect yield impacts.

Three main types of analyses are included:
- Estimation of HP effects on grapevine yields. Script "Causal_estimate_SP_gbm_dep_vines.R"
- Estimation of HP effects on downy mildew (DM) incidence and severity. Script "dr_pests_binomial_glmm.R"
- Yield modeling using DM and HP as predictors, with model validation. Script "tab_coef_figure_rmse_mff.R"

# Scripts

# "Causal_estimate_SP_gbm_dep_vines.R"

Purpose: Estimates the causal effects of high precipitation events on grapevine yields using double robust (DR) methods.

Methods: Combines outcome models (GLMM or GBM) with treatment models (logistic regression or GBM), adjusting for covariates and regional effects.

Inputs: vines_combined_2000_2023_treat_SP.rda

Outputs: Monthly text files with point estimates and bootstrapped confidence intervals.

Note: Run this script separately for each month from April to September by updating the treatment variable suffix (e.g., "apr", "may") and the corresponding climate variable index (e.g., _4, _5).


# "dr_pests_binomial_glmm.R"

Purpose: Estimates the causal effect of HP on downy mildew incidence and severity using DR models.

Methods: Uses binomial GLMM and GBM models; includes bootstrap resampling (200 iterations) for confidence interval estimation.

Inputs: output_files/vines_combined_2000_2023_treat_SP_pests_selected_departements.rda

Outputs: Monthly results saved in text files under results_DR_pests/.

Note: Run this script separately for each month from April to September.

MFF = downy mildew incidence (DI)

MFI = downy mildew severity (DS)

To switch between MFF and MFI, replace the relevant variable references in the script.


# "tab_coef_figure_rmse_mff.R"

Purpose: Models the effect of downy mildew incidence (MFF = DI) and HP (May) on grapevine yield using both statistical and machine learning approaches.

Methods: Predicts MFF using GLMER and GBM models. Fits three yield models (LM1, LMd1, LMd2) with both LMER and GBM. Performs leave-one-year-out cross-validation and calculates RMSE.

Inputs: output_files/vines_combined_2000_2023_treat_SP_pests_selected_departements.rda

Outputs: Tables and PNG plots for model coefficients, predicted vs. observed values, and RMSE values.

MFF = downy mildew incidence (DI)

MFI = downy mildew severity (DS)

To switch between MFF and MFI, replace the relevant variable references in the script.


# Notes on January–March Analyses

Separate scripts were used for the months January to March. These early-season analyses excluded the NDT30 variable, which is only relevant during the active growing season (April to
September). These scripts are not included in this repository.
