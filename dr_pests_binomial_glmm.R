###############################
### DR ANALYSIS FOR MFI DATA #
###      BINOMIAL GLMM       #
###############################

library(lme4)
library(gbm)
library(dplyr)
library(tidyr)

set.seed(123)

# Load input dataset
df <- readRDS("output_files/vines_combined_2000_2023_treat_SP_pests_selected_departements.rda")

results_list <- list()

# -------------------- Run DR estimation for one month --------------------
month <- "sep"
treat_col <- "treat2_sep"
cov_suffix <- "9"
covariates <- paste0("MM_", c("ST", "RG", "ETP", "SWI", "NJT30"), "_", cov_suffix)

# Check that required columns are present
if (!treat_col %in% colnames(df)) {
  cat("Skipping month", month, "- treatment column not found\n")
} else if (any(!covariates %in% colnames(df))) {
  missing_covs <- covariates[!covariates %in% colnames(df)]
  cat("Skipping month", month, "- missing covariates:", paste(missing_covs, collapse = ", "), "\n")
} else {
  
  # Prepare data subset for the selected month
  df$TREAT <- as.numeric(as.character(df[[treat_col]]))
  df_month <- df[!is.na(df$TREAT), ]
  
  if (nrow(df_month) == 0) {
    cat("Skipping month", month, "- no non-NA treatment rows\n")
  } else {
    
    # Standardize covariates
    for (cov in covariates) {
      df_month[[cov]] <- scale(df_month[[cov]])
    }
    
    dr_lm_boot <- c()
    dr_gbm_boot <- c()
    
    # -------------------- Bootstrap loop --------------------
    for (b in 1:200) {
      # Bootstrap resampling by department
      Tab_boot <- do.call(rbind, lapply(
        split(df_month, df_month$Departement),
        function(subdf) subdf[sample(nrow(subdf), replace = TRUE), ]
      ))
      
      # Clean bootstrap data
      Tab_boot <- Tab_boot[!is.na(Tab_boot$Departement), ]
      Tab_boot$Departement <- droplevels(as.factor(Tab_boot$Departement))
      if (nlevels(Tab_boot$Departement) < 2) next
      
      # Prepare binomial counts for MFI response
      Tab_boot <- Tab_boot[!is.na(Tab_boot$MFI), ]
      Tab_boot$MFI <- round(Tab_boot$MFI)
      Tab_boot$MFI[Tab_boot$MFI < 0] <- 0
      Tab_boot$MFI[Tab_boot$MFI > 100] <- 100
      Tab_boot$MFI <- as.integer(Tab_boot$MFI)
      Tab_boot$Failures <- 100L - Tab_boot$MFI
      
      form_covs <- paste(covariates, collapse = " + ")
      
      # Fit GLMM for treatment assignment (propensity score)
      fit_success <- tryCatch({
        Proba_glmer <- glmer(
          formula = as.formula(paste0("cbind(TREAT, 1 - TREAT) ~ ", form_covs, " + (1|Departement)")),
          family = binomial, data = Tab_boot
        )
        TRUE
      }, error = function(e) FALSE)
      
      if (!fit_success) next
      
      # Fit binomial GLMM for MFI outcome
      fit_success2 <- tryCatch({
        MFI_glmer <- glmer(
          formula = as.formula(paste0("cbind(MFI, Failures) ~ TREAT + ", form_covs, " + (1|Departement)")),
          family = binomial(link = "logit"),
          data = Tab_boot
        )
        cat("Model family:", family(MFI_glmer)$family, "- Link:", family(MFI_glmer)$link, "\n")
        TRUE
      }, error = function(e) {
        message("Model fitting failed: ", e$message)
        FALSE
      })
      
      if (!fit_success2) next
      
      # DR estimation using GLMM
      Score_glmer <- predict(Proba_glmer, type = "response")
      TAB <- data.frame(X = Tab_boot$TREAT, Y = Tab_boot$MFI / 100, Score = Score_glmer)
      TAB <- TAB[Score_glmer > 0.1 & Score_glmer < 0.9, ]
      
      TABPred <- Tab_boot[Score_glmer > 0.1 & Score_glmer < 0.9, ]
      TABPred$TREAT <- 1
      Pred_1 <- predict(MFI_glmer, newdata = TABPred, type = "response", allow.new.levels = TRUE)
      
      TABPred$TREAT <- 0
      Pred_0 <- predict(MFI_glmer, newdata = TABPred, type = "response", allow.new.levels = TRUE)
      
      cat("Bootstrap", b, "- Pred_1 range:", range(Pred_1), "\n")
      cat("Bootstrap", b, "- Pred_0 range:", range(Pred_0), "\n")
      
      # Compute double robust estimator (GLMM)
      DR1 <- mean(Pred_1 + (TAB$X / TAB$Score) * (TAB$Y - Pred_1))
      DR0 <- mean(Pred_0 + ((1 - TAB$X) / (1 - TAB$Score)) * (TAB$Y - Pred_0))
      dr_lm_boot <- c(dr_lm_boot, DR1 - DR0)
      
      # -------------------- DR estimation using GBM --------------------
      Proba_GBM <- gbm(as.formula(paste0("TREAT ~ ", form_covs, " + Departement")),
                       data = Tab_boot, distribution = "bernoulli", n.trees = 100)
      MFI_GBM <- gbm(as.formula(paste0("MFI / 100 ~ TREAT + ", form_covs, " + Departement")),
                     data = Tab_boot, distribution = "gaussian", n.trees = 100)
      Score_GBM <- predict(Proba_GBM, type = "response", n.trees = 100)
      
      TAB2 <- data.frame(X = Tab_boot$TREAT, Y = Tab_boot$MFI / 100, Score = Score_GBM)
      TAB2 <- TAB2[Score_GBM > 0.1 & Score_GBM < 0.9, ]
      
      TABPred <- Tab_boot[Score_GBM > 0.1 & Score_GBM < 0.9, ]
      TABPred$TREAT <- 1
      Pred_1 <- predict(MFI_GBM, newdata = TABPred, n.trees = 100)
      
      TABPred$TREAT <- 0
      Pred_0 <- predict(MFI_GBM, newdata = TABPred, n.trees = 100)
      
      # Compute double robust estimator (GBM)
      DR1 <- mean(Pred_1 + (TAB2$X / TAB2$Score) * (TAB2$Y - Pred_1))
      DR0 <- mean(Pred_0 + ((1 - TAB2$X) / (1 - TAB2$Score)) * (TAB2$Y - Pred_0))
      dr_gbm_boot <- c(dr_gbm_boot, DR1 - DR0)
      
      cat("Month:", month, "Bootstrap:", b, "\n")
    }
    
    # Save results for the current month
    month_result <- data.frame(
      Month = month,
      dr_lm = dr_lm_boot,
      dr_gbm = dr_gbm_boot
    )
    results_list[[month]] <- month_result
  }
}

# Combine and write all results to file
final_results <- bind_rows(results_list)
write.table(final_results, file = "results_DR_pests/results_DR_MFI_binomial_sep_200boot.txt", row.names = FALSE)
