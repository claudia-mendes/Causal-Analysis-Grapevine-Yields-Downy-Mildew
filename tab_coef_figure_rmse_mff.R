# =========================
# Script 3 (Final - without model_labels)
# =========================

library(pROC)
library(gbm)
library(lme4)
library(lmerTest)
library(Matching)
library(dplyr)
library(tibble)
library(readr)

# --- Data Setup ---
set.seed(3)
df <- readRDS("output_files/vines_combined_2000_2023_treat_SP_pests_selected_departements.rda")

df$TREAT <- as.numeric(as.character(df$treat2_may))
df$MM_ST_5 <- scale(df$MM_ST_5)
df$MM_SWI_5 <- scale(df$MM_SWI_5)
df$MM_RG_5 <- scale(df$MM_RG_5)
df$MM_ETP_5 <- scale(df$MM_ETP_5)
df$MM_NJT30_5 <- scale(df$MM_NJT30_5)
df$MM_SP_5 <- scale(df$MM_SP_5)
df$yield <- log(df$yield)
df$Departement <- as.factor(df$Departement)

dir.create("model_plots_DF", showWarnings = FALSE)

# --- Model Formulas ---
model_forms_lmer <- list(
  LMd1 = yield ~ predicted_MFF + MM_ST_5 + MM_NJT30_5 + MM_ETP_5 + MM_SWI_5 + MM_RG_5 + (1 | Departement),
  LMd2 = yield ~ predicted_MFF + TREAT + MM_ST_5 + MM_NJT30_5 + MM_ETP_5 + MM_SWI_5 + MM_RG_5 + (1 | Departement),
  LM1 = yield ~ TREAT + MM_ST_5 + MM_NJT30_5 + MM_ETP_5 + MM_SWI_5 + MM_RG_5 + (1 | Departement)
)

model_forms_gbm <- list(
  LMd1 = yield ~ predicted_MFF_gbm + MM_ST_5 + MM_NJT30_5 + MM_ETP_5 + MM_SWI_5 + MM_RG_5 + Departement,
  LMd2 = yield ~ predicted_MFF_gbm + TREAT + MM_ST_5 + MM_NJT30_5 + MM_ETP_5 + MM_SWI_5 + MM_RG_5 + Departement,
  LM1 = yield ~ TREAT + MM_ST_5 + MM_NJT30_5 + MM_ETP_5 + MM_SWI_5 + MM_RG_5 + Departement
)

# --- Outputs ---
yield_preds <- list()
for (name in names(model_forms_lmer)) yield_preds[[name]] <- list(lmer = c(), gbm = c())
rmse_table <- tibble(Model = character(), RMSE = numeric(), Type = character())
List_Year <- unique(df$YEAR)
Obs <- c()
ClassT <- c()

# --- Cross-validation loop ---
for (Year_i in List_Year) {
  df_mi <- df[df$YEAR != Year_i, ]
  df_i  <- df[df$YEAR == Year_i, ]
  
  Proba_glmer_MFF <- glmer(
    cbind(MFF, 100 - MFF) ~ TREAT + MM_ST_5 + MM_NJT30_5 + MM_ETP_5 + MM_SWI_5 + MM_RG_5 + (1 | Departement),
    family = binomial,
    data = df_mi
  )
  
  df_mi$predicted_MFF <- predict(Proba_glmer_MFF, type = "response")
  df_i$predicted_MFF <- predict(Proba_glmer_MFF, type = "response", newdata = df_i, allow.new.levels = TRUE)
  
  Proba_GBM <- gbm(
    formula = MFF / 100 ~ TREAT + MM_ST_5 + MM_NJT30_5 + MM_ETP_5 + MM_SWI_5 + MM_RG_5 + Departement,
    data = df_mi,
    distribution = "bernoulli",
    n.trees = 100
  )
  
  df_mi$predicted_MFF_gbm <- predict(Proba_GBM, newdata = df_mi, n.trees = 100, type = "response")
  df_i$predicted_MFF_gbm <- predict(Proba_GBM, newdata = df_i, n.trees = 100, type = "response")
  
  for (model_name in names(model_forms_lmer)) {
    lmer_model <- lmer(model_forms_lmer[[model_name]], data = df_mi)
    pred_lmer <- predict(lmer_model, newdata = df_i, allow.new.levels = TRUE)
    
    gbm_model <- gbm(
      formula = model_forms_gbm[[model_name]],
      data = df_mi,
      distribution = "gaussian",
      n.trees = 100
    )
    pred_gbm <- predict(gbm_model, newdata = df_i, n.trees = 100)
    
    yield_preds[[model_name]]$lmer <- c(yield_preds[[model_name]]$lmer, pred_lmer)
    yield_preds[[model_name]]$gbm <- c(yield_preds[[model_name]]$gbm, pred_gbm)
  }
  
  Obs <- c(Obs, df_i$yield)
  ClassT <- c(ClassT, df_i$TREAT)
}

# --- Fit full GLMER on filtered data ---
df_mff <- df %>% filter(!is.na(MFF))
Proba_glmer_MFF_full <- glmer(
  cbind(MFF, 100 - MFF) ~ TREAT + MM_ST_5 + MM_NJT30_5 + MM_ETP_5 + MM_SWI_5 + MM_RG_5 + (1 | Departement),
  family = binomial,
  data = df_mff
)
df_mff$predicted_MFF <- predict(Proba_glmer_MFF_full, type = "response")

# --- Extract and format full model table ---
coef_df <- summary(Proba_glmer_MFF_full)$coefficients
proba_glmer_table <- tibble(
  Variable = rownames(coef_df),
  Estimate = round(coef_df[, 1], 2),
  Std_Error = round(coef_df[, 2], 2),
  P_Value = round(coef_df[, 4], 4)
) %>%
  mutate(P_Value = ifelse(P_Value < 0.001, "< 0.001***", as.character(P_Value)))

# --- LMER yield summaries ---
yield_models_table <- tibble()
for (model_name in names(model_forms_lmer)) {
  fit <- lmer(model_forms_lmer[[model_name]], data = df_mff)
  smry <- summary(fit)$coefficients
  table <- tibble(
    Model = model_name,
    Variable = rownames(smry),
    Estimate = round(smry[, "Estimate"], 2),
    Std_Error = round(smry[, "Std. Error"], 2),
    P_Value = formatC(smry[, "Pr(>|t|)"], format = "e", digits = 2)
  )
  yield_models_table <- bind_rows(yield_models_table, table)
}

# --- Plot generation with accurate RMSE ---
for (model_name in names(model_forms_lmer)) {
  pred_lmer <- yield_preds[[model_name]]$lmer
  pred_gbm  <- yield_preds[[model_name]]$gbm
  
  rmse_lmer <- sqrt(mean((exp(Obs) - exp(pred_lmer))^2))
  rmse_gbm  <- sqrt(mean((exp(Obs) - exp(pred_gbm))^2))
  
  rmse_table <- bind_rows(
    rmse_table,
    tibble(Model = model_name, RMSE = as.numeric(sprintf("%.2f", rmse_lmer)), Type = "LMER"),
    tibble(Model = model_name, RMSE = as.numeric(sprintf("%.2f", rmse_gbm)), Type = "GBM")
  )
  
  filename <- model_name
  png(filename = paste0("model_plots/", filename, "_mff.png"), width = 960, height = 480)
  par(mfrow = c(1, 2), oma = c(2, 2, 2, 2))
  
  x_range <- range(c(exp(pred_lmer), exp(pred_gbm)), na.rm = TRUE)
  y_range <- range(exp(Obs), na.rm = TRUE)
  
  plot(exp(pred_lmer), exp(Obs), main = "LMER",
       xlab = "Predicted Yield", ylab = "Observed Yield",
       pch = 20, xlim = x_range, ylim = y_range)
  abline(0, 1, col = "blue")
  text(min(x_range), max(y_range)*0.95, paste("RMSE:", sprintf("%.2f", rmse_lmer)), adj = c(0, 1))
  mtext("A", side = 3, line = -1.5, adj = 0, outer = TRUE, cex = 1.5, font = 2)
  
  plot(exp(pred_gbm), exp(Obs), main = "GBM",
       xlab = "Predicted Yield", ylab = "Observed Yield",
       pch = 20, xlim = x_range, ylim = y_range)
  abline(0, 1, col = "blue")
  text(min(x_range), max(y_range)*0.95, paste("RMSE:", sprintf("%.2f", rmse_gbm)), adj = c(0, 1))
  mtext("B", side = 3, line = -1.5, adj = 0.5, outer = TRUE, cex = 1.5, font = 2)
  
  dev.off()
}

# --- Save tables ---
write_csv(proba_glmer_table, "model_plots/proba_glmer_mff.csv")
write_csv(yield_models_table, "model_plots/yield_models_mff.csv")
write_csv(rmse_table, "model_plots/rmse_all_models_mff.csv")
print(rmse_table)




# --- Print coefficients of MM_ST_5 and MM_SWI_5 in LMd2 model (coefficients appeared as 0 in the table) ---
fit_LMd2 <- lmer(model_forms_lmer[["LMd2"]], data = df_mff)
coef_LMd2 <- summary(fit_LMd2)$coefficients

# Print just MM_ST_5 and MM_RG_5 rows
print(coef_LMd2[rownames(coef_LMd2) %in% c("MM_ST_5", "MM_SWI_5"), ])

