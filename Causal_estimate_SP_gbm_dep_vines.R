##################################################
### SCRIPT 1 - Superposed Epoch Analysis (SEA) ###
##################################################

# Load necessary libraries
library(dplyr)
library(ggplot2)
library(readr)
library(tidyr)
library(purrr)
library(openxlsx)

# Load the dataset
df <- readRDS("output_files/vines_combined_2000_2023_treat_SP.rda")

# Display column names
print(colnames(df))

df$TREAT <- as.numeric(as.character(df$treat2_apr))

# Filter for sp between 2002 and 2021, and rename columns appropriately
sp <- df %>%
  filter(TREAT == 1, YEAR >= 2002, YEAR <= 2021) %>%
  rename(Department = Departement, sp_Year = YEAR) %>%
  mutate(Occurrence_id = paste(Department, sp_Year, mar = "_"))  # Assign a unique ID to each sp occurrence

# Generate co-occurring sp years (+/- 2 years around each sp year)
sp_expanded <- sp %>%
  group_by(Department) %>%
  arrange(sp_Year) %>%
  mutate(Cooccurring = map(sp_Year, ~ seq(.x - 2, .x + 2))) %>%
  unnest(cols = c(Cooccurring))

# Filter co-occurring sp and create Excluded_Year
cooccurring_sp <- sp_expanded %>%
  filter(Cooccurring %in% sp_Year & Cooccurring != sp_Year) %>%
  mutate(Excluded_Year = Cooccurring)

# Select relevant columns using base R indexing
cooccurring_sp <- cooccurring_sp[, c("Department", "sp_Year", "Excluded_Year", "Occurrence_id")]

# Rename the necessary columns in df
df_renamed <- df %>%
  rename(Department = Departement, Yield = yield)

# Join sp data with yield data
sp_yields <- sp %>%
  inner_join(df_renamed, by = c("Department", "sp_Year" = "YEAR"))

# Calculate the number of sp occurrences analyzed
number_of_sp_occurrences <- sp_yields %>%
  distinct(Department, sp_Year, Occurrence_id) %>%
  nrow()

# Print the number of sp occurrences analyzed
print(paste("Number of sp occurrences analyzed:", number_of_sp_occurrences))

# List of departments with sp occurrences
departments_with_sp <- unique(sp_yields$Department)

# Ensure that yield data for simulation is only from these departments
yield_data <- df_renamed %>%
  filter(Department %in% departments_with_sp)

# Generate occurrence ID and calculate the yield window (2 years before and 2 years after the event)
yield_window <- sp_yields %>%
  rowwise() %>%
  do({
    data.frame(Occurrence_id = .$Occurrence_id,
               Event_year = .$sp_Year,
               Department = .$Department,
               Centered_Year = -2:2,
               Year_yield = .$sp_Year + -2:2)
  }) %>%
  ungroup()

# Merge with yield data and check if Yield column exists
yield_window <- yield_window %>%
  left_join(yield_data, by = c("Department", "Year_yield" = "YEAR"))

# Ensure the Yield column exists
if (!"Yield" %in% colnames(yield_window)) {
  stop("The Yield column was not found after the join. Check the join operation.")
}

# Add a flag for co-occurring sp by joining with cooccurring_sp
yield_window <- yield_window %>%
  left_join(cooccurring_sp, by = c("Occurrence_id", "Department")) %>%
  mutate(Excluded_sp_Year = ifelse(!is.na(Excluded_Year) & Year_yield == Excluded_Year, "Yes", "No"))

# Add a flag for co-occurring sp
yield_window <- yield_window %>%
  mutate(Exclude = Excluded_sp_Year == "Yes")

# Calculate baseline excluding sp year and co-occurring sp years
baseline <- yield_window %>%
  filter(Centered_Year %in% c(-2, -1, 1, 2) & !Exclude) %>%
  group_by(Occurrence_id, Department, Event_year) %>%
  summarise(Baseline_mean = mean(Yield, na.rm = TRUE), .groups = 'drop')

# Normalize yields by baseline
normalized_yields <- yield_window %>%
  left_join(baseline, by = c("Occurrence_id", "Department", "Event_year")) %>%
  mutate(Normalized_Yield = ifelse(!is.na(Baseline_mean) & !is.na(Yield), Yield / Baseline_mean, NA))

# Identify and print points of years where co-occurring sp occurred (excluding Centered_Year == 0)
excluded_points <- normalized_yields %>%
  filter(Exclude & Centered_Year != 0)

# Select relevant columns manually using base R
excluded_points <- excluded_points[, c("Department", "Year_yield", "Centered_Year")]

print("Excluded points due to co-occurring sp:")
print(head(excluded_points))

# Remove points of years where co-occurring sp occurred
normalized_yields <- normalized_yields %>%
  filter(!(Exclude & Centered_Year != 0))

# Create a composite time series
composite_time_series <- normalized_yields %>%
  group_by(Centered_Year) %>%
  summarise(Global_Mean_Yield = mean(Normalized_Yield, na.rm = TRUE), .groups = 'drop')

# Plot the composite time series
ggplot(composite_time_series, aes(x = Centered_Year, y = Global_Mean_Yield)) +
  geom_line() +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_hline(yintercept = 1, linetype = "solid", color = "black") +
  scale_x_continuous(breaks = -2:2) +
  ylim(0.8, 1.2) +  # Setting the y-axis limits
  labs(title = "Composite Time Series of Global Mean Normalized Yield in France",
       x = "Year from Event",
       y = "Global Mean of Normalized Yield") +
  theme_minimal()



##################################
### SCRIPT 2 - CAUSAL ANALYSIS ###
##################################

library(lme4)
library(gbm)
library(corrplot)
library(Matching)

set.seed(0123456789)

######
# Data #
######

df <- readRDS("output_files/vines_combined_2000_2023_treat_SP.rda")
# Choice of treatment
df$TREAT <- as.numeric(as.character(df$treat2_apr))

# Count the number of unique "Departement" in the dataset
num_departements <- length(unique(df$Departement))
print(paste("Number of unique Departements:", num_departements))

# Adjust left margin to add space for y-axis labels
par(mfrow = c(2, 2), mar = c(5, 6, 4, 2))  # Increase the left margin (second value)

# Plot each boxplot with adjusted margins and labels for each panel
boxplot(df$MM_ETP_4 ~ df$TREAT, 
        xlab = "low precipitation occurrence in May", 
        ylab = "Average ETP in May (mm)", 
        range = 0, 
        cex.lab = 1.2,  # Increase font size of axis labels
        cex.axis = 1.2) # Increase font size of axis tick numbers
mtext("A", side = 3, adj = -0.1, line = 2, font = 2, cex = 1.5)  # Add label A

boxplot(df$MM_SWI_4 ~ df$TREAT, 
        xlab = "low precipitation occurrence in May", 
        ylab = "Soil humidity index in May (%)",  # Split y-axis label into two lines
        range = 0, 
        cex.lab = 1.2, 
        cex.axis = 1.2)
mtext("B", side = 3, adj = -0.1, line = 2, font = 2, cex = 1.5)  # Add label B

boxplot(df$MM_RG_4 ~ df$TREAT, 
        xlab = "low precipitation occurrence in May", 
        ylab = expression("Average radiation in May (J/cm"^2*")"), 
        range = 0, 
        cex.lab = 1.2, 
        cex.axis = 1.2)
mtext("C", side = 3, adj = -0.1, line = 2, font = 2, cex = 1.5)  # Add label C

boxplot(df$MM_ST_4 ~ df$TREAT, 
        xlab = "low precipitation occurrence in May", 
        ylab = "Average maximum\ntemperature in May (°C)",  # Split y-axis label into two lines
        range = 0, 
        cex.lab = 1.2, 
        cex.axis = 1.2)
mtext("D", side = 3, adj = -0.1, line = 2, font = 2, cex = 1.5)  # Add label D



# Open a new graphical window with specified dimensions
dev.new(width = 10, height = 5)

# Set up the layout for two plots side by side with adjustable width and spacing
layout(matrix(c(1, 2), 1, 2), widths = c(1, 1))

Vec_clim <- data.frame(Yield=log(df$yield), ETP=df$MM_ETP_4, SP=df$MM_SP_4, RG=df$MM_RG_4, ST=df$MM_ST_4, SWI=df$MM_SWI_4)
par(fig = c(0, 0.48, 0, 1), new = TRUE)  # Left half with a small right margin for spacing
corrplot(cor(Vec_clim), tl.col = "black", tl.cex = 1.5, number.cex = 1.5)
par(fig = c(0.52, 1, 0, 0.8), new = TRUE)  # Right side with a small left margin for spacing
boxplot(log(df$yield) ~ df$TREAT, ylab="Log Yield (t ha-1)", xlab="low precipitation occurrence in May", range=0)

#########
#Scaling#
#########

df$MM_ST_4 <- scale(df$MM_ST_4)
df$MM_RG_4 <- scale(df$MM_RG_4)
df$MM_ETP_4 <- scale(df$MM_ETP_4)
df$MM_SWI_4 <- scale(df$MM_SWI_4)
df$MM_NJT30_4 <- scale(df$MM_NJT30_4)
df$yield <- log(df$yield)
Weight <- df$TREAT
Weight[df$TREAT==1] <- 1
Weight[df$TREAT==0] <- 1
df$Weight <- Weight

################
#Naive estimate#
################

AVG_naive <- mean(df$yield[df$TREAT==1]) - mean(df$yield[df$TREAT==0])
AVG_naive

############
#Models lme#
############

Proba_glmer <- glmer(cbind(TREAT, 1-TREAT) ~ MM_ST_4 + MM_RG_4 + MM_ETP_4 + MM_SWI_4 + MM_NJT30_4 +  (1|Departement), 
                     family=binomial, data=df)
summary(Proba_glmer)

Yield_lmer <- lmer(yield ~ TREAT + MM_ST_4 + MM_RG_4 + MM_ETP_4 + MM_SWI_4 +  MM_NJT30_4 + (1|Departement), data=df)
summary(Yield_lmer)

Score_glmer <- predict(Proba_glmer, type="response")

############
#Models GBM#
############

Proba_GBM <- gbm(TREAT ~ MM_ST_4 + MM_RG_4 + MM_ETP_4 + MM_SWI_4 +  MM_NJT30_4 + Departement, data=df, distribution="bernoulli")
Proba_GBM

Yield_GBM <- gbm(yield ~ TREAT + MM_ST_4 + MM_RG_4 + MM_ETP_4 + MM_SWI_4 +  MM_NJT30_4 + Departement, data=df)
Yield_GBM

Score_GBM <- predict(Proba_GBM, type="response")

#######################
#Estimate based on glm#
#######################

Mod_Yield <- Yield_lmer
Score <- Score_glmer

#Inverse probability weighting#
TAB <- data.frame(X=df$TREAT, Y=df$yield, Score)
#Remove extreme scores
TAB <- TAB[(Score > 0.1 & Score < 0.9),]

E1 <- mean(TAB$X * TAB$Y / TAB$Score)
E0 <- mean((1 - TAB$X) * TAB$Y / (1 - TAB$Score))

AVG_ipw_lm <- E1 - E0
AVG_ipw_lm

#Standardization

TABPred <- df
TABPred <- TABPred[(Score > 0.1 & Score < 0.9),]

TABPred$TREAT <- rep(1, nrow(TABPred))
Pred_4 <- predict(Mod_Yield, newdata = TABPred)
#summary(Pred_4)

TABPred$TREAT <- rep(0, nrow(TABPred))
Pred_0 <- predict(Mod_Yield, newdata=TABPred)
#summary(Pred_0)

AVG_sdz_lm <- mean(Pred_4) - mean(Pred_0)
AVG_sdz_lm

#Double robust

DR1 <- mean(Pred_4 + (TAB$X / TAB$Score) * (TAB$Y - Pred_4))
DR0 <- mean(Pred_0 + ((1 - TAB$X) / (1 - TAB$Score)) * (TAB$Y - Pred_0))

AVG_dr_lm <- DR1 - DR0
AVG_dr_lm

#Matching

Appariement_lm <- Match(Y = TAB$Y, 
                        Tr = TAB$X, 
                        X = TAB$Score, 
                        estimand = "ATE",
                        M = 1,
                        caliper = 0.20,
                        replace = TRUE,
                        ties = FALSE)

Match_lm <- Appariement_lm$est

#######################
#Estimates based on RF#
#######################

Mod_Yield <- Yield_GBM
Score <- Score_GBM

#Inverse probability weighting#
TAB <- data.frame(X=df$TREAT, Y=df$yield, Score)
#Remove extreme scores
TAB <- TAB[(Score > 0.1 & Score < 0.9),]

E1 <- mean(TAB$X * TAB$Y / TAB$Score)
E0 <- mean((1 - TAB$X) * TAB$Y / (1 - TAB$Score))

AVG_ipw_gbm <- E1 - E0
AVG_ipw_gbm

#Standardization

TABPred <- df
TABPred <- TABPred[(Score > 0.1 & Score < 0.9),]

TABPred$TREAT <- rep(1, nrow(TABPred))
Pred_4 <- predict(Mod_Yield, newdata = TABPred)

TABPred$TREAT <- rep(0, nrow(TABPred))
Pred_0 <- predict(Mod_Yield, newdata=TABPred)

AVG_sdz_gbm <- mean(Pred_4) - mean(Pred_0)
AVG_sdz_gbm

#Double robust

DR1 <- mean(Pred_4 + (TAB$X / TAB$Score) * (TAB$Y - Pred_4))
DR0 <- mean(Pred_0 + ((1 - TAB$X) / (1 - TAB$Score)) * (TAB$Y - Pred_0))

AVG_dr_gbm <- DR1 - DR0
AVG_dr_gbm

#Matching

Appariement_gbm <- Match(Y = TAB$Y, 
                         Tr = TAB$X, 
                         X = TAB$Score, 
                         estimand = "ATE",
                         M = 1,
                         caliper = 0.20,
                         replace = TRUE,
                         ties = FALSE)

Match_gbm <- Appariement_gbm$est

###############
##All results##
###############

Results_all <- data.frame(naive=AVG_naive, ipw_lm=AVG_ipw_lm, sdz_lm=AVG_sdz_lm,
                          dr_lm=AVG_dr_lm, app_lm=Match_lm,
                          ipw_gbm=AVG_ipw_gbm, sdz_gbm=AVG_sdz_gbm, dr_gbm=AVG_dr_gbm, 
                          app_gbm=Match_gbm)

write.table(Results_all, file="results_DR_yield/Results_est_vines_SP_treat2_apr.txt")

# Combine SEA results with Results_all
SEA_results <- composite_time_series %>%
  filter(Centered_Year == 0) %>%
  summarize(Global_Mean_Yield = mean(Global_Mean_Yield, na.rm = TRUE))

Results_all$SEA <- SEA_results$Global_Mean_Yield

write.table(Results_all, file="results_DR_yield/Results_est_vines_SP_treat2_apr.txt")


#### SCRIPT 3 - BOOTSTRAPPING ANALYSIS ####

library(Matching)
library(gbm)
library(lme4)
library(dplyr)
library(tidyr)

######
#Data#
######
set.seed(3)
df <- readRDS("output_files/vines_combined_2000_2023_treat_SP.rda")
summary(df)
#Choice of treatment
df$TREAT <- as.numeric(as.character(df$treat2_apr))
#summary(df)
sum(df$TREAT[df$TREAT == 1])/nrow(df)

# Retain original yield for SEA analysis
df$original_yield <- df$yield

#########
#Scaling#
#########

df$MM_ST_4 <- scale(df$MM_ST_4)
df$MM_RG_4 <- scale(df$MM_RG_4)
df$MM_ETP_4 <- scale(df$MM_ETP_4)
df$MM_SWI_4 <- scale(df$MM_SWI_4)
df$MM_NJT30_4 <- scale(df$MM_NJT30_4)
df$yield <- log(df$yield)

##########
#Sampling#
##########
naive <- c()

ipw_lm <- c()
sdz_lm <- c()
dr_lm <- c()
app_lm <- c()

ipw_gbm <- c()
sdz_gbm <- c()
dr_gbm <- c()
app_gbm <- c()

sea_results <- c()

###############################
##Loop over bootstrap samples##
###############################

for (k in 1:200) {
  
  List_DPT <- unique(df$Departement)
  Tab_boot <- c()
  
  #Loop over departments for sampling by dpt.
  for (i in 1:length(unique(df$Departement))) {
    DPT_i <- List_DPT[i]
    df_i <- df[df$Departement == List_DPT[i],]
    Num_row_i <- 1:nrow(df_i)
    Rows_select_i <- sample(Num_row_i, replace = TRUE)
    boot_i <- df_i[Rows_select_i,]
    Tab_boot <- rbind(Tab_boot, boot_i)
  }
  
  ################
  #Naive estimate#
  ################
  
  AVG_naive <- mean(Tab_boot$yield[Tab_boot$TREAT == 1]) - mean(Tab_boot$yield[Tab_boot$TREAT == 0])
  naive <- c(naive, AVG_naive)
  
  ############
  #Models lme#
  ############
  
  Proba_glmer <- glmer(cbind(TREAT, 1 - TREAT) ~ MM_ST_4  + MM_RG_4 + MM_ETP_4 + MM_SWI_4 + MM_NJT30_4 +  (1 | Departement), 
                       family = binomial, data = Tab_boot)
  Yield_lmer <- lmer(yield ~ TREAT + MM_ST_4  + MM_RG_4 + MM_ETP_4 + MM_SWI_4 + MM_NJT30_4 +  (1 | Departement), data = Tab_boot)
  Score_glmer <- predict(Proba_glmer, type = "response")
  
  ###########
  #Models gbm#
  ###########
  
  Proba_GBM <- gbm(TREAT ~ MM_ST_4 + MM_RG_4 + MM_ETP_4 + MM_SWI_4 + MM_NJT30_4 +  Departement, data = Tab_boot, distribution = "bernoulli")
  Yield_GBM <- gbm(yield ~ TREAT + MM_ST_4 + MM_RG_4 + MM_ETP_4 + MM_SWI_4 + MM_NJT30_4 +  Departement, data = Tab_boot)
  Score_GBM <- predict(Proba_GBM, type = "response")
  
  #######################
  #Estimate based on glm#
  #######################
  
  Mod_Yield <- Yield_lmer
  Score <- Score_glmer
  
  #Inverse probability weighting#
  TAB <- data.frame(X = Tab_boot$TREAT, Y = Tab_boot$yield, Score)
  #Remove extreme scores
  TAB <- TAB[(Score > 0.1 & Score < 0.9),]
  
  E1 <- mean(TAB$X * TAB$Y / TAB$Score)
  E0 <- mean((1 - TAB$X) * TAB$Y / (1 - TAB$Score))
  
  AVG_ipw_lm <- E1 - E0
  ipw_lm <- c(ipw_lm, AVG_ipw_lm)
  
  #Standardization
  TABPred <- Tab_boot
  TABPred <- TABPred[(Score > 0.1 & Score < 0.9),]
  
  TABPred$TREAT <- rep(1, nrow(TABPred))
  Pred_4 <- predict(Mod_Yield, newdata = TABPred)
  
  TABPred$TREAT <- rep(0, nrow(TABPred))
  Pred_0 <- predict(Mod_Yield, newdata = TABPred)
  
  AVG_sdz_lm <- mean(Pred_4) - mean(Pred_0)
  sdz_lm <- c(sdz_lm, AVG_sdz_lm)
  
  #Double robust
  DR1 <- mean(Pred_4 + (TAB$X / TAB$Score) * (TAB$Y - Pred_4))
  DR0 <- mean(Pred_0 + ((1 - TAB$X) / (1 - TAB$Score)) * (TAB$Y - Pred_0))
  
  AVG_dr_lm <- DR1 - DR0
  dr_lm <- c(dr_lm, AVG_dr_lm)
  
  #Matching
  Appariement_lm <- Match(Y = TAB$Y, 
                          Tr = TAB$X, 
                          X = TAB$Score, 
                          estimand = "ATE",
                          M = 1,
                          caliper = 0.20,
                          replace = TRUE,
                          ties = FALSE)
  app_lm <- c(app_lm, Appariement_lm$est)
  
  #######################
  #Estimates based on GBM#
  #######################
  
  Mod_Yield <- Yield_GBM
  Score <- Score_GBM
  
  #Inverse probability weighting#
  TAB <- data.frame(X = Tab_boot$TREAT, Y = Tab_boot$yield, Score)
  #Remove extreme scores
  TAB <- TAB[(Score > 0.1 & Score < 0.9),]
  
  E1 <- mean(TAB$X * TAB$Y / TAB$Score)
  E0 <- mean((1 - TAB$X) * TAB$Y / (1 - TAB$Score))
  
  AVG_ipw_gbm <- E1 - E0
  ipw_gbm <- c(ipw_gbm, AVG_ipw_gbm)
  
  #Standardization
  TABPred <- Tab_boot
  TABPred <- TABPred[(Score > 0.1 & Score < 0.9),]
  
  TABPred$TREAT <- rep(1, nrow(TABPred))
  Pred_4 <- predict(Mod_Yield, newdata = TABPred)
  
  TABPred$TREAT <- rep(0, nrow(TABPred))
  Pred_0 <- predict(Mod_Yield, newdata = TABPred)
  
  AVG_sdz_gbm <- mean(Pred_4) - mean(Pred_0)
  sdz_gbm <- c(sdz_gbm, AVG_sdz_gbm)
  
  #Double robust
  DR1 <- mean(Pred_4 + (TAB$X / TAB$Score) * (TAB$Y - Pred_4))
  DR0 <- mean(Pred_0 + ((1 - TAB$X) / (1 - TAB$Score)) * (TAB$Y - Pred_0))
  
  AVG_dr_gbm <- DR1 - DR0
  dr_gbm <- c(dr_gbm, AVG_dr_gbm)
  
  #Matching
  Appariement_gbm <- Match(Y = TAB$Y, 
                           Tr = TAB$X, 
                           X = TAB$Score, 
                           estimand = "ATE",
                           M = 1,
                           caliper = 0.20,
                           replace = TRUE,
                           ties = FALSE)
  app_gbm <- c(app_gbm, Appariement_gbm$est)
  
  #######################
  # SEA Analysis
  #######################
  sp <- Tab_boot %>%
    filter(TREAT == 1, YEAR >= 2002, YEAR <= 2021) %>%
    rename(Department = Departement, sp_Year = YEAR) %>%
    mutate(Occurrence_id = paste(Department, sp_Year, mar = "_"))  # Assign a unique ID to each sp occurrence
  
  sp_expanded <- sp %>%
    group_by(Department) %>%
    arrange(sp_Year) %>%
    mutate(Cooccurring = map(sp_Year, ~ seq(.x - 2, .x + 2))) %>%
    unnest(cols = c(Cooccurring))
  
  cooccurring_sp <- sp_expanded %>%
    filter(Cooccurring %in% sp_Year & Cooccurring != sp_Year) %>%
    mutate(Excluded_Year = Cooccurring)
  
  cooccurring_sp <- cooccurring_sp[, c("Department", "sp_Year", "Excluded_Year", "Occurrence_id")]
  
  df_renamed <- Tab_boot %>%
    rename(Department = Departement, Yield = original_yield)
  
  sp_yields <- sp %>%
    inner_join(df_renamed, by = c("Department", "sp_Year" = "YEAR"))
  
  yield_data <- df_renamed %>%
    filter(Department %in% unique(sp_yields$Department))
  
  yield_window <- sp_yields %>%
    rowwise() %>%
    do({
      data.frame(Occurrence_id = .$Occurrence_id,
                 Event_year = .$sp_Year,
                 Department = .$Department,
                 Centered_Year = -2:2,
                 Year_yield = .$sp_Year + -2:2)
    }) %>%
    ungroup()
  
  yield_window <- yield_window %>%
    left_join(yield_data, by = c("Department", "Year_yield" = "YEAR"))
  
  yield_window <- yield_window %>%
    left_join(cooccurring_sp, by = c("Occurrence_id", "Department")) %>%
    mutate(Excluded_sp_Year = ifelse(!is.na(Excluded_Year) & Year_yield == Excluded_Year, "Yes", "No"))
  
  yield_window <- yield_window %>%
    mutate(Exclude = Excluded_sp_Year == "Yes")
  
  baseline <- yield_window %>%
    filter(Centered_Year %in% c(-2, -1, 1, 2) & !Exclude) %>%
    group_by(Occurrence_id, Department, Event_year) %>%
    summarise(Baseline_mean = mean(Yield, na.rm = TRUE), .groups = 'drop')
  
  normalized_yields <- yield_window %>%
    left_join(baseline, by = c("Occurrence_id", "Department", "Event_year")) %>%
    mutate(Normalized_Yield = ifelse(!is.na(Baseline_mean) & !is.na(Yield), Yield / Baseline_mean, NA))
  
  normalized_yields <- normalized_yields %>%
    filter(!(Exclude & Centered_Year != 0))
  
  composite_time_series <- normalized_yields %>%
    group_by(Centered_Year) %>%
    summarise(Global_Mean_Yield = mean(Normalized_Yield, na.rm = TRUE), .groups = 'drop')
  
  SEA_result <- composite_time_series %>%
    filter(Centered_Year == 0) %>%
    summarise(Global_Mean_Yield = mean(Global_Mean_Yield, na.rm = TRUE)) %>%
    pull(Global_Mean_Yield)
  
  SEA_effect_size <- SEA_result - 1
  sea_results <- c(sea_results, SEA_effect_size)
  
  print(k)
}

####Table of all results ####
MAT_result <- data.frame(naive, ipw_lm, sdz_lm, dr_lm, app_lm, ipw_gbm,
                         sdz_gbm, dr_gbm, app_gbm, sea_results)
apply(MAT_result, 2, quantile, probs = c(0.025, 0.975))
write.table(MAT_result, file = "results_DR_yield/Results_boot_vines_SP_treat2_apr.txt")



########################
## SCRIPT 4 - Dotchart##
########################

# Read point estimates and bootstrapped CIs
Point_est <- read.table("results_DR_yield/Results_est_vines_SP_treat2_apr.txt", header = TRUE)
CI_boot_tab <- read.table("results_DR_yield/Results_boot_vines_SP_treat2_apr.txt", header = TRUE)
CI <- apply(CI_boot_tab, 2, quantile, probs = c(0.025, 0.975))

# Calculate SEA effect size
SEA_effect_size <- Point_est$SEA - 1

# Add SEA effect size to Point_est as a new column
Point_est$SEA_effect_size <- SEA_effect_size

# Add SEA CI to CI
SEA_CI <- quantile(CI_boot_tab$sea_results, probs = c(0.025, 0.975))
CI <- cbind(CI[, -which(colnames(CI) == "sea_results")], SEA_effect_size = SEA_CI)

# Ensure labels and Vec are in the correct order
labels <- c("naive", "SEA_effect_size", "ipw_lm", "sdz_lm", "dr_lm", 
            "app_lm", "ipw_gbm", "sdz_gbm", "dr_gbm", 
            "app_gbm")
Vec <- Point_est[, labels]

# Rename "SEA_effect_size" to "SEA"
labels <- sub("SEA_effect_size", "SEA", labels)

# Create TAB_log as a data frame
TAB_log <- data.frame(
  Method = labels,
  Estimate = as.numeric(Point_est[, gsub("SEA", "SEA_effect_size", labels)]),
  LB = as.numeric(CI[1, gsub("SEA", "SEA_effect_size", labels)]),
  UB = as.numeric(CI[2, gsub("SEA", "SEA_effect_size", labels)])
)

# Convert TAB_log values to percentages for non-SEA methods
TAB_percentage <- TAB_log
TAB_percentage[TAB_percentage$Method != "SEA", "Estimate"] <- 100 * (exp(TAB_log[TAB_log$Method != "SEA", "Estimate"]) - 1)
TAB_percentage[TAB_percentage$Method != "SEA", "LB"] <- 100 * (exp(TAB_log[TAB_log$Method != "SEA", "LB"]) - 1)
TAB_percentage[TAB_percentage$Method != "SEA", "UB"] <- 100 * (exp(TAB_log[TAB_log$Method != "SEA", "UB"]) - 1)

# Convert SEA values to percentages directly
TAB_percentage[TAB_percentage$Method == "SEA", c("Estimate", "LB", "UB")] <- 100 * TAB_log[TAB_log$Method == "SEA", c("Estimate", "LB", "UB")]

# Plot the dotchart for TAB_log
par(mfrow = c(1, 1))
dotchart(rev(TAB_log$Estimate), labels = rev(TAB_log$Method), pch = 19, xlim = range(TAB_log$LB, TAB_log$UB), xlab = "Estimated relative effect (log scale)")
abline(v = 0, col = "blue", lty = 2)
for (i in 1:nrow(TAB_log)) {
  lines(c(TAB_log$LB[i], TAB_log$UB[i]), c(nrow(TAB_log)-i+1, nrow(TAB_log)-i+1))
}

# Plot the dotchart for TAB_percentage
dotchart(rev(TAB_percentage$Estimate), labels = rev(TAB_percentage$Method), pch = 19, xlim = c(-30, 30), xlab = "Estimated relative effect (%)")
abline(v = 0, col = "blue", lty = 2)
for (i in 1:nrow(TAB_percentage)) {
  lines(c(TAB_percentage$LB[i], TAB_percentage$UB[i]), c(nrow(TAB_percentage)-i+1, nrow(TAB_percentage)-i+1))
}


# Print dr_lm values from TAB_percentage
dr_lm_value <- round(TAB_percentage[TAB_percentage$Method == "dr_lm", "Estimate"], 2)
dr_lm_LB <- round(TAB_percentage[TAB_percentage$Method == "dr_lm", "LB"], 2)
dr_lm_UB <- round(TAB_percentage[TAB_percentage$Method == "dr_lm", "UB"], 2)

formatted_dr_lm <- paste0(dr_lm_value, "% [", dr_lm_LB, "%, ", dr_lm_UB, "%]")
print(paste("dr_lm result:", formatted_dr_lm))

# Print dr_gbm values from TAB_percentage
dr_gbm_value <- round(TAB_percentage[TAB_percentage$Method == "dr_gbm", "Estimate"], 2)
dr_gbm_LB <- round(TAB_percentage[TAB_percentage$Method == "dr_gbm", "LB"], 2)
dr_gbm_UB <- round(TAB_percentage[TAB_percentage$Method == "dr_gbm", "UB"], 2)

formatted_dr_gbm <- paste0(dr_gbm_value, "% [", dr_gbm_LB, "%, ", dr_gbm_UB, "%]")
print(paste("dr_gbm result:", formatted_dr_gbm))
