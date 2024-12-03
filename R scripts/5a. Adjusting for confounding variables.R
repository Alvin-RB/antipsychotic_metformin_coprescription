# Metformin and SGA co-prescrption #####
# Last run: 02/12/2024

# Confounder adjustment

# Clear memory
#rm(list = ls())

library(tidyverse)    
library(mice)         
library(gtsummary)

# Load data
load("tte_imputed_long.rdata")

# Imputation for missing values
imputed_data <- tte_imputed_long

fields <- TTEcohort_final %>%
  select(patid, diff, metformin_flag, pcos_flag) %>%
  distinct()

imputed_data <- imputed_data %>%
  left_join(fields, by = "patid")

metformin <- imputed_data %>%
  mutate(keep = case_when(is.na(diff) | diff >= -30 & diff <= 90 ~ 1,
                          TRUE ~ 0)) %>%
  filter(keep == 1) %>%
  select(-keep)

# Grouping together mixed and other ethnicities in imputed_data
metformin <- metformin %>%
  mutate(ethnicity_cat_cprdhes = case_when(ethnicity_cat_cprdhes == "mixed" ~ "mixed/other",
                                           ethnicity_cat_cprdhes == "other" ~ "mixed/other",
                                           TRUE ~ ethnicity_cat_cprdhes))

imputed_data$ethnicity_cat_cprdhes <- factor(imputed_data$ethnicity_cat_cprdhes, 
                                             levels = c("white", "black", "asian", "mixed/other"))

#Female sub-group for analysis with multiple covariates 
metformin_female <- metformin %>%
  filter(gender=="Female") %>%
  mice::as.mids()

# Model for 2-year weight (Females)
model2_2y_female <- with(metformin_female, 
                        lm(twoy_weightkg ~ metformin_flag + baseline_weightkg + I(baseline_weightkg^2) + 
                             trt_group + age_atprevcohortentry + I(age_atprevcohortentry^2) + 
                             ethnicity_cat_cprdhes + pat_2019imd_quintile + priordiabetes + pcos_flag))

#summary(model_2y_female)
tbl_regression(model2_2y_female, estimate_fun = function(x) sprintf("%.2f", x))
pooled_model2_2y_female <- pool(model2_2y_female)

# Model for 1-year weight (Females)
model2_1y_female <- with(metformin_female, 
                        lm(oney_weightkg ~ metformin_flag + baseline_weightkg + I(baseline_weightkg^2) + 
                             trt_group + age_atprevcohortentry + I(age_atprevcohortentry^2) + 
                             ethnicity_cat_cprdhes + pat_2019imd_quintile + priordiabetes + pcos_flag))

tbl_regression(model2_1y_female, estimate_fun = function(x) sprintf("%.2f", x))
pooled_model2_1y_female <- pool(model2_1y_female)

# Model for 6m weight (Females)
model2_6m_female <- with(metformin_female, 
                        lm(sixm_weightkg ~ metformin_flag + baseline_weightkg + I(baseline_weightkg^2) + 
                             trt_group + age_atprevcohortentry + I(age_atprevcohortentry^2) + 
                             ethnicity_cat_cprdhes + pat_2019imd_quintile + priordiabetes + pcos_flag))

tbl_regression(model2_6m_female, estimate_fun = function(x) sprintf("%.2f", x))
pooled_model2_6m_female <- pool(model2_6m_female)

# Generate forest plots for Female Models
# Extract coefficients and CI for each model
summary_model2_6m_female <- summary(pooled_model2_6m_female, conf.int = TRUE)
summary_model2_1y_female <- summary(pooled_model2_1y_female, conf.int = TRUE)
summary_model2_2y_female <- summary(pooled_model2_2y_female, conf.int = TRUE)

# Create data frames for the forest plot (Female)
model2_coefficients_df_female <- data.frame(
  Time = c("6m", "1y", "2y"),
  Estimate = c(summary_model2_6m_female$estimate[2], summary_model2_1y_female$estimate[2], summary_model2_2y_female$estimate[2]),
  Lower_CI = c(summary_model2_6m_female$`2.5 %`[2], summary_model2_1y_female$`2.5 %`[2], summary_model2_2y_female$`2.5 %`[2]),
  Upper_CI = c(summary_model2_6m_female$`97.5 %`[2], summary_model2_1y_female$`97.5 %`[2], summary_model2_2y_female$`97.5 %`[2]))

rm(summary_model2_6m_female, summary_model2_1y_female, summary_model2_2y_female)

# Modify the Time variable to be a factor with the desired order
model2_coefficients_df_female$Time <- factor(model2_coefficients_df_female$Time, levels = c("2y", "1y", "6m"))

# Create the forest plot (Female)
model2_forest_plot_female <- ggplot(model2_coefficients_df_female, aes(x = Time, y = Estimate)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.2, size = 1) +
  coord_flip() +  # Flip the axes for a horizontal forest plot
  labs(title = "Effect of Metformin on Weight Change Over Time,Adjusting for multiple confounding variables (Females)",
       x = "Time",
       y = "Change in weight in kg (95% CI)") +
  theme_minimal() +
  theme(text = element_text(family = "Times"))

# Display the plot for Females
print(model2_forest_plot_female)

#Female sub-group for analysis with baseline weight covariate

# Model for 2-year weight (Females)
model1_2y_female <- with(metformin_female, 
                        lm(twoy_weightkg ~ metformin_flag + baseline_weightkg + I(baseline_weightkg^2)))

#summary(model_2y_female)
tbl_regression(model1_2y_female, estimate_fun = function(x) sprintf("%.2f", x))
pooled_model1_2y_female <- pool(model1_2y_female)

# Model for 1-year weight (Females)
model1_1y_female <- with(metformin_female, 
                        lm(oney_weightkg ~ metformin_flag + baseline_weightkg + I(baseline_weightkg^2)))

tbl_regression(model1_1y_female, estimate_fun = function(x) sprintf("%.2f", x))
pooled_model1_1y_female <- pool(model1_1y_female)

# Model for 6m weight (Females)
model1_6m_female <- with(metformin_female, 
                        lm(sixm_weightkg ~ metformin_flag + baseline_weightkg + I(baseline_weightkg^2)))

tbl_regression(model1_6m_female, estimate_fun = function(x) sprintf("%.2f", x))
pooled_model1_6m_female <- pool(model1_6m_female)

# Generate forest plots for Female Models
# Extract coefficients and CI for each model (Female)
summary_model1_6m_female <- summary(pooled_model1_6m_female, conf.int = TRUE)
summary_model1_1y_female <- summary(pooled_model1_1y_female, conf.int = TRUE)
summary_model1_2y_female <- summary(pooled_model1_2y_female, conf.int = TRUE)

# Create data frames for the forest plot (Female)
model1_coefficients_df_female <- data.frame(
  Time = c("6m", "1y", "2y"),
  Estimate = c(summary_model1_6m_female$estimate[2], summary_model1_1y_female$estimate[2], summary_model1_2y_female$estimate[2]),
  Lower_CI = c(summary_model1_6m_female$`2.5 %`[2], summary_model1_1y_female$`2.5 %`[2], summary_model1_2y_female$`2.5 %`[2]),
  Upper_CI = c(summary_model1_6m_female$`97.5 %`[2], summary_model1_1y_female$`97.5 %`[2], summary_model1_2y_female$`97.5 %`[2]))

rm(summary_model1_6m_female, summary_model1_1y_female, summary_model1_2y_female)

# Modify the Time variable to be a factor with the desired order
model1_coefficients_df_female$Time <- factor(model1_coefficients_df_female$Time, levels = c("2y", "1y", "6m"))

# Create the forest plot (Female)
model1_forest_plot_female <- ggplot(model1_coefficients_df_female, aes(x = Time, y = Estimate)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.2, size = 1) +
  coord_flip() +  # Flip the axes for a horizontal forest plot
  labs(title = "Effect of Metformin on Weight Change Over Time,Adjusting for baseline weight (Females)",
       x = "Time",
       y = "Change in weight in kg (95% CI)") +
  theme_minimal() +
  theme(text = element_text(family = "Times"))

# Display the plot for Females
print(model1_forest_plot_female)

#Male sub-group for analysis for multiple covariates
metformin_male <- metformin %>%
  filter(gender=="Male") %>%
  mice::as.mids()

# Model for 2-year weight (Males)
model2_2y_male <- with(metformin_male, 
                      lm(twoy_weightkg ~ metformin_flag + baseline_weightkg + I(baseline_weightkg^2) + 
                           trt_group + age_atprevcohortentry + I(age_atprevcohortentry^2) + 
                           ethnicity_cat_cprdhes + pat_2019imd_quintile + priordiabetes))
#summary(model_2y_male)
pooled_model2_2y_male <- pool(model2_2y_male)
tbl_regression(model2_2y_male, estimate_fun = function(x) sprintf("%.2f", x))

# Model for 1-year weight (Males)
model2_1y_male <- with(metformin_male, 
                      lm(oney_weightkg ~ metformin_flag + baseline_weightkg + I(baseline_weightkg^2) + 
                           trt_group + age_atprevcohortentry + I(age_atprevcohortentry^2) + 
                           ethnicity_cat_cprdhes + pat_2019imd_quintile + priordiabetes))

#summary(model_1y_male)
pooled_model2_1y_male <- pool(model2_1y_male)
tbl_regression(model2_1y_male, estimate_fun = function(x) sprintf("%.2f", x))

# Model for 6-month weight (Males)
model2_6m_male <- with(metformin_male, 
                      lm(sixm_weightkg ~ metformin_flag + baseline_weightkg + I(baseline_weightkg^2) + 
                           trt_group + age_atprevcohortentry + I(age_atprevcohortentry^2) + 
                           ethnicity_cat_cprdhes + pat_2019imd_quintile + priordiabetes))

#summary(model_6m_male)
pooled_model2_6m_male <- pool(model2_6m_male)
tbl_regression(model2_6m_male, estimate_fun = function(x) sprintf("%.2f", x))

# Extract coefficients and CI for each model (Male)
summary_model2_6m_male <- summary(pooled_model2_6m_male, conf.int = TRUE)
summary_model2_1y_male <- summary(pooled_model2_1y_male, conf.int = TRUE)
summary_model2_2y_male <- summary(pooled_model2_2y_male, conf.int = TRUE)

# Create data frame for the forest plot (Male)
model2_coefficients_df_male <- data.frame(
  Time = c("6m", "1y", "2y"),
  Estimate = c(summary_model2_6m_male$estimate[2], summary_model2_1y_male$estimate[2], summary_model2_2y_male$estimate[2]),
  Lower_CI = c(summary_model2_6m_male$`2.5 %`[2], summary_model2_1y_male$`2.5 %`[2], summary_model2_2y_male$`2.5 %`[2]),
  Upper_CI = c(summary_model2_6m_male$`97.5 %`[2], summary_model2_1y_male$`97.5 %`[2], summary_model2_2y_male$`97.5 %`[2]))

rm(summary_model2_6m_male, summary_model2_1y_male, summary_model2_2y_male)

# Modify the Time variable to be a factor with the desired order
model2_coefficients_df_male$Time <- factor(model2_coefficients_df_male$Time, levels = c( "2y","1y","6m"))

# Create the forest plot (Male)
model2_forest_plot_male <- ggplot(model2_coefficients_df_male, aes(x = Time, y = Estimate)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.2, size = 1) +
  coord_flip() +  # Flip the axes for a horizontal forest plot
  labs(title = "Effect of Metformin on Weight Change Over Time, Adjusting for multiple confounding variables (Males)",
       x = "Time",
       y = "Change in weight in kg (95% CI)") +
  theme_minimal() +
  theme(text = element_text(family = "Times"))

# Display the plot for Males
print(model2_forest_plot_male)

#Male sub-group for analysis for baseline weight covariate. 

# Model for 2-year weight (Males)
model1_2y_male <- with(metformin_male, 
                      lm(twoy_weightkg ~ metformin_flag + baseline_weightkg + I(baseline_weightkg^2)))

#summary(model_2y_male)
pooled_model1_2y_male <- pool(model1_2y_male)
tbl_regression(model1_2y_male, estimate_fun = function(x) sprintf("%.2f", x))

# Model for 1-year weight (Males)
model1_1y_male <- with(metformin_male, 
                      lm(oney_weightkg ~ metformin_flag + baseline_weightkg + I(baseline_weightkg^2)))

#summary(model_1y_male)
pooled_model1_1y_male <- pool(model1_1y_male)
tbl_regression(model1_1y_male, estimate_fun = function(x) sprintf("%.2f", x))

# Model for 6-month weight (Males)
model1_6m_male <- with(metformin_male, 
                      lm(sixm_weightkg ~ metformin_flag + baseline_weightkg + I(baseline_weightkg^2)))

#summary(model_6m_male)
pooled_model1_6m_male <- pool(model1_6m_male)
tbl_regression(model1_6m_male, estimate_fun = function(x) sprintf("%.2f", x))

# Extract coefficients and CI for each model (Male)
summary_model1_6m_male <- summary(pooled_model1_6m_male, conf.int = TRUE)
summary_model1_1y_male <- summary(pooled_model1_1y_male, conf.int = TRUE)
summary_model1_2y_male <- summary(pooled_model1_2y_male, conf.int = TRUE)

# Create data frame for the forest plot (Male)
model1_coefficients_df_male <- data.frame(
  Time = c("6m", "1y", "2y"),
  Estimate = c(summary_model1_6m_male$estimate[2], summary_model1_1y_male$estimate[2], summary_model1_2y_male$estimate[2]),
  Lower_CI = c(summary_model1_6m_male$`2.5 %`[2], summary_model1_1y_male$`2.5 %`[2], summary_model1_2y_male$`2.5 %`[2]),
  Upper_CI = c(summary_model1_6m_male$`97.5 %`[2], summary_model1_1y_male$`97.5 %`[2], summary_model1_2y_male$`97.5 %`[2]))

rm(summary_model1_6m_male, summary_model1_1y_male, summary_model1_2y_male)

# Modify the Time variable to be a factor with the desired order
model1_coefficients_df_male$Time <- factor(model1_coefficients_df_male$Time, levels = c( "2y","1y","6m"))

# Create the forest plot (Male)
model1_forest_plot_male <- ggplot(model1_coefficients_df_male, aes(x = Time, y = Estimate)) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = Lower_CI, ymax = Upper_CI), width = 0.2, size = 1) +
  coord_flip() +  # Flip the axes for a horizontal forest plot
  labs(title = "Effect of Metformin on Weight Change Over Time, Adjusting for baseline weight (Males)",
       x = "Time",
       y = "Change in weight in kg (95% CI)") +
  theme_minimal() +
  theme(text = element_text(family = "Times"))

# Display the plot for Males
print(model1_forest_plot_male)

# Combine all coefficient data frames into one for unified plotting
combined_coefficients <- data.frame(
  Time = rep(c("6 months", "1 year", "2 years"), 4),
  Estimate = c(
    model2_coefficients_df_female$Estimate,  # Female (Multiple Covariates)
    model1_coefficients_df_female$Estimate,  # Female (Baseline Weight)
    model2_coefficients_df_male$Estimate,    # Male (Multiple Covariates)
    model1_coefficients_df_male$Estimate     # Male (Baseline Weight)
  ),
  Lower_CI = c(
    model2_coefficients_df_female$Lower_CI,
    model1_coefficients_df_female$Lower_CI,
    model2_coefficients_df_male$Lower_CI,
    model1_coefficients_df_male$Lower_CI),
  Upper_CI = c(
    model2_coefficients_df_female$Upper_CI,
    model1_coefficients_df_female$Upper_CI,
    model2_coefficients_df_male$Upper_CI,
    model1_coefficients_df_male$Upper_CI),
  Group = rep(
    c("Female (Multiple Covariates)", "Female (Baseline Weight)",
      "Male (Multiple Covariates)", "Male (Baseline Weight)"),
    each = 3))

# Format the Time variable for desired y-axis order
combined_coefficients$Time <- factor(
  combined_coefficients$Time,
  levels = c("2 years", "1 year", "6 months")  # Ensure correct y-axis order
)


# Update Group labels to "Model 1" and "Model 2" and reorder them
combined_coefficients$Group <- factor(
  combined_coefficients$Group,
  levels = c(
    "Female (Multiple Covariates)", 
    "Male (Multiple Covariates)", 
    "Female (Baseline Weight)", 
    "Male (Baseline Weight)"
  ),
  labels = c(
    "Female (Model 2)", 
    "Male (Model 2)", 
    "Female (Model 1)", 
    "Male (Model 1)"))

# Unified Forest Plot with Updated Aesthetics
unified_forest_plot <- ggplot(combined_coefficients, aes(x = Time, y = Estimate, color = Group)) +
  # Assign shapes and colours
  geom_point(
    aes(shape = Group),
    size = 4,
    position = position_dodge(width = 0.7)
  ) +
  geom_errorbar(
    aes(ymin = Lower_CI, ymax = Upper_CI),
    width = 0.2,
    size = 0.8,
    position = position_dodge(width = 0.7)
  ) +
  coord_flip() +  # Flip axes for horizontal plot
  labs(
    x = "Time-point",                    # Updated Y-axis title
    y = "Change in Weight (kg, 95% CI)"
  ) +
  theme_minimal(base_family = "Times") +  # Times font
  # Custom colour-blind-friendly colours: orange and blue
  scale_color_manual(values = c(
    "Female (Model 1)" = "#E69F00",      # Orange
    "Female (Model 2)" = "#0072B2",      # Blue
    "Male (Model 1)" = "#E69F00",        # Orange
    "Male (Model 2)" = "#0072B2"         # Blue
  )) +
  # Custom shapes: triangles for females, circles for males
  scale_shape_manual(values = c(
    "Female (Model 1)" = 17,  # Triangle
    "Female (Model 2)" = 17,  # Triangle
    "Male (Model 1)" = 16,    # Circle
    "Male (Model 2)" = 16     # Circle
  )) +
  theme(
    legend.position = "bottom",          # Legend at the bottom
    legend.title = element_blank(),      # Remove legend title
    panel.grid = element_blank(),        # Remove grid lines
    axis.line = element_line(size = 0.5) # Add axis lines for clarity
  )

# Display the updated plot
print(unified_forest_plot)

# Save the updated plot
ggsave(
  filename = "unified_forest_plot_updated.pdf",  # Name of the output file
  plot = unified_forest_plot,                    # The plot object
  device = "pdf",                                # File format
  width = 10,                                    # Width of the PDF in inches
  height = 8                                     # Height of the PDF in inches
)

# Package versions used for this script
# dplyr version: 1.1.4 
# mice version: 3.16.0 
# gtsummary version: 2.0.2 
