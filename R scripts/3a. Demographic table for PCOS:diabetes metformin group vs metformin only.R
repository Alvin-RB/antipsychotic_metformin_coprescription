# Metformin and SGA co-prescrption #####
# Last run: 02/12/2024

# Demographics of indication vs. no indication

# Load necessary libraries
library(dplyr)
library(gtsummary)
library(openxlsx)

# Define Group 3 and Group 4 based on the updated criteria
TTEcohort_indication <- TTEcohort_final %>%
  mutate(Group = case_when(metformin_flag == 1 & group == "Study cohort" & !(pcos_flag == 1 | priordiabetes == 1) ~ "Group 3",
                           metformin_flag == 1 & group == "Study cohort" & (pcos_flag == 1 | priordiabetes == 1) ~ "Group 4",
                           TRUE ~ NA_character_)) %>%
  filter(!is.na(Group))

#Create table
table_summary2 <- TTEcohort_indication %>%
  select(
    Group, gender, age_atfirstap, age_atprevcohortentry, age_atfirstdiag, ethnicity_cat_cprdhes,
    last_smi_diag, region, patprac_2019imd_quintile, baseline_bmi_cat,
    priorcerebrovasculardisease, priormyocardialinfarction, priorliverdisease, priorrenaldisease,
    priorhypertension, priordyslipidaemia, priordiabetes, pcos_flag,
    antidepressant_prior2years, lipiddrugs_prior2years, insulin, prioralcoholabuse, priorsubstanceabuse,
    smoking_status_cat, baseline_hba1c, baseline_glucose
  ) %>%
  tbl_summary(
    by = Group,
    missing = "ifany",  # Include missing values without affecting percentages
    type = list(
      all_categorical() ~ "categorical",  # Specify categorical variables
      all_continuous() ~ "continuous"     # Specify continuous variables
    ),
    statistic = list(
      all_categorical() ~ "{n} ({p}%)",     # Format as value (percentage)
      all_continuous() ~ "{mean} ({sd})"    # Format as mean (SD)
    ),
    digits = list(all_continuous() ~ 2, all_categorical() ~ c(0, 2))  # Set results to 2 decimal places
  ) %>%
  modify_header(label = "**Variable**")

# Manually order rows as specified
table_summary2 <- table_summary2 %>%
  modify_table_body(
    ~ .x %>% arrange(
      factor(label, levels = c(
        "Gender", "Mean Age (SD)", "Age at First AP Prescription", "Age At Cohort Entry", 
        "Age At First SMI Diagnosis", "Ethnicity, N (%)", "White", "Black", "Asian", 
        "Mixed/ Other", "SMI Diagnosis, N (%)", "Bipolar", "Schizophrenia", 
        "Other, non-organic Psychosis", "Geographical Region, N (%)", 
        "Socioeconomic Deprivation, N (%)", "Most Deprived", "BMI, N (%)", "Underweight", 
        "Healthy", "Obese", "Overweight", "Co-Morbidities, N (%)", 
        "Prior Cerebrovascular Disease", "Prior Myocardial Infarction", "Prior Liver Disease", 
        "Prior Renal Disease", "Prior Hypertension", "Prior Dyslipidaemia", "Prior Diabetes", 
        "Prior PCOS", "Other Medications Prescribed in the last 2 Years, N (%)", 
        "Prior Anti-Depressants", "Prior Lipid Drugs", "Prior Insulin", 
        "Previous Alcohol, Smoking or Substance Abuse, N (%)", "Prior Alcohol Abuse", 
        "Prior Substance Abuse", "Ex-Smoker", "Current Smoker", "Mean Biochemical Markers (SD)", 
        "Baseline HbA1c", "Baseline Glucose"))))

# Print the ordered table
print(table_summary2)

# Step 4: Convert the summary table to a tibble (data frame) and export to Excel
table_summary_df2 <- as_tibble(table_summary2)

# Additional variables
# Calculate mean and standard deviation of baseline_weightkg for each group
baseline_weight_stats <- TTEcohort_indication %>%
  group_by(Group) %>%  # Group by Group 3 and Group 4
  summarise(
    Mean_Baseline_Weight = mean(baseline_weightkg, na.rm = TRUE),
    SD_Baseline_Weight = sd(baseline_weightkg, na.rm = TRUE),
    .groups = "drop")  # Ensure the output is ungrouped

# Print the results
cat("Mean and Standard Deviation of Baseline Weight (kg) by Group:\n")
print(baseline_weight_stats)

# Packages verions used in this script
# dplyr version: 1.1.4 
# gtsummary version: 2.0.2 
# openxlsx version: 4.2.7.1 
