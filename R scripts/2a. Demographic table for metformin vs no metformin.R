# Metformin and SGA co-prescrption #####
# Last run: 02/12/2024

# Demographics of SGAs vs SGA_metformin groups

# Clear memory
rm(list = ls())

# Packages
library(dplyr)
library(gtsummary)
library(ggplot2)
library(tidyr)
library(lubridate)
library(tableone)
library(readr)
library(survival)
library(survminer)
library(gt)
library(mice)
library(tidyselect)
library(broom.mixed)

load("SMIantidiabetics_C.Rdata")
load("tte_prevalent_cohort.rdata")
load("SMIPCOS_C.Rdata")
load ("tte_imputed_long.rdata")
load("TTEcohortwithmetforminflag_andissuedate.RData")

#Table 1 - demographic characteristics of metformin naive vs metformin group

# Organise data
# Derive PCOS covariate
pcos <- TTEcohortwithmetforminflag_andissuedate %>%
  select(patid, cohortentrydate) %>%
  left_join(SMIPCOS_C, by = "patid") %>%
  filter(!is.na(eventdate) & eventdate <= cohortentrydate) %>%  # Ensure eventdate exists and is on/before cohortentrydate
  select(patid) %>%
  distinct() %>%
  mutate(pcos_flag = 1)

# Step 2: Add pcos_flag to TTEcohortwithmetforminflag_andissuedate
TTEcohortwithmetforminflag_andissuedate <- TTEcohortwithmetforminflag_andissuedate %>%
  left_join(pcos, by = "patid") %>%
  mutate(pcos_flag = replace_na(pcos_flag, 0))  # Set NA pcos_flag to 0

# Add Insulin
SMIantidiabetics_C_unique <- SMIantidiabetics_C %>%
  inner_join(TTEcohortwithmetforminflag_andissuedate %>% select(patid, cohortentrydate), by = "patid") %>%
  group_by(patid) %>%
  summarise(insulin = case_when(
    any(grepl("Insulin", Antidiabetic) & issuedate <= cohortentrydate) ~ 1,
    all(!grepl("Insulin", Antidiabetic)) ~ 0,
    TRUE ~ 0)) %>%
  ungroup()

TTEcohort_final <- TTEcohortwithmetforminflag_andissuedate %>%
  left_join(SMIantidiabetics_C_unique %>% select(patid, insulin), by = "patid") %>%
  mutate(insulin = replace_na(insulin, 0)) %>%  # Replace NA insulin with 0
  distinct(patid, .keep_all = TRUE)

# Creating a table
# Define Group 1 and Group 2 based on the updated criteria
group1 <- TTEcohort_final %>% 
  filter(metformin_flag == 0)

group2 <- TTEcohort_final %>%
  filter(metformin_flag == 1 & group == "Study cohort")

# Add Group identifier to TTEcohort_final
TTEcohort_final <- TTEcohort_final %>%
  mutate(Group = case_when(metformin_flag == 0 ~ "Group 1",
                           metformin_flag == 1 & group == "Study cohort" ~ "Group 2",
                           TRUE ~ NA_character_)) %>%
  filter(!is.na(Group))

# Create a summary table with the specified order, formatting, and exclusions
table_summary <- TTEcohort_final %>%
  select(Group, gender, age_atfirstap, age_atprevcohortentry, age_atfirstdiag, ethnicity_cat_cprdhes,
    last_smi_diag, region, patprac_2019imd_quintile, baseline_bmi_cat,
    priorcerebrovasculardisease, priormyocardialinfarction, priorliverdisease, priorrenaldisease,
    priorhypertension, priordyslipidaemia, priordiabetes, pcos_flag,
    antidepressant_prior2years, lipiddrugs_prior2years, insulin, prioralcoholabuse, priorsubstanceabuse,
    smoking_status_cat, baseline_hba1c, baseline_glucose) %>%
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
table_summary <- table_summary %>%
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
print(table_summary)

# Convert table_summary to a tibble (data frame) for export
table_summary_df <- as_tibble(table_summary)

# Additional variables
# Baseline weight
baseline_weight_stats <- TTEcohort_final %>%
  group_by(Group) %>%  # Group by Group 1 and Group 2
  summarise(Mean_Baseline_Weight = round(mean(baseline_weightkg, na.rm = TRUE), 5),
            SD_Baseline_Weight = round(sd(baseline_weightkg, na.rm = TRUE), 5),
            .groups = "drop")  # Ensure the output is ungrouped

# Print the results ensuring three decimal places in display
cat("Mean and Standard Deviation of Baseline Weight (kg) by Group:\n")
print(as.data.frame(baseline_weight_stats), digits = 5, nsmall = 5)

# Package version used for this script
# dplyr version: 1.1.4 
# gtsummary version: 2.0.2 
# ggplot2 version: 3.5.1 
# tidyr version: 1.3.1 
# lubridate version: 1.9.3 
# tableone version: 0.13.2 
# readr version: 2.1.5 
# survival version: 3.7.0 
# survminer version: 0.4.9 
# gt version: 0.11.0 
# mice version: 3.16.0 
# tidyselect version: 1.2.1 
# broom.mixed version: 0.2.9.6 
