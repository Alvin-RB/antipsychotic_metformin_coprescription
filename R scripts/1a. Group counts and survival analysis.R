# Metformin and SGA co-prescrption #####
# Last run: 02/12/2024

# Group counts and survival analysis

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

# Creating SMI group with Metformin data
metformin <- SMIantidiabetics_C %>%
  filter(Antidiabetic == "Metformin" & issuedate < '2020-01-01')

ids <- tte_prevalent_cohort %>%
  select(patid)

merge <- ids %>%
  left_join(metformin, by = "patid") %>%
  select(patid, Antidiabetic) %>%
  distinct()

merge <- merge %>%
  mutate(metformin_flag = ifelse(Antidiabetic == "Metformin", 1, 0))

merge2 <- merge %>%
  select(metformin_flag) %>%
  tbl_summary()

TTEcohortwithmetforminflag <- merge %>%
  left_join(tte_prevalent_cohort, by = "patid")

TTEcohortwithmetforminflag <- TTEcohortwithmetforminflag %>%
  mutate(metformin_flag = if_else(is.na(metformin_flag), 0, metformin_flag))

# Add earliest issue date of Metformin to TTEcohortwithmetformin flag. 

# Step 1: Identify the earliest issuedate for each patid in SMIantidiabetics_C
earliest_issuedate <- SMIantidiabetics_C %>%
  filter(Antidiabetic == "Metformin" & issuedate < '2020-01-01') %>%
  group_by(patid) %>%
  summarise(earliest_issuedate = min(issuedate, na.rm = TRUE))  # Get the earliest issuedate per patid

# Step 2: Perform a left join with TTEcohortwithmetforminflag
merged_data <- TTEcohortwithmetforminflag %>%
  left_join(earliest_issuedate, by = "patid")

# Step 3: Reorder columns so that 'earliest_issuedate' is next to 'metformin_flag'
TTEcohortwithmetforminflag_andissuedate <- merged_data %>%
  select(patid, Antidiabetic, metformin_flag, earliest_issuedate, everything())  # Use everything() and then reorder

# check if there are any metformin flags without earliest dates (and vice versa): there are none.
check1 <- TTEcohortwithmetforminflag_andissuedate %>%
  filter(!is.na(metformin_flag) & is.na(earliest_issuedate) | !is.na(earliest_issuedate) & is.na(metformin_flag))
rm(check1)

# View the final result
print(TTEcohortwithmetforminflag_andissuedate)

table(TTEcohortwithmetforminflag_andissuedate$metformin_flag, useNA = "always")


df <- TTEcohortwithmetforminflag_andissuedate

df <- df %>%
  mutate(earliest_issuedate = as.Date(earliest_issuedate),
         cohortentrydate = as.Date(cohortentrydate))

# Sorting into groups
TTEcohortwithmetforminflag_andissuedate <- TTEcohortwithmetforminflag_andissuedate %>%
  mutate(diff = earliest_issuedate - cohortentrydate) %>%
  mutate(group = case_when(!is.na(earliest_issuedate) & diff >= -30.25 & diff <= 730.5 ~ "Study cohort",
                           !is.na(earliest_issuedate) & diff < 30.25 ~ "Before",
                           !is.na(earliest_issuedate) & diff > 730.5 ~ "After",
                           TRUE ~ NA))

table(TTEcohortwithmetforminflag_andissuedate$group, useNA = "always")

# Save as an .RData file
save(TTEcohortwithmetforminflag_andissuedate, file = "TTEcohortwithmetforminflag_andissuedate.RData")

# Prepare the survival data
survivaldata <- TTEcohortwithmetforminflag_andissuedate %>%
  select(patid, earliest_issuedate, cohortentrydate, regenddate, lcd, deathdate, metformin_flag, group) %>% 
  mutate(event = case_when(is.na(metformin_flag) ~ 0, TRUE ~ metformin_flag)) %>%
  mutate(event = case_when(group == "After" ~ 0,
                           TRUE ~ event)) %>% # censor events at 2y
  mutate(end_date = pmin(earliest_issuedate, regenddate, lcd, deathdate, na.rm = TRUE)) %>%
  mutate(time = as.numeric(end_date - cohortentrydate)) %>%
  mutate(timenew = case_when(earliest_issuedate == cohortentrydate ~ 0.5, # assign those starting on cohortentrydate with minimal f-up time
                             time < -30.25 ~ 0, # set f-up time to 0 if started more than 1 month before
                             time >= -30.25 & time <= -1 ~ 0.5, # assign those starting within 1m prior to cohortentrydate with minimal f-up time
                             time > 730.5 ~ 730, # censor at 2y
                             TRUE ~ time)) %>%
  filter(timenew > 0)

# Fit survival model
survobj <- Surv(time = survivaldata$timenew / 365, event = survivaldata$event)
km_fit <- survfit(survobj ~ 1)

# Plot the survival curves using ggsurvplot
figure <- ggsurvplot(
  km_fit,                          # the fitted survival object
  data = survivaldata,              # dataset used for the survival fit
  fun = "event",                    # flip the plot to show incidence (cumulative events)
  conf.int = FALSE,                 # remove confidence intervals
  xlab = "Time (years)",            # x-axis label
  ylab = "Cumulative Incidence of Metformin Prescription",    # y-axis label to reflect incidence
  xlim = c(0, 2),                   # limit x-axis to 0 to 2 years
  ylim = c(0, 0.05),                # set y-axis limit to 5% (0.05)
  break.time.by = 0.5,              # set x-axis breaks by 0.5-year intervals
  risk.table = "nrisk_cumevents",   # show risk table with cumulative events
  risk.table.col = "strata",        # color by strata
  legend.title = "Metformin Prescriptions", # legend characteristics
  font.legend = 10,                 # set smaller font size for legend
  font.tickslab = 8,                # reduce font size of tick labels
  font.x = 10,                      # adjust x-axis label font size
  font.y = 10,                      # adjust y-axis label font size
  palette = "Dark2",                # use Dark2 color palette
  risk.table.height = 0.2,          # adjust the height of the risk table
  lwd = 0.2,                        # set line width to 0.2 for a finer line
  ggtheme = theme_minimal(base_size = 9) + # use a minimal theme with a base font size
    theme(
      axis.text = element_text(size = 8),  # smaller axis text
      axis.title = element_text(size = 9), # smaller axis titles
      legend.text = element_text(size = 8),# smaller legend text
      plot.title = element_text(size = 11, face = "bold") # adjust plot title font
    ),
  tables.theme = theme_cleantable() +      # cleaner table theme for risk table
    theme(
      text = element_text(size = 7)))        # set smaller font size for risk table

# Display the plot
print(figure)

tbl_survfit(km_fit,
            times = c(0, 0.5, 1, 1.5, 2), type = "risk")


# Package versions used for this script
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
