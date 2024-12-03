# Metformin and SGA co-prescrption #####
# Last run: 02/12/2024

# Annual period prevalences

# Clear memory
#rm(list = ls())

# Load required libraries
library(dplyr)
library(ggplot2)
library(ggthemes)

active_by_year <- TTEcohort_final %>%
  select(patid, cohortentrydate, regstartdate, regenddate, lcd, deathdate) %>%
  mutate(endfup_date = cohortentrydate + 730,
         end = pmin(regenddate, lcd, deathdate, endfup_date, na.rm = TRUE)) %>%
  mutate(start_year = year(cohortentrydate),
         end_year = year(end)) %>% 
  select(patid, start_year, end_year)

add_year_indicators <- function(df, start_year, end_year) {
  years <- 2005:2017 # Create a vector of study years
  for (year in years) {
    df <- df %>%
      mutate(!!paste0("year_", year) := as.integer(year >= start_year & year <= end_year))} # Add column for each year with active status
  return(df)}

# Call the function for each data frame
active_by_year <- add_year_indicators(active_by_year, start_year, end_year)
active_by_year <- tidyr::gather(active_by_year, key = "year", value = "active", starts_with("year_"))
active_by_year$year <- as.numeric(sub("year_", "", active_by_year$year))

aggregate_count <- function(df) { 
  df %>%
    group_by(year) %>%
    summarize(count_active = sum(as.numeric(active), na.rm = TRUE))}

active_counts <- active_by_year %>%
  aggregate_count()

rm(list = ls(pattern = "active_by_year"))

overall_presc_rate_counts <- SMIantidiabetics_C %>% 
  filter(Antidiabetic == "Metformin" & issuedate < '2020-01-01') %>%
  select(patid, issuedate, Antidiabetic) %>%
  mutate(issue_year = year(issuedate)) %>%
  select(-issuedate) %>%
  distinct() %>%
  inner_join(select(TTEcohort_final, patid, cohortentrydate), by = "patid", relationship = "many-to-many") %>%
  mutate(entry_year = year(cohortentrydate)) %>%
  select(-cohortentrydate) %>%
  distinct() %>%
  filter(issue_year >= entry_year) %>% # remove prescriptions that are before cohortentryyear, as we are only monitoring from cohort entry
  select(-entry_year) %>%
  group_by(issue_year) %>% 
  summarise(n_exp = n_distinct(patid)) %>% 
  ungroup() %>%
  rename(year = issue_year) %>%
  left_join(active_counts, by = "year")

# Standardised to 1,000 and add CIs
calculate_CI <- function(exposed, total, confidence_level = 0.95) {
  # Calculate sample rate (per 1000)
  rate <- exposed / total * 1000
  
  # Calculate standard error
  SE <- sqrt(rate * (1000 - rate) / total)
  
  # Calculate z value based on the confidence level
  z <- qnorm((1 + confidence_level) / 2)
  
  # Calculate margin of error
  MOE <- z * SE
  
  # Calculate confidence interval
  lower_bound <- rate - MOE
  upper_bound <- rate + MOE
  
  # Return CI as a character string
  return(paste(round(lower_bound, 2), ",", round(upper_bound, 2)))}

# Mutate to add CI to the data frame
overall_presc_rate_counts <- overall_presc_rate_counts %>%
  mutate(rate_obs_1000 = n_exp / count_active * 1000,
         CI = calculate_CI(n_exp, count_active)) %>%
  separate(CI, into = c("lower_ci", "upper_ci"), sep = " , ") %>%
  mutate_at(vars(4:6), as.numeric) %>%
  mutate_at(vars(4:6), ~ round(., digits = 1))

overall_presc_rate_counts %>%
  gt() %>%
  cols_label( # tidy column labels
    n_exp = "Number prescribed AP",
    count_active = "Number active",
    year = "Year",
    rate_obs_1000 = "Rate/1000 (95% CI)") %>%
  cols_merge(columns = c("rate_obs_1000", "lower_ci", "upper_ci"), pattern = "{1} ({2}, {3})")

# Create the plot
period_prev <- ggplot(overall_presc_rate_counts, aes(x = year, y = rate_obs_1000)) +
  # Confidence interval ribbon
  geom_ribbon(aes(ymin = lower_ci, ymax = upper_ci), 
              alpha = 0.2, fill = "#1E90FF") +
  
  # Line plot
  geom_line(color = "#1E90FF", size = 1.2) +
  
  # Points at each data point
  geom_point(color = "#1E90FF", size = 3) +
  
  # Formatting for publication-ready appearance
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5)
  ) +
  
  # Labels
  labs(
    x = "Year",
    y = "Period prevalence (per 1,000 patients)"
  ) +
  
  # Set x-axis to specific years
  scale_x_continuous(
    breaks = c(2005, 2007, 2009, 2011, 2013, 2015, 2017),
    limits = c(2005, 2017)
  ) +
  
  # Set y-axis limit
  scale_y_continuous(
    expand = c(0, 0), 
    limits = c(0, 1000))

print(period_prev)

ggsave(period_prev,
       filename = "period_prev.png", 
       dpi = 300, width = 12, height = 9, bg = "white")

# Package version used in this script
# dplyr version: 1.1.4 
# ggplot2 version: 3.5.1