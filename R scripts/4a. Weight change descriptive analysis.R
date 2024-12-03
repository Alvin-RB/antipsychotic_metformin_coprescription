# Metformin and SGA co-prescrption #####
# Last run: 02/12/2024

# Weight change over time (descriptive)

# Clear memory
rm(list = ls())

library(dplyr)
library(tidyr)
library(ggplot2)
library(gridExtra)
library(grid)
library(purrr)
library(janitor)
library(patchwork)

# Unadjusted percentage change in weight
# Using imputed data
pct_change_imp <- tte_imputed_long %>%
  select(.imp, .id, patid, trt_group, baseline_weightkg, sixm_weightkg, oney_weightkg, twoy_weightkg) %>%
  left_join(select(TTEcohort_final, patid, diff, metformin_flag, pcos_flag), by = "patid") %>%
  filter(.imp > 0 & !is.na(metformin_flag)) %>%
  dplyr::mutate(pct_change_6m = ((sixm_weightkg - baseline_weightkg) / baseline_weightkg) * 100,
                pct_change_1y = ((oney_weightkg - baseline_weightkg) / baseline_weightkg) * 100,
                pct_change_2y = ((twoy_weightkg - baseline_weightkg) / baseline_weightkg) * 100) %>%
  mutate(keep = case_when(is.na(diff) | diff >= -30 & diff <= 90 ~ 1,
                          TRUE ~ 0)) %>%
  filter(keep == 1)

# First calculate means and variances for each imputed dataset
imp_stats <- pct_change_imp %>%
  group_by(.imp, metformin_flag) %>%
  summarize(n = n(),
            mean_6m = mean(pct_change_6m, na.rm = TRUE),
            var_6m = var(pct_change_6m, na.rm = TRUE)/n,  # Using SE^2
            mean_1y = mean(pct_change_1y, na.rm = TRUE),
            var_1y = var(pct_change_1y, na.rm = TRUE)/n,
            mean_2y = mean(pct_change_2y, na.rm = TRUE),
            var_2y = var(pct_change_2y, na.rm = TRUE)/n)

# Apply Rubin's rules
pool_estimates <- function(means, variances) {
  m <- length(means)  # number of imputations
  
  # Overall estimate (Q-bar)
  q_bar <- mean(means)
  
  # Within-imputation variance (U-bar)
  u_bar <- mean(variances)
  
  # Between-imputation variance (B)
  b <- sum((means - q_bar)^2) / (m - 1)
  
  # Total variance (T)
  t <- u_bar + (1 + 1/m) * b
  
  # Degrees of freedom using Barnard-Rubin adjustment
  lambda <- (b * (1 + 1/m)) / t
  df_old <- (m - 1) / lambda^2
  df_obs <- (1/lambda - 1) / (1 - lambda)  # approximate df from complete data
  df <- 1 / (1/df_old + 1/df_obs)
  
  # Calculate 95% CI
  ci_lower <- q_bar - qt(0.975, df) * sqrt(t)
  ci_upper <- q_bar + qt(0.975, df) * sqrt(t)
  
  return(c(mean = q_bar, 
           ci_lower = ci_lower, 
           ci_upper = ci_upper,
           se = sqrt(t)))}  # also returning SE for potential use

# Apply Rubin's rules for each metformin group and timepoint
results <- imp_stats %>%
  group_by(metformin_flag) %>%
  summarize(
    # 6 months
    res_6m = list(pool_estimates(mean_6m, var_6m)),
    # 1 year
    res_1y = list(pool_estimates(mean_1y, var_1y)),
    # 2 years
    res_2y = list(pool_estimates(mean_2y, var_2y))) %>%
  mutate(
    # 6 months
    mean_6m = map_dbl(res_6m, ~.x["mean"]),
    ci_lower_6m = map_dbl(res_6m, ~.x["ci_lower"]),
    ci_upper_6m = map_dbl(res_6m, ~.x["ci_upper"]),
    # 1 year
    mean_1y = map_dbl(res_1y, ~.x["mean"]),
    ci_lower_1y = map_dbl(res_1y, ~.x["ci_lower"]),
    ci_upper_1y = map_dbl(res_1y, ~.x["ci_upper"]),
    # 2 years
    mean_2y = map_dbl(res_2y, ~.x["mean"]),
    ci_lower_2y = map_dbl(res_2y, ~.x["ci_lower"]),
    ci_upper_2y = map_dbl(res_2y, ~.x["ci_upper"])) %>%
  select(-res_6m, -res_1y, -res_2y)

final_results_long <- results %>%
  pivot_longer(cols = -metformin_flag,
               names_to = c(".value", "timepoint"),
               names_pattern = "(.+)_(.+)") %>%
         mutate(timepoint = factor(timepoint, levels = c("Baseline", "6m", "1y", "2y")),
         formatted = sprintf("%.2f (%.2f, %.2f)", mean, ci_lower, ci_upper),
         metformin_flag = factor(metformin_flag, levels = c(0, 1), labels = c("SGA only", "SGA+Metformin"))) %>%
  bind_rows(
    tibble(metformin_flag = factor(c("SGA only", "SGA+Metformin"), levels = c("SGA only", "SGA+Metformin")), 
           timepoint = "Baseline", 
           mean = 0, 
           ci_lower = 0, 
           ci_upper = 0, 
           formatted = "0.00 (0.00, 0.00)")) %>%
  mutate(timepoint = factor(timepoint, levels = c("Baseline", "6m", "1y", "2y"))) %>%
  arrange(timepoint)

# Plot
graph_imp <- ggplot(final_results_long, aes(x = timepoint, y = mean, color = metformin_flag, group = metformin_flag)) + 
  geom_line(size = 1.2) + 
  geom_point(size = 3) + 
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.1) + 
  labs( 
    y = "Change in weight (%)", 
    x = "Time-point") + 
  scale_color_manual(values = c("SGA only" = "#1f77b4", "SGA+Metformin" = "#ff7f0e")) + 
  scale_y_continuous(breaks = seq(-5, 10, by = 2.5)) +  # Add y-axis labels every 2.5 units
  theme_minimal(base_size = 14) + 
  theme( 
    legend.title = element_text(size = 16), 
    legend.text = element_text(size = 14), 
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 14), 
    plot.title = element_blank(), 
    panel.grid.major = element_line(size = 0.5, color = "gray85"), 
    panel.grid.minor = element_blank(), 
    panel.border = element_rect(color = "dark grey", fill = NA, size = 2), 
    legend.position = "bottom") + 
  guides(color = guide_legend(title = "Group"))

# Complete case

cc <- tte_imputed_long %>%
  select(.imp, .id, patid, trt_group, baseline_weightkg, sixm_weightkg, oney_weightkg, twoy_weightkg) %>%
  left_join(select(TTEcohort_final, patid, diff, metformin_flag, pcos_flag), by = "patid") %>%
  filter(.imp == 0 & !is.na(metformin_flag)) %>%
  dplyr::mutate(pct_change_6m = ((sixm_weightkg - baseline_weightkg) / baseline_weightkg) * 100,
                pct_change_1y = ((oney_weightkg - baseline_weightkg) / baseline_weightkg) * 100,
                pct_change_2y = ((twoy_weightkg - baseline_weightkg) / baseline_weightkg) * 100) %>%
  mutate(keep = case_when(is.na(diff) | diff >= -30 & diff <= 90 ~ 1,
                          TRUE ~ 0)) %>%
  filter(keep == 1)

cc_f <- cc %>%
  pivot_longer(cols = starts_with("pct_change"), 
               names_to = "timepoint", 
               values_to = "pct_change") %>%
  group_by(metformin_flag, timepoint) %>%
  summarise(
    mean = mean(pct_change, na.rm = TRUE),
    se = sd(pct_change, na.rm = TRUE) / sqrt(sum(!is.na(pct_change))),  # Standard error
    ci_lower = mean - qt(0.975, df = sum(!is.na(pct_change)) - 1) * se,  # 95% CI lower bound
    ci_upper = mean + qt(0.975, df = sum(!is.na(pct_change)) - 1) * se,  # 95% CI upper bound
    n_obs = sum(!is.na(pct_change)),  # Number of non-missing observations
    .groups = "drop"
  ) %>%
  mutate(timepoint = case_when(timepoint == "pct_change_6m" ~ "6m",
                               timepoint == "pct_change_1y" ~ "1y",
                               timepoint == "pct_change_2y" ~ "2y"),
         timepoint = factor(timepoint, levels = c("Baseline", "6m", "1y", "2y")),
         formatted = sprintf("%.2f (%.2f, %.2f)", mean, ci_lower, ci_upper),
         metformin_flag = factor(metformin_flag, levels = c(0, 1), labels = c("SGA only", "SGA+Metformin"))) %>%
  bind_rows(
    tibble(metformin_flag = factor(c("SGA only", "SGA+Metformin"), levels = c("SGA only", "SGA+Metformin")), 
           timepoint = "Baseline", 
           mean = 0, 
           ci_lower = 0, 
           ci_upper = 0, 
           formatted = "0.00 (0.00, 0.00)", 
           n_obs = 0)) %>%
  mutate(timepoint = factor(timepoint, levels = c("Baseline", "6m", "1y", "2y"))) %>%
  arrange(timepoint)

graph_cc <- ggplot(cc_f, aes(x = timepoint, y = mean, color = metformin_flag, group = metformin_flag)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.1) +
  labs(
    y = "Change in weight (%)",
    x = "Time-point") +
  scale_color_manual(values = c("SGA only" = "#1f77b4", "SGA+Metformin" = "#ff7f0e")) +  # Custom colors
  theme_minimal(base_size = 14) +  # Larger base font size for publication
  theme(
    legend.title = element_text(size = 16),  # Change the legend title size
    legend.text = element_text(size = 14),   # Change the legend text size
    axis.title = element_text(size = 16),    # Axis titles size
    axis.text = element_text(size = 14),     # Axis text size
    plot.title = element_blank(),  # Remove the plot title
    panel.grid.major = element_line(size = 0.5, color = "gray85"),    # Lighter grid lines
    panel.grid.minor = element_blank(),     # No minor grid lines
    panel.border = element_rect(color = "dark grey", fill = NA, size = 2), 
    legend.position = "bottom") +  # Move the legend to the bottom
  guides(color = guide_legend(title = "Group"))  # Change legend title to "Group"


# Individual SGAs

analyze_weight_change <- function(data, tte_cohort, treatment_group) {
  # Imputed dataset analysis
  pct_change_imp <- data %>%
    select(.imp, .id, patid, trt_group, baseline_weightkg, sixm_weightkg, oney_weightkg, twoy_weightkg) %>%
    filter(trt_group == treatment_group) %>%
    left_join(select(tte_cohort, patid, diff, metformin_flag, pcos_flag), by = "patid") %>%
    filter(.imp > 0 & !is.na(metformin_flag)) %>%
    dplyr::mutate(
      pct_change_6m = ((sixm_weightkg - baseline_weightkg) / baseline_weightkg) * 100,
      pct_change_1y = ((oney_weightkg - baseline_weightkg) / baseline_weightkg) * 100,
      pct_change_2y = ((twoy_weightkg - baseline_weightkg) / baseline_weightkg) * 100
    ) %>%
    mutate(keep = case_when(is.na(diff) | diff >= -30 & diff <= 90 ~ 1, TRUE ~ 0)) %>%
    filter(keep == 1)
  
  # Imputed stats calculation (same as original code)
  imp_stats <- pct_change_imp %>%
    group_by(.imp, metformin_flag) %>%
    summarize(
      n = n(),
      mean_6m = mean(pct_change_6m, na.rm = TRUE),
      var_6m = var(pct_change_6m, na.rm = TRUE)/n,
      mean_1y = mean(pct_change_1y, na.rm = TRUE),
      var_1y = var(pct_change_1y, na.rm = TRUE)/n,
      mean_2y = mean(pct_change_2y, na.rm = TRUE),
      var_2y = var(pct_change_2y, na.rm = TRUE)/n
    )
  
  # Pooling function (same as original)
  pool_estimates <- function(means, variances) {
    m <- length(means)
    q_bar <- mean(means)
    u_bar <- mean(variances)
    b <- sum((means - q_bar)^2) / (m - 1)
    t <- u_bar + (1 + 1/m) * b
    lambda <- (b * (1 + 1/m)) / t
    df_old <- (m - 1) / lambda^2
    df_obs <- (1/lambda - 1) / (1 - lambda)
    df <- 1 / (1/df_old + 1/df_obs)
    
    ci_lower <- q_bar - qt(0.975, df) * sqrt(t)
    ci_upper <- q_bar + qt(0.975, df) * sqrt(t)
    
    return(c(mean = q_bar, ci_lower = ci_lower, ci_upper = ci_upper, se = sqrt(t)))
  }
  
  # Apply Rubin's rules
  results <- imp_stats %>%
    group_by(metformin_flag) %>%
    summarize(
      res_6m = list(pool_estimates(mean_6m, var_6m)),
      res_1y = list(pool_estimates(mean_1y, var_1y)),
      res_2y = list(pool_estimates(mean_2y, var_2y))
    ) %>%
    mutate(
      mean_6m = map_dbl(res_6m, ~.x["mean"]),
      ci_lower_6m = map_dbl(res_6m, ~.x["ci_lower"]),
      ci_upper_6m = map_dbl(res_6m, ~.x["ci_upper"]),
      mean_1y = map_dbl(res_1y, ~.x["mean"]),
      ci_lower_1y = map_dbl(res_1y, ~.x["ci_lower"]),
      ci_upper_1y = map_dbl(res_1y, ~.x["ci_upper"]),
      mean_2y = map_dbl(res_2y, ~.x["mean"]),
      ci_lower_2y = map_dbl(res_2y, ~.x["ci_lower"]),
      ci_upper_2y = map_dbl(res_2y, ~.x["ci_upper"])
    ) %>%
    select(-res_6m, -res_1y, -res_2y)
  
  # Long format results for plotting
  final_results_long <- results %>%
    pivot_longer(cols = -metformin_flag,
                 names_to = c(".value", "timepoint"),
                 names_pattern = "(.+)_(.+)") %>%
    mutate(
      timepoint = factor(timepoint, levels = c("baseline", "6m", "1y", "2y")),
      formatted = sprintf("%.2f (%.2f, %.2f)", mean, ci_lower, ci_upper),
      metformin_flag = factor(metformin_flag, levels = c(0, 1), labels = c("SGA only", "SGA+Metformin"))
    ) %>%
    bind_rows(
      tibble(
        metformin_flag = factor(c("SGA only", "SGA+Metformin"), levels = c("SGA only", "SGA+Metformin")), 
        timepoint = "baseline", 
        mean = 0, 
        ci_lower = 0, 
        ci_upper = 0, 
        formatted = "0.00 (0.00, 0.00)")) %>%
    mutate(timepoint = factor(timepoint, levels = c("baseline", "6m", "1y", "2y"))) %>%
    arrange(timepoint)
  
  # Imputed dataset plot
  graph_imp <- ggplot(final_results_long, aes(x = timepoint, y = mean, color = metformin_flag, group = metformin_flag)) +
    geom_line(size = 1.2) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.1) +
    labs(y = "Change in weight (%)", x = "Time-point") +
    scale_color_manual(values = c("SGA only" = "#1f77b4", "SGA+Metformin" = "#ff7f0e")) +
    theme_minimal(base_size = 14) +
    theme(
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14),
      plot.title = element_blank(),
      panel.grid.major = element_line(size = 0.5, color = "gray85"),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      legend.position = "bottom"
    ) +
    guides(color = guide_legend(title = "Group"))
  
  # Complete case analysis
  cc <- data %>%
    select(.imp, .id, patid, trt_group, baseline_weightkg, sixm_weightkg, oney_weightkg, twoy_weightkg) %>%
    filter(trt_group == treatment_group) %>%
    left_join(select(tte_cohort, patid, diff, metformin_flag, pcos_flag), by = "patid") %>%
    filter(.imp == 0 & !is.na(metformin_flag)) %>%
    dplyr::mutate(
      pct_change_6m = ((sixm_weightkg - baseline_weightkg) / baseline_weightkg) * 100,
      pct_change_1y = ((oney_weightkg - baseline_weightkg) / baseline_weightkg) * 100,
      pct_change_2y = ((twoy_weightkg - baseline_weightkg) / baseline_weightkg) * 100
    ) %>%
    mutate(keep = case_when(is.na(diff) | diff >= -30 & diff <= 90 ~ 1, TRUE ~ 0)) %>%
    filter(keep == 1)
  
  cc_f <- cc %>%
    pivot_longer(cols = starts_with("pct_change"), 
                 names_to = "timepoint", 
                 values_to = "pct_change") %>%
    group_by(metformin_flag, timepoint) %>%
    summarise(
      mean = mean(pct_change, na.rm = TRUE),
      se = sd(pct_change, na.rm = TRUE) / sqrt(sum(!is.na(pct_change))),
      ci_lower = mean - qt(0.975, df = sum(!is.na(pct_change)) - 1) * se,
      ci_upper = mean + qt(0.975, df = sum(!is.na(pct_change)) - 1) * se,
      n_obs = sum(!is.na(pct_change)),
      .groups = "drop"
    ) %>%
    mutate(
      timepoint = case_when(
        timepoint == "pct_change_6m" ~ "6m",
        timepoint == "pct_change_1y" ~ "1y",
        timepoint == "pct_change_2y" ~ "2y"
      ),
      timepoint = factor(timepoint, levels = c("baseline", "6m", "1y", "2y")),
      formatted = sprintf("%.2f (%.2f, %.2f)", mean, ci_lower, ci_upper),
      metformin_flag = factor(metformin_flag, levels = c(0, 1), labels = c("SGA only", "SGA+Metformin"))
    ) %>%
    bind_rows(
      tibble(
        metformin_flag = factor(c("SGA only", "SGA+Metformin"), levels = c("SGA only", "SGA+Metformin")), 
        timepoint = "baseline", 
        mean = 0, 
        ci_lower = 0, 
        ci_upper = 0, 
        formatted = "0.00 (0.00, 0.00)", 
        n_obs = 0
      )
    ) %>%
    mutate(timepoint = factor(timepoint, levels = c("baseline", "6m", "1y", "2y"))) %>%
    arrange(timepoint)
  
  # Complete case plot
  graph_cc <- ggplot(cc_f, aes(x = timepoint, y = mean, color = metformin_flag, group = metformin_flag)) +
    geom_line(size = 1.2) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.1) +
    labs(y = "Change in weight (%)", x = "Time-point") +
    scale_color_manual(values = c("SGA only" = "#1f77b4", "SGA+Metformin" = "#ff7f0e")) +
    theme_minimal(base_size = 14) +
    theme(
      legend.title = element_text(size = 16),
      legend.text = element_text(size = 14),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 14),
      plot.title = element_blank(),
      panel.grid.major = element_line(size = 0.5, color = "gray85"),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      legend.position = "bottom"
    ) +
    guides(color = guide_legend(title = "Group"))
  
  # Return results
  return(list(
    imputed_results = results,
    imputed_plot = graph_imp,
    complete_case_results = cc_f,
    complete_case_plot = graph_cc
  ))}

ari_results <- analyze_weight_change(tte_imputed_long, TTEcohort_final, "Aripiprazole")
ola_results <- analyze_weight_change(tte_imputed_long, TTEcohort_final, "Olanzapine")
que_results <- analyze_weight_change(tte_imputed_long, TTEcohort_final, "Quetiapine")
ris_results <-  analyze_weight_change(tte_imputed_long, TTEcohort_final, "Risperidone")

# Add titles to each plot with axis label removal
# Function to create consistent plot processing
process_plots <- function(plot_list, is_imputed = TRUE, ylim_imputed = c(-10, 10.5), ylim_cc = c(-6, 6)) {
  lapply(names(plot_list), function(med) {
    plot <- plot_list[[med]] + 
      labs(title = med)
    
    # Remove x-axis label for Aripiprazole and Olanzapine
    if (med %in% c("Aripiprazole", "Olanzapine")) {
      plot <- plot + 
        theme(axis.title.x = element_blank()) +
        theme(legend.position = "none")  # Remove legend for these plots
    }
    
    # Remove y-axis label for Olanzapine and Risperidone
    if (med %in% c("Olanzapine", "Risperidone")) {
      plot <- plot + theme(axis.title.y = element_blank())
    }
    
    plot + 
      theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
      coord_cartesian(ylim = if (is_imputed) ylim_imputed else ylim_cc)
  })
}

# Process imputed and complete case plots
medication_plots <- list(
  Aripiprazole = ari_results$imputed_plot,
  Olanzapine = ola_results$imputed_plot,
  Quetiapine = que_results$imputed_plot,
  Risperidone = ris_results$imputed_plot)


pct_imputed_grid <- wrap_plots(process_plots(medication_plots, is_imputed = TRUE), ncol = 2) + 
  plot_annotation(
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
  ) & 
  theme(legend.position = "bottom")

cc_medication_plots <- list(
  Aripiprazole = ari_results$complete_case_plot,
  Olanzapine = ola_results$complete_case_plot,
  Quetiapine = que_results$complete_case_plot,
  Risperidone = ris_results$complete_case_plot)

pct_cc_grid <- wrap_plots(process_plots(cc_medication_plots, is_imputed = FALSE), ncol = 2) + 
  plot_annotation(
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
  ) & 
  theme(legend.position = "bottom")

# Display or save
pct_imputed_grid
pct_cc_grid

graph_imp
graph_cc

# Absolute change

abs_change_imp <- tte_imputed_long %>%
  select(.imp, .id, patid, trt_group, baseline_weightkg, sixm_weightkg, oney_weightkg, twoy_weightkg) %>%
  left_join(select(TTEcohort_final, patid, diff, metformin_flag, pcos_flag), by = "patid") %>%
  filter(.imp > 0 & !is.na(metformin_flag)) %>%
  mutate(keep = case_when(is.na(diff) | diff >= -30 & diff <= 90 ~ 1,
                          TRUE ~ 0)) %>%
  filter(keep == 1)

# First calculate means and variances for each imputed dataset
imp_stats <- abs_change_imp %>%
  group_by(.imp, metformin_flag) %>%
  summarize(n = n(),
            mean_baseline = mean(baseline_weightkg, na.rm = TRUE),
            var_baseline =  var(baseline_weightkg, na.rm = TRUE)/n,
            mean_6m = mean(sixm_weightkg, na.rm = TRUE),
            var_6m = var(sixm_weightkg, na.rm = TRUE)/n,  # Using SE^2
            mean_1y = mean(oney_weightkg, na.rm = TRUE),
            var_1y = var(oney_weightkg, na.rm = TRUE)/n,
            mean_2y = mean(twoy_weightkg, na.rm = TRUE),
            var_2y = var(twoy_weightkg, na.rm = TRUE)/n)

# Revised function to apply Rubin's rules
# Following formulas from:
# https://www.stata-journal.com/article.html?article=st0139
pool_estimates <- function(means, variances) {
  m <- length(means)  # number of imputations
  
  # Overall estimate (Q-bar)
  q_bar <- mean(means)
  
  # Within-imputation variance (U-bar)
  u_bar <- mean(variances)
  
  # Between-imputation variance (B)
  b <- sum((means - q_bar)^2) / (m - 1)
  
  # Total variance (T)
  t <- u_bar + (1 + 1/m) * b
  
  # Degrees of freedom using Barnard-Rubin adjustment
  lambda <- (b * (1 + 1/m)) / t
  df_old <- (m - 1) / lambda^2
  df_obs <- (1/lambda - 1) / (1 - lambda)  # approximate df from complete data
  df <- 1 / (1/df_old + 1/df_obs)
  
  # Calculate 95% CI
  ci_lower <- q_bar - qt(0.975, df) * sqrt(t)
  ci_upper <- q_bar + qt(0.975, df) * sqrt(t)
  
  return(c(mean = q_bar, 
           ci_lower = ci_lower, 
           ci_upper = ci_upper,
           se = sqrt(t)))}  # also returning SE for potential use

# Apply Rubin's rules for each metformin group and timepoint
results <- imp_stats %>%
  group_by(metformin_flag) %>%
  summarize(
    # baseline
    res_baseline = list(pool_estimates(mean_baseline, var_baseline)),
    # 6 months
    res_6m = list(pool_estimates(mean_6m, var_6m)),
    # 1 year
    res_1y = list(pool_estimates(mean_1y, var_1y)),
    # 2 years
    res_2y = list(pool_estimates(mean_2y, var_2y))) %>%
  mutate(
    # Baseline
    mean_baseline = map_dbl(res_baseline, ~.x["mean"]),
    ci_lower_baseline = map_dbl(res_baseline, ~.x["ci_lower"]),
    ci_upper_baseline = map_dbl(res_baseline, ~.x["ci_upper"]),
    # 6 months
    mean_6m = map_dbl(res_6m, ~.x["mean"]),
    ci_lower_6m = map_dbl(res_6m, ~.x["ci_lower"]),
    ci_upper_6m = map_dbl(res_6m, ~.x["ci_upper"]),
    # 1 year
    mean_1y = map_dbl(res_1y, ~.x["mean"]),
    ci_lower_1y = map_dbl(res_1y, ~.x["ci_lower"]),
    ci_upper_1y = map_dbl(res_1y, ~.x["ci_upper"]),
    # 2 years
    mean_2y = map_dbl(res_2y, ~.x["mean"]),
    ci_lower_2y = map_dbl(res_2y, ~.x["ci_lower"]),
    ci_upper_2y = map_dbl(res_2y, ~.x["ci_upper"])) %>%
  select(-res_baseline, -res_6m, -res_1y, -res_2y)

final_results_long_abs <- results %>%
  pivot_longer(cols = -metformin_flag,
               names_to = c(".value", "timepoint"),
               names_pattern = "(.+)_(.+)") %>%
  mutate(timepoint = factor(timepoint, levels = c("baseline", "6m", "1y", "2y")),
         formatted = sprintf("%.2f (%.2f, %.2f)", mean, ci_lower, ci_upper),
         metformin_flag = factor(metformin_flag, levels = c(0, 1), labels = c("SGA only", "SGA+Metformin"))) %>%
  arrange(timepoint)

# Plot
abs_graph_imp <- ggplot(final_results_long_abs, aes(x = timepoint, y = mean, color = metformin_flag, group = metformin_flag)) + 
  geom_line(size = 1.2) + 
  geom_point(size = 3) + 
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.1) + 
  labs( 
    y = "Body weight (kg)", 
    x = "Time-point") + 
  scale_color_manual(values = c("SGA only" = "#1f77b4", "SGA+Metformin" = "#ff7f0e")) + 
  theme_minimal(base_size = 14) + 
  theme( 
    legend.title = element_text(size = 16), 
    legend.text = element_text(size = 14), 
    axis.title = element_text(size = 16), 
    axis.text = element_text(size = 14), 
    plot.title = element_blank(), 
    panel.grid.major = element_line(size = 0.5, color = "gray85"), 
    panel.grid.minor = element_blank(), 
    panel.border = element_rect(color = "dark grey", fill = NA, size = 2), 
    legend.position = "bottom") + 
  guides(color = guide_legend(title = "Group"))

# Complete case

cc <- tte_imputed_long %>%
  select(.imp, .id, patid, trt_group, baseline_weightkg, sixm_weightkg, oney_weightkg, twoy_weightkg) %>%
  left_join(select(TTEcohort_final, patid, diff, metformin_flag, pcos_flag), by = "patid") %>%
  filter(.imp == 0 & !is.na(metformin_flag)) %>%
  mutate(keep = case_when(is.na(diff) | diff >= -30 & diff <= 90 ~ 1,
                          TRUE ~ 0)) %>%
  filter(keep == 1)

cc_f <- cc %>%
  pivot_longer(cols = ends_with("weightkg"), 
               names_to = "timepoint", 
               values_to = "weightkg") %>%
  group_by(metformin_flag, timepoint) %>%
  summarise(mean = mean(weightkg, na.rm = TRUE),
            se = sd(weightkg, na.rm = TRUE) / sqrt(sum(!is.na(weightkg))),  # Standard error
            ci_lower = mean - qt(0.975, df = sum(!is.na(weightkg)) - 1) * se,  # 95% CI lower bound
            ci_upper = mean + qt(0.975, df = sum(!is.na(weightkg)) - 1) * se,  # 95% CI upper bound
            n_obs = sum(!is.na(weightkg)),  # Number of non-missing observations
            .groups = "drop") %>%
  mutate(timepoint = case_when(timepoint == "baseline_weightkg" ~ "Baseline",
                               timepoint == "sixm_weightkg" ~ "6m",
                               timepoint == "oney_weightkg" ~ "1y",
                               timepoint == "twoy_weightkg" ~ "2y"),
         timepoint = factor(timepoint, levels = c("Baseline", "6m", "1y", "2y")),
         formatted = sprintf("%.2f (%.2f, %.2f)", mean, ci_lower, ci_upper),
         metformin_flag = factor(metformin_flag, levels = c(0, 1), labels = c("SGA only", "SGA+Metformin"))) %>%
  arrange(timepoint)

abs_graph_cc <- ggplot(cc_f, aes(x = timepoint, y = mean, color = metformin_flag, group = metformin_flag)) +
  geom_line(size = 1.2) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.1) +
  labs(
    y = "Body weight (kg)",
    x = "Time-point") +
  scale_color_manual(values = c("SGA only" = "#1f77b4", "SGA+Metformin" = "#ff7f0e")) +  # Custom colors
  theme_minimal(base_size = 14) +  # Larger base font size for publication
  theme(
    legend.title = element_text(size = 16),  # Change the legend title size
    legend.text = element_text(size = 14),   # Change the legend text size
    axis.title = element_text(size = 16),    # Axis titles size
    axis.text = element_text(size = 14),     # Axis text size
    plot.title = element_blank(),  # Remove the plot title
    panel.grid.major = element_line(size = 0.5, color = "gray85"),    # Lighter grid lines
    panel.grid.minor = element_blank(),     # No minor grid lines
    panel.border = element_rect(color = "dark grey", fill = NA, size = 2), 
    legend.position = "bottom") +  # Move the legend to the bottom
  guides(color = guide_legend(title = "Group"))  # Change legend title to "Group"

# View plots
abs_graph_imp
abs_graph_cc

# Individual SGAs

analyse_abs_weight_change <- function(data, tte_cohort, treatment_group) {
  # Imputed dataset analysis
  
  abs_change_imp <- tte_imputed_long %>%
    select(.imp, .id, patid, trt_group, baseline_weightkg, sixm_weightkg, oney_weightkg, twoy_weightkg) %>%
    filter(trt_group == treatment_group) %>%
    left_join(select(TTEcohort_final, patid, diff, metformin_flag, pcos_flag), by = "patid") %>%
    filter(.imp > 0 & !is.na(metformin_flag)) %>%
    mutate(keep = case_when(is.na(diff) | diff >= -30 & diff <= 90 ~ 1,
                            TRUE ~ 0)) %>%
    filter(keep == 1)
  
  # First calculate means and variances for each imputed dataset
  imp_stats <- abs_change_imp %>%
    group_by(.imp, metformin_flag) %>%
    summarize(n = n(),
              mean_baseline = mean(baseline_weightkg, na.rm = TRUE),
              var_baseline =  var(baseline_weightkg, na.rm = TRUE)/n,
              mean_6m = mean(sixm_weightkg, na.rm = TRUE),
              var_6m = var(sixm_weightkg, na.rm = TRUE)/n,  # Using SE^2
              mean_1y = mean(oney_weightkg, na.rm = TRUE),
              var_1y = var(oney_weightkg, na.rm = TRUE)/n,
              mean_2y = mean(twoy_weightkg, na.rm = TRUE),
              var_2y = var(twoy_weightkg, na.rm = TRUE)/n)
  
  # Revised function to apply Rubin's rules
  pool_estimates <- function(means, variances) {
    m <- length(means)  # number of imputations
    
    # Overall estimate (Q-bar)
    q_bar <- mean(means)
    
    # Within-imputation variance (U-bar)
    u_bar <- mean(variances)
    
    # Between-imputation variance (B)
    b <- sum((means - q_bar)^2) / (m - 1)
    
    # Total variance (T)
    t <- u_bar + (1 + 1/m) * b
    
    # Degrees of freedom using Barnard-Rubin adjustment
    lambda <- (b * (1 + 1/m)) / t
    df_old <- (m - 1) / lambda^2
    df_obs <- (1/lambda - 1) / (1 - lambda)  # approximate df from complete data
    df <- 1 / (1/df_old + 1/df_obs)
    
    # Calculate 95% CI
    ci_lower <- q_bar - qt(0.975, df) * sqrt(t)
    ci_upper <- q_bar + qt(0.975, df) * sqrt(t)
    
    return(c(mean = q_bar, 
             ci_lower = ci_lower, 
             ci_upper = ci_upper,
             se = sqrt(t)))}  # also returning SE for potential use
  
  # Apply Rubin's rules for each metformin group and timepoint
  results <- imp_stats %>%
    group_by(metformin_flag) %>%
    summarize(
      # baseline
      res_baseline = list(pool_estimates(mean_baseline, var_baseline)),
      # 6 months
      res_6m = list(pool_estimates(mean_6m, var_6m)),
      # 1 year
      res_1y = list(pool_estimates(mean_1y, var_1y)),
      # 2 years
      res_2y = list(pool_estimates(mean_2y, var_2y))) %>%
    mutate(
      # Baseline
      mean_baseline = map_dbl(res_baseline, ~.x["mean"]),
      ci_lower_baseline = map_dbl(res_baseline, ~.x["ci_lower"]),
      ci_upper_baseline = map_dbl(res_baseline, ~.x["ci_upper"]),
      # 6 months
      mean_6m = map_dbl(res_6m, ~.x["mean"]),
      ci_lower_6m = map_dbl(res_6m, ~.x["ci_lower"]),
      ci_upper_6m = map_dbl(res_6m, ~.x["ci_upper"]),
      # 1 year
      mean_1y = map_dbl(res_1y, ~.x["mean"]),
      ci_lower_1y = map_dbl(res_1y, ~.x["ci_lower"]),
      ci_upper_1y = map_dbl(res_1y, ~.x["ci_upper"]),
      # 2 years
      mean_2y = map_dbl(res_2y, ~.x["mean"]),
      ci_lower_2y = map_dbl(res_2y, ~.x["ci_lower"]),
      ci_upper_2y = map_dbl(res_2y, ~.x["ci_upper"])) %>%
    select(-res_baseline, -res_6m, -res_1y, -res_2y)
  
  final_results_long <- results %>%
    pivot_longer(cols = -metformin_flag,
                 names_to = c(".value", "timepoint"),
                 names_pattern = "(.+)_(.+)") %>%
    mutate(timepoint = factor(timepoint, levels = c("baseline", "6m", "1y", "2y")),
           formatted = sprintf("%.2f (%.2f, %.2f)", mean, ci_lower, ci_upper),
           metformin_flag = factor(metformin_flag, levels = c(0, 1), labels = c("SGA only", "SGA+Metformin"))) %>%
    arrange(timepoint)
  
  # Plot
  abs_graph_imp <- ggplot(final_results_long, aes(x = timepoint, y = mean, color = metformin_flag, group = metformin_flag)) + 
    geom_line(size = 1.2) + 
    geom_point(size = 3) + 
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.1) + 
    labs( 
      y = "Body weight (kg)", 
      x = "Time-point") + 
    scale_color_manual(values = c("SGA only" = "#1f77b4", "SGA+Metformin" = "#ff7f0e")) + 
    theme_minimal(base_size = 14) + 
    theme( 
      legend.title = element_text(size = 16), 
      legend.text = element_text(size = 14), 
      axis.title = element_text(size = 16), 
      axis.text = element_text(size = 14), 
      plot.title = element_blank(), 
      panel.grid.major = element_line(size = 0.5, color = "gray85"), 
      panel.grid.minor = element_blank(), 
      panel.border = element_rect(color = "dark grey", fill = NA, size = 2), 
      legend.position = "bottom") + 
    guides(color = guide_legend(title = "Group"))
  
  # Complete case
  
  abs_cc <- tte_imputed_long %>%
    select(.imp, .id, patid, trt_group, baseline_weightkg, sixm_weightkg, oney_weightkg, twoy_weightkg) %>%
    filter(trt_group == treatment_group) %>%
    left_join(select(TTEcohort_final, patid, diff, metformin_flag, pcos_flag), by = "patid") %>%
    filter(.imp == 0 & !is.na(metformin_flag)) %>%
    mutate(keep = case_when(is.na(diff) | diff >= -30 & diff <= 90 ~ 1,
                            TRUE ~ 0)) %>%
    filter(keep == 1)
  
  abs_cc_f <- abs_cc %>%
    pivot_longer(cols = ends_with("weightkg"), 
                 names_to = "timepoint", 
                 values_to = "weightkg") %>%
    group_by(metformin_flag, timepoint) %>%
    summarise(mean = mean(weightkg, na.rm = TRUE),
              se = sd(weightkg, na.rm = TRUE) / sqrt(sum(!is.na(weightkg))),  # Standard error
              ci_lower = mean - qt(0.975, df = sum(!is.na(weightkg)) - 1) * se,  # 95% CI lower bound
              ci_upper = mean + qt(0.975, df = sum(!is.na(weightkg)) - 1) * se,  # 95% CI upper bound
              n_obs = sum(!is.na(weightkg)),  # Number of non-missing observations
              .groups = "drop") %>%
    mutate(timepoint = case_when(timepoint == "baseline_weightkg" ~ "Baseline",
                                 timepoint == "sixm_weightkg" ~ "6m",
                                 timepoint == "oney_weightkg" ~ "1y",
                                 timepoint == "twoy_weightkg" ~ "2y"),
           timepoint = factor(timepoint, levels = c("Baseline", "6m", "1y", "2y")),
           formatted = sprintf("%.2f (%.2f, %.2f)", mean, ci_lower, ci_upper),
           metformin_flag = factor(metformin_flag, levels = c(0, 1), labels = c("SGA only", "SGA+Metformin"))) %>%
    arrange(timepoint)
  
  abs_graph_cc <- ggplot(abs_cc_f, aes(x = timepoint, y = mean, color = metformin_flag, group = metformin_flag)) +
    geom_line(size = 1.2) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0.1) +
    labs(
      y = "Body weight (kg)",
      x = "Time-point") +
    scale_color_manual(values = c("SGA only" = "#1f77b4", "SGA+Metformin" = "#ff7f0e")) +  # Custom colors
    theme_minimal(base_size = 14) +  # Larger base font size for publication
    theme(
      legend.title = element_text(size = 16),  # Change the legend title size
      legend.text = element_text(size = 14),   # Change the legend text size
      axis.title = element_text(size = 16),    # Axis titles size
      axis.text = element_text(size = 14),     # Axis text size
      plot.title = element_blank(),  # Remove the plot title
      panel.grid.major = element_line(size = 0.5, color = "gray85"),    # Lighter grid lines
      panel.grid.minor = element_blank(),     # No minor grid lines
      panel.border = element_rect(color = "dark grey", fill = NA, size = 2), 
      legend.position = "bottom") +  # Move the legend to the bottom
    guides(color = guide_legend(title = "Group"))  # Change legend title to "Group"
  
  # Return results
  return(list(
    imputed_results = results,
    imputed_plot = abs_graph_imp,
    complete_case_results = abs_cc_f,
    complete_case_plot = abs_graph_cc
  ))}

ari_abs_results <- analyse_abs_weight_change(tte_imputed_long, TTEcohort_final, "Aripiprazole")
ola_abs_results <- analyse_abs_weight_change(tte_imputed_long, TTEcohort_final, "Olanzapine")
que_abs_results <- analyse_abs_weight_change(tte_imputed_long, TTEcohort_final, "Quetiapine")
ris_abs_results <-  analyse_abs_weight_change(tte_imputed_long, TTEcohort_final, "Risperidone")

# Add titles to each plot with axis label removal
# Function to create consistent plot processing
process_plots <- function(plot_list, is_imputed = TRUE, ylim_imputed = c(70, 103), ylim_cc = c(70, 110)) {
  lapply(names(plot_list), function(med) {
    plot <- plot_list[[med]] + 
      labs(title = med)
    
    # Remove x-axis label for Aripiprazole and Olanzapine
    if (med %in% c("Aripiprazole", "Olanzapine")) {
      plot <- plot + 
        theme(axis.title.x = element_blank()) +
        theme(legend.position = "none")  # Remove legend for these plots
    }
    
    # Remove y-axis label for Olanzapine and Risperidone
    if (med %in% c("Olanzapine", "Risperidone")) {
      plot <- plot + theme(axis.title.y = element_blank())
    }
    
    plot + 
      theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
      coord_cartesian(ylim = if (is_imputed) ylim_imputed else ylim_cc)
  })
}

# Process imputed and complete case plots
abs_medication_plots <- list(
  Aripiprazole = ari_abs_results$imputed_plot,
  Olanzapine = ola_abs_results$imputed_plot,
  Quetiapine = que_abs_results$imputed_plot,
  Risperidone = ris_abs_results$imputed_plot)

# For imputed plots
abs_imputed_grid <- wrap_plots(process_plots(abs_medication_plots, is_imputed = TRUE), ncol = 2) + 
  plot_annotation(
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
  ) & 
  theme(legend.position = "bottom")

abs_cc_medication_plots <- list(
  Aripiprazole = ari_abs_results$complete_case_plot,
  Olanzapine = ola_abs_results$complete_case_plot,
  Quetiapine = que_abs_results$complete_case_plot,
  Risperidone = ris_abs_results$complete_case_plot)

abs_cc_grid <- wrap_plots(process_plots(abs_cc_medication_plots, is_imputed = FALSE), ncol = 2) + 
  plot_annotation(
    theme = theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 16))
  ) & 
  theme(legend.position = "bottom")


# Save plots
ggsave("Organised R scripts/pct_change_imp.tiff", graph_imp, width = 10, height = 8, background = "white")
ggsave("Organised R scripts/pct_change_imp.pdf", graph_imp, width = 10, height = 8)

ggsave("Organised R scripts/pct_change_cc.tiff", graph_cc, width = 10, height = 8, background = "white")
ggsave("Organised R scripts/pct_change_cc.pdf", graph_cc, width = 10, height = 8)

ggsave("Organised R scripts/pct_change_stratified_imp.tiff", pct_imputed_grid, width = 12, height = 9, background = "white")
ggsave("Organised R scripts/pct_change_stratified_imp.pdf", pct_imputed_grid, width = 12, height = 9)

ggsave("Organised R scripts/pct_change_stratified_cc.tiff", pct_cc_grid, width = 12, height = 9, background = "white")
ggsave("Organised R scripts/pct_change_stratified_cc.pdf", pct_cc_grid, width = 12, height = 9)

ggsave("Organised R scripts/abs_change_imp.tiff", abs_graph_imp, width = 10, height = 8, background = "white")
ggsave("Organised R scripts/abs_change_imp.pdf", abs_graph_imp, width = 10, height = 8)

ggsave("Organised R scripts/abs_change_cc.tiff", abs_graph_cc, width = 10, height = 8, background = "white")
ggsave("Organised R scripts/abs_change_cc.pdf", abs_graph_cc, width = 10, height = 8)

ggsave("Organised R scripts/abs_change_stratified_imp.tiff", abs_imputed_grid, width = 12, height = 9, background = "white")
ggsave("Organised R scripts/abs_change_stratified_imp.pdf", abs_imputed_grid, width = 12, height = 9)

ggsave("Organised R scripts/abs_change_stratified_cc.tiff", abs_cc_grid, width = 12, height = 9, background = "white")
ggsave("Organised R scripts/abs_change_stratified_cc.pdf", abs_cc_grid, width = 12, height = 9)

overall_combined <- abs_graph_imp + graph_imp
ggsave("Organised R scripts/pctandabs_change_combined_imp.tiff", overall_combined, width = 12, height = 7.5, background = "white")
ggsave("Organised R scripts/pctandabs_change_combined_imp.pdf", overall_combined, width = 12, height = 7.5)

overall_combined_2 <- abs_graph_cc + graph_cc
ggsave("Organised R scripts/pctandabs_change_combined_cc.tiff", overall_combined_2, width = 12, height = 7.5, background = "white")
ggsave("Organised R scripts/pctandabs_change_combined_cc.pdf", overall_combined_2, width = 12, height = 7.5)



# Table
table1 <- final_results_long_abs %>%
  select(metformin_flag, timepoint, abs = formatted)
table2 <- final_results_long %>%
  select(metformin_flag, timepoint, pct = formatted)

table <- table1 %>%
  left_join(table2, by = c("metformin_flag", "timepoint")) %>%
  pivot_wider(
    names_from = metformin_flag,
    values_from = c(abs, pct),
    names_sep = "_") %>%
  clean_names() %>%
  select(timepoint, abs_sga_only, pct_sga_only, abs_sga_metformin, pct_sga_metformin) %>%
  gt(rowname_col = "timepoint") %>%
  cols_label(
    abs_sga_only = "Body weight (kg)",
    abs_sga_metformin = "Body weight (kg)",
    pct_sga_only = "Change (%)",
    pct_sga_metformin = "Change (%)"
  ) %>%
  cols_align(
    align = "center", # center allign columns
    columns = everything()) %>%
  tab_spanner(
    label = "SGA Only",
    columns = c(abs_sga_only, pct_sga_only)
  ) %>%
  tab_spanner(
    label = "SGA+Metformin",
    columns = c(abs_sga_metformin, pct_sga_metformin)
  ) %>%
  sub_missing(missing_text = "-") %>% # use - for empty cells
  tab_style(
    style = cell_borders(sides = "bottom", color = "black", weight = px(1)),
    locations = cells_body(columns = everything(), rows = everything())
  ) %>%
  tab_options(
    table.font.size = "small",
    table.border.top.color = "black",
    table.border.bottom.color = "black",
    heading.align = "center"
  ) %>%
  tab_footnote("Values are mean with 95% confidence intervals. Missing data was imputed using multiple imputation. Estimates were pooled according to Rubin's rules. Change is the percentage change from baseline.")

table %>%
  gtsave(filename = paste0("pctandabs_weight_results_table_", today(), ".docx"), 
         path = "Organised R scripts/")

# Packages version used in this script
# dplyr version: 1.1.4 
# tidyr version: 1.3.1 
# ggplot2 version: 3.5.1 
# gridExtra version: 2.3 
# grid version: 4.4.1 
# purrr version: 1.0.2
# patchwork version: 1.3.0
# janitor version: 2.2.0

