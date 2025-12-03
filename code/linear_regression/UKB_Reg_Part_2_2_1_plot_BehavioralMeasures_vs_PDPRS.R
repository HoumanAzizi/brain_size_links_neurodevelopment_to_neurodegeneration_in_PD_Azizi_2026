# Script to create visualization of behavioral measures regression results
library(tidyverse)
library(ggplot2)
library(reshape2)
library(viridis)
library(scales)

rm(list = ls())
cat("\014")

# Set working directory
setwd("/Users/houmanazizi/Library/Mobile\ Documents/com~apple~CloudDocs/Education/GitHub/UKB_Analysis/")











########################## HEATMAP PLOTTING ##########################
### Analysis
# Define the sex groups to process
sex_groups <- c("all", "male", "female")

# Define the base path for input files and output path
plot_type <- 'imaging' # 'fullUKB' or 'imaging'
base_path <- paste0("Outputs/PD_PRS_Regression_Results/with_all_pathwayPRS_Behavioral_", plot_type, "/using_Linear_regression/Behavioral_Measures_All_PRS_Regression_Results_")
output_dir <- paste0("Outputs/PD_PRS_Regression_Results/plots/with_all_pathwayPRS_Behavioral_", plot_type, "/using_Linear_regression/")
if (!dir.exists(output_dir)) { dir.create(output_dir, recursive = TRUE) }

# Loop through each sex group
for(sex in sex_groups) {
  # Construct the full file path
  file_path <- paste0(base_path, sex, ".csv")
  
  # Read the data
  behavioral_results <- read.csv(file_path)
  
  # Keep important PRS values only
  behavioral_results <- behavioral_results %>% filter(str_starts(PRS, "PD_PRScs_GP2"))
  
  ### Prepare data for plotting
  plot_data <- behavioral_results %>%
    mutate(
      # Clean up PRS names by removing the PD_PRScs_GP2 prefix
      PRS = str_replace(PRS, "PD_PRScs_GP2_", ""),
      # Clean up Measure names by replacing underscores with spaces
      Measure = str_replace_all(Measure, "_", " ")) %>% 
    filter(!(PRS %in% c('Autophagy_direct_genes_4_exclude', 'Lysosomal_direct_genes_2_exclude', 'Mitochondrial_genes_1_exclude')))
  
  ### Manually rename values
  # Manual renaming of PRS values
  prs_mapping <- c(
    "PD_PRScs_GP2" = "all genes (PD-PRS)",
    "Autophagy_direct_genes_4" = "autophagy pathway",
    "Lysosomal_direct_genes_2" = "lysosomal pathway",
    "Mitochondrial_genes_1" = "mitochondrial pathway",
    "Mitochondrial_Lysosomal_Autophagy_direct_genes_9" = "three pathways combined",
    "Mitochondrial_Lysosomal_Autophagy_direct_genes_9_exclude" = "non-pathways genes"
  )
  
  # Manual renaming of Measure values
  measure_mapping <- c(
    "Coffee intake" = "Coffee intake",
    "Fluid intelligence score" = "Fluid Intelligence",
    "BMI" = "BMI",
    "Sleep_duration" = "Sleep duration",
    "Multiple deprivation index" = "Multiple Deprivation Index",
    "Smoking PackYears" = "Smoking (Pack-Years)",
    "Household income" = "Household Income",
    "Alcohol usage total" = "Alcohol Usage",
    "Education level" = "Education Level",
    "Addiction score total" = "Addiction Score",
    "Constipation" = 'Constipation', 
    "Orthostatic hypotension" = 'Orthostatic Hypotension', 
    "Urinary incontinence" = "Urinary Incontinence", 
    "Hand grip strength" = "Hand Grip Strength", 
    "Apathy past2weeks" = "Apathy (Last 2 Weeks)", 
    "Depressed mood past2weeks" = 'Depressed Mood (Last 2 Weeks)'
  )
  
  # Apply the mappings
  plot_data <- plot_data %>%
    mutate(
      PRS = recode(PRS, !!!prs_mapping),
      Measure = recode(Measure, !!!measure_mapping),
      # Create a significance indicator
      sig_indicator = case_when(
        p_FDR < 0.001 ~ "***",
        p_FDR < 0.01 ~ "**",
        p_FDR < 0.05 ~ "*",
        TRUE ~ ""))
  
  # Define the order for PRS (Y-axis)
  prs_order <- c(
    "all genes (PD-PRS)",
    "non-pathways genes",
    "three pathways combined",
    "mitochondrial pathway",
    "lysosomal pathway",
    "autophagy pathway"
  )
  
  # Define the order for Measures (X-axis)
  measure_order <- c(
    "Education Level",
    "Household Income",
    "Fluid Intelligence",
    "Coffee intake",
    "Sleep duration",
    "BMI",
    "Hand Grip Strength",
    "Multiple Deprivation Index",
    "Smoking (Pack-Years)",
    "Alcohol Usage",
    "Addiction Score",
    "Constipation",
    "Orthostatic Hypotension",
    "Urinary Incontinence",
    "Apathy (Last 2 Weeks)",
    "Depressed Mood (Last 2 Weeks)"
  )
  
  # Convert PRS and Measure to factors with the specified order
  plot_data <- plot_data %>%
    mutate(
      PRS = factor(PRS, levels = rev(prs_order)),
      Measure = factor(Measure, levels = measure_order)
    )
  
  ### Plot as heatmap
  # Find the maximum absolute effect size and round up to the next convenient value
  max_abs_effect <- max(abs(plot_data$effect_size))
  # Round up to the next 0.01 increment
  max_limit <- ceiling(max_abs_effect * 100) / 100
  # OR set manually if needed (uncomment to use)
  # max_limit <- 0.10
  
  # Create the heatmap
  p <- ggplot(plot_data, aes(x = Measure, y = PRS)) +
    # Use effect_size for the color mapping
    geom_tile(aes(fill = effect_size), color = "white", size = 0.5) +
    # Add significance indicators
    geom_text(aes(label = sig_indicator), 
              color = "black", size = 4) +
    # Use a diverging color scale similar to Python's RdYlBu but with white at zero
    scale_fill_gradientn(
      colors = c("#313695", "#4575b4", "#74add1", "#abd9e9", "#e0f3f8", 
                 "#FFFBF3",  # White at center
                 "#fee090", "#fdae61", "#f46d43", "#d73027", "#a50026"),
      values = scales::rescale(c(-max_limit, -max_limit*0.8, -max_limit*0.6, -max_limit*0.4, -max_limit*0.2, 
                                 0,  # Center point
                                 max_limit*0.2, max_limit*0.4, max_limit*0.6, max_limit*0.8, max_limit)),
      name = "Effect Size",
      limits = c(-max_limit, max_limit)
    ) +
    # Customize theme
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
      axis.text.y = element_text(size = 10),
      axis.title = element_blank(),
      panel.grid = element_blank(),
      legend.position = "right",
      plot.title = element_text(size = 12, face = "bold"),
      plot.subtitle = element_text(size = 10)
    ) +
    # Add informative title
    labs(
      title = paste0("Regression Results of PD PRS on Behavioral Measures (", sex, ")"),
      subtitle = "Effect sizes with significance indicators (* p<0.05, ** p<0.01, *** p<0.001 - FDR)"
    )
  
  ### Saving
  # Save as SVG with sex in filename
  svg_filename <- paste0(output_dir, "Behavioral_Measures_vs_PD_PRS_Heatmap_", sex, ".svg")
  ggsave(svg_filename, plot = p, width = 13, height = 6)
  
  # Save as PNG with sex in filename
  png_filename <- paste0(output_dir, "Behavioral_Measures_vs_PD_PRS_Heatmap_", sex, ".png")
  ggsave(png_filename, plot = p, width = 13, height = 6, dpi = 300)
}














########################## HORIZONTAL BAR PLOTTING FOR BEHAVIORAL MEASURES ##########################
# Define the sex groups to process
sex_groups <- c("all", "male", "female")

# Define the base path for input files and output path
#plot_type <- 'imaging' # 'fullUKB' or 'imaging'
base_path <- paste0("Outputs/PD_PRS_Regression_Results/with_all_pathwayPRS_Behavioral_", plot_type, "/using_Linear_regression/Behavioral_Measures_All_PRS_Regression_Results_")
output_dir <- paste0("Outputs/PD_PRS_Regression_Results/plots/with_all_pathwayPRS_Behavioral_", plot_type, "/using_Linear_regression/")
if (!dir.exists(output_dir)) { dir.create(output_dir, recursive = TRUE) }

# Loop through each sex group
for(sex in sex_groups) {
  # Construct the full file path
  file_path <- paste0(base_path, sex, ".csv")
  
  # Read the data
  behavioral_results <- read.csv(file_path)
  
  # Filter for only PD_PRScs_GP2 (exact match)
  plot_data <- behavioral_results %>% 
    filter(PRS == "PD_PRScs_GP2")
  
  # Clean up Measure names by replacing underscores with spaces
  plot_data <- plot_data %>%
    mutate(Measure = str_replace_all(Measure, "_", " "))
  
  # Manual renaming of Measure values
  measure_mapping <- c(
    "Coffee intake" = "Coffee intake",
    "Fluid intelligence score" = "Fluid Intelligence",
    "BMI" = "BMI",
    "Sleep_duration" = "Sleep duration",
    "Multiple deprivation index" = "Multiple Deprivation Index",
    "Smoking PackYears" = "Smoking (Pack-Years)",
    "Household income" = "Household Income",
    "Alcohol usage total" = "Alcohol Usage",
    "Education level" = "Education Level",
    "Addiction score total" = "Addiction Score",
    "Constipation" = 'Constipation', 
    "Orthostatic hypotension" = 'Orthostatic Hypotension', 
    "Urinary incontinence" = "Urinary Incontinence", 
    "Hand grip strength" = "Hand Grip Strength", 
    "Apathy past2weeks" = "Apathy (Last 2 Weeks)", 
    "Depressed mood past2weeks" = 'Depressed Mood (Last 2 Weeks)'
  )
  
  # Apply the measure mapping
  plot_data <- plot_data %>%
    mutate(Measure = recode(Measure, !!!measure_mapping))
  
  # Define the order for Measures
  measure_order <- c(
    "Education Level",
    "Household Income",
    "Fluid Intelligence",
    "Coffee intake",
    "Sleep duration",
    "BMI",
    "Hand Grip Strength",
    "Multiple Deprivation Index",
    "Smoking (Pack-Years)",
    "Alcohol Usage",
    "Addiction Score",
    "Constipation",
    "Orthostatic Hypotension",
    "Urinary Incontinence",
    "Apathy (Last 2 Weeks)",
    "Depressed Mood (Last 2 Weeks)"
  )
  
  # Convert Measure to factor with the specified order
  plot_data <- plot_data %>%
    mutate(Measure = factor(Measure, levels = rev(measure_order)))
  
  # Determine if the confidence interval includes zero
  plot_data <- plot_data %>%
    mutate(
      significant = ((CI_low > 0 & effect_size > 0) | (CI_high < 0 & effect_size < 0)) & (p_FDR <= 0.05)
    )
  
  # Create a color vector based on significance and direction
  plot_data <- plot_data %>%
    mutate(
      color = case_when(
        significant & effect_size > 0 ~ "#d62727",  # Red for significant positive
        significant & effect_size < 0 ~ "#2077b4",  # Blue for significant negative
        TRUE ~ "gray60"                            # Gray for non-significant
      )
    )
  
  # Find the maximum absolute value considering CI bounds
  max_abs_value <- max(abs(c(plot_data$CI_low, plot_data$CI_high)))
  
  # Round up to the next 0.01 increment
  max_limit <- ceiling(max_abs_value * 100) / 100
  
  # Set manual limit for consistency across sexes
  max_limit <- 0.06
  
  # Create the horizontal bar plot
  p <- ggplot(plot_data, aes(x = effect_size, y = Measure)) +
    # Add a vertical line at x=0
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5) +
    # Add the bars
    geom_bar(stat = "identity", aes(fill = color), width = 0.6) +
    # Add error bars using the actual CI bounds
    geom_errorbarh(aes(xmin = CI_low, xmax = CI_high), 
                   height = 0.15, color = "black", linewidth = 0.5) +
    # Set the color scale
    scale_fill_identity() +
    # Set x-axis limits symmetrically
    xlim(-max_limit, max_limit) +
    # Customize theme
    theme_minimal() +
    theme(
      axis.text.y = element_text(size = 10),
      axis.title.y = element_blank(),
      axis.title.x = element_text(size = 12),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "none",
      plot.title = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 10),
      # Add a black box around the plot panel
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      # Remove the default gray panel background
      panel.background = element_blank()
    ) +
    # Add labels
    labs(
      x = "Effect Size (β with 95% CI)",
      title = paste0("Effect of PD PRS on Behavioral Measures (", sex, ")"),
      subtitle = paste0("PD_PRScs_GP2 regression coefficients with 95% confidence intervals")
    )
  
  ### Saving
  # Save as SVG with sex in filename
  svg_filename <- paste0(output_dir, "PD_PRScs_Effect_on_Behavioral_Measures_", sex, ".svg")
  ggsave(svg_filename, plot = p, width = 7, height = 8)
  
  # Save as PNG with sex in filename
  png_filename <- paste0(output_dir, "PD_PRScs_Effect_on_Behavioral_Measures_", sex, ".png")
  ggsave(png_filename, plot = p, width = 7, height = 8, dpi = 300)
}




















########################## COMPARISON PLOT: fullUKB vs imaging BETA VALUES ##########################
# Define the sex groups to process
sex_groups <- c("all", "male", "female")

# Define the output path (save in fullUKB folder)
output_dir <- "Outputs/PD_PRS_Regression_Results/plots/with_all_pathwayPRS_Behavioral_fullUKB/using_Linear_regression/"
if (!dir.exists(output_dir)) { dir.create(output_dir, recursive = TRUE) }

# Loop through each sex group
for(sex in sex_groups) {
  # Read fullUKB data
  fullUKB_path <- paste0("Outputs/PD_PRS_Regression_Results/with_all_pathwayPRS_Behavioral_fullUKB/using_Linear_regression/Behavioral_Measures_All_PRS_Regression_Results_", sex, ".csv")
  fullUKB_data <- read.csv(fullUKB_path)
  
  # Read imaging data
  imaging_path <- paste0("Outputs/PD_PRS_Regression_Results/with_all_pathwayPRS_Behavioral_imaging/using_Linear_regression/Behavioral_Measures_All_PRS_Regression_Results_", sex, ".csv")
  imaging_data <- read.csv(imaging_path)
  
  # Filter for only PD_PRScs_GP2 and clean measure names
  fullUKB_filtered <- fullUKB_data %>% 
    filter(PRS == "PD_PRScs_GP2") %>%
    mutate(Measure = str_replace_all(Measure, "_", " ")) %>%
    select(Measure, effect_size, p_FDR, CI_low, CI_high) %>%
    rename(beta_fullUKB = effect_size, 
           p_FDR_fullUKB = p_FDR,
           CI_low_fullUKB = CI_low,
           CI_high_fullUKB = CI_high)
  
  imaging_filtered <- imaging_data %>% 
    filter(PRS == "PD_PRScs_GP2") %>%
    mutate(Measure = str_replace_all(Measure, "_", " ")) %>%
    select(Measure, effect_size, p_FDR, CI_low, CI_high) %>%
    rename(beta_imaging = effect_size, 
           p_FDR_imaging = p_FDR,
           CI_low_imaging = CI_low,
           CI_high_imaging = CI_high)
  
  # Manual renaming of Measure values (same as before)
  measure_mapping <- c(
    "Coffee intake" = "Coffee intake",
    "Fluid intelligence score" = "Fluid Intelligence",
    "BMI" = "BMI",
    "Sleep_duration" = "Sleep duration",
    "Multiple deprivation index" = "Multiple Deprivation Index",
    "Smoking PackYears" = "Smoking (Pack-Years)",
    "Household income" = "Household Income",
    "Alcohol usage total" = "Alcohol Usage",
    "Education level" = "Education Level",
    "Addiction score total" = "Addiction Score",
    "Constipation" = 'Constipation', 
    "Orthostatic hypotension" = 'Orthostatic Hypotension', 
    "Urinary incontinence" = "Urinary Incontinence", 
    "Hand grip strength" = "Hand Grip Strength", 
    "Apathy past2weeks" = "Apathy (Last 2 Weeks)", 
    "Depressed mood past2weeks" = 'Depressed Mood (Last 2 Weeks)'
  )
  
  # Apply measure mapping to both datasets
  fullUKB_filtered <- fullUKB_filtered %>%
    mutate(Measure = recode(Measure, !!!measure_mapping))
  
  imaging_filtered <- imaging_filtered %>%
    mutate(Measure = recode(Measure, !!!measure_mapping))
  
  # Merge the datasets
  comparison_data <- inner_join(fullUKB_filtered, imaging_filtered, by = "Measure")
  
  # Create significance categories for color coding
  comparison_data <- comparison_data %>%
    mutate(
      significance_category = case_when(
        p_FDR_fullUKB < 0.05 & p_FDR_imaging < 0.05 ~ "Both Significant",
        p_FDR_fullUKB < 0.05 & p_FDR_imaging >= 0.05 ~ "FullUKB Only",
        p_FDR_fullUKB >= 0.05 & p_FDR_imaging < 0.05 ~ "Imaging Only",
        TRUE ~ "Neither Significant"
      )
    )
  
  # Define colors for significance categories
  sig_colors <- c(
    "Both Significant" = "#d62727",      # Red
    "FullUKB Only" = "#ff7f0e",          # Orange
    "Imaging Only" = "#2ca02c",          # Green
    "Neither Significant" = "gray60"      # Gray
  )
  
  # Calculate correlation and linear regression
  correlation <- cor(comparison_data$beta_fullUKB, comparison_data$beta_imaging, use = "complete.obs")
  lm_model <- lm(beta_imaging ~ beta_fullUKB, data = comparison_data)
  
  # Determine axis limits (symmetric around 0)
  max_abs_beta <- max(abs(c(comparison_data$beta_fullUKB, comparison_data$beta_imaging)), na.rm = TRUE)
  axis_limit <- ceiling(max_abs_beta * 100) / 100
  
  # Create the scatter plot
  p <- ggplot(comparison_data, aes(x = beta_fullUKB, y = beta_imaging)) +
    # Add diagonal reference line (y = x)
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5) +
    # Add vertical and horizontal lines at 0
    geom_vline(xintercept = 0, linetype = "dotted", color = "gray40", linewidth = 0.3) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "gray40", linewidth = 0.3) +
    # Add error bars for fullUKB (horizontal)
    geom_errorbarh(aes(xmin = CI_low_fullUKB, xmax = CI_high_fullUKB, color = significance_category), 
                   height = 0, alpha = 0.6, linewidth = 0.5) +
    # Add error bars for imaging (vertical)
    geom_errorbar(aes(ymin = CI_low_imaging, ymax = CI_high_imaging, color = significance_category), 
                  width = 0, alpha = 0.6, linewidth = 0.5) +
    # Add points
    geom_point(aes(color = significance_category), size = 3, alpha = 0.8) +
    # Add line of best fit
    geom_smooth(method = "lm", se = TRUE, color = "black", linewidth = 1, alpha = 0.3) +
    # Add measure labels
    geom_text(aes(label = Measure), size = 2.5, nudge_y = axis_limit * 0.03, check_overlap = TRUE) +
    # Set colors
    scale_color_manual(values = sig_colors, name = "Significance\n(FDR < 0.05)") +
    # Set axis limits
    xlim(-axis_limit, axis_limit) +
    ylim(-axis_limit, axis_limit) +
    # Customize theme
    theme_minimal() +
    theme(
      axis.title = element_text(size = 12),
      legend.position = "right",
      plot.title = element_text(size = 13, face = "bold"),
      plot.subtitle = element_text(size = 10),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      panel.background = element_blank(),
      legend.title = element_text(size = 10),
      legend.text = element_text(size = 9)
    ) +
    # Add labels
    labs(
      x = "Effect Size (β) - Full UKB",
      y = "Effect Size (β) - Imaging Subset",
      title = paste0("Comparison of PD PRS Effects: Full UKB vs Imaging Subset (", sex, ")"),
      subtitle = paste0("PD_PRScs_GP2 regression coefficients with 95% CI | r = ", 
                        round(correlation, 3), " | R² = ", round(summary(lm_model)$r.squared, 3))
    )
  
  ### Saving
  # Save as SVG
  svg_filename <- paste0(output_dir, "FullUKB_vs_Imaging_Beta_Comparison_", sex, ".svg")
  ggsave(svg_filename, plot = p, width = 10, height = 8)
  
  # Save as PNG
  png_filename <- paste0(output_dir, "FullUKB_vs_Imaging_Beta_Comparison_", sex, ".png")
  ggsave(png_filename, plot = p, width = 10, height = 8, dpi = 300)
  
  # Print correlation and R-squared for this sex group
  cat("Sex group:", sex, "\n")
  cat("Correlation:", round(correlation, 3), "\n")
  cat("R-squared:", round(summary(lm_model)$r.squared, 3), "\n\n")
}
