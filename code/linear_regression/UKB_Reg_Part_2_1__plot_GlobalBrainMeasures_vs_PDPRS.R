# Script to create visualization of global measures regression results
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

# Define the base path for input files
base_path <- "Outputs/PD_PRS_Regression_Results/with_all_pathwayPRS_Brain_GlobalMeasures/with_all_pathwayPRS_Brain_GlobalMeasures_withMotion_"

# Loop through each sex group
for(sex in sex_groups) {
  # Construct the full file path
  file_path <- paste0(base_path, sex, "/Global_Measures_All_PRS_Regression_Results.csv")
  
  # Read the data
  global_results <- read.csv(file_path)
  
  
  # Keep important PRS values only
  global_results <- global_results %>% filter(str_starts(PRS, "PD_PRScs_GP2")) %>% 
    filter(Measure %in% c('Total_SA_DKT', 'Total_CT_DKT', 'Total_Subcortical_Volume', 'Total_FA', 'Total_MD', 'Total_SWI'))
  
  ### Prepare data for plotting
  plot_data <- global_results %>%
    mutate(
      # Clean up PRS names by removing the PD_PRScs_GP2 prefix
      PRS = str_replace(PRS, "PD_PRScs_GP2_", ""),
      # Clean up Measure names by removing 'Total_' and all underscores
      Measure = str_replace(Measure, "Total_", ""),
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
    "SA DKT" = "Surface Area",
    "CT DKT" = "Cortical Thickness",
    "Subcortical Volume" = "Subcortical Volume",
    "FA" = "Fractional Anisotropy (FA)",
    "MD" = "Mean Diffusivity (MD)",
    "SWI" = "SWI"
  )
  
  # Apply the mappings
  plot_data <- plot_data %>%
    mutate(
      PRS = recode(PRS, !!!prs_mapping),
      Measure = recode(Measure, !!!measure_mapping),
      # Create a significance indicator
      sig_indicator = case_when(
        p_FDR < 0.001 ~ "***", # or p_value
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
    "Subcortical Volume",
    "Surface Area",
    "Cortical Thickness",
    "Fractional Anisotropy (FA)",
    "Mean Diffusivity (MD)",
    "SWI"
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
  # OR set manually
  max_limit <- 0.10
  
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
    # # Use a diverging color scale similar to Python's RdYlBu_r
    # scale_fill_gradientn(
    #   colors = c("#313695", "#4575b4", "#74add1", "#abd9e9", "#e0f3f8", 
    #              "#ffffbf", "#fee090", "#fdae61", "#f46d43", "#d73027", "#a50026"),
    #   name = "Effect Size",
    #   limits = c(-max_limit, max_limit)
    # ) +
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
      title = paste0("Regression Results of PD PRS on Global Brain Measures (", sex, ")"),
      subtitle = "Effect sizes with FDR significance indicators (* p<0.05, ** p<0.01, *** p<0.001)" # with or without FDR
    )
  
  
  
  ### Saving
  # Create output directory if it doesn't exist
  output_dir <- "Outputs/PD_PRS_Regression_Results/plots/with_all_pathwayPRS_Brain_GlobalMeasures/"
  if (!dir.exists(output_dir)) { dir.create(output_dir, recursive = TRUE) }
  
  # Save as SVG with sex in filename
  svg_filename <- paste0(output_dir, "Global_Measures_vs_PD_PRS_Heatmap_", sex, ".svg")
  ggsave(svg_filename, plot = p, width = 7, height = 6)
  
  # Save as PNG with sex in filename
  png_filename <- paste0(output_dir, "Global_Measures_vs_PD_PRS_Heatmap_", sex, ".png")
  ggsave(png_filename, plot = p, width = 7, height = 6, dpi = 300)
    
}














########################## HORIZONTAL BAR PLOTTING ##########################
# Define the sex groups to process
sex_groups <- c("all", "male", "female")

# Define the base path for input files
base_path <- "Outputs/PD_PRS_Regression_Results/with_all_pathwayPRS_Brain_GlobalMeasures/with_all_pathwayPRS_Brain_GlobalMeasures_withMotion_"

# Loop through each sex group
for(sex in sex_groups) {
  # Construct the full file path
  file_path <- paste0(base_path, sex, "/Global_Measures_All_PRS_Regression_Results.csv")
  
  # Read the data
  global_results <- read.csv(file_path)
  
  # Filter for only PD_PRScs_GP2 (exact match) and the specific measures
  plot_data <- global_results %>% 
    filter(PRS == "PD_PRScs_GP2") %>% 
    filter(Measure %in% c('Total_SA_DKT', 'Total_CT_DKT', 'Total_Subcortical_Volume', 'Total_FA', 'Total_MD', 'Total_SWI'))
  
  # Clean up Measure names by removing 'Total_' and all underscores
  plot_data <- plot_data %>%
    mutate(Measure = str_replace(Measure, "Total_", ""),
           Measure = str_replace_all(Measure, "_", " "))
  
  # Manual renaming of Measure values
  measure_mapping <- c(
    "SA DKT" = "Surface Area",
    "CT DKT" = "Cortical Thickness",
    "Subcortical Volume" = "Subcortical Volume",
    "FA" = "Fractional Anisotropy (FA)",
    "MD" = "Mean Diffusivity (MD)",
    "SWI" = "SWI"
  )
  
  # Apply the measure mapping
  plot_data <- plot_data %>%
    mutate(Measure = recode(Measure, !!!measure_mapping))
  
  # Define the order for Measures
  measure_order <- c(
    "Subcortical Volume",
    "Surface Area",
    "Cortical Thickness",
    "Fractional Anisotropy (FA)",
    "Mean Diffusivity (MD)",
    "SWI"
  )
  
  # Convert Measure to factor with the specified order
  plot_data <- plot_data %>%
    mutate(Measure = factor(Measure, levels = rev(measure_order)))
  
  # Determine if the confidence interval includes zero
  plot_data <- plot_data %>%
    mutate(
      significant = (CI_low > 0 & effect_size > 0) | (CI_high < 0 & effect_size < 0)
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
  
  # Create the horizontal bar plot
  p <- ggplot(plot_data, aes(x = effect_size, y = Measure)) +
    # Add a vertical line at x=0
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray40", linewidth = 0.5) +
    # Add the bars
    geom_bar(stat = "identity", aes(fill = color), width = 0.6) +
    # Add error bars using the actual CI bounds with shorter cap width (reduced from 0.3 to 0.15)
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
      title = paste0("Effect of PD PRS on Global Brain Measures (", sex, ")"),
      subtitle = paste0("PD_PRScs_GP2 regression coefficients with 95% confidence intervals")
    )
  
  ### Saving
  # Create output directory if it doesn't exist
  output_dir <- "Outputs/PD_PRS_Regression_Results/plots/with_all_pathwayPRS_Brain_GlobalMeasures/"
  if (!dir.exists(output_dir)) { dir.create(output_dir, recursive = TRUE) }
  
  # Save as SVG with sex in filename
  svg_filename <- paste0(output_dir, "PD_PRScs_Effect_on_Global_Measures_", sex, ".svg")
  ggsave(svg_filename, plot = p, width = 7, height = 4.5)
  
  # Save as PNG with sex in filename
  png_filename <- paste0(output_dir, "PD_PRScs_Effect_on_Global_Measures_", sex, ".png")
  ggsave(png_filename, plot = p, width = 7, height = 4.5, dpi = 300)
}



