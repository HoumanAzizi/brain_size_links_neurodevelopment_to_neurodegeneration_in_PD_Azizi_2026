# Script to create heatmap visualization of ENIGMA PD correlations
library(tidyverse)
library(ggplot2)
library(reshape2)
library(viridis)
library(scales)

rm(list = ls())
cat("\014")

# Set working directory
setwd("/Users/houmanazizi/Library/Mobile\ Documents/com~apple~CloudDocs/Education/GitHub/ENIGMA_PD_Correlation/")

########################## HEATMAP VISUALIZATION ##########################
# Define PRS measures to process
prs_measures <- c("SA", "CT")

# Process each PRS measure
for(prs_measure in prs_measures) {
  # Define file paths for this measure
  enigma_sa_file <- paste0('Outputs/', prs_measure, '_PRS_Results/SpinTest_ENIGMA_SA_vs_', prs_measure, '_PRS_Results.csv')
  enigma_ct_file <- paste0('Outputs/', prs_measure, '_PRS_Results/SpinTest_ENIGMA_CT_vs_', prs_measure, '_PRS_Results.csv')
  
  # Create output directory for this measure
  output_dir <- paste0("Plots/heatmaps/", prs_measure, "_PRS/")
  if (!dir.exists(output_dir)) { dir.create(output_dir, recursive = TRUE) }
  
  # Process each ENIGMA measure type (SA and CT)
  for(enigma_measure in c("SA", "CT")) {
    # Select the appropriate file based on ENIGMA measure type
    file_path <- if(enigma_measure == "SA") enigma_sa_file else enigma_ct_file
    
    # Read the ENIGMA results
    enigma_data <- read.csv(file_path)
    
    # Keep only spearman correlations
    enigma_data <- enigma_data %>% filter(Correlation == 'spearman') %>% select(-Correlation)
    
    ########################## PREPROCESSING ##########################
    # Filter to keep only specific PRS types that match the original analysis
    # We'll filter and rename in one step using case_when
    plot_data <- enigma_data %>%
      mutate(
        # Extract the core PRS type for filtering
        PRS_type = case_when(
          str_detect(PD_PRS, "PD_PRScs_GP2$") ~ "all_genes",
          str_detect(PD_PRS, "PD_PRScs_GP2_Autophagy_direct_genes_4$") ~ "autophagy",
          str_detect(PD_PRS, "PD_PRScs_GP2_Lysosomal_direct_genes_2$") ~ "lysosomal",
          str_detect(PD_PRS, "PD_PRScs_GP2_Mitochondrial_genes_1$") ~ "mitochondrial",
          str_detect(PD_PRS, "PD_PRScs_GP2_Mitochondrial_Lysosomal_Autophagy_direct_genes_9$") ~ "three_pathways",
          str_detect(PD_PRS, "PD_PRScs_GP2_Mitochondrial_Lysosomal_Autophagy_direct_genes_9_exclude$") ~ "non_pathways",
          TRUE ~ NA_character_
        )
      ) %>%
      filter(!is.na(PRS_type))
    
    # Manual renaming of PRS values
    prs_mapping <- c(
      "all_genes" = "all genes (PD-PRS)",
      "autophagy" = "autophagy pathway",
      "lysosomal" = "lysosomal pathway",
      "mitochondrial" = "mitochondrial pathway",
      "three_pathways" = "three pathways combined",
      "non_pathways" = "non-pathways genes"
    )
    
    # Manual renaming of Map values
    map_mapping <- c(
      "group_avg" = "Group Average",
      "hy_1" = "HY Stage 1",
      "hy_2" = "HY Stage 2",
      "hy_3" = "HY Stage 3",
      "hy_4" = "HY Stage 4"
    )
    
    # Apply the mappings
    plot_data <- plot_data %>%
      mutate(
        PRS = prs_mapping[PRS_type],
        Map = map_mapping[Map],
        # Create a significance indicator based on spin-test p-values
        sig_indicator = case_when(
          Spin_P < 0.001 ~ "***",
          Spin_P < 0.01 ~ "**",
          Spin_P < 0.05 ~ "*",
          TRUE ~ "")
      )
    
    # Define the order for PRS (Y-axis)
    prs_order <- c(
      "all genes (PD-PRS)",
      "non-pathways genes",
      "three pathways combined",
      "mitochondrial pathway",
      "lysosomal pathway",
      "autophagy pathway"
    )
    
    # Define the order for Map values (X-axis)
    map_order <- c(
      "Group Average",
      "HY Stage 1",
      "HY Stage 2",
      "HY Stage 3",
      "HY Stage 4"
    )
    
    # Convert PRS and Map to factors with the specified order
    plot_data <- plot_data %>%
      mutate(
        PRS = factor(PRS, levels = rev(prs_order)),
        Map = factor(Map, levels = map_order)
      )
    
    ### Actual visualization
    # Loop through each sex group
    for(sex in c("all", "male", "female")) {
      # Filter data for the current sex
      sex_data <- plot_data %>% filter(Sex == sex)
      
      # Set fixed limit for consistency across plots
      max_limit <- 0.6
      
      # Create the heatmap
      p <- ggplot(sex_data, aes(x = Map, y = PRS)) +
        # Use Rho for the color mapping
        geom_tile(aes(fill = Rho), color = "white", linewidth = 0.5) +
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
          name = "Correlation (ρ)",
          limits = c(-max_limit, max_limit),
          breaks = seq(-max_limit, max_limit, by = 0.3),  # Show more tick marks on legend
          labels = function(x) sprintf("%.1f", x)  # Format to 1 decimal place
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
        # Add informative title with your requested changes
        labs(
          title = paste0("Correlation between ", prs_measure, "-PRS & ENIGMA PD ", enigma_measure, " Maps (", 
                         ifelse(sex == "all", "All Participants", 
                                ifelse(sex == "male", "Males", "Females")), ")"),
          subtitle = "Spearman spatial correlations (* p<0.05, ** p<0.01, *** p<0.001)"
        )
      
      # Save as SVG with measure types and sex in filename
      svg_filename <- paste0(output_dir, "ENIGMA_PD_", enigma_measure, "_vs_", prs_measure, "_PRS_Correlation_Heatmap_", sex, ".svg")
      ggsave(svg_filename, plot = p, width = 7.5, height = 6)
      
      # Save as PNG with measure types and sex in filename
      png_filename <- paste0(output_dir, "ENIGMA_PD_", enigma_measure, "_vs_", prs_measure, "_PRS_Correlation_Heatmap_", sex, ".png")
      ggsave(png_filename, plot = p, width = 7.5, height = 6, dpi = 300)
      
      cat("Created", enigma_measure, "vs", prs_measure, "heatmap for", sex, "participants\n")
    }
  }
}

####################### BOX PLOT VISUALIZATION #######################
# Process each PRS measure
for(prs_measure in c("SA", "CT")) {
  # Create output directory for box plots for this measure
  output_dir <- paste0("Plots/boxplots/", prs_measure, "_PRS/")
  if (!dir.exists(output_dir)) { dir.create(output_dir, recursive = TRUE) }
  
  # Process each ENIGMA measure type (SA and CT)
  for(enigma_measure in c("SA", "CT")) {
    # Define file paths for results and null distributions
    result_file <- paste0('Outputs/', prs_measure, '_PRS_Results/SpinTest_ENIGMA_', enigma_measure, '_vs_', prs_measure, '_PRS_Results.csv')
    null_file <- paste0('Outputs/', prs_measure, '_PRS_Results/SpinTest_ENIGMA_', enigma_measure, '_vs_', prs_measure, '_PRS_Nulls_wide.rds')
    
    # Read the results data
    results_data <- read.csv(result_file) %>%
      filter(Correlation == 'spearman') %>%
      # Filter for only all_genes (PD_PRScs_GP2)
      filter(str_detect(PD_PRS, "PD_PRScs_GP2$"))
    
    # Read the null distribution data
    nulls_data <- readRDS(null_file) %>%
      filter(Correlation == 'spearman') %>%
      # Filter for only all_genes (PD_PRScs_GP2)
      filter(str_detect(PD_PRS, "PD_PRScs_GP2$"))
    
    # Get the null spin columns
    spin_cols <- grep("Spin_[0-9]+", colnames(nulls_data), value = TRUE)
    
    # Map names for plots
    map_mapping <- c(
      "group_avg" = "Group Average",
      "hy_1" = "HY Stage 1",
      "hy_2" = "HY Stage 2",
      "hy_3" = "HY Stage 3",
      "hy_4" = "HY Stage 4"
    )
    
    # Define the order for Map values - include group_avg
    map_order <- c("group_avg", "hy_1", "hy_2", "hy_3", "hy_4")
    map_labels <- c("Group Average", "HY Stage 1", "HY Stage 2", "HY Stage 3", "HY Stage 4")
    
    # Loop through each sex group
    for(sex in c("all", "male", "female")) {
      # Filter data for the current sex
      sex_results <- results_data %>% filter(Sex == sex) %>% filter(Map %in% map_order)
      sex_nulls <- nulls_data %>% filter(Sex == sex) %>% filter(Map %in% map_order)
      
      # Prepare data for boxplot
      plot_data <- NULL
      actual_values <- NULL
      
      for(map_val in map_order) {
        # Get the actual correlation value for this map
        actual_corr <- sex_results %>% 
          filter(Map == map_val) %>% 
          pull(Rho)
        
        # Get the p-value for significance
        p_value <- sex_results %>% 
          filter(Map == map_val) %>% 
          pull(Spin_P)
        
        # Get the null distribution values for this map
        null_row <- sex_nulls %>% filter(Map == map_val)
        # Extract the spin values and convert to a proper vector
        map_nulls <- as.numeric(as.matrix(null_row[, spin_cols]))
        
        # Create a data frame for this map
        map_data <- data.frame(
          Map = rep(map_mapping[map_val], length(map_nulls)),
          Value = map_nulls,
          Type = rep("Null Distribution", length(map_nulls)),
          Significant = rep(p_value < 0.05, length(map_nulls)),
          Direction = rep(ifelse(actual_corr > 0, "positive", "negative"), length(map_nulls)),
          CorrType = rep(
            ifelse(p_value >= 0.05, "Non-significant",
                   ifelse(actual_corr > 0, "Positive (p_spin < 0.05)", "Negative (p_spin < 0.05)")),
            length(map_nulls))
        )
        
        # Add the actual correlation as a separate row
        actual_data <- data.frame(
          Map = map_mapping[map_val],
          Value = actual_corr,
          Type = "Actual Correlation",
          Significant = p_value < 0.05,
          Direction = ifelse(actual_corr > 0, "positive", "negative"),
          CorrType = ifelse(p_value >= 0.05, "Non-significant",
                            ifelse(actual_corr > 0, "Positive (p_spin < 0.05)", "Negative (p_spin < 0.05)"))
        )
        
        # Combine with the main data frame
        if(is.null(plot_data)) {
          plot_data <- map_data
          actual_values <- actual_data
        } else {
          plot_data <- rbind(plot_data, map_data)
          actual_values <- rbind(actual_values, actual_data)
        }
      }
      
      # Convert Map to a factor with the specified order
      plot_data$Map <- factor(plot_data$Map, levels = map_labels)
      actual_values$Map <- factor(actual_values$Map, levels = map_labels)
      
      # Create a legend dataframe with exactly one entry per category
      legend_data <- data.frame(
        CorrType = c("Positive (p_spin < 0.05)", "Negative (p_spin < 0.05)", "Non-significant"),
        stringsAsFactors = FALSE
      )
      # Make sure it's in the right order
      legend_data$CorrType <- factor(legend_data$CorrType, 
                                     levels = c("Positive (p_spin < 0.05)", 
                                                "Negative (p_spin < 0.05)", 
                                                "Non-significant"))
      
      # Create the boxplot with enhanced visuals
      p <- ggplot() +
        # Add jittered points for null distribution (with transparency)
        # Non-significant points are grey
        geom_jitter(
          data = plot_data %>% filter(CorrType == "Non-significant"),
          aes(x = Map, y = Value),
          width = 0.2, alpha = 0.1, size = 1, color = "grey70"
        ) +
        # Significant negative points are blue
        geom_jitter(
          data = plot_data %>% filter(CorrType == "Negative (p_spin < 0.05)"),
          aes(x = Map, y = Value),
          width = 0.2, alpha = 0.1, size = 1, color = "#4575b4"
        ) +
        # Significant positive points are red
        geom_jitter(
          data = plot_data %>% filter(CorrType == "Positive (p_spin < 0.05)"),
          aes(x = Map, y = Value),
          width = 0.2, alpha = 0.1, size = 1, color = "#d73027"
        ) +
        
        # Add boxplots with fill according to correlation type
        geom_boxplot(
          data = plot_data %>% filter(CorrType == "Non-significant"),
          aes(x = Map, y = Value),
          fill = "grey90", 
          alpha = 0.7,
          color = "grey30",
          width = 0.5,
          outlier.shape = NA
        ) +
        
        geom_boxplot(
          data = plot_data %>% filter(CorrType == "Negative (p_spin < 0.05)"),
          aes(x = Map, y = Value),
          fill = "#4575b4", 
          alpha = 0.7,
          color = "grey30",
          width = 0.5,
          outlier.shape = NA
        ) +
        
        geom_boxplot(
          data = plot_data %>% filter(CorrType == "Positive (p_spin < 0.05)"),
          aes(x = Map, y = Value),
          fill = "#d73027", 
          alpha = 0.7,
          color = "grey30",
          width = 0.5,
          outlier.shape = NA
        ) +
        
        # Add actual correlation points
        geom_point(
          data = actual_values %>% filter(CorrType == "Non-significant"),
          aes(x = Map, y = Value),
          shape = 21,
          fill = "grey90",
          size = 4.5,
          color = "black",
          stroke = 1
        ) +
        
        geom_point(
          data = actual_values %>% filter(CorrType == "Negative (p_spin < 0.05)"),
          aes(x = Map, y = Value),
          shape = 23,
          fill = "#4575b4",
          size = 4.5,
          color = "black",
          stroke = 1
        ) +
        
        geom_point(
          data = actual_values %>% filter(CorrType == "Positive (p_spin < 0.05)"),
          aes(x = Map, y = Value),
          shape = 23,
          fill = "#d73027",
          size = 4.5,
          color = "black",
          stroke = 1
        ) +
        
        # Add horizontal line at y=0 for reference
        geom_hline(yintercept = 0, linetype = "dotted", color = "gray50", size = 0.5) +
        
        # Add x axis order
        scale_x_discrete(limits = map_labels) +
        
        # Customize theme
        theme_minimal() +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 11, face = "bold"),
          axis.text.y = element_text(size = 10),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 12, face = "bold"),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.major.y = element_line(linetype = "dotted", color = "gray80"),
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          plot.subtitle = element_text(size = 11, hjust = 0.5, margin = margin(b = 15)),
          legend.position = "right",
          legend.title = element_text(face = "bold"),
          plot.background = element_rect(fill = "white", color = NA),
          panel.background = element_rect(fill = "white", color = NA),
          plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
        ) +
        # Add informative title
        labs(
          title = paste0("Full PD-PRS & ENIGMA PD Maps Correlation"),
          subtitle = paste0(prs_measure, "-PRS vs ENIGMA ", enigma_measure, " null distributions (", 
                            ifelse(sex == "all", "All Participants", 
                                   ifelse(sex == "male", "Males", "Females")), ")"),
          y = "Correlation (ρ)"
        ) +
        # Set y-axis limits to be symmetric and appropriate for correlation values
        scale_y_continuous(
          limits = c(-0.6, 0.6), 
          breaks = seq(-0.6, 0.6, by = 0.2),
          labels = function(x) sprintf("%.1f", x)  # Format to 1 decimal place
        )
      
      # Manually add a custom legend
      # Create a blank annotation layer for the legend
      p <- p + 
        annotation_custom(
          grob = grid::rectGrob(gp = grid::gpar(fill = NA, col = NA)),
          xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf
        )
      
      # Create empty data for the legend
      legend_data <- data.frame(
        x = c(1, 1, 1),
        y = c(1, 1, 1),
        CorrType = factor(c("Positive (p_spin < 0.05)", 
                            "Negative (p_spin < 0.05)", 
                            "Non-significant"),
                          levels = c("Positive (p_spin < 0.05)", 
                                     "Negative (p_spin < 0.05)", 
                                     "Non-significant"))
      )
      
      # Add the legend manually
      p <- p + 
        # This adds an invisible geom to generate the legend
        geom_point(
          data = legend_data,
          aes(x = x, y = y, fill = CorrType, shape = CorrType),
          alpha = 0
        ) +
        scale_shape_manual(
          name = "Correlation",
          values = c(
            "Positive (p_spin < 0.05)" = 23,
            "Negative (p_spin < 0.05)" = 23,
            "Non-significant" = 21
          )
        ) +
        scale_fill_manual(
          name = "Correlation",
          values = c(
            "Positive (p_spin < 0.05)" = "#d73027",
            "Negative (p_spin < 0.05)" = "#4575b4",
            "Non-significant" = "grey90"
          )
        ) +
        guides(fill = guide_legend(override.aes = list(
          alpha = 1,
          shape = c(23, 23, 21),
          size = 4
        )))
      
      # Save as SVG
      svg_filename <- paste0(output_dir, "ENIGMA_PD_", enigma_measure, "_vs_", prs_measure, "_PRS_Null_Distribution_", sex, ".svg")
      ggsave(svg_filename, plot = p, width = 8, height = 10)
      
      # Save as PNG
      png_filename <- paste0(output_dir, "ENIGMA_PD_", enigma_measure, "_vs_", prs_measure, "_PRS_Null_Distribution_", sex, ".png")
      ggsave(png_filename, plot = p, width = 8, height = 10, dpi = 300)
      
      cat("Created enhanced", enigma_measure, "vs", prs_measure, "null distribution boxplot for", sex, "participants\n")
    }
  }
}
