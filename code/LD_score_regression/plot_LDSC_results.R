# Load required libraries
library(dplyr)
library(data.table)
library(ggplot2)
library(reshape2)
library(viridis)
library(scales)
library(stringr)

rm(list = ls())
cat("\014")

setwd("/Users/houmanazizi/Library/Mobile\ Documents/com~apple~CloudDocs/Education/GitHub/LD_Score_Regression/")

# Load the data
data <- read.table('./results_by_Lang/LDSC_correlation_noGlobalCorrection.tsv', 
                   sep = '\t', header = TRUE, stringsAsFactors = FALSE)

# Get unique traits + keep important ones
all_traits <- unique(c(data$Trait1, data$Trait2))

traits_PD <- c('PD.2025', 'Mitochondrial_Lysosomal_Autophagy_direct_genes_9_excluded_PD_GWAS', 'Mitochondrial_Lysosomal_Autophagy_direct_genes_9_PD_GWAS', 
               'Autophagy_direct_genes_4_PD_GWAS', 'Mitochondrial_genes_1_PD_GWAS', 'Lysosomal_direct_genes_2_PD_GWAS')
traits_brain <- c('Accumbens_GarciaMarin', 'Amygdala_GarciaMarin', 'Brainstem_GarciaMarin', 'Caudate_GarciaMarin', 'Hippocampus_GarciaMarin', 'Pallidum_GarciaMarin', 'Putamen_GarciaMarin', 'Thalamus_GarciaMarin', all_traits[grepl("^ENIGMA3_", all_traits)])

# Filter data
filtered_data <- data[data$Trait1 %in% traits_PD & data$Trait2 %in% traits_brain, ]
filtered_data[(filtered_data$Correlation > 1 | filtered_data$Correlation < -1) & !is.na(filtered_data$Correlation), 3:5] <- NA

# Rename PD and subcortical names
filtered_data$Trait1[filtered_data$Trait1 == 'PD.2025'] <- 'PD'
filtered_data$Trait2 <- tolower(gsub('_GarciaMarin', '', filtered_data$Trait2))


# Apply FDR correction
filtered_data$p_fdr <- p.adjust(filtered_data$p, method = "fdr")

########################## PREPROCESSING ##########################

# Function to format brain region names properly
format_brain_name <- function(name) {
  # Handle special cases first
  special_cases <- c(
    "bankssts" = "bankssts",
    "frontalpole" = "frontalpole",
    "temporalpole" = "temporalpole",
    "superiorfrontal" = "superior frontal",
    "rostralmiddlefrontal" = "rostral middle frontal",
    "caudalmiddlefrontal" = "caudal middle frontal",
    "parsopercularis" = "pars opercularis",
    "parstriangularis" = "pars triangularis",
    "parsorbitalis" = "pars orbitalis",
    "lateralorbitofrontal" = "lateral orbitofrontal",
    "medialorbitofrontal" = "medial orbitofrontal",
    "precentral" = "precentral",
    "postcentral" = "postcentral",
    "supramarginal" = "supramarginal",
    "superiorparietal" = "superior parietal",
    "inferiorparietal" = "inferior parietal",
    "precuneus" = "precuneus",
    "cuneus" = "cuneus",
    "pericalcarine" = "pericalcarine",
    "lateraloccipital" = "lateral occipital",
    "lingual" = "lingual",
    "fusiform" = "fusiform",
    "parahippocampal" = "parahippocampal",
    "entorhinal" = "entorhinal",
    "superiortemporal" = "superior temporal",
    "middletemporal" = "middle temporal",
    "inferiortemporal" = "inferior temporal",
    "transversetemporal" = "transverse temporal",
    "insula" = "insula",
    "isthmuscingulate" = "isthmus cingulate",
    "posteriorcingulate" = "posterior cingulate",
    "rostralanteriorcingulate" = "rostral anterior cingulate",
    "caudalanteriorcingulate" = "caudal anterior cingulate"
  )
  
  if (name %in% names(special_cases)) {
    return(special_cases[name])
  }
  
  # For names not in special cases, try to format automatically
  # Insert spaces before capital letters (for camelCase)
  formatted <- gsub("([a-z])([A-Z])", "\\1 \\2", name)
  
  # Convert to lowercase
  formatted <- tolower(formatted)
  
  return(formatted)
}

# Updated function to categorize and rename brain traits
categorize_brain_traits <- function(trait_name) {
  if (trait_name %in% c('accumbens', 'amygdala', 'brainstem', 'caudate', 'pallidum', 'putamen', 'thalamus', 'hippocampus')) {
    return(list(category = "subcortical volumes", display_name = format_brain_name(trait_name)))
  } else if (str_detect(trait_name, "surfavg")) {
    # Extract region name from ENIGMA3 surface area measures
    region <- str_extract(trait_name, "(?<=enigma3_mixed_se_wo_mean_)([^_]+)(?=_surfavg)")
    return(list(category = "surface area", display_name = format_brain_name(region)))
  } else if (str_detect(trait_name, "thickavg")) {
    # Extract region name from ENIGMA3 thickness measures
    region <- str_extract(trait_name, "(?<=enigma3_mixed_se_wo_mean_)([^_]+)(?=_thickavg)")
    return(list(category = "thickness", display_name = format_brain_name(region)))
  } else if (str_detect(trait_name, 'full_surfarea')) {
    return(list(category = "surface area", display_name = 'total cortical surface area'))
  } else if (str_detect(trait_name, 'full_thickness')) {
    return(list(category = "thickness", display_name = 'total cortical thickness'))   
  } else {
    return(list(category = "other", display_name = format_brain_name(trait_name)))
  }
}

# Function to rename PD traits
rename_pd_traits <- function(trait_name) {
  pd_mapping <- c(
    "PD" = "Parkinson's GWAS (all genes)",
    "Mitochondrial_Lysosomal_Autophagy_direct_genes_9_excluded_PD_GWAS" = "non-pathways genes",
    "Mitochondrial_Lysosomal_Autophagy_direct_genes_9_PD_GWAS" = "three pathways combined",
    "Autophagy_direct_genes_4_PD_GWAS" = "autophagy pathway",
    "Mitochondrial_genes_1_PD_GWAS" = "mitochondrial pathway",
    "Lysosomal_direct_genes_2_PD_GWAS" = "lysosomal pathway"
  )
  return(pd_mapping[trait_name])
}

# Process the data
plot_data <- filtered_data %>%
  rowwise() %>%
  mutate(
    brain_info = list(categorize_brain_traits(Trait2)),
    brain_category = brain_info$category,
    brain_display_name = brain_info$display_name,
    pd_display_name = rename_pd_traits(Trait1),
    # Create significance indicator based on FDR-corrected p-values
    sig_indicator = case_when(
      p_fdr < 0.001 ~ "***",
      p_fdr < 0.01 ~ "**",
      p_fdr < 0.05 ~ "*",
      TRUE ~ ""
    )
  ) %>%
  select(-brain_info) %>%
  filter(!is.na(Correlation))

# remove unwanted regions
unwanted_regions <- c("bankssts", "frontalpole", "temporalpole")
plot_data <- plot_data %>% filter(!brain_display_name %in% unwanted_regions)

# Define the order for PD traits (Y-axis)
pd_order <- c(
  "Parkinson's GWAS (all genes)",
  "non-pathways genes",
  "three pathways combined",
  "mitochondrial pathway",
  "lysosomal pathway",
  "autophagy pathway"
)

# Define category order for better visualization
category_order <- c(
  "subcortical volumes",
  "surface area",
  "thickness"
)

# Get unique brain traits and order them by category, then alphabetically within category
brain_order_df <- plot_data %>%
  select(brain_display_name, brain_category) %>%
  distinct() %>%
  # Create custom ordering within each category
  arrange(
    match(brain_category, category_order),
    case_when(
      # For surface area: put 'total cortical surface area' first, then alphabetical
      brain_category == "surface area" & brain_display_name == "total cortical surface area" ~ "0",
      brain_category == "surface area" ~ paste("1", brain_display_name),
      # For thickness: put 'total cortical thickness' first, then alphabetical  
      brain_category == "thickness" & brain_display_name == "total cortical thickness" ~ "0",
      brain_category == "thickness" ~ paste("1", brain_display_name),
      # For subcortical volumes: just alphabetical
      TRUE ~ brain_display_name
    )
  )

# Create position index with spacing between categories
spacing_width <- 2  # Width of white space between categories

brain_order_df <- brain_order_df %>%
  group_by(brain_category) %>%
  mutate(position_within_category = row_number()) %>%
  ungroup() %>%
  mutate(
    category_num = as.numeric(factor(brain_category, levels = category_order)),
    # Add spacing between categories
    position = row_number() + (category_num - 1) * spacing_width
  )


# Update plot_data with new positions
plot_data <- plot_data %>%
  left_join(brain_order_df %>% select(brain_display_name, position, brain_category), by = c("brain_display_name", 'brain_category')) %>%
  mutate(
    pd_display_name = factor(pd_display_name, levels = rev(pd_order))
  ) %>%
  filter(!is.na(position))

# Create category breaks for labels
category_breaks <- brain_order_df %>%
  group_by(brain_category) %>%
  summarise(
    start_pos = min(position),
    end_pos = max(position),
    mid_pos = (start_pos + end_pos) / 2,
    .groups = 'drop'
  )

# Set fixed limit for consistency
max_limit <- max(abs(plot_data$Correlation), na.rm = TRUE)
max_limit <- ceiling(max_limit * 10) / 10  # Round up to nearest 0.1

# Create directory if it doesn't exist
if (!dir.exists('./Figures')) { dir.create('./Figures') }

########################## VISUALIZATION ##########################

# Create the heatmap
p <- ggplot(plot_data, aes(x = position, y = pd_display_name)) +
  # Use Correlation for the color mapping
  geom_tile(aes(fill = Correlation), color = "white", linewidth = 0.5, width = 0.9, height = 0.9) +
  # Add significance indicators
  geom_text(aes(label = sig_indicator),
            color = "black", size = 4, fontface = "bold") +
  # Use a diverging color scale
  scale_fill_gradientn(
    colors = c("#313695", "#4575b4", "#74add1", "#abd9e9", "#e0f3f8",
               "#FFFBF3",  # White at center
               "#fee090", "#fdae61", "#f46d43", "#d73027", "#a50026"),
    values = scales::rescale(c(-max_limit, -max_limit*0.8, -max_limit*0.6, -max_limit*0.4, -max_limit*0.2,
                               0,  # Center point
                               max_limit*0.2, max_limit*0.4, max_limit*0.6, max_limit*0.8, max_limit)),
    name = "genetic correlation (rg)",
    limits = c(-max_limit, max_limit),
    breaks = seq(-max_limit, max_limit, by = 0.2),
    labels = function(x) sprintf("%.1f", x),
    na.value = "grey90"
  ) +
  # Set x-axis breaks and labels
  scale_x_continuous(
    breaks = brain_order_df$position,
    labels = brain_order_df$brain_display_name,
    expand = c(0.01, 0.01)
  ) +
  # Customize theme
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 16),
    axis.text.y = element_text(size = 16), # originally 11
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11),
    plot.margin = margin(t = 60, r = 20, b = 80, l = 20)  # Extra top margin for category labels
  ) +
  # Add informative title
  labs(
    title = "genetic correlations: pd traits vs brain traits",
    subtitle = "fdr-corrected significance: * p<0.05, ** p<0.01, *** p<0.001"
  )

# Add category labels at the top
for(i in 1:nrow(category_breaks)) {
  p <- p + annotate("text",
                    x = category_breaks$mid_pos[i],
                    y = length(pd_order) + 1.2,  # Position above the plot
                    label = category_breaks$brain_category[i],
                    size = 4,
                    fontface = "bold",
                    color = "gray30")
}

# Extend plot area to accommodate category labels
max_pos <- max(brain_order_df$position)
p <- p + coord_cartesian(ylim = c(0.5, length(pd_order) + 1.5),
                         xlim = c(0.5, max_pos + 0.5),
                         clip = "off")

# Save as SVG
svg_filename <- './Figures/PD_Brain_genetic_Correlation_improved.svg'
ggsave(svg_filename, plot = p, width = 45, height = 8.5)

# Save as PNG
png_filename <- './Figures/PD_Brain_genetic_Correlation_improved.png'
ggsave(png_filename, plot = p, width = 45, height = 8.5, dpi = 300, bg = 'white')

cat("Improved plots saved to ./Figures/\n")
cat(paste("Number of significant correlations (FDR < 0.05):", sum(plot_data$p_fdr < 0.05, na.rm = TRUE), "\n"))
cat(paste("Number of significant correlations (uncorrected p < 0.05):", sum(plot_data$p < 0.05, na.rm = TRUE), "\n"))

# Print summary statistics
cat("\nSummary:\n")
cat(paste("Total correlations plotted:", nrow(plot_data), "\n"))
cat(paste("Range of correlations:", round(min(plot_data$Correlation, na.rm = TRUE), 3), 
          "to", round(max(plot_data$Correlation, na.rm = TRUE), 3), "\n"))
cat(paste("Categories:", paste(unique(plot_data$brain_category), collapse = ", "), "\n"))


# Save the result file
final_summary <- plot_data %>% select("Trait1", "Trait2", "Correlation", "p", "SE", "p_fdr", "brain_category", "brain_display_name", "pd_display_name") %>% arrange(Trait1, brain_category)
final_summary_PD <- final_summary %>% filter(Trait1 == 'PD')

write.csv(final_summary, './Figures/correlation_CSVs/final_summary_LDSC.csv', row.names = FALSE)
write.csv(final_summary_PD, './Figures/correlation_CSVs/final_summary_LDSC_PD.csv', row.names = FALSE)
