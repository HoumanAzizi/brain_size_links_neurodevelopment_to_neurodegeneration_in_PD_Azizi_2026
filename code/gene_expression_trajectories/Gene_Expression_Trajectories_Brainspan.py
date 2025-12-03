
'''
Script to visualize gene expression trajectories across developmental stages using BrainSpan dataset (by Moohebat)
'''

### Load required packages
import numpy as np
import pandas as pd
from scipy.stats import zscore
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import colormaps as cmaps
import seaborn as sns
import statsmodels
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches  # For creating custom legend handles


### Define path and load preprocessed data
data_path = '/Users/houmanazizi/Library/Mobile Documents/com~apple~CloudDocs/Education/GitHub/Gene_Developmental_Expression/Data'

## old QC
# bs_exp_old = pd.read_csv(f'{data_path}/bs_expression_cleaned.csv').drop('Unnamed: 0', axis=1)
# sample_info_old = pd.read_csv(f'{data_path}/sample_data_cleaned.csv').drop('Unnamed: 0', axis=1)

## new QC
bs_exp = pd.read_csv(f'{data_path}/brainspan_rna_exp_qc.csv').drop('Unnamed: 0', axis=1)
bs_exp = bs_exp.set_index('gene_symbol').T.reset_index().drop(columns='index')
bs_exp = np.log(bs_exp + 1) # to log transform

sample_info = pd.read_csv(f'{data_path}/sample_data_rna_qc.csv').drop('Unnamed: 0', axis=1)
sample_info['age_weeks'] = sample_info['age_days'] / 7
sample_info['log_age_weeks'] = np.log(sample_info['age_weeks'])
sample_info['age_group2'] = sample_info['age_group2'].str.replace('_', ' ')


########################
#### Cross-Brain normalization using 75th percentile
### Step 1: Merge sample_info and bs_exp
# Add donor information to bs_exp by merging with sample_info
merged_df = pd.concat([sample_info[['donor_name']], bs_exp], axis=1)




### Step 2: Group by donor_name and normalize
## v1: using for loop
# Group the data by donor_name
grouped = merged_df.groupby('donor_name')
# Initialize list to store normalized dataframes
normalized_dfs = []
# Process each group
for name, group in grouped:
    # Calculate 75th percentile across all genes/region/values for the current group
    #percentile_75 = group.iloc[:, 1:].quantile(0.75).mean(axis=0) # normalize over columns (i.e. 1 value per gene)
    #percentile_75 = group.iloc[:, 1:].quantile(0.75, axis=1).mean(axis=0) # normalize over rows (i.e. 1 value per brain region)
    percentile_75 = np.percentile(group.iloc[:, 1:].values.flatten(), 75) # flatten all values and normalize that way
    # Normalize gene expression values by the 75th percentile
    normalized_group = group.iloc[:, 1:] / percentile_75
    # Add back the donor_name column
    normalized_group.insert(0, 'donor_name', group['donor_name'].values)
    # Append to list
    normalized_dfs.append(normalized_group)
# Combine all normalized dataframes
normalized_data = pd.concat(normalized_dfs)




# Step 3: Drop the donor_name column and convert back to a format similar to bs_exp
bs_exp = normalized_data.drop('donor_name', axis=1).copy()

# create dataframe for gene set plots
bs_df = bs_exp.T.reset_index()
bs_df = bs_df.rename(columns={'index': 'gene_symbol'})

'''
## Log transform
# Get expression columns (all except gene_symbol)
expression_cols = bs_df.columns[1:]
bs_df_transformed = bs_df.copy()
bs_df_transformed[expression_cols] = np.log2(bs_df[expression_cols] + 1)
bs_df = bs_df_transformed.copy()
'''



############################################### Gene Sets ###############################################
# Load the pathway genes file
#pathway_genes = pd.read_csv(f'{data_path}/pathway_genes_in_Nalls_PD_GWAS.csv') # Nalls version
rank_limit = '3EricYu'
pathway_genes = pd.read_csv(f'{data_path}/pathway_genes_in_GP2_PD_GWAS_rank{rank_limit}.csv') # GP2 version

# Iterate through each row (gene list) in the pathway genes file
gene_list = {}
for index, row in pathway_genes.iterrows():
    gene_list_name = row['gene_list_name']
    
    # Get the genes that have a value of 1 in this row
    # First, we exclude the 'gene_list_name' column
    gene_indicators = row.drop('gene_list_name')
    
    # Find column names (genes) where value is 1
    genes_in_list = gene_indicators[gene_indicators == 1].index.tolist()
    # Find column names (genes) where value is 0
    genes_not_in_list = gene_indicators[gene_indicators == 0].index.tolist()
    
    # Add to gene_list dictionary
    # Add both lists to gene_list dictionary
    gene_list[gene_list_name] = genes_in_list
    gene_list[f'GP2_except_{gene_list_name}'] = genes_not_in_list # or GP2 except
    

# Create dictionaries for gene set expressions
bs_exp_energy = {} # contains expression data for each pathway's genes
bs_mean_energy = {} # contains the mean expression across all genes in each pathway
for key, value in gene_list.items():
    # For each pathway, get expression data for its genes
    bs_exp_energy[key] = bs_df[bs_df.gene_symbol.isin(value)]
    numeric_data = bs_exp_energy[key].drop('gene_symbol', axis=1) # Drop the 'gene_symbol' column before calculating mean
    bs_exp_energy[key] = numeric_data
    # Calculate mean expression across all genes in the pathway
    bs_mean_energy[key] = np.mean(bs_exp_energy[key], axis=0)


# Converting mean expression dict to dataframe
bs_mean_energy_df = (pd.DataFrame.from_dict(bs_mean_energy).reset_index().drop('index', axis=1))
# Combine sample information with mean expression data
df = pd.concat([sample_info, bs_mean_energy_df], axis=1)

# Optional: removing outlier samples
remove_outliers = 0 # Flag to control outlier removal
if remove_outliers:
    df = df[df['age'] != '2 yrs']
    df = df[df['age'] != '18 yrs']





## Option 3.2: using polynomial regression for plot with 4 specific gene sets
# Define the specific gene sets we want to plot
specific_gene_sets = [
    'Mitochondrial_genes_1',
    'Lysosomal_direct_genes_2',
    'Autophagy_direct_genes_4',
    'GP2_except_Mitochondrial_Lysosomal_Autophagy_direct_genes_9'
]

# Define colors for each gene set
colors = {
    'Mitochondrial_genes_1': 'blue', #blue
    'Lysosomal_direct_genes_2': 'green', #green
    'Autophagy_direct_genes_4': 'purple', #purple
    'GP2_except_Mitochondrial_Lysosomal_Autophagy_direct_genes_9': 'coral' #coral
}

# Define labels for legend (shorter names)
legend_labels = {
    'Mitochondrial_genes_1': 'Mitochondrial',
    'Lysosomal_direct_genes_2': 'Lysosomal',
    'Autophagy_direct_genes_4': 'Autophagy',
    'GP2_except_Mitochondrial_Lysosomal_Autophagy_direct_genes_9': 'Other PD risk genes'
}

# Create a figure with three subplots - HEIGHT SET TO 8
fig = plt.figure(figsize=(15, 8))
gs = gridspec.GridSpec(1, 3, width_ratios=[1, 1, 1])
ax1 = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1])
ax3 = plt.subplot(gs[2])

#### First subplot with log_age_weeks for all 4 gene sets
for gene_set in specific_gene_sets:
    sns.regplot(
        data=df[df['ctx'] == 'cortex'],
        x='log_age_weeks',
        y=gene_set,
        order=3,  # Cubic polynomial fit
        scatter=False,  # Do not plot scatter points
        line_kws={'lw': 1.5, 'color': colors[gene_set]},  # Customize regression line appearance
        ax=ax1
    )

# Create a custom legend for the first subplot
legend_handles = [mpatches.Patch(color=colors[gene_set], label=legend_labels[gene_set]) for gene_set in specific_gene_sets]
ax1.legend(
    handles=legend_handles,
    title='Gene Sets',
    bbox_to_anchor=(0.5, -0.35),  # Adjusted position for figure height of 8
    loc='upper center',
    frameon=True
)

ax1.set_xlabel('Log age weeks', fontsize=11)  # Increased font size
ax1.set_ylabel('Normalized expression', fontsize=11)  # Increased font size
ax1.set_title('Gene Set Comparison', fontsize=12)  # Increased font size


#### Second subplot with log_age_weeks and age group markers
for gene_set in specific_gene_sets:
    sns.regplot(
        data=df[df['ctx'] == 'cortex'],
        x='log_age_weeks',
        y=gene_set,
        order=3,
        scatter=False,
        line_kws={'lw': 1.5, 'color': colors[gene_set]},
        ax=ax2
    )

# Get the current y-axis limits
y_min, y_max = ax2.get_ylim()
# Calculate position for text
text_y_position = y_min - (y_max - y_min) * 0.13  # Adjusted for figure height of 8

# Add vertical lines for age groups using max values
for group in df['age_group2'].unique():
    max_age = df[df['age_group2'] == group]['log_age_weeks'].max()
    ax2.axvline(x=max_age, color='gray', linestyle='--', alpha=0.3)
    ax2.text(max_age, text_y_position, group, 
            rotation=45, ha='right', va='top', fontsize=9,
            color='dimgray')
    
# Add birth line
birth_x = np.log(((401.72 - (4*30)) / 7))
ax2.axvline(x=birth_x, color='red', linestyle='--', alpha=0.3)

# Set xlabel with adjusted position using labelpad - INCREASED LABELPAD FROM 45 TO 60
ax2.set_xlabel('Developmental stage (log age weeks)', labelpad=60, fontsize=11)

# Legend
legend_handles = [mpatches.Patch(color=colors[gene_set], label=legend_labels[gene_set]) for gene_set in specific_gene_sets]
legend_handles.append(Line2D([0], [0], color='red', linestyle='--', alpha=0.3, label='Birth'))  # Added birth line to legend

ax2.legend(
    handles=legend_handles,
    title='Gene Sets',
    bbox_to_anchor=(0.5, -0.5),  # Adjusted for figure height of 8
    loc='upper center',
    frameon=True
)

ax2.set_ylabel('Normalized expression', fontsize=11)  # Increased font size
ax2.set_title('Gene Set Comparison', fontsize=12)  # Increased font size


#### Third subplot (box plots with individual points)
# For each gene set, add box plot, strip plot, and median line
for gene_set in specific_gene_sets:
    # Box plot
    sns.boxplot(
        data=df[df['ctx'] == 'cortex'],
        x='age_group2',
        y=gene_set,
        boxprops={'facecolor': 'none', 'edgecolor': colors[gene_set]},
        whiskerprops={'color': colors[gene_set]},
        capprops={'visible': False},
        medianprops={'color': colors[gene_set]},
        showfliers=False,
        linewidth=0.7,
        width=0.4,
        ax=ax3
    )

    # Strip plot
    sns.stripplot(
        data=df[df['ctx'] == 'cortex'],
        x='age_group2',
        y=gene_set,
        color=colors[gene_set],
        size=3,
        alpha=0.2,
        jitter=True,
        zorder=-5,
        ax=ax3
    )

    # Median line
    sns.lineplot(
        data=df[df['ctx'] == 'cortex'],
        x='age_group2',
        y=gene_set,
        estimator=np.median,
        ci=None,
        sort=False,
        color=colors[gene_set],
        marker='o',
        markersize=3,
        markeredgewidth=0,
        ax=ax3,
        linewidth=0.7,
    )

# Customize third subplot
ax3.set_xticklabels(ax3.get_xticklabels(), rotation=45, ha='right', fontsize=10)  # Increased font size
ax3.set_xlabel('Developmental stage', fontsize=11)  # Increased font size
ax3.set_ylabel('Normalized expression', fontsize=11)  # Increased font size
ax3.set_title('Detailed Expression Pattern', fontsize=12)  # Increased font size

# Add legend to third plot
legend_handles = [mpatches.Patch(color=colors[gene_set], label=legend_labels[gene_set]) for gene_set in specific_gene_sets]
ax3.legend(
    handles=legend_handles,
    title='Gene Sets',
    bbox_to_anchor=(0.5, -0.5),  # Adjusted for figure height of 8
    loc='upper center',
    frameon=True
)

# Finalize and save
plt.tight_layout()

# Save the figure
save_path = f'{data_path}/../Outputs/gene_trajectories_pathways_in_GP2/Combined_PD_pathways_trajectory_rank{rank_limit}.svg'
plt.savefig(save_path, format='svg', dpi=300, bbox_inches='tight')

plt.show()

# Close the figure to free up memory
plt.close(fig)













## Option 3.3: Box plots comparing prenatal vs postnatal expression for each gene set
# Define the specific gene sets we want to plot
specific_gene_sets = [
    'Mitochondrial_genes_1',
    'Lysosomal_direct_genes_2',
    'Autophagy_direct_genes_4',
    'GP2_except_Mitochondrial_Lysosomal_Autophagy_direct_genes_9'
]

# Define labels for legend (shorter names)
legend_labels = {
    'Mitochondrial_genes_1': 'Mitochondrial',
    'Lysosomal_direct_genes_2': 'Lysosomal',
    'Autophagy_direct_genes_4': 'Autophagy',
    'GP2_except_Mitochondrial_Lysosomal_Autophagy_direct_genes_9': 'Other PD risk genes'
}

# Manually define colors for each gene set and developmental period
# Format: (gene_set_label, period): color in HEX or RGB
color_mapping = {
    ('Mitochondrial', 'Prenatal'): '#98BEF7',  # Light blue
    ('Mitochondrial', 'Postnatal'): '#0F52BA',  # Dark blue
    ('Lysosomal', 'Prenatal'): '#92E6A7',  # Light green
    ('Lysosomal', 'Postnatal'): '#155D27',  # Dark green
    ('Autophagy', 'Prenatal'): '#B884EC',  # Light purple
    ('Autophagy', 'Postnatal'): '#5A189A',  # Dark purple
    ('Other PD risk genes', 'Prenatal'): '#F19962',  # Light salmon
    ('Other PD risk genes', 'Postnatal'): '#E35E14'  # Dark orange/coral
}

# Create a copy of the dataframe to add prenatal/postnatal classification
df_boxplot = df.copy()

# Define the cutoff for prenatal vs postnatal (as in your code)
prenatal_cutoff = (401.72 - (4*30))

# Add a new column for developmental period
df_boxplot['developmental_period'] = np.where(df_boxplot['age_days'] < prenatal_cutoff, 'Prenatal', 'Postnatal')

# Create a new column combining gene set and developmental period for plotting
plot_data = pd.DataFrame()
for gene_set in specific_gene_sets:
    temp_df = df_boxplot[df_boxplot['ctx'] == 'cortex'][['developmental_period', gene_set]].copy()
    temp_df['gene_set'] = legend_labels[gene_set]
    temp_df = temp_df.rename(columns={gene_set: 'expression'})
    plot_data = pd.concat([plot_data, temp_df])

# Create a figure with two subplots - reduced width by another 5%
# Previous width was 10.69, now reducing by 5% more: 10.69 * 0.95 = 10.16
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(8, 5))

# Create patches for the legend using gray colors
import matplotlib.patches as mpatches
legend_patches = []
for period, color in zip(['Prenatal', 'Postnatal'], ['lightgray', 'darkgray']):
    patch = mpatches.Patch(
        color=color, 
        label=period
    )
    legend_patches.append(patch)

# First subplot: Box plots with outliers using a different approach
# Group by gene_set and developmental_period to plot each box individually
gene_sets = plot_data['gene_set'].unique()
dev_periods = ['Prenatal', 'Postnatal']

# Calculate positions for boxes - making paired boxes closer to each other
# Reduced the offset from 0.4 to 0.3 to bring pairs closer together
positions = []
labels = []
colors_list = []

for i, gene in enumerate(gene_sets):
    for j, period in enumerate(dev_periods):
        # Position: base position for gene + offset for period
        # Reduced offset from 0.4 to 0.3 to make pairs closer
        pos = i + (j * 0.3 - 0.15)
        positions.append(pos)
        labels.append(f"{gene}_{period}")
        colors_list.append(color_mapping[(gene, period)])

# Create custom box plot for first subplot
for i, gene in enumerate(gene_sets):
    for j, period in enumerate(dev_periods):
        # Position with reduced offset
        pos = i + (j * 0.3 - 0.15)
        data = plot_data[(plot_data['gene_set'] == gene) & 
                        (plot_data['developmental_period'] == period)]['expression']
        
        # Use the color from our mapping
        color = color_mapping[(gene, period)]
        
        # Draw the box plot - increased linewidth from default to 1.5
        ax1.boxplot(data, positions=[pos], widths=0.3, 
                   patch_artist=True, 
                   boxprops=dict(facecolor=color, color='black', alpha=0.7, linewidth=1.5),
                   medianprops=dict(color='black', linewidth=1.5),
                   whiskerprops=dict(color='black', linewidth=1.5),
                   capprops=dict(color='black', linewidth=1.5),
                   flierprops=dict(marker='o', markerfacecolor=color, markersize=3, 
                                  markeredgecolor='black', alpha=0.5),
                   capwidths=0)  # Set capwidths to 0 to remove horizontal lines

# Set x-ticks and labels for first subplot
ax1.set_xticks(range(len(gene_sets)))
ax1.set_xticklabels(gene_sets, rotation=45, ha='right', fontsize=9)
ax1.set_xlim(-0.5, len(gene_sets) - 0.5)

# Add legend to top left with gray colors
ax1.legend(handles=legend_patches, title='Developmental period', loc='upper left',
           prop={'size': 8},  # Smaller font size for legend text
           title_fontsize=9   # Smaller font size for legend title)
)

# Customize first subplot
ax1.set_xlabel('Gene set', fontsize=10)
ax1.set_ylabel('Normalized expression', fontsize=10)
ax1.set_title('Prenatal vs Postnatal Expression Comparison', fontsize=11)

# Second subplot: Box plots with all data points
# Create custom box plot for second subplot (same as first but without fliers)
for i, gene in enumerate(gene_sets):
    for j, period in enumerate(dev_periods):
        # Position with reduced offset
        pos = i + (j * 0.3 - 0.15)
        data = plot_data[(plot_data['gene_set'] == gene) & 
                        (plot_data['developmental_period'] == period)]['expression']
        
        # Use the color from our mapping
        color = color_mapping[(gene, period)]
        
        # Draw the box plot without fliers and with thicker lines
        ax2.boxplot(data, positions=[pos], widths=0.3, 
                   showfliers=False,
                   patch_artist=True, 
                   boxprops=dict(facecolor=color, color='black', alpha=0.7, linewidth=1.5),
                   medianprops=dict(color='black', linewidth=1.5),
                   whiskerprops=dict(color='black', linewidth=1.5),
                   capprops=dict(color='black', linewidth=1.5),
                   capwidths=0)  # Set capwidths to 0 to remove horizontal lines
        
        # Add individual data points
        y = data.values
        x = np.random.normal(pos, 0.05, size=len(y))
        ax2.scatter(x, y, s=10, alpha=0.4, c='black', edgecolor=None)

# Set x-ticks and labels for second subplot
ax2.set_xticks(range(len(gene_sets)))
ax2.set_xticklabels(gene_sets, rotation=45, ha='right', fontsize=9)
ax2.set_xlim(-0.5, len(gene_sets) - 0.5)

# Add legend to top left with gray colors
ax2.legend(handles=legend_patches, title='Developmental period', loc='upper left', 
           prop={'size': 8},  # Smaller font size for legend text
           title_fontsize=9   # Smaller font size for legend title)
)

# Customize second subplot
ax2.set_xlabel('Gene set', fontsize=10)
ax2.set_ylabel('Normalized expression', fontsize=10)
ax2.set_title('Prenatal vs Postnatal Expression (All Data Points)', fontsize=11)

# Finalize and save
plt.tight_layout()

# Save the figure
save_path = f'{data_path}/../Outputs/gene_trajectories_pathways_in_GP2/Prenatal_Postnatal_Comparison_rank{rank_limit}.svg'
plt.savefig(save_path, format='svg', dpi=300, bbox_inches='tight')

plt.show()

# Close the figure to free up memory
plt.close(fig)


