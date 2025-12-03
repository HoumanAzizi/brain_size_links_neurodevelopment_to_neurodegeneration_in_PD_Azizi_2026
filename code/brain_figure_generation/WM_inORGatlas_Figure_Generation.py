import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import nibabel as nib
from nilearn import plotting
from matplotlib import cm
import os
import glob
from surfplot import Plot
from brainspace.datasets import load_parcellation
from neuromaps.datasets import fetch_fslr
from brainspace.datasets import load_conte69
import matplotlib.colors as mcolors

### Define Regression Output Directory
regression_results_path = "/Users/houmanazizi/Library/Mobile Documents/com~apple~CloudDocs/Education/GitHub/UKB_Analysis/Outputs/PD_PRS_Regression_Results/with_all_pathwayPRS_Brain/"
sex_specific_folders = [
    'with_all_pathwayPRS_withMotion_all',
    'with_all_pathwayPRS_withMotion_female',
    'with_all_pathwayPRS_withMotion_male']
figure_output_path = "/Users/houmanazizi/Library/Mobile Documents/com~apple~CloudDocs/Education/GitHub/Brain_Figure_Generation_Regionwise/Outputs/UKB_PhD_Project/Paper1/Regression_Figures_all_pathwayPRS/"

def process_tract_data(df):
    """Process tract data to create paired left/right hemisphere data"""
    
    # Create a copy and clean the data
    data = df[['Tract_file_name', 'Tract_full_name', 'Tract_category', 't_value', 'p_FDR']].copy()
    
    # Extract hemisphere information and base tract name
    data['hemisphere'] = data['Tract_full_name'].str.extract(r'\(([LR])\)$')
    data['base_tract_name'] = data['Tract_full_name'].str.replace(r'\s*\([LR]\)$', '', regex=True)
    
    # Handle tracts without hemisphere specification (like CC tracts)
    data['hemisphere'] = data['hemisphere'].fillna('bilateral')
    data.loc[data['hemisphere'] == 'bilateral', 'base_tract_name'] = data.loc[data['hemisphere'] == 'bilateral', 'Tract_full_name']
    
    # Create paired data structure
    paired_data = []
    
    # Get unique base tract names
    unique_tracts = data['base_tract_name'].unique()
    
    for tract in unique_tracts:
        tract_subset = data[data['base_tract_name'] == tract]
        
        # Get category (should be same for all hemispheres of a tract)
        category = tract_subset['Tract_category'].iloc[0]
        
        # Initialize values
        left_t = right_t = 0
        left_p = right_p = 1
        is_central = False
        
        if len(tract_subset) == 1:
            # Bilateral tract (like CC) - use same value for both sides
            if tract_subset['hemisphere'].iloc[0] == 'bilateral':
                left_t = right_t = tract_subset['t_value'].iloc[0]
                left_p = right_p = tract_subset['p_FDR'].iloc[0]
                is_central = True
            else:
                # Single hemisphere tract
                if tract_subset['hemisphere'].iloc[0] == 'L':
                    left_t = tract_subset['t_value'].iloc[0]
                    left_p = tract_subset['p_FDR'].iloc[0]
                else:
                    right_t = tract_subset['t_value'].iloc[0]
                    right_p = tract_subset['p_FDR'].iloc[0]
        else:
            # Paired hemispheres
            for _, row in tract_subset.iterrows():
                if row['hemisphere'] == 'L':
                    left_t = row['t_value']
                    left_p = row['p_FDR']
                elif row['hemisphere'] == 'R':
                    right_t = row['t_value']
                    right_p = row['p_FDR']
        
        # Add (C) suffix for central tracts
        display_name = f"{tract} (C)" if is_central else tract
        
        paired_data.append({
            'base_tract_name': display_name,
            'category': category,
            'left_t': left_t,
            'right_t': right_t,
            'left_p': left_p,
            'right_p': right_p,
            'is_central': is_central
        })
    
    return pd.DataFrame(paired_data)

def create_spaced_layout(df):
    """Create layout with spacing between categories"""
    categories = df['category'].unique()
    spaced_data = []
    current_y = 0
    
    for cat in categories:
        cat_data = df[df['category'] == cat].copy()
        cat_data['y_position'] = range(current_y, current_y + len(cat_data))
        cat_data['category_start'] = current_y
        cat_data['category_end'] = current_y + len(cat_data) - 1
        cat_data['category_middle'] = (current_y + current_y + len(cat_data) - 1) / 2
        
        spaced_data.append(cat_data)
        current_y += len(cat_data) + 1
    
    return pd.concat(spaced_data, ignore_index=True)

def plot_horizontal_barplot(df, measure, threshold, output_path):
    """Create horizontal bar plot for tract data"""
    
    # Process the data
    paired_df = process_tract_data(df)
    
    # Sort by category and tract name
    paired_df = paired_df.sort_values(['category', 'base_tract_name']).reset_index(drop=True)
    
    # Create spaced layout
    spaced_df = create_spaced_layout(paired_df)
    
    # Define colors for categories
    Color_set_6 = ['#357266', '#C2D693', '#386FA4', '#E3C0FF', '#99DFE2', '#FC9E4F']
    category_colors = dict(zip(paired_df['category'].unique(), Color_set_6))
    
    # Create the plot
    total_height = spaced_df['y_position'].max() + 1
    fig, ax = plt.subplots(figsize=(9.1, total_height * 0.4))
    
    # Remove box around plot
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    
    # Plot parameters
    bar_height = 0.5
    bar_spacing = 0.02
    left_bar_offset = -bar_height/2 - bar_spacing/2
    right_bar_offset = bar_spacing/2
    
    # Plot bars
    for _, row in spaced_df.iterrows():
        category = row['category']
        base_color = category_colors[category]
        y_pos = row['y_position']
        
        # Determine colors based on significance
        left_color = base_color if row['left_p'] <= 0.05 else 'lightgray'
        right_color = base_color if row['right_p'] <= 0.05 else 'lightgray'
        
        if row['is_central']:
            # For central tracts (bilateral), show only one bar in the middle
            ax.barh(y_pos, row['right_t'], left=0, height=0.35, 
                    color=right_color, alpha=0.8, edgecolor='white', linewidth=1)
        else:
            # For paired tracts, show both bars
            ax.barh(y_pos + left_bar_offset, row['left_t'], left=0, height=bar_height/2.2, 
                    color=left_color, alpha=0.8, edgecolor='white', linewidth=1)
            ax.barh(y_pos + right_bar_offset, row['right_t'], left=0, height=bar_height/2.2, 
                    color=right_color, alpha=0.8, edgecolor='white', linewidth=1)
    
    # Add vertical lines at x=0 for each category
    categories = spaced_df['category'].unique()
    for category in categories:
        cat_data = spaced_df[spaced_df['category'] == category]
        cat_start = cat_data['category_start'].iloc[0]
        cat_end = cat_data['category_end'].iloc[0]
        
        ylim = ax.get_ylim()
        y_range = ylim[0] - ylim[1]
        
        ymin_relative = ((cat_end + 0.5 - ylim[1]) / y_range)
        ymax_relative = ((cat_start - 0.5 - ylim[1]) / y_range)
        ymin_relative, ymax_relative = ymax_relative, ymin_relative
        
        ax.axvline(x=0, ymin=ymin_relative, ymax=ymax_relative, 
                  color='black', linewidth=1.5, clip_on=False)
    
    # Add threshold line
    ax.axvline(x=threshold, color='red', linestyle='--', linewidth=2, alpha=0.35)
    
    # Set axis limits
    if measure == 'MD':
        min_t_value = min(spaced_df['left_t'].min(), spaced_df['right_t'].min())
        max_t_value = max(spaced_df['left_t'].max(), spaced_df['right_t'].max())
        ax.set_xlim(min_t_value - 0.5, max_t_value + 0.5)
    else:
        max_t_value = max(spaced_df['left_t'].max(), spaced_df['right_t'].max())
        ax.set_xlim(-1, max_t_value + 0.5)
    
    # Customize the plot
    ax.set_yticks(spaced_df['y_position'])
    ax.set_yticklabels(spaced_df['base_tract_name'], fontsize=10)
    ax.set_xlabel('t-value', fontsize=12, fontweight='bold')
    
    # Add category labels
    for category in categories:
        cat_data = spaced_df[spaced_df['category'] == category]
        cat_start = cat_data['category_start'].iloc[0]
        
        category_name = category.replace(' Tracts', '') + ' Tracts'
        ax.text(-1.2, cat_start - 0.55, category_name, ha='right', va='bottom',
               fontsize=11, fontweight='bold', color=category_colors[category])
    
    # Invert y-axis
    ax.invert_yaxis()
    
    # Add legend
    from matplotlib.patches import Rectangle
    from matplotlib.lines import Line2D
    
    legend_elements = []
    legend_elements.append(Line2D([0], [0], color='red', linestyle='--', linewidth=2, alpha=0.5, 
                                 label=f'Significance threshold\n(t = ±{threshold:.2f}, p_FDR < 0.05)'))
    legend_elements.append(Rectangle((0,0), 1, 1, facecolor='gray', alpha=0.6, 
                                    label='Upper bar = Left hemisphere'))
    legend_elements.append(Rectangle((0,0), 1, 1, facecolor='gray', alpha=0.6, 
                                    label='Lower bar = Right hemisphere'))
    
    legend = ax.legend(handles=legend_elements, loc='upper right', bbox_to_anchor=(1.25, 1), 
                      fontsize=9, title='Legend')
    legend.get_title().set_fontweight('bold')
    
    # Adjust layout
    plt.tight_layout()
    plt.subplots_adjust(right=0.8)
    
    # Save the figure
    plt.savefig(output_path, format='svg', bbox_inches='tight')
    plt.close()

def create_category_brainmaps(df, measure, current_prs, output_dir):
    """Create brain maps for each tract category"""
    
    # Get the atlas template
    atlas_template_image = nib.load('/Users/houmanazizi/Library/Mobile Documents/com~apple~CloudDocs/Education/GitHub/ORG_Tracts_to_Nifti_Volumes/Data/bilateral_nifti_tracts_ICBM/T_AF_left.nii')
    atlas_template = atlas_template_image.get_fdata()
    
    # Get unique tract categories
    tract_categories = df['Tract_category'].unique()
    
    # Create brain maps for each category
    for category in tract_categories:
        print(f"Processing {measure} category: {category}")
        
        # Filter data for current category
        category_data = df[df['Tract_category'] == category].copy()
        
        # Prepare the data
        tracts_for_brainmap = category_data[['Tract_file_name', 'p_FDR', 't_value']].copy()
        tracts_for_brainmap['Tract_file_name'] = tracts_for_brainmap['Tract_file_name'].apply(lambda x: x.replace("I&P", "IP"))
        
        # Threshold t-values (set non-significant to 0)
        tracts_for_brainmap.loc[tracts_for_brainmap['p_FDR'] > 0.05, 't_value'] = 0
        
        # Create brain map for this category
        results_tmp = np.zeros_like(atlas_template)
        for idx, row in tracts_for_brainmap.iterrows():
            # Read the atlas related to current tract
            current_Tract_filename = row['Tract_file_name']
            current_tract_image = nib.load(f'/Users/houmanazizi/Library/Mobile Documents/com~apple~CloudDocs/Education/GitHub/ORG_Tracts_to_Nifti_Volumes/Data/bilateral_nifti_tracts_ICBM/{current_Tract_filename}')
            current_tract = current_tract_image.get_fdata()
            
            # Add t-value for those voxels (only if significant)
            current_tract_tf = (current_tract == 1)
            if row['t_value'] != 0:
                results_tmp[current_tract_tf] = (results_tmp[current_tract_tf] + row['t_value'])/2
            del current_tract_tf
        
        # Create NIfTI image
        results_image = nib.Nifti1Image(results_tmp, atlas_template_image.affine, atlas_template_image.header)
        
        # Clean category name for filename
        category_clean = category.replace(' ', '_').replace('/', '_')
        
        # Define colormaps based on category
        if category == 'Association Tracts':
            color_m = 'Greens'
        elif category == 'Projection Tracts':
            RdPu_light = mcolors.LinearSegmentedColormap.from_list('RdPu_light', plt.cm.RdPu(np.linspace(0, 0.7, 256)))
            try:
                plt.colormaps.unregister('RdPu_light')
            except:
                pass
            plt.colormaps.register(RdPu_light, name='RdPu_light')
            color_m = 'RdPu_light'
        elif category == 'Superficial Tracts':
            oranges_light = mcolors.LinearSegmentedColormap.from_list('oranges_light', plt.cm.Oranges(np.linspace(0, 0.6, 256)))
            try:
                plt.colormaps.unregister('oranges_light')
            except:
                pass
            plt.colormaps.register(oranges_light, name='oranges_light')
            color_m = 'oranges_light'
        elif category == 'Cerebellar Tracts':
            greens_light = mcolors.LinearSegmentedColormap.from_list('greens_light', plt.cm.Greens(np.linspace(0, 0.4, 256)))
            try:
                plt.colormaps.unregister('greens_light')
            except:
                pass
            plt.colormaps.register(greens_light, name='greens_light')
            color_m = 'greens_light'
        elif category == 'Striatal Tracts':
            blues_light = mcolors.LinearSegmentedColormap.from_list('blues_light', plt.cm.Blues(np.linspace(0, 0.6, 256)))
            try:
                plt.colormaps.unregister('blues_light')
            except:
                pass
            plt.colormaps.register(blues_light, name='blues_light')
            color_m = 'blues_light'
        elif category == 'Commissural Tracts':
            color_m = 'Blues'
        else:
            color_m = 'viridis'  # Default colormap
        
        # Plot glass brain for this category
        display = plotting.plot_glass_brain(results_image, 
                                            display_mode='ortho', 
                                            title=f'{category} - {measure} t-values of FDR Significant Tracts', 
                                            colorbar=True,
                                            alpha=0.7,
                                            plot_abs=False,
                                            annotate=False,
                                            cmap=color_m,
                                            output_file=f'{output_dir}/{measure}_vs_{current_prs}_{category_clean}_Brainmap.svg')
        
        plt.close()

### Main Loop
for sex_folder in sex_specific_folders:
    sex = sex_folder.split('_')[-1] if sex_folder != 'with_all_pathwayPRS_withMotion_all' else 'all'
    prs_folder_path = os.path.join(regression_results_path, sex_folder)
    prs_folders = [f for f in os.listdir(prs_folder_path) if os.path.isdir(os.path.join(prs_folder_path, f))]
    
    for current_prs in prs_folders:
        print(f"Processing {current_prs} for {sex}")
        
        # Create output directory
        output_dir = f'{figure_output_path}/{sex_folder}/{current_prs}'
        os.makedirs(output_dir, exist_ok=True)
        
        ########## Processing FA ##########
        prs_folder_full = os.path.join(prs_folder_path, current_prs)
        csv_files = glob.glob(os.path.join(prs_folder_full, 'FA_Tracts_*.csv'))
        prs_csv = csv_files[0]
        FA_ORG = pd.read_csv(prs_csv)
        
        # Read category names and fix the names to match FA_ORG
        tract_realnames_category = pd.read_csv('/Users/houmanazizi/Library/Mobile Documents/com~apple~CloudDocs/Education/GitHub/UKB_PLS/Data/White_Matter/ORG_Tract_Names_and_Categories_forPLS.csv', header=None)
        tract_realnames_category[0] = tract_realnames_category[0].apply(lambda x: x.replace(".", "-"))
        tract_realnames_category[0] = tract_realnames_category[0].apply(lambda x: x.replace("-nii", ".nii"))
        tract_realnames_category[0] = tract_realnames_category[0].apply(lambda x: x.replace("I-P", "I&P"))
        tract_realnames_category.columns = ['Tract_file_name', 'Tract_full_name', 'Tract_category']
        
        # Merge tract names and category with FA_ORG
        FA_ORG = pd.merge(FA_ORG, tract_realnames_category, on='Tract_file_name', how='left')
        
        # Calculate threshold
        fa_threshold = (FA_ORG[FA_ORG['p_FDR'] > 0.05]['t_value'].max() + 
                       FA_ORG[FA_ORG['p_FDR'] <= 0.05]['t_value'].min()) / 2
        
        # Create horizontal bar plot
        plot_horizontal_barplot(FA_ORG, 'FA', fa_threshold, 
                               f'{output_dir}/FA_vs_{current_prs}_Regression_Figure_v2.svg')
        
        # Create category-specific brain maps
        create_category_brainmaps(FA_ORG, 'FA', current_prs, output_dir)
        
        ########## Processing MD ##########
        csv_files = glob.glob(os.path.join(prs_folder_full, 'MD_Tracts_*.csv'))
        prs_csv = csv_files[0]
        MD_ORG = pd.read_csv(prs_csv)
        
        # Merge tract names and category with MD_ORG
        MD_ORG = pd.merge(MD_ORG, tract_realnames_category, on='Tract_file_name', how='left')
        
        # For MD, we don't need a threshold calculation since it can have negative values
        md_threshold = 0  # Or calculate similar to FA if needed
        
        # Create horizontal bar plot
        plot_horizontal_barplot(MD_ORG, 'MD', md_threshold, 
                               f'{output_dir}/MD_vs_{current_prs}_Regression_Figure_v2.svg')
        
        # Create category-specific brain maps
        create_category_brainmaps(MD_ORG, 'MD', current_prs, output_dir)
