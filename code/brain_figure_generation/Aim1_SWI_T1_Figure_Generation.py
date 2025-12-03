import nilearn.plotting as plotting
import matplotlib.pyplot as plt
import numpy as np
import nibabel as nib
import pandas as pd
import os
import glob

################### READING ATLASES ###################
atlas_img = "/Users/houmanazizi/Library/Mobile Documents/com~apple~CloudDocs/Education/GitHub/Brain_Figure_Generation_Regionwise/Required_Files/Schaefer2018_200Parcels_7Networks_Xiao_2019_SubCorSeg.nii.gz"
template_img = "/Users/houmanazizi/Library/Mobile Documents/com~apple~CloudDocs/Education/GitHub/Brain_Figure_Generation_Regionwise/Required_Files/mni_icbm152_t1_tal_nlin_asym_09c.nii.gz"
template_mask = "/Users/houmanazizi/Library/Mobile Documents/com~apple~CloudDocs/Education/GitHub/Brain_Figure_Generation_Regionwise/Required_Files/mni_icbm152_t1_tal_nlin_asym_09c_mask.nii.gz"
atlas_data = nib.load(atlas_img).get_fdata()
atlas_labels = pd.read_csv('/Users/houmanazizi/Library/Mobile Documents/com~apple~CloudDocs/Education/GitHub/Brain_Figure_Generation_Regionwise/Required_Files/Labels/subcortical-labels.csv')

# Load and prepare masked template
mask_data = nib.load(template_mask).get_fdata()
template_data = nib.load(template_img).get_fdata()
template_data_masked = template_data.copy()
template_data_masked[mask_data == 0] = 0
template_img_masked = nib.Nifti1Image(template_data_masked, nib.load(template_img).affine)

# Get atlas affine transformation
atlas_affine = nib.load(atlas_img).affine

################### DEFINE PATHS ###################
regression_results_path = "/Users/houmanazizi/Library/Mobile Documents/com~apple~CloudDocs/Education/GitHub/UKB_Analysis/Outputs/PD_PRS_Regression_Results/with_all_pathwayPRS_Brain/"
figure_output_path = "/Users/houmanazizi/Library/Mobile Documents/com~apple~CloudDocs/Education/GitHub/Brain_Figure_Generation_Regionwise/Outputs/UKB_PhD_Project/Paper1/Regression_Figures_all_pathwayPRS/"

# Define sex-specific folders
sex_specific_folders = [
    'with_all_pathwayPRS_withMotion_all',
    'with_all_pathwayPRS_withMotion_female',
    'with_all_pathwayPRS_withMotion_male'
]

################### REGION MAPPING FUNCTION ###################
def create_region_mapping():
    """Create mapping from region names to atlas values"""
    region_mapping = {}
    
    def map_region(region_name):
        # Remove prefix and get the key part of the region name
        region_clean = region_name.replace('Median_MagneticSusceptibility_', '')
        
        # Map to atlas labels
        if 'Accumbens_L' in region_clean:
            return 219  # left nucleus accumbens (19 + 200)
        elif 'Accumbens_R' in region_clean:
            return 220  # right nucleus accumbens (20 + 200)
        elif 'Amygdala_L' in region_clean:
            return 221  # left amygdala (21 + 200)
        elif 'Amygdala_R' in region_clean:
            return 222  # right amygdala (22 + 200)
        elif 'Caudate_L' in region_clean:
            return 207  # left caudate (7 + 200)
        elif 'Caudate_R' in region_clean:
            return 208  # right caudate (8 + 200)
        elif 'Hippocampus_L' in region_clean:
            return 217  # left hippocampus (17 + 200)
        elif 'Hippocampus_R' in region_clean:
            return 218  # right hippocampus (18 + 200)
        elif 'Pallidum_L' in region_clean:
            return 211  # left globus pallidus externa (11 + 200)
        elif 'Pallidum_R' in region_clean:
            return 212  # right globus pallidus externa (12 + 200)
        elif 'Putamen_L' in region_clean:
            return 209  # left putamen (9 + 200)
        elif 'Putamen_R' in region_clean:
            return 210  # right putamen (10 + 200)
        elif 'SubstantiaNigra_L' in region_clean:
            return 203  # left substantia nigra (3 + 200)
        elif 'SubstantiaNigra_R' in region_clean:
            return 204  # right substantia nigra (4 + 200)
        elif 'Thalamus_L' in region_clean:
            return 215  # left thalamus (15 + 200)
        elif 'Thalamus_R' in region_clean:
            return 216  # right thalamus (16 + 200)
        else:
            return None
    
    return map_region

################### CREATE BRAIN MAPS ###################
def create_brain_map(results_df, value_column, atlas_data, atlas_affine):
    """Create a brain map with specified values for regions"""
    brain_map = np.zeros_like(atlas_data)
    
    for idx, row in results_df.iterrows():
        if not pd.isna(row['atlas_region']):
            atlas_value = int(row['atlas_region'])
            t_value = row[value_column]
            
            # Set the t-value for all voxels in this region
            brain_map[atlas_data == atlas_value] = t_value
    
    return brain_map

################### PLOTTING FUNCTION ###################
def create_plots(results_subset, output_path, plot_name):
    """Create and save all 4 plots for a given result file"""
    
    # 1. All regions with t-values
    all_regions_map = create_brain_map(results_subset, 't_value', atlas_data, atlas_affine)
    all_regions_img = nib.Nifti1Image(all_regions_map, atlas_affine)

    # 2. Only FDR significant regions (p_FDR < 0.05)
    significant_results = results_subset[results_subset['p_FDR'] < 0.05].copy()
    significant_map = create_brain_map(significant_results, 't_value', atlas_data, atlas_affine)
    significant_img = nib.Nifti1Image(significant_map, atlas_affine)
    
    # Determine symmetric color range for t-values
    all_t_values = results_subset['t_value'].values
    max_abs_t = np.max(np.abs(all_t_values))
    vmin, vmax = -max_abs_t, max_abs_t

    # Plot 1: All regions with skull
    display1 = plotting.plot_stat_map(
        stat_map_img=all_regions_img,
        bg_img=template_img,
        display_mode='mosaic',
        cut_coords=25,
        threshold=0.01,
        colorbar=True,
        title=f'{plot_name} - All Regions T-values (With Skull)',
        cmap='RdYlBu_r',
        annotate=False,
        draw_cross=False,
        vmin=vmin,
        vmax=vmax
    )
    display1.savefig(os.path.join(output_path, f'{plot_name}_Figure_withskull_noThreshold.svg'), dpi=300)
    display1.close()

    # Plot 2: All regions without skull
    display2 = plotting.plot_stat_map(
        stat_map_img=all_regions_img,
        bg_img=template_img_masked,
        display_mode='mosaic',
        cut_coords=25,
        threshold=0.01,
        colorbar=True,
        title=f'{plot_name} - All Regions T-values (Skull Removed)',
        cmap='RdYlBu_r',
        annotate=False,
        draw_cross=False,
        vmin=vmin,
        vmax=vmax
    )
    display2.savefig(os.path.join(output_path, f'{plot_name}_Figure_masked_noThreshold.svg'), dpi=300)
    display2.close()

    # Version 2: Only FDR significant regions
    if len(significant_results) > -1: # or 0 to not plot empty brins
        sig_t_values = significant_results['t_value'].values
        # max_abs_sig_t = np.max(np.abs(sig_t_values))
        if len(sig_t_values) == 0: # for empty plots
            max_abs_sig_t = 999 
        else:
            max_abs_sig_t = np.max(np.abs(sig_t_values))
        sig_vmin, sig_vmax = -max_abs_sig_t, max_abs_sig_t
        
        # Plot 3: Significant regions with skull
        display3 = plotting.plot_stat_map(
            stat_map_img=significant_img,
            bg_img=template_img,
            display_mode='mosaic',
            cut_coords=25,
            threshold=0.01,
            colorbar=True,
            title=f'{plot_name} - FDR Significant Regions T-values (With Skull)',
            cmap='RdYlBu_r',
            annotate=False,
            draw_cross=False,
            vmin=sig_vmin,
            vmax=sig_vmax
        )
        display3.savefig(os.path.join(output_path, f'{plot_name}_Figure_withskull.svg'), dpi=300)
        display3.close()
        
        # Plot 4: Significant regions without skull
        display4 = plotting.plot_stat_map(
            stat_map_img=significant_img,
            bg_img=template_img_masked,
            display_mode='mosaic',
            cut_coords=25,
            threshold=0.01,
            colorbar=True,
            title=f'{plot_name} - FDR Significant Regions T-values (Skull Removed)',
            cmap='RdYlBu_r',
            annotate=False,
            draw_cross=False,
            vmin=sig_vmin,
            vmax=sig_vmax
        )
        display4.savefig(os.path.join(output_path, f'{plot_name}__Figure_masked.svg'), dpi=300)
        display4.close()
        
        plt.close('all')  # Close any remaining matplotlib figures

        return len(significant_results)
    else:
        print(f"No FDR significant regions found for {plot_name}")
        return 0

################### MAIN LOOP ###################
# Initialize region mapping function
map_region = create_region_mapping()

# Loop through sex-specific folders
for sex_folder in sex_specific_folders:
    sex = sex_folder.split('_')[-1] if sex_folder != 'with_all_pathwayPRS_withMotion_all' else 'all'
    prs_folder_path = os.path.join(regression_results_path, sex_folder)
    
    # Get all PRS folders
    prs_folders = [f for f in os.listdir(prs_folder_path) if os.path.isdir(os.path.join(prs_folder_path, f))]
    
    for current_prs in prs_folders:
        print(f"Processing SWI results for {current_prs} in {sex_folder}")
        
        # Find SWI CSV files
        prs_folder_full = os.path.join(prs_folder_path, current_prs)
        csv_files = glob.glob(os.path.join(prs_folder_full, 'SWI_vs_*.csv'))
        
        if not csv_files:
            print(f"No SWI CSV files found in {prs_folder_full}")
            continue
            
        for csv_file in csv_files:
            # Read the results
            current_results = pd.read_csv(csv_file)
            
            # Extract relevant columns from results
            results_subset = current_results[['Region', 't_value', 'p_FDR']].copy()
            
            # Create region mapping
            results_subset['atlas_region'] = results_subset['Region'].apply(map_region)
            
            # Remove rows where mapping failed
            # results_subset = results_subset.dropna(subset=['atlas_region'])
            
            '''if len(results_subset) == 0:
                print(f"No regions could be mapped for {csv_file}")
                continue'''
            
            # Create output directory
            output_dir = os.path.join(figure_output_path, sex_folder, current_prs)
            os.makedirs(output_dir, exist_ok=True)
            
            # Create plot name from filename
            base_filename = os.path.basename(csv_file).replace('.csv', '')
            base_filename = base_filename[:-8]
            
            # Generate plots
            num_significant = create_plots(results_subset, output_dir, base_filename)
            
            print(f"Generated 4 brain plots for {base_filename}")
            print(f"Number of significant regions: {num_significant}")
            print(f"Total regions mapped: {len(results_subset)}")
            print(f"Plots saved in: {output_dir}")
            print("-" * 50)

print("All plots generated successfully!")
