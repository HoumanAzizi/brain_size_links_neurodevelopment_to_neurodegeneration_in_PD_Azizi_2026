import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import matplotlib.cm as cm
import matplotlib.colors as mcolors
from enigmatoolbox.plotting import plot_subcortical

### Create cmap with 0 being white
def create_precise_white_centered_rdylbu():
    """
    Create a custom colormap that closely resembles RdYlBu_r but with pure white at center.
    This preserves the yellow/orange/cyan transitions of the original colormap.
    """
    # Get the original RdYlBu_r colormap
    original_cmap = plt.cm.RdYlBu_r
    
    # Extract the colors from the original colormap
    original_colors = original_cmap(np.linspace(0, 1, 256))
    
    # Find the middle index
    mid_idx = 128
    
    # Create a smooth transition to white at the center
    # We'll modify colors near the center to ensure a smooth gradient to white
    modified_colors = original_colors.copy()
    
    # Set the exact center to pure white
    modified_colors[mid_idx] = [1.0, 1.0, 1.0, 1.0]
    
    # Create a smooth transition to white (10 steps on each side)
    transition_steps = 10
    for i in range(1, transition_steps + 1):
        # Weight for blending: closer to center = more white
        weight = i / (transition_steps + 1)
        
        # Modify colors before the center (approaching from blue side)
        idx_before = mid_idx - i
        modified_colors[idx_before] = (1-weight) * original_colors[idx_before] + weight * np.array([1.0, 1.0, 1.0, 1.0])
        
        # Modify colors after the center (approaching from red side)
        idx_after = mid_idx + i
        modified_colors[idx_after] = (1-weight) * original_colors[idx_after] + weight * np.array([1.0, 1.0, 1.0, 1.0])
    
    # Create a new colormap from the modified colors
    white_centered_cmap = mcolors.LinearSegmentedColormap.from_list(
        'precise_white_centered_rdylbu', 
        modified_colors
    )
    
    # Register the colormap with matplotlib
    plt.cm.register_cmap(name='precise_white_centered_rdylbu', cmap=white_centered_cmap)
    
    return 'precise_white_centered_rdylbu'
# Create and register the colormap
custom_cmap_name = create_precise_white_centered_rdylbu()


### Define Regression Output Directory
regression_results_path = "/Users/houmanazizi/Library/Mobile Documents/com~apple~CloudDocs/Education/GitHub/UKB_Analysis/Outputs/PD_PRS_Regression_Results/with_all_pathwayPRS_Brain/"
sex_specific_folders = [
    'with_all_pathwayPRS_withMotion_all',
    'with_all_pathwayPRS_withMotion_female',
    'with_all_pathwayPRS_withMotion_male']
figure_output_path = "/Users/houmanazizi/Library/Mobile Documents/com~apple~CloudDocs/Education/GitHub/Brain_Figure_Generation_Regionwise/Outputs/UKB_PhD_Project/Paper1/Regression_Figures_all_pathwayPRS/"





########## Plotting Subcortical Volumes - Enigma Version ##########
### Loop over each PRS and plot subcortical and WM measures
for sex_folder in sex_specific_folders:
    sex = sex_folder.split('_')[-1] if sex_folder != 'with_all_pathwayPRS_withMotion_all' else 'all'
    prs_folder_path = os.path.join(regression_results_path, sex_folder)
    prs_folders = [f for f in os.listdir(prs_folder_path) if os.path.isdir(os.path.join(prs_folder_path, f))]

    for current_prs in prs_folders:
        print(f'{current_prs} for {sex_folder}')


        # Read SV value
        prs_folder_full = os.path.join(prs_folder_path, current_prs)
        csv_files = glob.glob(os.path.join(prs_folder_full, 'SubcorticaVolume_vs_*.csv'))
        prs_csv = csv_files[0]
        SV_HarvardOxford = pd.read_csv(prs_csv)

        # Sort regions in the correct order for ENIGMA
        # Set region names in this order: alphabetical AND L before R
        # note: final order should be {'Accumbens_L', 'Amygdala_L', 'Caudate_L', 'Hippocampus_L', 'Pallidum_L', 'Putamen_L', 'Thalamus_L', 'Accumbens_R', 'Amygdala_R', 'Caudate_R', 'Hippocampus_R', 'Pallidum_R', 'Putamen_R', 'Thalamus_R'}
        region_names_L = sorted([name for name in SV_HarvardOxford['Region'] if name.endswith('_L')])
        region_names_R = sorted([name for name in SV_HarvardOxford['Region'] if name.endswith('_R')])
        region_names_L = sorted(region_names_L)
        region_names_R = sorted(region_names_R)
        region_names_sorted = region_names_L + region_names_R
        # Apply the sorting to results
        SV_HarvardOxford = SV_HarvardOxford.set_index('Region').loc[region_names_sorted].reset_index()



        ###### Plotting only significant regions
        # Get significant t_values and set the rest to 0
        SV_t_value_inOrder = SV_HarvardOxford.apply(lambda row: 0 if pd.isnull(row['isSignificant']) else row['t_value'], axis=1)
        SV_t_value_inOrder = SV_t_value_inOrder.to_numpy()

        # Generate the plot with the data
        fig_title = 'HarvardOxford SV ~ PD-PRS Regression - Significant t_Values'
        max_val = np.ceil(np.max(np.abs(SV_t_value_inOrder)) * 2) / 2
        if max_val == 0:
            max_val = 1  # Default if all values are zero

        # Plot subcortical data and save
        color_range = [((-max_val, max_val))]
        # interactive version
        '''label_text = {'bottom': ['t-value']}'''
        '''plot_subcortical(array_name=SV_t_value_inOrder, 
                        ventricles=False, 
                        cmap='RdYlBu_r',  # Using the built-in RdYlBu_r colormap
                        color_bar=True,
                        color_range=color_range,
                        label_text=label_text,
                        size=(1600, 800))'''

        # plot as high res png
        plot_subcortical(array_name=SV_t_value_inOrder, 
                    ventricles=False, 
                    cmap=custom_cmap_name, #'RdYlBu_r',
                    color_bar=True,
                    color_range=color_range,
                    size=(1600, 800),
                    screenshot=True,
                    filename=f'{figure_output_path}/{sex_folder}/{current_prs}/SV_vs_{current_prs}_Regression_Figure.png',
                    scale=5)

        # plot as high res svg
        plot_subcortical(array_name=SV_t_value_inOrder, 
                    ventricles=False, 
                    cmap=custom_cmap_name, #'RdYlBu_r',
                    color_bar=True,
                    color_range=color_range,
                    size=(1600, 800),
                    screenshot=True,
                    transparent_bg=True,
                    filename=f'{figure_output_path}/{sex_folder}/{current_prs}/SV_vs_{current_prs}_Regression_Figure.svg',
                    scale=5)
        



        ###### Plotting ALL regions
        # Get ALL t_values
        SV_t_value_inOrder = SV_HarvardOxford['t_value']
        SV_t_value_inOrder = SV_t_value_inOrder.to_numpy()

        # Generate the plot with the data
        fig_title = 'HarvardOxford SV ~ PD-PRS Regression - All t_Values'
        max_val = np.ceil(np.max(np.abs(SV_t_value_inOrder)) * 2) / 2
        if max_val == 0:
            max_val = 1  # Default if all values are zero

        # Plot subcortical data and save
        color_range = [((-max_val, max_val))]

        # plot as high res png
        plot_subcortical(array_name=SV_t_value_inOrder, 
                    ventricles=False, 
                    cmap=custom_cmap_name, #'RdYlBu_r',
                    color_bar=True,
                    color_range=color_range,
                    size=(1600, 800),
                    screenshot=True,
                    filename=f'{figure_output_path}/{sex_folder}/{current_prs}/SV_vs_{current_prs}_Regression_Figure_noThreshold.png',
                    scale=5)

        # plot as high res svg
        plot_subcortical(array_name=SV_t_value_inOrder, 
                    ventricles=False, 
                    cmap=custom_cmap_name, #'RdYlBu_r',
                    color_bar=True,
                    color_range=color_range,
                    size=(1600, 800),
                    screenshot=True,
                    transparent_bg=True,
                    filename=f'{figure_output_path}/{sex_folder}/{current_prs}/SV_vs_{current_prs}_Regression_Figure_noThreshold.svg',
                    scale=5)