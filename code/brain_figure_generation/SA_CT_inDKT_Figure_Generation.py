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
from neuromaps.datasets import fetch_civet
from neuromaps.datasets import fetch_fslr
from neuromaps.datasets import fetch_fsaverage
from brainspace.datasets import load_conte69
from neuromaps.transforms import civet_to_fslr



### Define Regression Output Directory
regression_results_path = "/Users/houmanazizi/Library/Mobile Documents/com~apple~CloudDocs/Education/GitHub/UKB_Analysis/Outputs/PD_PRS_Regression_Results/with_all_pathwayPRS_Brain/"
sex_specific_folders = [
    'with_all_pathwayPRS_withMotion_all',
    'with_all_pathwayPRS_withMotion_female',
    'with_all_pathwayPRS_withMotion_male']
figure_output_path = "/Users/houmanazizi/Library/Mobile Documents/com~apple~CloudDocs/Education/GitHub/Brain_Figure_Generation_Regionwise/Outputs/UKB_PhD_Project/Paper1/Regression_Figures_all_pathwayPRS/"


### Loop through measures, sexes, and PRS values to plot
Measure_list = ["SA_DKT", "CT_DKT"]

for current_measure in Measure_list:
    
    for sex_folder in sex_specific_folders:
        sex = sex_folder.split('_')[-1] if sex_folder != 'with_all_pathwayPRS_withMotion_all' else 'all'
        prs_folder_path = os.path.join(regression_results_path, sex_folder)
        prs_folders = [f for f in os.listdir(prs_folder_path) if os.path.isdir(os.path.join(prs_folder_path, f))]

        
        for current_prs in prs_folders:
            print(current_measure, current_prs)
            
            ### Read SA/CT value
            prs_folder_full = os.path.join(prs_folder_path, current_prs)
            csv_files = glob.glob(os.path.join(prs_folder_full, f'{current_measure}_*.csv'))
            prs_csv = csv_files[0]
            model_result_DKT = pd.read_csv(prs_csv, index_col=0)

            ### Create DKT atlas and Surfaces
            ## Option1: CIVET type image
            # DKT_parcellation = pd.read_csv('/Users/houmanazizi/Library/Mobile Documents/com~apple~CloudDocs/Education/GitHub/Brain_Figure_Generation_Regionwise/Required_Files/DKT_atlas/DKT_atlas_in_CIVET_vertexwise.csv', header=None)
            # surfaces = fetch_civet() # or fetch_fslr()
            # lh, rh = surfaces['inflated']
            ## Option2: conte69 type image
            DKT_parcellation = pd.read_csv('/Users/houmanazizi/Library/Mobile Documents/com~apple~CloudDocs/Education/GitHub/Brain_Figure_Generation_Regionwise/Required_Files/DKT_atlas/DKT_atlas_in_fslr_conte69_32k_vertexwise.csv', header=None)
            lh, rh = load_conte69()


            ### Separate DKT atlas into left and right hemispheres
            DKT_parcellation = DKT_parcellation.to_numpy().squeeze()
            half_length = DKT_parcellation.shape[0] // 2
            lh_parc = DKT_parcellation[:half_length].astype(np.int32)
            rh_parc = DKT_parcellation[half_length:].astype(np.int32)


            ##### Plotting SIGNIFICANT t values
            ### Prepare data to plot significant regions
            # Find significant region numbers
            region_numbers = list(model_result_DKT.loc[model_result_DKT['isSignificant'].notnull(), 'Region_Number'])
            # Set values not in region_numbers to 0
            lh_parc_sig = np.where(np.isin(lh_parc, region_numbers), lh_parc.astype(float), 0)
            rh_parc_sig = np.where(np.isin(rh_parc, region_numbers), rh_parc.astype(float), 0)
            # Set values in the atlas
            for index, row in model_result_DKT.iterrows():
                region_number = row['Region_Number']
                if pd.isnull(row['isSignificant']):
                    value = 0
                else:
                    value = row['t_value']
                
                # Update the corresponding values in lh_parc_sig and rh_parc_sig
                lh_parc_sig[lh_parc == region_number] = value
                rh_parc_sig[rh_parc == region_number] = value

            # Get max value for plot range
            max_abs_value = np.max(np.abs(np.concatenate((lh_parc_sig, rh_parc_sig))))
            max_abs_value = np.ceil(max_abs_value)

            ### Plotting significant t-values
            p = Plot(lh, rh, zoom = 1.1, size = (1500, 1500))
            p.add_layer({'left': lh_parc_sig, 'right': rh_parc_sig}, cmap='RdYlBu_r', cbar=True, color_range=[-max_abs_value, max_abs_value])
            #p.add_layer({'left': lh_parc_sig, 'right': rh_parc_sig}, cmap='gray', as_outline=True, cbar=False)
            fig = p.build()
            # Add title to the figure
            fig.suptitle(f'{current_measure} ~ PD-PRS Regression - Significant t_Values', fontsize=20)
            # Save as svg
            os.makedirs(f'{figure_output_path}/{sex_folder}/{current_prs}', exist_ok=True)
            fig.savefig(f'{figure_output_path}/{sex_folder}/{current_prs}/{current_measure}_vs_{current_prs}_Regression_Figure.svg', format='svg')
            plt.close(fig)  # Close the figure to free up memory



            ##### Plotting ALL t values
            ### Prepare data to plot all t-values
            lh_parc_all = np.zeros_like(lh_parc.astype(float))
            rh_parc_all = np.zeros_like(rh_parc.astype(float))

            for index, row in model_result_DKT.iterrows():
                region_number = row['Region_Number']
                value = row['t_value'] if not pd.isnull(row['t_value']) else 0

                lh_parc_all[lh_parc == region_number] = value
                rh_parc_all[rh_parc == region_number] = value

            # Plotting all t-values
            max_abs_value_all = np.max(np.abs(np.concatenate((lh_parc_all, rh_parc_all))))
            max_abs_value_all = np.ceil(max_abs_value_all)

            p_no_threshold = Plot(lh, rh, zoom=1.1, size=(1500, 1500))
            p_no_threshold.add_layer({'left': lh_parc_all, 'right': rh_parc_all}, cmap='RdYlBu_r', cbar=True, color_range=[-max_abs_value_all, max_abs_value_all])
            fig_no_threshold = p_no_threshold.build()
            fig_no_threshold.suptitle(f'{current_measure} ~ PD-PRS Regression - All t_Values', fontsize=20)
            fig_no_threshold.savefig(f'{figure_output_path}/{sex_folder}/{current_prs}/{current_measure}_vs_{current_prs}_Regression_Figure_noThreshold.svg', format='svg')
            plt.close(fig_no_threshold)  # Close the figure to free up memory


