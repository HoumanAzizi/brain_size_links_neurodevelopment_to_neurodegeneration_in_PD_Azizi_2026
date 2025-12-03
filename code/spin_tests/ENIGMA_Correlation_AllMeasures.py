"""
A script to relate PRS*SA/CT and ENIGMA-PD maps with spatial null permutation testing.
@author: Andrew Vo (Original version) & Houman Azizi (Final version)
"""

import os
import glob
import numpy as np
import pandas as pd
from netneurotools import freesurfer, stats
from scipy.stats import pearsonr, spearmanr
import pyreadr




# ------------- CONFIGURATION -------------

base_sa_path = '/Users/houmanazizi/Library/Mobile Documents/com~apple~CloudDocs/Education/GitHub/PD_Map_Correlations_Spin_Test/'
base_prs_path = '/Users/houmanazizi/Library/Mobile Documents/com~apple~CloudDocs/Education/GitHub/UKB_Analysis/Outputs/PD_PRS_Regression_Results/with_all_pathwayPRS_Brain/'
sex_specific_folders = [
    'with_all_pathwayPRS_withMotion_all',
    'with_all_pathwayPRS_withMotion_female',
    'with_all_pathwayPRS_withMotion_male']

maps = ['group_avg', 'hy_1', 'hy_2', 'hy_3', 'hy_4']

n_spins = 10000

annot_lh = base_sa_path + 'Data/DKT_annotations/lh.DKT.annot'
annot_rh = base_sa_path + 'Data/DKT_annotations/rh.DKT.annot'

# Load ENIGMA-PD atrophy maps
enigma_sa = pd.read_csv(base_sa_path + 'Data/ENIGMA_maps/enigma_sa.csv')
enigma_th = pd.read_csv(base_sa_path + 'Data/ENIGMA_maps/enigma_th.csv')

# Find centroids for spin test
parcel_centroids, parcel_hemi = freesurfer.find_parcel_centroids(
    lhannot=annot_lh, rhannot=annot_rh, version='fsaverage', surf='sphere', method='surface'
)
parcel_hemi = parcel_hemi == 1
# Create spins
spins = stats.gen_spinsamples(parcel_centroids, parcel_hemi, method='vasa', n_rotate=n_spins, seed=1234)




# ------------- SPIN TEST FUNCTION -------------

def corr_spin(x, y, spins, nspins, method):
    null = np.zeros((nspins,))

    if method == 'pearson':
        rho, pval = pearsonr(x, y)
        for i in range(nspins):
            null[i], _ = pearsonr(x, y[spins[:, i]])
    elif method == 'spearman':
        rho, pval = spearmanr(x, y)
        for i in range(nspins):
            null[i], _ = spearmanr(x, y[spins[:, i]])

    pval = (1 + sum(abs((null - np.mean(null))) > abs((rho - np.mean(null))))) / (nspins + 1)
    return rho, pval, null




# ------------- SAVE RESULTS FUNCTIONS -------------

def save_nulls_wide(results, outfile):
    # Prepare identifier columns
    keys = ['Sex', 'PD_PRS', 'Map', 'Correlation']
    # Find number of spins (assuming all nulls are the same length)
    n_spins = len(results[0]['Null_Distribution']) if results else 0
    spin_cols = [f'Spin_{i}' for i in range(n_spins)]
    # Build rows
    records = []
    for res in results:
        row = {k: res[k] for k in keys}
        nulls = res['Null_Distribution']
        for i, val in enumerate(nulls):
            row[f'Spin_{i}'] = val
        records.append(row)
    # Create DataFrame
    df_nulls = pd.DataFrame(records)
    df_nulls.to_csv(outfile, index=False)

def save_nulls_rds(results, outfile):
    # Prepare identifier columns
    keys = ['Sex', 'PD_PRS', 'Map', 'Correlation']
    n_spins = len(results[0]['Null_Distribution']) if results else 0
    spin_cols = [f'Spin_{i}' for i in range(n_spins)]
    records = []
    for res in results:
        row = {k: res[k] for k in keys}
        nulls = res['Null_Distribution']
        for i, val in enumerate(nulls):
            row[f'Spin_{i}'] = val
        records.append(row)
    df_nulls = pd.DataFrame(records)
    # Save as RDS
    pyreadr.write_rds(outfile, df_nulls)




# ------------- MAIN ANALYSIS -------------

# Create base output directory
base_output_dir = os.path.join(os.getcwd(), 'Outputs')
os.makedirs(base_output_dir, exist_ok=True)

# Define file patterns to process
file_patterns = ['SA_DKT_*.csv', 'CT_DKT_*.csv']

for file_pattern in file_patterns:
    # Determine if we're working with SA or CT PRS data
    prs_measure = file_pattern.split('_')[0]  # 'SA' or 'CT'
    
    # Create measure-specific output directory
    output_dir = os.path.join(base_output_dir, f"{prs_measure}_PRS_Results")
    os.makedirs(output_dir, exist_ok=True)
    
    # Initialize results containers
    results_sa = []  # ENIGMA SA vs PRS
    results_th = []  # ENIGMA CT vs PRS
    
    for sex_folder in sex_specific_folders:
        sex = sex_folder.split('_')[-1] if sex_folder != 'with_all_pathwayPRS_withMotion_all' else 'all'
        prs_folder_path = os.path.join(base_prs_path, sex_folder)
        prs_folders = [f for f in os.listdir(prs_folder_path) if os.path.isdir(os.path.join(prs_folder_path, f))]
        
        for prs_folder in prs_folders:
            PRS_name = prs_folder

            print(f"---Analyzing {PRS_name} ---For {sex} ---Using {prs_measure} PRS data")

            prs_folder_full = os.path.join(prs_folder_path, prs_folder)
            csv_files = glob.glob(os.path.join(prs_folder_full, file_pattern))
            if not csv_files:
                print(f'No {file_pattern} found in {prs_folder_full}')
                continue
            prs_csv = csv_files[0]
            
            # Load and process PRS data
            prs_df = pd.read_csv(prs_csv)
            if 'Region' in prs_df.columns:
                prs_df['Region'] = prs_df['Region'].apply(lambda x: x.split('_')[-1] + '_' + '_'.join(x.split('_')[:-1]))
                prs_df = prs_df.rename(columns={'Region': 'region'})
            prs_df = prs_df.set_index('region').reindex(enigma_sa['region']).reset_index()
            prs_vals = prs_df['effect_size'].values
            
            # Surface Area - ENIGMA SA vs PRS
            for mapname in maps:
                enigma_map = enigma_sa[mapname].values
                for corr_method in ['pearson', 'spearman']:
                    rho, pval, null = corr_spin(prs_vals, enigma_map, spins, n_spins, method=corr_method)
                    results_sa.append({
                        'Sex': sex,
                        'PD_PRS': PRS_name,
                        'Map': mapname,
                        'Correlation': corr_method,
                        'Rho': rho,
                        'Spin_P': pval,
                        'Null_Distribution': null
                    })
            
            # Thickness - ENIGMA CT vs PRS
            for mapname in maps:
                enigma_map = enigma_th[mapname].values
                for corr_method in ['pearson', 'spearman']:
                    rho, pval, null = corr_spin(prs_vals, enigma_map, spins, n_spins, method=corr_method)
                    results_th.append({
                        'Sex': sex,
                        'PD_PRS': PRS_name,
                        'Map': mapname,
                        'Correlation': corr_method,
                        'Rho': rho,
                        'Spin_P': pval,
                        'Null_Distribution': null
                    })
    
    # Save results with appropriate naming
    # Convert results to dataframes
    df_sa = pd.DataFrame(results_sa)
    df_th = pd.DataFrame(results_th)
    
    # Save summary results (without nulls)
    sa_results_file = os.path.join(output_dir, f'SpinTest_ENIGMA_SA_vs_{prs_measure}_PRS_Results.csv')
    th_results_file = os.path.join(output_dir, f'SpinTest_ENIGMA_CT_vs_{prs_measure}_PRS_Results.csv')
    
    df_sa.drop('Null_Distribution', axis=1).to_csv(sa_results_file, index=False)
    df_th.drop('Null_Distribution', axis=1).to_csv(th_results_file, index=False)
    
    # Save null distributions in wide format
    sa_nulls_wide_file = os.path.join(output_dir, f'SpinTest_ENIGMA_SA_vs_{prs_measure}_PRS_Nulls_wide.csv')
    th_nulls_wide_file = os.path.join(output_dir, f'SpinTest_ENIGMA_CT_vs_{prs_measure}_PRS_Nulls_wide.csv')
    
    save_nulls_wide(results_sa, sa_nulls_wide_file)
    save_nulls_wide(results_th, th_nulls_wide_file)
    
    # Save null distributions in rds format for R
    sa_nulls_rds_file = os.path.join(output_dir, f'SpinTest_ENIGMA_SA_vs_{prs_measure}_PRS_Nulls_wide.rds')
    th_nulls_rds_file = os.path.join(output_dir, f'SpinTest_ENIGMA_CT_vs_{prs_measure}_PRS_Nulls_wide.rds')
    
    save_nulls_rds(results_sa, sa_nulls_rds_file)
    save_nulls_rds(results_th, th_nulls_rds_file)
    
    print(f"Completed analysis for {prs_measure} PRS data")

print("All analyses completed!")
