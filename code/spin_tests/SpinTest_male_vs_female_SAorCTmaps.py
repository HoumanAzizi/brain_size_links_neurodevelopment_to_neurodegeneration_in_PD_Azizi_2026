import os
import numpy as np
import pandas as pd
from netneurotools import freesurfer, stats
from scipy.stats import pearsonr, spearmanr
import pyreadr


# ------------- CONFIGURATION -------------

base_path = '/Users/houmanazizi/Library/Mobile Documents/com~apple~CloudDocs/Education/GitHub/UKB_Analysis/Outputs/PD_PRS_Regression_Results/with_all_pathwayPRS_Brain/'
spin_test_base = '/Users/houmanazizi/Library/Mobile Documents/com~apple~CloudDocs/Education/GitHub/PD_Map_Correlations_Spin_Test/'

# Input file
input_csv = os.path.join(base_path, 'other/male_female_CT_beta_maps.csv') # SA or CT

# Annotation files for spin test
annot_lh = os.path.join(spin_test_base, 'Data/DKT_annotations/lh.DKT.annot')
annot_rh = os.path.join(spin_test_base, 'Data/DKT_annotations/rh.DKT.annot')

# Reference file to get correct region ordering
enigma_sa = pd.read_csv(os.path.join(spin_test_base, 'Data/ENIGMA_maps/enigma_sa.csv'))

# Spin test parameters
n_spins = 10000


# ------------- LOAD DATA -------------

# Load male vs female effect sizes
df = pd.read_csv(input_csv)
print(f"Loaded data with {len(df)} regions")
print(f"Columns: {df.columns.tolist()}\n")

# Reformat region names to match ENIGMA format (if needed)
# Your data has format like "caudalmiddlefrontal_L", need to convert to match enigma_sa
if 'Region' in df.columns:
    # Convert from "regionname_L/R" to "L/R_regionname" format
    df['region'] = df['Region'].apply(lambda x: x.split('_')[-1] + '_' + '_'.join(x.split('_')[:-1]))
else:
    print("Warning: 'Region' column not found, assuming 'region' column exists")

# Reorder data to match the order in enigma_sa (which matches annotation file order)
print(f"Enigma SA has {len(enigma_sa)} regions")
print(f"Region column in enigma_sa: {enigma_sa['region'].iloc[:5].tolist()}")
print(f"Region column in your data: {df['region'].iloc[:5].tolist()}\n")

df_reordered = df.set_index('region').reindex(enigma_sa['region']).reset_index()

# Check for missing regions
missing_regions = df_reordered[df_reordered['effect_size_males'].isna()]
if len(missing_regions) > 0:
    print(f"Warning: {len(missing_regions)} regions missing from your data:")
    print(missing_regions['region'].tolist())
    print("\nDropping these regions from analysis...")
    df_reordered = df_reordered.dropna(subset=['effect_size_males', 'effect_size_females'])

print(f"Final data has {len(df_reordered)} regions for analysis\n")

# Extract effect sizes in correct order
male_effects = df_reordered['effect_size_males'].values
female_effects = df_reordered['effect_size_females'].values


# ------------- GENERATE SPINS -------------

print("Generating spin samples...")
# Find centroids for spin test
parcel_centroids, parcel_hemi = freesurfer.find_parcel_centroids(
    lhannot=annot_lh, 
    rhannot=annot_rh, 
    version='fsaverage', 
    surf='sphere', 
    method='surface'
)
parcel_hemi = parcel_hemi == 1

# If we had to drop regions, we need to adjust centroids and spins accordingly
if len(df_reordered) < len(enigma_sa):
    # Find indices of kept regions
    kept_indices = df_reordered.index.tolist()
    parcel_centroids = parcel_centroids[kept_indices]
    parcel_hemi = parcel_hemi[kept_indices]

# Create spins
spins = stats.gen_spinsamples(
    parcel_centroids, 
    parcel_hemi, 
    method='vasa', 
    n_rotate=n_spins, 
    seed=1234
)
print(f"Generated {n_spins} spin samples\n")


# ------------- SPIN TEST FUNCTION -------------

def corr_spin(x, y, spins, nspins, method):
    """
    Perform spin test for spatial correlation.
    
    Parameters:
    -----------
    x : array
        First variable (e.g., male effect sizes)
    y : array
        Second variable (e.g., female effect sizes)
    spins : array
        Spin permutations
    nspins : int
        Number of spins
    method : str
        'pearson' or 'spearman'
    
    Returns:
    --------
    rho : float
        Correlation coefficient
    pval : float
        Spin-test p-value
    null : array
        Null distribution
    """
    null = np.zeros((nspins,))

    if method == 'pearson':
        rho, _ = pearsonr(x, y)
        for i in range(nspins):
            null[i], _ = pearsonr(x, y[spins[:, i]])
    elif method == 'spearman':
        rho, _ = spearmanr(x, y)
        for i in range(nspins):
            null[i], _ = spearmanr(x, y[spins[:, i]])

    pval = (1 + sum(abs((null - np.mean(null))) > abs((rho - np.mean(null))))) / (nspins + 1)
    
    return rho, pval, null


# ------------- PERFORM ANALYSIS -------------

print("Performing spin test analysis...")
results = []

for corr_method in ['pearson', 'spearman']:
    print(f"  Computing {corr_method} correlation...")
    rho, pval, null = corr_spin(male_effects, female_effects, spins, n_spins, method=corr_method)
    
    results.append({
        'Comparison': 'Male_vs_Female',
        'Measure': 'Surface_Area',
        'Correlation_Method': corr_method,
        'Rho': rho,
        'Spin_P': pval,
        'Null_Distribution': null
    })
    
    print(f"    {corr_method.capitalize()} r = {rho:.4f}, p = {pval:.4f}")


# ------------- SAVE RESULTS -------------

# Create output directory
output_dir = os.path.join(base_path, 'other/spin_test_results')
os.makedirs(output_dir, exist_ok=True)

# Save summary results
df_results = pd.DataFrame(results)
summary_file = os.path.join(output_dir, 'Male_vs_Female_SA_SpinTest_Results.csv')
df_results.drop('Null_Distribution', axis=1).to_csv(summary_file, index=False)
print(f"\nSaved summary results to: {summary_file}")

# Save null distributions in wide format (CSV)
nulls_wide = []
for res in results:
    row = {
        'Comparison': res['Comparison'],
        'Measure': res['Measure'],
        'Correlation_Method': res['Correlation_Method']
    }
    for i, val in enumerate(res['Null_Distribution']):
        row[f'Spin_{i}'] = val
    nulls_wide.append(row)

df_nulls_wide = pd.DataFrame(nulls_wide)
nulls_wide_file = os.path.join(output_dir, 'Male_vs_Female_SA_SpinTest_Nulls_wide.csv')
df_nulls_wide.to_csv(nulls_wide_file, index=False)
print(f"Saved null distributions (wide format) to: {nulls_wide_file}")

# Save null distributions as RDS for R
nulls_rds_file = os.path.join(output_dir, 'Male_vs_Female_SA_SpinTest_Nulls_wide.rds')
pyreadr.write_rds(nulls_rds_file, df_nulls_wide)
print(f"Saved null distributions (RDS format) to: {nulls_rds_file}")

# Also save reordered input data with results merged
df_with_results = df_reordered.copy()
for res in results:
    method = res['Correlation_Method']
    df_with_results[f'{method}_rho'] = res['Rho']
    df_with_results[f'{method}_spin_p'] = res['Spin_P']

merged_file = os.path.join(output_dir, 'Male_vs_Female_SA_with_SpinTest.csv')
df_with_results.to_csv(merged_file, index=False)
print(f"Saved merged data with results to: {merged_file}")

print("\nAnalysis completed successfully!")
