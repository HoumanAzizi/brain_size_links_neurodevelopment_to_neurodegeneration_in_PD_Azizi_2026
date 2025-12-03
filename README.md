# Brain size links neurodevelopment to neurodegeneration in Parkinson's disease (Azizi et al. 2026)

This repository contains the code used for creating the results presented in "Brain size links neurodevelopment to neurodegeneration in Parkinson's disease".

# Codes

### linear_regression/
`UKB_Reg_Part_0_UKB_Wrangling_and_Exclusion.R`: 
Performs comprehensive data wrangling and quality control for the UK Biobank dataset, including exclusion criteria application, PCA calculation, and preparation of brain imaging and behavioral measures for downstream analyses.

`UKB_Reg_Part_1_1_AllMeasures_vs_PDPRS_Regression.R`: 
Conducts region-wise linear regression analyses across multiple brain imaging modalities (SA, CT, subcortical volume, FA, MD, SWI, T2*) against PD-PRSs with sex and motion confound stratification.

`UKB_Reg_Part_1_2_GlobalBrainMeasures_vs_PDPRS_Regression.R`: 
Performs linear regression of global (whole-brain averaged) imaging measures against PD-PRSs, with separate analyses by sex and motion confound status.

`UKB_Reg_Part_1_3_1_BehavioralMeasures_vs_PDPRS_Regression_newMeasures_imaging.R`: 
Analyzes associations between behavioral and cognitive measures and PD-PRS in the imaging subset of UK Biobank.

`UKB_Reg_Part_1_4_1_BehavioralMeasures_vs_PDPRS_Regression_newMeasures_fullUKB.R`: 
Analyzes associations between behavioral and cognitive measures and PD-PRS in the full UK Biobank.

`UKB_Reg_Part_2_1__plot_GlobalBrainMeasures_vs_PDPRS.R`: 
Generates heatmaps and horizontal bar plots summarizing global brain measure associations with PD-PRSs, stratified by sex with effect sizes and confidence intervals.

`UKB_Reg_Part_2_2_1_plot_BehavioralMeasures_vs_PDPRS.R`: 
Creates heatmaps and bar plots visualizing behavioral measure associations with PD-PRS, including direct comparison between full UK Biobank and imaging-subset cohorts.

`UKB_Reg_Part_4_create_Tables.R`: 
Compiles regression results from all brain imaging measures into a formatted Excel workbook with multiple sheets organized by PRS variant, imaging measure, and sex, with FDR-significant results highlighted.

### brain_figure_generation/
`WM_inORGatlas_Figure_Generation.py`:
Generates visualizations of white matter tract associations with PD-PRS, including paired left/right hemisphere horizontal bar plots and category-specific brain maps for FA and MD measures.

`SWI_T1_Figure_Generation.py`:
Creates mosaic-style brain visualizations displaying regional associations between susceptibility-weighted imaging (SWI) and T1 measurements with PD-PRS, with options for displaying all regions or only FDR-significant findings.

`SubcorticalVolume_ENIGMA_Figure_Generation.py`: 
Generates ENIGMA-style subcortical volume visualizations with a custom white-centered colormap, displaying both FDR-significant and all regional t-values for subcortical structures in relation to PD-PRS.

`SA_CT_inDKT_Figure_Generation.py`: 
Produces surface-based brain visualizations of cortical surface area and thickness associations with PD-PRS using the DKT atlas.

### mendelian_randomization/
`MR_Part_0_add_N_to_GP2_GWAS.R`: 
Adds sample size information (N_case, N_proxy, N_case_proxy, N_control) to the GP2 PD GWAS summary statistics file for downstream MR analyses.

`MR_Part_1_1_SmithGWAS_vsPD_MR.R`: 
Performs two-sample Mendelian randomization analysis testing whether brain imaging measures from Smith et al. (surface area, cortical thickness, subcortical volume, white matter metrics) have a causal effect on Parkinson's disease risk using the GP2 PD GWAS as the outcome.

`MR_Part_1_2_WarrierGWAS_vsPD_MR.R`: 
Conducts two-sample MR analysis examining causal relationships between Warrier et al. brain imaging GWAS (surface area and cortical thickness) and PD risk.

`MR_Part_1_4_ZhaoGWAS_vsPD_MR.R`: 
Performs MR analysis testing whether white matter diffusion metrics (FA and MD) from Zhao et al. GWAS have causal effects on PD susceptibility.

`MR_Part_2_reverse_MR_PD_vsBrain.R`: 
Conducts reverse Mendelian randomization to test whether PD has a causal effect on various brain imaging measures, using PD GWAS as exposure and multiple brain imaging GWASes as outcomes.

`MR_Part_3_process_MR_results.R`: 
Aggregates and processes MR results across all brain measures, calculates FDR-corrected p-values, computes statistical power, and generates summary tables of significant findings.

### spin_tests/
`ENIGMA_Correlation_AllMeasures.py`: 
Performs spatial null permutation testing (spin tests) to assess correlations between polygenic risk scores and ENIGMA-PD atrophy maps while controlling for spatial autocorrelation across brain regions.

`SpinTest_male_vs_female_SAorCTmaps.py`: 
Conducts spin test analysis comparing sex differences in brain imaging effect sizes to determine whether male-female PD risk-associated patterns are statistically spatially similar.

`ENIGMA_Correlation_Plotting.R`: 
Generates heatmaps and boxplots visualizing spin test results showing correlations between pathway-specific PRS and ENIGMA-PD brain maps across disease stages and sex groups.

### pathway_gene_processing/
`prepare_Pathway_PRS_Gene_List.R`: 
Curates and organizes gene lists for mitochondrial, lysosomal, and autophagy pathways from multiple sources (MitoCarta, HADb, GO annotations) and maps them to genomic positions for PRS construction.

`match_PathwayPRS_genes_with_BrainSpan_and_GP2PD.R`: 
Matches pathway-specific genes with PD GWAS loci and BrainSpan developmental expression data, creating comprehensive annotations linking PD genetic variants to specific biological pathways.

`subset_GWAS.R`: 
Partitions the GP2 PD GWAS into pathway-specific and non-pathway SNP subsets for pathway-stratified polygenic risk score analyses.

### gene_expression_trajectories/
`Gene_Expression_Trajectories_Brainspan.py`: 
Visualizes developmental gene expression trajectories across prenatal and postnatal stages for pathway-specific gene sets from the BrainSpan dataset, including polynomial regression fits and prenatal/postnatal comparisons.

### LD_score_regression/
`plot_LDSC_results.R`: 
Creates heatmaps displaying genetic correlations between pathway-specific PD-PRS and brain imaging traits, with FDR-corrected significance indicators.