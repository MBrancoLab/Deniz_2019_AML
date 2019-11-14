# Deniz_2019_AML
Scripts and files from Deniz et al 2019 bioRxiv preprint:

Deniz O, Ahmed M, Todd CD, Rio-Machin A, Dawson MA, Branco MR (2019)
"Endogenous retroviruses are a source of oncogenic enhancers in acute myeloid leukemia"
bioRxiv, 772954; doi: https://doi.org/10.1101/772954


## DHS analysis

*1. Permutation test*

A permutation test was used to identify TE families enriched for DNAse hypersensitive sites (DHS). The main script is permutation_test.sh, which is a HPC cluster array job script that runs the permutation_pipeline.R script for all of the peak files listed in dhs_list.txt. Requires bedtools. See permutation_pipeline.R for more details about other files.


*2. TE family level analysis*

The family_level_analysis.R script takes the permutation test results from above and selects significantly DHS-enriched TE families based on the criteria described in the paper. It produces the heatmaps shown in Figure 1B and Supplementary Figure 1A.  
HOMER was used to generate the data for the DHS profiles in Figure 1C (DHS_profiles.sh in Shell_scripts folder). The draw_heatmaps.R script makes the respective heatmaps.
The average trend plots in Supplementary Figure 1C were made with HOMER (DHS_trends.sh in Shell_scripts folder) and the draw_trends.R script.

