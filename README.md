# Deniz_2019_AML
Scripts and files from Deniz et al 2019 bioRxiv preprint:

Deniz O, Ahmed M, Todd CD, Rio-Machin A, Dawson MA, Branco MR (2019)
"Endogenous retroviruses are a source of oncogenic enhancers in acute myeloid leukemia"
bioRxiv, 772954; doi: https://doi.org/10.1101/772954


## DHS analysis

*1. Permutation test*

A permutation test was used to identify TE families enriched for DNAse hypersensitive sites (DHS). The main script is permutation_test.sh, which is a HPC cluster array job script that runs the permutation_pipeline.R script for all of the peak files listed in dhs_list.txt (peak files not included in the repository due to space restrictions). See permutation_pipeline.R for more details about other files used.

*2. TE family level analysis*

The family_level_analysis.R script takes the permutation test results from above and selects significantly DHS-enriched TE families based on the criteria described in the paper. It produces the heatmaps shown in Figure 1B and Supplementary Figure 1A. HOMER was used to generate the data for the DHS profiles in Figure 1C (DHS_profiles.sh in Shell_scripts folder). The draw_heatmaps.R script makes the respective heatmaps. The average trend plots in Supplementary Figure 1C were made with HOMER (DHS_trends.sh in Shell_scripts folder) and the draw_trends.R script.

*3. Element level analysis*

To analyse DHS patterns at individual elements from the selected TE families, a matrix of DHS-TE overlaps was made using the get_overlaps.R script. The element_level_analysis.R script takes these overlaps and, amongst other plots, checks for an association with the mutational profiles of the samples (Supplementary Figure 2). The plot_genotypes.R function is used to represent the mutational profiles as in Supplementary Figure 2A.

## Expression

*1. Generate expression table*

A table of FPKM values for selected Blueprint samples was generated using the make_expression_table.R script. This downloads the relevant files from Blueprint, merges them, and adds annotation for the nearest LTR (from selected families).

*2. Plot expression depending on DHS status*

The LTR_nearest_expression.R script merges the expression table above with DHS overlap info, and then generates the plots displayed in Figures 1D and 1E of the paper.

## Chimeric transcripts

*1. Get TSSs from spliced transcripts*

De novo transcriptomes from Blueprint's AML samples were generated using stringtie (stringtie.sh in Shell_scripts folder) on STAR-aligned RNA-seq data. Single exon transcripts were removed (remove_single_exon.R) and the coordinates of multi-exonic transcripts extracted (get_TSSs.R).

*2. Overlap multi-exon TSSs with LTRs

The find_LTR_TSSs.R script finds overlap between the TSSs determined above and the LTR families of interest. It generates the plot in Figure 2A and a table of AML-associated transcripts with LTR-coupled TSSs.
