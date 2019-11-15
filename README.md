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

*2. Overlap multi-exon TSSs with LTRs*

The find_LTR_TSSs.R script finds overlap between the TSSs determined above and the LTR families of interest. It generates the plot in Figure 2A and a table of AML-associated transcripts with LTR-coupled TSSs.

## ChIP-seq

*1. Histone ChIP-seq at LTRs*

Histone ChIP-seq peak files from relevant Blueprint samples were downoaded Download with download_peak_files.R, and the overlaps with LTR families of interest determined (get_histone_overlaps.R). These overlap tables were used to plot the proportion of elements from each family that is marked by a particular histone modification (plot_proportions.R), as shown in Supplementary Figure 3A.

*2. Cluster analysis*

The define_clusters.R script was used to performd k-means clustering of the LTR-histone overlaps defined above, and to display clusters  in heatmaps, as shown in Figure 2C and Supplementary Figure S3B. HOMER was used to generate data for trend plots across LTR2B elements belonging to each of the defined clusters (ChIPseq_trends.sh in Shell_scripts folder), and the cluster_trends.R script used to draw these trend plots (Figure 2D).

*3. Cell line data*

The cell_line_H3K27ac.R script, first finds overlaps between LTR families of interest and H3K27ac ChIP-seq data generated in this study. It then finds how many of these peaks are also seen in Blueprint data (from get_histone_overlaps.R ). It also looks at overlaps with a K562 ChromHMM annotation. Finally, it generates the plot in Figure 2E.

## K562 ENCODE TF data

*1. Get LTR-TF overlaps*

Select (get_peak_files.R) and download relevant hg19 data (get_K562_peaks.sh), get overlaps with real and shuffled LTRs (upload count_overlaps.R). LTR hg19 annotations as zip file. Shuffled versions generated with bedtools shuffle.

*2. Find enriched TFs*

The find_enriched_TFs.R script compares real and shuffled overlaps to identify enriched TFs. It then generates a table with average enrichment values for each significant TF, and expression metrics for that TF extracted from Blueprint AML data. A heatmap with enrichment values of a few selected TFs is plotted (Figure 3A).

*3. TF ChIP-seq profiles*

Heatmap data from HOMER (TF_trends_heatmaps.sh) and plots made with make_heatmaps.R. Heatmap data included as a zip file.

## TF motifs

*1. Enriched motifs against random*

xx

*2. Enriched motifs in DHS+ LTRs*

Not used.

*3. Motif frequency*

FIMO

