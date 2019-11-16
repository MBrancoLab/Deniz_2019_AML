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

ChIP-seq peak files from K562 ENCODE data were selected (get_peak_files.R) and downloaded (get_K562_peaks.sh). The count_overlaps.R then counts the number of LTR-TF overlaps for each family of interest. Because hg19 peak files were used, hg19 annotations of the LTRs are included (as zip file). The same script was used to find overlaps with shuffled versions of the LTR annotations (generated using bedtools shuffle; included as a zip file).

*2. Find enriched TFs*

The find_enriched_TFs.R script compares the TF overlaps with real and shuffled version of the LTRs to identify enriched TFs. It then generates a table with average enrichment values for each significant TF, and expression metrics for that TF extracted from Blueprint AML data. A heatmap with enrichment values of a few selected TFs is plotted (Figure 3A).

*3. TF ChIP-seq profiles*

HOMER was used to generate TF ChIP-seq profile data (TF_trends_heatmaps.sh), and the respective heatmaps generated using make_heatmaps.R (Figure 3B). The heatmap data are included as a zip file.

## TF motifs

*1. Motif frequency*

The FIMO tool from the MEME suite was used to find TF motifs at all LTRs of interest (fimo.sh in Shell_scripts folder). The plot_motif_frequency.R script calculates the percentage of elements within each family that have the motif, and uses this to plot the frequency of selected motifs (Figure 3C) and how these compare between DHS+ and DHS- copies (Figure 3E). The output of FIMO is included as a zip file.

*2. Enriched motifs*

The AME tool from the MEME suite was used to identify TF motifs enriched in each LTR family over shuffled sequences (ame_shuffled.sh in Shell_scripts folder). A comparison between DHS+ and DHS- copies was also performed (ame_dhs.sh in Shell_scripts folder), but not included in the paper. Fasta files of DHS+ and DHS- copies (generated with split_fa_files.R), as well as shuffled versions of the LTRs (generated using the fasta-dinucleotide-shuffle tool of MEME) are included in the repository. A summary of the two analyses was generated using merge_ame_results.R.

## CRISPRi

*1. dCas9 ChIP-seq*

ChIP-seq profiles for dCas9 at LTR2B and LTR2 elements (Figure 5A) were generated using HOMER (as in several instances above) and the LTR2_heatmaps.R script. The annotate_dCas9_peaks.R takes the dCas9 peaks from MACS2 and annotates them with LTRs and genic features (Figure 5B). It also produces a list of genes lying nearest to each dCas9 peak, which is used for the RNA-seq analysis.

*2. H3K27ac/H3K9me3 ChIP-seq*

Quantification of the ChIP-seq signal at dCas9 peaks was done using Seqmonk. This quantitation is then used by ChIPseq_quantification.R to plot the fold change in signal and highlight the LTR2/2B targets (Figure 5C and Supplementary Figure 5C).

*3. RNA-seq*

Differential expression analysis was done using DESeq2 (RNA_DESeq.R). This script also generates normalised gene expression values, and includes a function to make expression barplots for specific genes, as in Figure 5G. The expression_analysis.R script plots all gene expression values and highlights genes close to dCas9 peaks, distinguishing those at LTR2/2B elements from other targets (Figure 5F and Supplementary Figure 5D).
