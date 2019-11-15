#!/bin/sh
#$ -cwd
#$ -V
#$ -l h_vmem=4G
#$ -t 1-176

module load use.own
module load homer

##get name of bigwig file to process (not included in this repository)
bw=$(sed -n "${SGE_TASK_ID}p" bw_list.txt)

##convert to bedGraph
bedgraph=$(echo $bw | sed s/bw$/bedGraph/)
./bigWigToBedGraph $bw $bedgraph

##generate data to draw trend plots across LTR2B elements (+/- 1.5kb from centre) belonging to different clusters
##clusters defined in define_clusters.R and respective annotation in Annotations folder
for k in {1..7}; do
	annotatePeaks.pl ../Annotations/LTR2B_cluster${k}_3kb.bed hg38 -size 3000 -hist 100 -bedGraph $bedgraph > LTR2B_cluster${k}_${bedgraph}_trend.txt
done
