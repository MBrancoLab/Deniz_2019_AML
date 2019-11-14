#!/bin/sh
#$ -cwd
#$ -V
#$ -l h_vmem=4G
#$ -t 1-37

module load use.own
module load homer

##get bigwig file name to process from bw_list.txt
bw=$(sed -n "${SGE_TASK_ID}p" bw_list.txt)

##convert to bedGraph
bedgraph=$(echo $bw | sed s/bw$/bedGraph/)
./bigWigToBedGraph $bw $bedgraph

##cycle through TE annotation files (+/- 1.5kb from TE centre) and get trend data
for te in LTR2B LTR2C LTR5B LTR5_Hs LTR12C LTR13A; do
	annotatePeaks.pl ${te}_hg38_3kb.bed hg38 -size 3000 -hist 10 -bedGraph $bedgraph > ${te}_${bedgraph}_trend.txt
done

