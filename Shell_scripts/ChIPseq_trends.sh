#!/bin/sh
#$ -cwd
#$ -V
#$ -l h_vmem=4G
#$ -t 1-176

module load use.own
module load homer

bw=$(sed -n "${SGE_TASK_ID}p" bw_list.txt)

bedgraph=$(echo $bw | sed s/bw$/bedGraph/)
./bigWigToBedGraph $bw $bedgraph

for k in {1..7}; do
	annotatePeaks.pl LTR2B_cluster${k}_3kb.bed hg38 -size 3000 -hist 100 -bedGraph $bedgraph > LTR2B_cluster${k}_${bedgraph}_trend.txt
done
