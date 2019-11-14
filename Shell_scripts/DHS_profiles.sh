#!/bin/sh
#$ -cwd
#$ -V
#$ -l h_vmem=4G
#$ -t 1-6

module load samtools
module load use.own
module load homer

##get bam file name to process from bam_list.txt (bam files not included in the repository)
bam=$(sed -n "${SGE_TASK_ID}p" bam_list.txt)

set=$(echo $bam | sed s/.bam//)

samtools view $bam > $set.sam

makeTagDirectory $set-tags/ $set.sam -format sam -unique

##cycle through TE annotation files (+/- 1.5kb from TE centre) and get heatmap data
for te in LTR2B LTR2C LTR5B LTR5_Hs LTR12C LTR13A; do
	annotatePeaks.pl ../Annotations/${te}_hg38_3kb.bed hg38 -size 3000 -hist 10 -ghist -d $set-tags/ > $te-$set-heatmap.txt
done
