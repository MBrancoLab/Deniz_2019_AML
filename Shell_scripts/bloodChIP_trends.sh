#!/bin/sh
#$ -cwd
#$ -V
#$ -l h_vmem=4G
#$ -t 1-8

module load samtools
module load use.own
module load homer

bam=$(sed -n "${SGE_TASK_ID}p" bam_list.txt)
set=$(echo $bam | sed s/.bam//)


#deduplicate

samtools collate -o $set-collated.bam $bam
samtools fixmate -m $set-collated.bam $set-fixmate.bam
samtools sort -o $set-sorted.bam $set-fixmate.bam
samtools markdup -r $set-sorted.bam $set-dedup.bam


#make trends

samtools view $set-dedup.bam > $set.sam

makeTagDirectory $set-tags/ $set.sam -format sam -unique

for te in LTR2B LTR2C LTR5B LTR5_Hs LTR12C LTR13A; do
	annotatePeaks.pl ${te}_hg38_3kb.bed hg38 -size 3000 -hist 10 -d $set-tags/ > $te-$set-trend.txt
done
