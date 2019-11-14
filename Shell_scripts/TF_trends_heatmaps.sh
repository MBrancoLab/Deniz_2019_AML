#!/bin/sh
#$ -cwd
#$ -V
#$ -l h_vmem=4G
#$ -m bea

module load use.own
module load homer


while read file; do

	tf=$(echo $file | sed s/_sets.txt//)


	##trend plot

	while read set; do
		makeTagDirectory $set-tags/ $set-hg38.sam.gz -format sam -keepAll
	done < $file

	for te in LTR2B LTR2C LTR5B LTR5_Hs LTR12C LTR13A; do
		annotatePeaks.pl ${te}_hg38_3kb.bed hg38 -size 3000 -hist 10 -d *-tags/ > $te-$tf-trend.txt
	done


	##heatmap

	while read set;	do
		rm -r $set-tags
		makeTagDirectory $set-tags/ $set-hg38.sam.gz -format sam -unique
		for te in LTR2B LTR2C LTR5B LTR5_Hs LTR12C LTR13A; do
			annotatePeaks.pl ${te}_hg38_3kb.bed hg38 -size 3000 -hist 10 -ghist -d $set-tags/ > $te-$set-heatmap.txt
		done
		rm -r $set-tags
	done < $file

done < set_lists.txt
