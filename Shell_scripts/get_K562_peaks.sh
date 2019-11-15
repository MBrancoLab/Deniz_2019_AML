#!/bin/sh
#$ -cwd
#$ -V
#$ -l h_vmem=1G
#$ -l h_rt=1:0:0
#$ -t 1-1943


##download selected peak files (lots!)
file=$(sed -n "${SGE_TASK_ID}p" K562_hg19_peaks.txt | cut -f 2)
wget $file

##rename files
enc_name=$(echo $file | rev | cut -d / -f 1 | rev)
new_name=$(sed -n "${SGE_TASK_ID}p" K562_hg19_peaks.txt | cut -f 1)
mv $enc_name $new_name-$enc_name

