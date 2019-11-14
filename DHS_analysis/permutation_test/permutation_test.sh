#!/bin/sh
#$ -cwd
#$ -V
#$ -l h_rt=12:0:0
#$ -l h_vmem=2G
#$ -t 1-43

module load bedtools/2.26.0
module load R

##dhs_list.txt contains list of DHS peak files
file=$(sed -n "${SGE_TASK_ID}p" dhs_list.txt)

Rscript permutation_pipeline.R $file
