#!/bin/sh
#$ -cwd
#$ -V
#$ -l h_rt=4:0:0
#$ -l h_vmem=4G
#$ -t 1-5

bam=$(sed -n "${SGE_TASK_ID}p" bam_list.txt)
out=$(echo $bam | sed s/.bam/.gtf/)

~/stringtie/stringtie $bam —rf -G ~/scratch/Genomes/Human/GRCh38/GRCh38_gencode_v26_CTAT_lib_July192017/ctat_genome_lib_build_dir/ref_annot.gtf -o $out
