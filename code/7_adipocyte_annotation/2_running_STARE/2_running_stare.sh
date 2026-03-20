#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --time=30

##############################
#Let's load the modules first#
##############################

module load anaconda3/5.3.1
module load bedtools/2.31.0

#conda create -n stare_mgu 
source activate stare_mgu
#conda install -c conda-forge stare-abc

####################
#Let's set the data#
####################

input=/maps/projects/kilpelainen-AUDIT/people/zlc436/IR_GSEM_2025/output/7_functional_annotation/2_enhancer_gene_links/1_stare/input
output=/maps/projects/kilpelainen-AUDIT/people/zlc436/IR_GSEM_2025/output/7_functional_annotation/2_enhancer_gene_links/1_stare/output
gene_annotation=/maps/projects/kilpelainen-AUDIT/data/gencode/gencode.v19.annotation.gtf
code=/maps/projects/kilpelainen-AUDIT/people/zlc436/IR_GSEM_2025/code/7_adipocyte_annotation/2_running_STARE/STARE-main/Code

#Let's run this for SGBS day 0:

STARE_ABCpp -b $input/sgbs_day_0_4_stare.bed -n 4 -a $gene_annotation -o $output/sgbs_day_0
STARE_ABCpp -b $input/sgbs_day_4_4_stare.bed -n 4 -a $gene_annotation -o $output/sgbs_day_4
STARE_ABCpp -b $input/sgbs_day_14_4_stare.bed -n 4 -a $gene_annotation -o $output/sgbs_day_14

source deactivate 



