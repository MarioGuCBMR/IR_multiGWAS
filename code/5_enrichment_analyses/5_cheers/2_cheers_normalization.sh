#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --time=10

module load anaconda3/5.3.1
source activate cheers_update

#########################################
#STEP 1: Let's get in the working folder#
#########################################

cd /projects/kilpelainen-AUDIT/people/zlc436/IR_GSEM_2025/output/5_cheers/1_cheers_input/peak_data

code=/maps/projects/kilpelainen-AUDIT/people/zlc436/IR_GSEM_2025/code/5_enrichment_analyses/5_cheers
input=/projects/kilpelainen-AUDIT/people/zlc436/IR_GSEM_2025/output/5_cheers/1_cheers_input/peak_data/combat_formatted_normalized_data/

##########################
#Let's normalize the data#
##########################

#This version of CHEERs has positional inputs:
#prefix output
#outdir
#prefix input

python $code/CHEERS_normalize.py SGBS $input ./* #this worked! It is not a prefix, but just selecting all the data
