#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --time=60
######################################
#Let's creat environment with python2#
######################################

module load anaconda3/5.3.1
#conda config --add channels conda-forge
#conda create -n goshifter_env python=2.7

################################
#Let's install the dependencies#
################################

source activate goshifter_env
#conda install numpy
#conda install -c bioconda bx-python=0.8.8 #last version for python2.7!

###############################################
#We need to gzip the bed files for some reason#
###############################################

annotation=/maps/projects/kilpelainen-AUDIT/people/zlc436/IR_GSEM_2025/output/3_go_shifter/curated_atac_seq

#all data should be gzipped!

##################################
#Let's run goshifter for our data#
##################################

#Set paths:

code=/maps/projects/kilpelainen-AUDIT/people/zlc436/IR_GSEM_2025/code/5_enrichment_analyses/3_go_shifter/3_running_go_shifter/goshifter-master

cd $code #just in case we require direct paths...

input=/maps/projects/kilpelainen-AUDIT/people/zlc436/IR_GSEM_2025/output/3_go_shifter
annotation=/maps/projects/kilpelainen-AUDIT/people/zlc436/IR_GSEM_2025/output/3_go_shifter/curated_atac_seq
output=/maps/projects/kilpelainen-AUDIT/people/zlc436/IR_GSEM_2025/output/3_go_shifter/output/bmi_neg

#Make output folders:

mkdir $output
mkdir $output/metsim
mkdir $output/roadmap_systematic

#Run analysis

$code/goshifter.py --snpmap $input/input_bmi_neg.txt --annotation $annotation/peaks_adipose_tissue.bed.gz --permute 10000 --proxies $input/ld_bmi_neg.txt --out $output/metsim/adipose #P=0.0012
$code/goshifter.py --snpmap $input/input_bmi_neg.txt --annotation $annotation/peaks_sgbs_day0.bed.gz --permute 10000 --proxies $input/ld_bmi_neg.txt --out $output/metsim/sgbs_day0 #P=0.1157
$code/goshifter.py --snpmap $input/input_bmi_neg.txt --annotation $annotation/peaks_sgbs_day4.bed.gz --permute 10000 --proxies $input/ld_bmi_neg.txt --out $output/metsim/sgbs_day4 #P=0.0013
$code/goshifter.py --snpmap $input/input_bmi_neg.txt --annotation $annotation/peaks_sgbs_day14.bed.gz --permute 10000 --proxies $input/ld_bmi_neg.txt --out $output/metsim/sgbs_day14 #P=0.003

#Let's run the analyses systematically for all roadmap:
#First save all files to work with in an array

files=()
for file in $annotation/E*.bed.gz; do
    [[ -e $file ]] || continue  # Skip if no matching files
    files+=("${file%.bed.gz}")  # Remove .bed.gz and add to array
done

#The files have the total path - let's use them as they are!
#And loop over them and run the analysis:

for file in "${files[@]}"; do

   filename=$(basename $file)
   $code/goshifter.py --snpmap $input/input_bmi_neg.txt --annotation $file.bed.gz --permute 10000 --proxies $input/ld_bmi_neg.txt --out $output/roadmap_systematic/$filename
done
