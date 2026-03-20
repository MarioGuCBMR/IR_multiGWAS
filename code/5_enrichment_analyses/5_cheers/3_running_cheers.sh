#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --time=30

###########################
#Load modules and the data#
###########################

module load anaconda3/5.3.1
source activate cheers_updated

#Also let's get this ready:

code=/maps/projects/kilpelainen-AUDIT/people/zlc436/IR_GSEM_2025/code/5_enrichment_analyses/5_cheers
input=/projects/kilpelainen-AUDIT/people/zlc436/IR_GSEM_2025/output/5_cheers/1_cheers_input/peak_data/combat_formatted_normalized_data/
proxy=/projects/kilpelainen-AUDIT/people/zlc436/IR_GSEM_2025/output/5_cheers/1_cheers_input/variant_data
output=/projects/kilpelainen-AUDIT/people/zlc436/IR_GSEM_2025/output/5_cheers/2_cheers_output/

mkdir $output

###############################
#Let's do enrichment analysis!#
###############################

#CHEERS_computeEnrichment.py
python $code/CHEERS_computeEnrichment.py proxies_all $output $input/SGBS_counts_normToMax_quantileNorm_euclideanNorm.txt $proxy/proxies_full.txt
python $code/CHEERS_computeEnrichment.py proxies_bmi_ns $output $input/SGBS_counts_normToMax_quantileNorm_euclideanNorm.txt $proxy/proxies_bmi_ns.txt
python $code/CHEERS_computeEnrichment.py proxies_bmi_neg $output $input/SGBS_counts_normToMax_quantileNorm_euclideanNorm.txt $proxy/proxies_bmi_neg.txt
python $code/CHEERS_computeEnrichment.py proxies_bmi_pos $output $input/SGBS_counts_normToMax_quantileNorm_euclideanNorm.txt $proxy/proxies_bmi_pos.txt
