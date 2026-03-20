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
input=/maps/projects/kilpelainen-AUDIT/people/zlc436/IR_GSEM_2025/output/5_cheers/1_cheers_input/peak_data/combat_formatted_normalized_data/
proxy=/maps/projects/kilpelainen-AUDIT/people/zlc436/IR_GSEM_2025/output/5_cheers/3_sensitivity_tests/1_cheers_input/variant_data
output=/maps/projects/kilpelainen-AUDIT/people/zlc436/IR_GSEM_2025/output/5_cheers/3_sensitivity_tests/2_cheers_output/

mkdir $output

###############################
#Let's do enrichment analysis!#
###############################

#CHEERS_computeEnrichment.py
python $code/CHEERS_computeEnrichment.py IR_proxies_cluster_lotta $output $input/SGBS_counts_normToMax_quantileNorm_euclideanNorm.txt $proxy/proxies_lotta.txt
python $code/CHEERS_computeEnrichment.py IR_proxies_cluster_lipodystrophy $output $input/SGBS_counts_normToMax_quantileNorm_euclideanNorm.txt $proxy/proxies_lipodystrophy.txt

#Let's do the rest of clusters from Suzuki:

python $code/CHEERS_computeEnrichment.py IR_proxies_cluster_beta_cell_neg $output $input/SGBS_counts_normToMax_quantileNorm_euclideanNorm.txt $proxy/proxies_beta_cell_neg.txt
python $code/CHEERS_computeEnrichment.py IR_proxies_cluster_beta_cell_pos $output $input/SGBS_counts_normToMax_quantileNorm_euclideanNorm.txt $proxy/proxies_beta_cell_pos.txt
python $code/CHEERS_computeEnrichment.py IR_proxies_cluster_residual $output $input/SGBS_counts_normToMax_quantileNorm_euclideanNorm.txt $proxy/proxies_residual.txt
python $code/CHEERS_computeEnrichment.py IR_proxies_cluster_obesity $output $input/SGBS_counts_normToMax_quantileNorm_euclideanNorm.txt $proxy/proxies_obesity.txt
python $code/CHEERS_computeEnrichment.py IR_proxies_cluster_metabolic_syndrome $output $input/SGBS_counts_normToMax_quantileNorm_euclideanNorm.txt $proxy/proxies_metabolic_syndrome.txt
python $code/CHEERS_computeEnrichment.py IR_proxies_cluster_body_fat $output $input/SGBS_counts_normToMax_quantileNorm_euclideanNorm.txt $proxy/proxies_body_fat.txt
python $code/CHEERS_computeEnrichment.py IR_proxies_cluster_liver $output $input/SGBS_counts_normToMax_quantileNorm_euclideanNorm.txt $proxy/proxies_liver.txt