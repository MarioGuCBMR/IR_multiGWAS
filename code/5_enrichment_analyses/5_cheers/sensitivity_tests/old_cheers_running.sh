#!/bin/bash
#SBATCH --cpus-per-task=1
#SBATCH --time=30

###########################
#Load modules and the data#
###########################

module load anaconda3/5.3.1
source activate cheers_updated

#Also let's get this ready:

code=/projects/kilpelainen-AUDIT/people/zlc436/IR_GSEM_functional_annotation/code/1_clustering/3_enrichment_analysis/2_cheers
input=/projects/kilpelainen-AUDIT/people/zlc436/IR_GSEM_functional_annotation/output/3_cheers/1_cheers_input/peak_data/combat_formatted_normalized_data
proxy=/projects/kilpelainen-AUDIT/people/zlc436/IR_GSEM_functional_annotation/output/3_cheers/1_cheers_input/variant_data
output=/projects/kilpelainen-AUDIT/people/zlc436/IR_GSEM_functional_annotation/output/3_cheers/2_cheers_output/

mkdir $output

###############################
#Let's do enrichment analysis!#
###############################

#CHEERS_computeEnrichment.py
python $code/CHEERS_computeEnrichment.py IR_proxies_all $output $input/SGBS_counts_normToMax_quantileNorm_euclideanNorm.txt $proxy/proxies_full.txt
python $code/CHEERS_computeEnrichment.py IR_proxies_cluster_1 $output $input/SGBS_counts_normToMax_quantileNorm_euclideanNorm.txt $proxy/proxies_cluster_1.txt
python $code/CHEERS_computeEnrichment.py IR_proxies_cluster_2 $output $input/SGBS_counts_normToMax_quantileNorm_euclideanNorm.txt $proxy/proxies_cluster_2.txt
python $code/CHEERS_computeEnrichment.py IR_proxies_cluster_3 $output $input/SGBS_counts_normToMax_quantileNorm_euclideanNorm.txt $proxy/proxies_cluster_3.txt
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
