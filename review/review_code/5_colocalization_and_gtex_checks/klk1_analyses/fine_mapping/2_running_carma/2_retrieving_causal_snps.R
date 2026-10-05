##############
#INTRODUCTION#
##############

#This code assesses whether any of the leads from the MR analyses are in fine-mapped credible sets in FIadjBMI by CARMA.

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)

################################################################################################
#Let's load the results first, cuz we need the indexes of the SNPs of interest to retrieve this#
################################################################################################
 
carma_res <- readRDS("/maps/projects/kilpelainen-AUDIT/people/zlc436/IR_GSEM_2025/review/review_output/4_colocalization_and_gtex_checks/klk1_analyses/all_loci/carma_res/19_48707206_49707206.RDS")
ss_df <- fread("/maps/projects/kilpelainen-AUDIT/people/zlc436/IR_GSEM_2025/review/review_output/4_colocalization_and_gtex_checks/klk1_analyses/all_loci/loci_ss_aligned/19_48707206_49707206.txt")

#So far there are not credible sets or models so..., things are not very optimistic.

#Let's get the SNPs that we are intersted in:

klk1_snps = c("rs5516", "rs601338", "rs2638281", "rs3212815", "rs4802745", "rs8106430", "rs60621230", "rs143820692")
klk1_snps = klk1_snps[which(klk1_snps%in%ss_df$variant)] #some might be out of the locus

#We have 3, let's get their indexes:

index_snps = which(ss_df$variant%in%klk1_snps)

#finally, let's get the PIP

pips <- as.numeric(unlist(carma_res[[1]][1]))[index_snps]

#pips
#3.357993e-05 1.693747e-03 6.538918e-04
