##############
#INTRODUCTION#
##############

#This is a code to curate age of menarche

#Let's go:

###########
#libraries#
###########

library(data.table)
library(tidyverse)

##################################
#Loading posterior_thight_fat data#
##################################

#We are gonna load the posterior_thight_fat from 2024. 

project_path <- "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/" #change it with your own path.

setwd(project_path)

posterior_thight_fat <- fread("raw_data/summary_statistics/GCST90267355.tsv.gz") #https://www.ebi.ac.uk/gwas/studies/GCST90267351

#Checked -they are all in build 37

######################
#STEP 1: GENERATE IDs#
######################

#We are gonna set the chr_pos:

posterior_thight_fat$chr_pos <- paste("chr", posterior_thight_fat$chromosome, ":", posterior_thight_fat$base_pair_location, sep = "")

#Remove alleles:

yes_vect <- c("A", "G", "C", "T")

posterior_thight_fat <- posterior_thight_fat[which(posterior_thight_fat$effect_allele%in%yes_vect & posterior_thight_fat$other_allele%in%yes_vect),]

###############################
#STEP 2: CHECK FOR WEIRD STUFF#
###############################

summary(as.numeric(posterior_thight_fat$chromosome)) #only autosomal

###########################
#STEP 3: REMOVE MAF < 0.01#
###########################

summary(as.numeric(posterior_thight_fat$effect_allele_frequency)) #No EAF!

posterior_thight_fat$effect_allele_frequency <- NA

posterior_thight_fat_corrected_eaf <- posterior_thight_fat
  
##############################
#CURATION OF THE DATA IS DONE#
##############################

#Let's check now if there are any in the MHC region.

posterior_thight_fat_mhc <- posterior_thight_fat_corrected_eaf[which(as.numeric(posterior_thight_fat_corrected_eaf$chromosome) == 6 & as.numeric(posterior_thight_fat_corrected_eaf$base_pair_location) >= 26000000 & as.numeric(posterior_thight_fat_corrected_eaf$base_pair_location) <= 34000000),]

summary(as.numeric(posterior_thight_fat_mhc$chromosome)) #perfect!!
summary(as.numeric(posterior_thight_fat_mhc$base_pair_location)) #perfect!!

posterior_thight_fat_end <- posterior_thight_fat_corrected_eaf[which(!(posterior_thight_fat_corrected_eaf$chr_pos%in%posterior_thight_fat_mhc$chr_pos)),]

#Let's adjust the data:

posterior_thight_fat_end <- posterior_thight_fat_end %>%
  select("variant_id", "chromosome", "base_pair_location", "effect_allele", "other_allele","effect_allele_frequency", "beta", "standard_error", "p_value", "n",  "chr_pos")

colnames(posterior_thight_fat_end) <- c("variant", "chromosome", "base_pair_location", "effect_allele", "other_allele","effect_allele_frequency", "beta", "standard_error", "p_value",  "sample_size", "chr_pos")

#########################
#We can save this data!!#
#########################

fwrite(posterior_thight_fat_end, "output/1_curated_gwas/posterior_thight_fat_curated.txt")
