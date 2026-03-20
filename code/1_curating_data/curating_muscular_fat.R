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
#Loading anterior_thight_fat data#
##################################

#We are gonna load the anterior_thight_fat from 2024. 

project_path <- "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/" #change it with your own path.

setwd(project_path)

anterior_thight_fat <- fread("raw_data/summary_statistics/GCST90267351.tsv.gz") #https://www.ebi.ac.uk/gwas/studies/GCST90267351

#Checked -they are all in build 37

######################
#STEP 1: GENERATE IDs#
######################

#We are gonna set the chr_pos:

anterior_thight_fat$chr_pos <- paste("chr", anterior_thight_fat$chromosome, ":", anterior_thight_fat$base_pair_location, sep = "")

#Remove alleles:

yes_vect <- c("A", "G", "C", "T")

anterior_thight_fat <- anterior_thight_fat[which(anterior_thight_fat$effect_allele%in%yes_vect & anterior_thight_fat$other_allele%in%yes_vect),]

###############################
#STEP 2: CHECK FOR WEIRD STUFF#
###############################

summary(as.numeric(anterior_thight_fat$chromosome)) #only autosomal

###########################
#STEP 3: REMOVE MAF < 0.01#
###########################

summary(as.numeric(anterior_thight_fat$effect_allele_frequency)) #No EAF!

anterior_thight_fat$effect_allele_frequency <- NA

anterior_thight_fat_corrected_eaf <- anterior_thight_fat
  
##############################
#CURATION OF THE DATA IS DONE#
##############################

#Let's check now if there are any in the MHC region.

anterior_thight_fat_mhc <- anterior_thight_fat_corrected_eaf[which(as.numeric(anterior_thight_fat_corrected_eaf$chromosome) == 6 & as.numeric(anterior_thight_fat_corrected_eaf$base_pair_location) >= 26000000 & as.numeric(anterior_thight_fat_corrected_eaf$base_pair_location) <= 34000000),]

summary(as.numeric(anterior_thight_fat_mhc$chromosome)) #perfect!!
summary(as.numeric(anterior_thight_fat_mhc$base_pair_location)) #perfect!!

anterior_thight_fat_end <- anterior_thight_fat_corrected_eaf[which(!(anterior_thight_fat_corrected_eaf$chr_pos%in%anterior_thight_fat_mhc$chr_pos)),]

#Let's adjust the data:

anterior_thight_fat_end <- anterior_thight_fat_end %>%
  select("variant_id", "chromosome", "base_pair_location", "effect_allele", "other_allele","effect_allele_frequency", "beta", "standard_error", "p_value", "n",  "chr_pos")

colnames(anterior_thight_fat_end) <- c("variant", "chromosome", "base_pair_location", "effect_allele", "other_allele","effect_allele_frequency", "beta", "standard_error", "p_value",  "sample_size", "chr_pos")

#########################
#We can save this data!!#
#########################

fwrite(anterior_thight_fat_end, "output/1_curated_gwas/anterior_thight_fat_curated.txt")
