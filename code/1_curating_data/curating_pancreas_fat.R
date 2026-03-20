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

#######################
#Loading pancreas_fat data#
#######################

#We are gonna load the pancreas_fat from 2024. 

project_path <- "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/" #change it with your own path.

setwd(project_path)

pancreas_fat <- fread("raw_data/summary_statistics/GCST90016675_buildGRCh37.tsv.gz") #https://www.ebi.ac.uk/gwas/studies/GCST90016673

######################
#STEP 1: GENERATE IDs#
######################

#We are gonna set the chr_pos:

pancreas_fat$chr_pos <- paste("chr", pancreas_fat$chromosome, ":", pancreas_fat$base_pair_location, sep = "")

#Remove alleles:

yes_vect <- c("A", "G", "C", "T")

pancreas_fat <- pancreas_fat[which(pancreas_fat$effect_allele%in%yes_vect & pancreas_fat$other_allele%in%yes_vect),]

###############################
#STEP 2: CHECK FOR WEIRD STUFF#
###############################

summary(as.numeric(pancreas_fat$chromosome)) #only autosomal

###########################
#STEP 3: REMOVE MAF < 0.01#
###########################

summary(as.numeric(pancreas_fat$effect_allele_frequency)) #We need to remove the NAs, and those with EAF > 0.99 or EAF < 0.01

pancreas_fat_eaf_OK <- pancreas_fat[which(pancreas_fat$effect_allele_frequency > 0.01),]
pancreas_fat_eaf_OK <- pancreas_fat_eaf_OK[which(pancreas_fat_eaf_OK$effect_allele_frequency < 0.99),]

summary(pancreas_fat_eaf_OK$effect_allele_frequency) #worked like a charm.

pancreas_fat_corrected_eaf <- pancreas_fat_eaf_OK
  
##############################
#CURATION OF THE DATA IS DONE#
##############################

#Let's check now if there are any in the MHC region.

pancreas_fat_mhc <- pancreas_fat_corrected_eaf[which(as.numeric(pancreas_fat_corrected_eaf$chromosome) == 6 & as.numeric(pancreas_fat_corrected_eaf$base_pair_location) >= 26000000 & as.numeric(pancreas_fat_corrected_eaf$base_pair_location) <= 34000000),]

summary(as.numeric(pancreas_fat_mhc$chromosome)) #perfect!!
summary(as.numeric(pancreas_fat_mhc$base_pair_location)) #perfect!!

pancreas_fat_end <- pancreas_fat_corrected_eaf[which(!(pancreas_fat_corrected_eaf$chr_pos%in%pancreas_fat_mhc$chr_pos)),]

#Let's adjust the data:

pancreas_fat_end <- pancreas_fat_end %>%
  select("variant_id", "chromosome", "base_pair_location", "effect_allele", "other_allele","effect_allele_frequency", "beta", "standard_error", "p_value",  "chr_pos")

colnames(pancreas_fat_end) <- c("variant", "chromosome", "base_pair_location", "effect_allele", "other_allele","effect_allele_frequency", "beta", "standard_error", "p_value",  "chr_pos")

pancreas_fat_end$sample_size <- 25617

#########################
#We can save this data!!#
#########################

fwrite(pancreas_fat_end, "output/1_curated_gwas/pancreas_fat_curated.txt")
