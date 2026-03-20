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
#Loading liver_fat data#
#######################

#We are gonna load the liver_fat from 2024. 

project_path <- "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/" #change it with your own path.

setwd(project_path)

liver_fat <- fread("raw_data/summary_statistics/GCST90016673_buildGRCh37.tsv.gz") #https://www.ebi.ac.uk/gwas/studies/GCST90016673

######################
#STEP 1: GENERATE IDs#
######################

#We are gonna set the chr_pos:

liver_fat$chr_pos <- paste("chr", liver_fat$chromosome, ":", liver_fat$base_pair_location, sep = "")

#Remove alleles:

yes_vect <- c("A", "G", "C", "T")

liver_fat <- liver_fat[which(liver_fat$effect_allele%in%yes_vect & liver_fat$other_allele%in%yes_vect),]

###############################
#STEP 2: CHECK FOR WEIRD STUFF#
###############################

summary(as.numeric(liver_fat$chromosome)) #only autosomal

###########################
#STEP 3: REMOVE MAF < 0.01#
###########################

summary(as.numeric(liver_fat$effect_allele_frequency)) #We need to remove the NAs, and those with EAF > 0.99 or EAF < 0.01

liver_fat_eaf_OK <- liver_fat[which(liver_fat$effect_allele_frequency > 0.01),]
liver_fat_eaf_OK <- liver_fat_eaf_OK[which(liver_fat_eaf_OK$effect_allele_frequency < 0.99),]

summary(liver_fat_eaf_OK$effect_allele_frequency) #worked like a charm.

liver_fat_corrected_eaf <- liver_fat_eaf_OK
  
##############################
#CURATION OF THE DATA IS DONE#
##############################

#Let's check now if there are any in the MHC region.

liver_fat_mhc <- liver_fat_corrected_eaf[which(as.numeric(liver_fat_corrected_eaf$chromosome) == 6 & as.numeric(liver_fat_corrected_eaf$base_pair_location) >= 26000000 & as.numeric(liver_fat_corrected_eaf$base_pair_location) <= 34000000),]

summary(as.numeric(liver_fat_mhc$chromosome)) #perfect!!
summary(as.numeric(liver_fat_mhc$base_pair_location)) #perfect!!

liver_fat_end <- liver_fat_corrected_eaf[which(!(liver_fat_corrected_eaf$chr_pos%in%liver_fat_mhc$chr_pos)),]

#Let's adjust the data:

liver_fat_end <- liver_fat_end %>%
  select("variant_id", "chromosome", "base_pair_location", "effect_allele", "other_allele","effect_allele_frequency", "beta", "standard_error", "p_value",  "chr_pos")

colnames(liver_fat_end) <- c("variant", "chromosome", "base_pair_location", "effect_allele", "other_allele","effect_allele_frequency", "beta", "standard_error", "p_value",  "chr_pos")

liver_fat_end$sample_size <- 32858

#########################
#We can save this data!!#
#########################

fwrite(liver_fat_end, "output/1_curated_gwas/liver_fat_curated.txt")
