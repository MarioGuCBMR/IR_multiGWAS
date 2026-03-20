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

###########################
#Loading liver_volume data#
###########################

#We are gonna load the liver_volume from 2024. 

project_path <- "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/" #change it with your own path.

setwd(project_path)

liver_volume <- fread("raw_data/summary_statistics/GCST90016666_buildGRCh37.tsv.gz") #https://www.ebi.ac.uk/gwas/studies/GCST90016673

######################
#STEP 1: GENERATE IDs#
######################

#We are gonna set the chr_pos:

liver_volume$chr_pos <- paste("chr", liver_volume$chromosome, ":", liver_volume$base_pair_location, sep = "")

#Remove alleles:

yes_vect <- c("A", "G", "C", "T")

liver_volume <- liver_volume[which(liver_volume$effect_allele%in%yes_vect & liver_volume$other_allele%in%yes_vect),]

###############################
#STEP 2: CHECK FOR WEIRD STUFF#
###############################

summary(as.numeric(liver_volume$chromosome)) #only autosomal

###########################
#STEP 3: REMOVE MAF < 0.01#
###########################

summary(as.numeric(liver_volume$effect_allele_frequency)) #We need to remove the NAs, and those with EAF > 0.99 or EAF < 0.01

liver_volume_eaf_OK <- liver_volume[which(liver_volume$effect_allele_frequency > 0.01),]
liver_volume_eaf_OK <- liver_volume_eaf_OK[which(liver_volume_eaf_OK$effect_allele_frequency < 0.99),]

summary(liver_volume_eaf_OK$effect_allele_frequency) #worked like a charm.

liver_volume_corrected_eaf <- liver_volume_eaf_OK
  
##############################
#CURATION OF THE DATA IS DONE#
##############################

#Let's check now if there are any in the MHC region.

liver_volume_mhc <- liver_volume_corrected_eaf[which(as.numeric(liver_volume_corrected_eaf$chromosome) == 6 & as.numeric(liver_volume_corrected_eaf$base_pair_location) >= 26000000 & as.numeric(liver_volume_corrected_eaf$base_pair_location) <= 34000000),]

summary(as.numeric(liver_volume_mhc$chromosome)) #perfect!!
summary(as.numeric(liver_volume_mhc$base_pair_location)) #perfect!!

liver_volume_end <- liver_volume_corrected_eaf[which(!(liver_volume_corrected_eaf$chr_pos%in%liver_volume_mhc$chr_pos)),]

#Let's adjust the data:

liver_volume_end <- liver_volume_end %>%
  select("variant_id", "chromosome", "base_pair_location", "effect_allele", "other_allele","effect_allele_frequency", "beta", "standard_error", "p_value",  "chr_pos")

colnames(liver_volume_end) <- c("variant", "chromosome", "base_pair_location", "effect_allele", "other_allele","effect_allele_frequency", "beta", "standard_error", "p_value",  "chr_pos")

liver_volume_end$sample_size <- 32860

#########################
#We can save this data!!#
#########################

fwrite(liver_volume_end, "output/1_curated_gwas/liver_volume_curated.txt")
