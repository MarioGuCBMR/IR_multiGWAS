##############
#INTRODUCTION#
##############

#This is a code to curate nafld data from Che et al 2021 (European data)
#Data was downloaded 03/08/2022

#Let's go:

###########
#libraries#
###########

library(data.table)
library(tidyverse)

#######################
#Loading nafld data#
#######################

#We are gonna load the nafld from 2024. 

project_path <- "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2023/" #change it with your own path.

setwd(project_path)

nafld <- fread("raw_data/GCST90091033_buildGRCh37.tsv.gz")

######################
#STEP 1: GENERATE IDs#
######################

#We are gonna set the chr_pos:

nafld$chr_pos <- paste("chr", nafld$chromosome, ":", nafld$base_pair_location, sep = "")

###############################
#STEP 2: CHECK FOR WEIRD STUFF#
###############################

summary(as.numeric(nafld$chromosome)) #no problems

###########################
#STEP 3: REMOVE MAF < 0.01#
###########################

summary(as.numeric(nafld$effect_allele_frequency)) #We need to remove the NAs, and those with EAF > 0.99 or EAF < 0.01

nafld_eaf_OK <- nafld[which(nafld$effect_allele_frequency > 0.01),]
nafld_eaf_OK <- nafld_eaf_OK[which(nafld_eaf_OK$effect_allele_frequency < 0.99),]

summary(nafld_eaf_OK$effect_allele_frequency) #worked like a charm.

nafld_corrected_eaf <- nafld_eaf_OK
  
##############################
#CURATION OF THE DATA IS DONE#
##############################

#Let's check now if there are any in the MHC region.

nafld_mhc <- nafld_corrected_eaf[which(as.numeric(nafld_corrected_eaf$chromosome) == 6 & as.numeric(nafld_corrected_eaf$base_pair_location) >= 26000000 & as.numeric(nafld_corrected_eaf$base_pair_location) <= 34000000),]

summary(as.numeric(nafld_mhc$chromosome)) #perfect!!
summary(as.numeric(nafld_mhc$base_pair_location)) #perfect!!

nafld_end <- nafld_corrected_eaf[which(!(nafld_corrected_eaf$chr_pos%in%nafld_mhc$chr_pos)),]

nafld_end <- nafld_end %>%
  select(variant_id, chromosome, base_pair_location, effect_allele, other_allele, effect_allele_frequency, beta, standard_error, p_value, chr_pos)

colnames(nafld_end) <- c("variant", "chromosome", "base_pair_location", "effect_allele", "other_allele", "effect_allele_frequency", "beta", "standard_error", "p_value", "chr_pos")

nafld_end$sample_size <-  480698
nafld_end$sample_cases <-  41395
nafld_end$sample_controls <- 439303

#########################
#We can save this data!!#
#########################

fwrite(nafld_end, "output/1_curated_data/nafld_latest_curated.txt")
