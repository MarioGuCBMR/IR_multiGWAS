##############
#INTRODUCTION#
##############

#This is a code to curate hypertension data from Che et al 2021 (European data)
#Data was downloaded 03/08/2022

#Let's go:

###########
#libraries#
###########

library(data.table)
library(tidyverse)

#######################
#Loading hypertension data#
#######################

#We are gonna load the hypertension from 2024. 

project_path <- "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2023/" #change it with your own path.

setwd(project_path)

hypertension <- fread("raw_data/ZhuZ_30940143_ukbb.bolt_460K_selfRepWhite.doctor_highbloodpressure.assoc.gz")

######################
#STEP 1: GENERATE IDs#
######################

#We are gonna set the chr_pos:

hypertension$chr_pos <- paste("chr", hypertension$CHR, ":", hypertension$BP, sep = "")

###############################
#STEP 2: CHECK FOR WEIRD STUFF#
###############################

summary(as.numeric(hypertension$Chr)) #no problems

###########################
#STEP 3: REMOVE MAF < 0.01#
###########################

summary(as.numeric(hypertension$MAF)) #We need to remove the NAs, and those with EAF > 0.99 or EAF < 0.01

hypertension_eaf_OK <- hypertension[which(as.numeric(hypertension$MAF) > 0.01),]

hypertension_corrected_eaf <- hypertension_eaf_OK
  
##############################
#CURATION OF THE DATA IS DONE#
##############################

#Let's check now if there are any in the MHC region.

hypertension_mhc <- hypertension_corrected_eaf[which(as.numeric(hypertension_corrected_eaf$CHR) == 6 & as.numeric(hypertension_corrected_eaf$BP) >= 26000000 & as.numeric(hypertension_corrected_eaf$BP) <= 34000000),]

summary(as.numeric(hypertension_mhc$CHR)) #perfect!!
summary(as.numeric(hypertension_mhc$BP)) #perfect!!

hypertension_end <- hypertension_corrected_eaf[which(!(hypertension_corrected_eaf$chr_pos%in%hypertension_mhc$chr_pos)),]

#################
#Change columns!#
#################

colnames(hypertension_end) <- c("variant", "chromosome", "base_pair_location", "effect_allele", "other_allele", "minimum_allele_frequency", "HWEP", "INFO", "beta", "standard_error", "p_value", "chr_pos")

hypertension_end <- hypertension_end %>%
  select(variant, chromosome, base_pair_location, effect_allele, other_allele, minimum_allele_frequency, beta, standard_error, p_value, chr_pos)

#########################
#We can save this data!!#
#########################

fwrite(hypertension_end, "output/1_curated_data/hypertension_latest_curated.txt")
