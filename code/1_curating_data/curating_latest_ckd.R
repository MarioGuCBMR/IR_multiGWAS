##############
#INTRODUCTION#
##############

#This is a code to curate ckd data from Che et al 2021 (European data)
#Data was downloaded 03/08/2022

#Let's go:

###########
#libraries#
###########

library(data.table)
library(tidyverse)

#######################
#Loading ckd data#
#######################

#We are gonna load the ckd from 2024. 

project_path <- "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2023/" #change it with your own path.

setwd(project_path)

ckd <- fread("raw_data/CKD_overall_EA_JW_20180223_nstud23.dbgap.txt.gz")

######################
#STEP 1: GENERATE IDs#
######################

#We are gonna set the chr_pos:

ckd$chr_pos <- paste("chr", ckd$Chr, ":", ckd$Pos_b37, sep = "")

#Let's try make the effect allele and other allele:

ckd$Allele1 <- toupper(ckd$Allele1)
ckd$Allele2 <- toupper(ckd$Allele2)

###############################
#STEP 2: CHECK FOR WEIRD STUFF#
###############################

summary(as.numeric(ckd$Chr)) #no problems

###########################
#STEP 3: REMOVE MAF < 0.01#
###########################

summary(as.numeric(ckd$Freq1)) #We need to remove the NAs, and those with EAF > 0.99 or EAF < 0.01

ckd_eaf_OK <- ckd[which(as.numeric(ckd$Freq1) > 0.01),]
ckd_eaf_OK <- ckd_eaf_OK[which(as.numeric(ckd_eaf_OK$Freq1) < 0.99),]

summary(ckd_eaf_OK$Freq1) #worked like a charm.

ckd_corrected_eaf <- ckd_eaf_OK
  
##############################
#CURATION OF THE DATA IS DONE#
##############################

#Let's check now if there are any in the MHC region.

ckd_mhc <- ckd_corrected_eaf[which(as.numeric(ckd_corrected_eaf$Chr) == 6 & as.numeric(ckd_corrected_eaf$Pos_b37) >= 26000000 & as.numeric(ckd_corrected_eaf$Pos_b37) <= 34000000),]

summary(as.numeric(ckd_mhc$Chr)) #perfect!!
summary(as.numeric(ckd_mhc$Pos_b37)) #perfect!!

ckd_end <- ckd_corrected_eaf[which(!(ckd_corrected_eaf$chr_pos%in%ckd_mhc$chr_pos)),]

#################
#Change columns!#
#################

ckd_end <- ckd_end %>%
  select(RSID, Chr, Pos_b37, Allele1, Allele2, Freq1, Effect, StdErr, "P-value", n_total_sum, chr_pos)

colnames(ckd_end) <- c("variant", "chromosome", "base_pair_location", "effect_allele", "other_allele", "effect_allele_frequency", "beta", "standard_error", "p_value", "sample_size", "chr_pos")

#########################
#We can save this data!!#
#########################

fwrite(ckd_end, "output/1_curated_data/ckd_latest_curated.txt")
