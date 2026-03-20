##############
#INTRODUCTION#
##############

#This is a code to curate cad data from Che et al 2021 (European data)
#Data was downloaded 03/08/2022

#Let's go:

###########
#libraries#
###########

library(data.table)
library(tidyverse)

#######################
#Loading cad data#
#######################

#We are gonna load the cad from 2024. 

project_path <- "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2023/" #change it with your own path.

setwd(project_path)

cad <- fread("raw_data/CAD_GWAS_primary_discovery_meta.tsv")

######################
#STEP 1: GENERATE IDs#
######################

#We are gonna set the chr_pos:

cad$chr_pos <- paste("chr", cad$CHR, ":", cad$BP, sep = "")

#Let's try make the effect allele and other allele:

cad$Allele1 <- toupper(cad$Allele1)
cad$Allele2 <- toupper(cad$Allele2)

###############################
#STEP 2: CHECK FOR WEIRD STUFF#
###############################

summary(as.numeric(cad$CHR)) #no problems

###########################
#STEP 3: REMOVE MAF < 0.01#
###########################

summary(as.numeric(cad$Freq1)) #We need to remove the NAs, and those with EAF > 0.99 or EAF < 0.01

cad_eaf_OK <- cad[which(cad$Freq1 > 0.01),]
cad_eaf_OK <- cad_eaf_OK[which(cad_eaf_OK$Freq1 < 0.99),]

summary(cad_eaf_OK$Freq1) #worked like a charm.

cad_corrected_eaf <- cad_eaf_OK
  
##############################
#CURATION OF THE DATA IS DONE#
##############################

#Let's check now if there are any in the MHC region.

cad_mhc <- cad_corrected_eaf[which(as.numeric(cad_corrected_eaf$CHR) == 6 & as.numeric(cad_corrected_eaf$BP) >= 26000000 & as.numeric(cad_corrected_eaf$BP) <= 34000000),]

summary(as.numeric(cad_mhc$CHR)) #perfect!!
summary(as.numeric(cad_mhc$BP)) #perfect!!

cad_end <- cad_corrected_eaf[which(!(cad_corrected_eaf$chr_pos%in%cad_mhc$chr_pos)),]

#################
#Change columns!#
#################

cad_end <- cad_end %>%
  select(MarkerName, CHR, BP, Allele1, Allele2, Freq1, Effect, StdErr, "P-value", N, chr_pos)

colnames(cad_end) <- c("variant", "chromosome", "base_pair_location", "effect_allele", "other_allele", "effect_allele_frequency", "beta", "standard_error", "p_value", "sample_size", "chr_pos")

#########################
#We can save this data!!#
#########################

fwrite(cad_end, "output/1_curated_data/cad_curated.txt")
