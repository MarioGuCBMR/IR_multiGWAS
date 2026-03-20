##############
#INTRODUCTION#
##############

#This is a code to curate tg_hdl_ratio data from Che et al 2021 (European data)
#Data was downloaded 03/08/2022

#Let's go:

###########
#libraries#
###########

library(data.table)
library(tidyverse)

#######################
#Loading tg_hdl_ratio data#
#######################

#We are gonna load the tg_hdl_ratio from 2024. 

project_path <- "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2023/" #change it with your own path.

setwd(project_path)

tg_hdl_ratio <- fread("raw_data/GCST90295950.tsv")

###################
#Curating the data#
###################

#Let's check that the one we have downloaded is the EUR with the effect sizes reported in the Supl.Table:

check <- tg_hdl_ratio[which(tg_hdl_ratio$rs_id == "rs72786786"),] #PERFECT MATCH!!! But we have to convert the freaking pvals...

tg_hdl_ratio$p_value <- 10^-(tg_hdl_ratio$neg_log_10_p_value) #checked and it matches with male!!

######################
#STEP 1: GENERATE IDs#
######################

#We are gonna set the chr_pos:

tg_hdl_ratio$chr_pos <- paste("chr", tg_hdl_ratio$chromosome, ":", tg_hdl_ratio$base_pair_location, sep = "")

#PERFECT!!! Those t

###############################
#STEP 2: CHECK FOR WEIRD STUFF#
###############################

summary(as.numeric(tg_hdl_ratio$chromosome)) #no problems

###########################
#STEP 3: REMOVE MAF < 0.01#
###########################

summary(as.numeric(tg_hdl_ratio$effect_allele_frequency)) #We need to remove the NAs, and those with EAF > 0.99 or EAF < 0.01

tg_hdl_ratio_eaf_OK <- tg_hdl_ratio[which(tg_hdl_ratio$effect_allele_frequency > 0.01),]
tg_hdl_ratio_eaf_OK <- tg_hdl_ratio_eaf_OK[which(tg_hdl_ratio_eaf_OK$effect_allele_frequency < 0.99),]

summary(tg_hdl_ratio_eaf_OK$effect_allele_frequency) #worked like a charm.

tg_hdl_ratio_corrected_eaf <- tg_hdl_ratio_eaf_OK
  
##############################
#CURATION OF THE DATA IS DONE#
##############################

#Let's check now if there are any in the MHC region.

tg_hdl_ratio_mhc <- tg_hdl_ratio_corrected_eaf[which(as.numeric(tg_hdl_ratio_corrected_eaf$chromosome) == 6 & as.numeric(tg_hdl_ratio_corrected_eaf$base_pair_location) >= 26000000 & as.numeric(tg_hdl_ratio_corrected_eaf$base_pair_location) <= 34000000),]

summary(as.numeric(tg_hdl_ratio_mhc$chromosome)) #perfect!!
summary(as.numeric(tg_hdl_ratio_mhc$base_pair_location)) #perfect!!

tg_hdl_ratio_end <- tg_hdl_ratio_corrected_eaf[which(!(tg_hdl_ratio_corrected_eaf$chr_pos%in%tg_hdl_ratio_mhc$chr_pos)),]

colnames(tg_hdl_ratio_end) <- c("chromosome",              "base_pair_location",      "effect_allele",           "other_allele",            "beta",                    "standard_error",          "effect_allele_frequency", "neg_log_10_p_value",      "variant_id",              "variant",                   "sample_size",                       "p_value",                 "chr_pos")

#########################
#We can save this data!!#
#########################

fwrite(tg_hdl_ratio_end, "output/1_curated_data/tg_hdl_ratio_male_curated.txt")
