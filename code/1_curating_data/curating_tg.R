##############
#INTRODUCTION#
##############

#This code is to curate lipid data from Graham 2022.

###########
#Libraries#
###########

library(data.table)
library(tidyverse)

memory.limit(size=80000000)

###################
#Loading functions#
###################

chr_parser <- function(chr_pos){
  
  #This function retrieves chromosomes in chr1 format.
  
  tmp <- strsplit(as.character(chr_pos), ":")[[1]][1]
  
  return(tmp)
  
}

pos_parser <- function(chr_pos){
  
  #This function retrieves chromosomes in chr1 format.
  
  tmp <- strsplit(as.character(chr_pos), ":")[[1]][2]
  
  return(tmp)
  
}

###############
#Loading files#
###############

project_path <- "J:/CBMR/SUN-CBMR-Kilpelainen-Group/insulin_resistance_variants_common_info/" #change it with your own path.

setwd(project_path)

tg <- fread("raw_data/logTG_INV_EUR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz")

#####################################################
#Let's check the dataframes and what do we have here#
#####################################################

head(tg)

#Here the insertions and deletions were change to D and I.
#Let's curate and do some checks in the end.

###################
#tg data analysis#
###################

################################################
#1. Now we are going to remove those SNPs that #
#Have low sample (from June Version is 10.000) #
################################################

summary(tg$N) #no NAs. Min = 797, Max = 1.2M. Works for me.

tg <- tg[which(tg$N > 10000),] #We go from 47M to 46M. 

######################################################
#2. Remove those that have a EAF > 0.99 or EAF < 0.01#
######################################################

#Let's check first for any NAs:

summary(tg$POOLED_ALT_AF) #no NAs. Interestingly, Max is 0.5, but it is not MAF, but Alternative allele frequency. This means they selected the allele with minimum allele frequency as the effect allele.

#We are gonna remove those with MAF < 0.01.

tg_eaf_OK <- tg[which(as.numeric(tg$POOLED_ALT_AF) > 0.01),] #We go from 46M to...9.9M
tg_eaf_OK <- tg_eaf_OK[which(as.numeric(tg_eaf_OK$POOLED_ALT_AF) < 0.99),] #This remained the same.

summary(tg_eaf_OK$POOLED_ALT_AF) #perfect

#From  here on we are going to have the chr:pos:

tg_eaf_OK$chr_pos <- paste("chr", tg_eaf_OK$CHROM, ":", tg_eaf_OK$POS_b37, sep="")

#########################################
#Let's check on insertions and deletions#
#########################################

#Let's load fiadjbmi variants and see what is going on:

tg_eaf_alleles <- tg_eaf_OK

#############################
#Let's remove the MHC region#
#############################

summary(tg_eaf_alleles$POS_b37) #all good, no issues here.

tg_eaf_alleles_mhc <- tg_eaf_alleles[which(as.numeric(tg_eaf_alleles$CHROM) == 6),]
tg_eaf_alleles_mhc <- tg_eaf_alleles_mhc[which(as.numeric(tg_eaf_alleles_mhc$POS_b37) >= 26000000),]
tg_eaf_alleles_mhc <- tg_eaf_alleles_mhc[which(as.numeric(tg_eaf_alleles_mhc$POS_b37) <= 34000000),]

head((tg_eaf_alleles_mhc$CHROM)) #perfect.
summary(as.numeric(tg_eaf_alleles_mhc$POS_b37)) #perfect.

tg_eaf_alleles_NO_mhc <- tg_eaf_alleles[which(!(tg_eaf_alleles$chr_pos%in%tg_eaf_alleles_mhc$chr_pos)),] #some RSIDs have NAs, so this way of matching is going to fail.

#Let's check if this is done properly...

length(tg_eaf_alleles$POS_b37)-length(tg_eaf_alleles_NO_mhc$POS_b37) #perfect matching.

###########################################################################################################
#Let's take the columns that we want and make them in the same format as fiadjbmi, which is honestly great#
###########################################################################################################

tg_end <- tg_eaf_alleles_NO_mhc %>%
  select(rsID, CHROM, POS_b37, ALT, REF, POOLED_ALT_AF, EFFECT_SIZE, SE, pvalue_GC, N,chr_pos)

colnames(tg_end) <- c("variant", "chromosome","base_pair_location", "effect_allele", "other_allele", "effect_allele_frequency", "beta", "standard_error", "p_value", "sample_size", "chr_pos")

fwrite(tg_end, "output/1_curated_data/tg_curated.txt")

#############
#WE ARE DONE#
#############