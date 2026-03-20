##############
#INTRODUCTION#
##############

#This code is to curate lipid data from Graham 2022.

###########
#Libraries#
###########

library(data.table)
library(tidyverse)

#memory.limit(size=80000000)

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

project_path <- "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/" #change it with your own path.

setwd(project_path)

hdl <- fread("raw_data/summary_statistics/HDL_INV_AFR_HRC_1KGP3_others_ALL.meta.singlevar.results.gz")

#####################################################
#Let's check the dataframes and what do we have here#
#####################################################

head(hdl)

#Here the insertions and deletions were change to D and I.
#Let's curate and do some checks in the end.

###################
#HDL data analysis#
###################

################################################
#1. Now we are going to remove those SNPs that #
#Have low sample (from June Version is 10.000) #
################################################

summary(hdl$N) #no NAs. Min = 1288, Max = . Works for me.

######################################################
#2. Remove those that have a EAF > 0.99 or EAF < 0.01#
######################################################

#Let's check first for any NAs:

summary(hdl$POOLED_ALT_AF) #no NAs. Interestingly, Max is 0.5, but it is not MAF, but Alternative allele frequency. This means they selected the allele with minimum allele frequency as the effect allele.

#We are gonna remove those with MAF < 0.01.

hdl_eaf_OK <- hdl[which(as.numeric(hdl$POOLED_ALT_AF) > 0.01),] #We go from 46M to...9.9M
hdl_eaf_OK <- hdl_eaf_OK[which(as.numeric(hdl_eaf_OK$POOLED_ALT_AF) < 0.99),] #This remained the same.

summary(hdl_eaf_OK$POOLED_ALT_AF) #perfect

#From  here on we are going to have the chr:pos:

hdl_eaf_OK$chr_pos <- paste("chr", hdl_eaf_OK$CHROM, ":", hdl_eaf_OK$POS_b37, sep="")

#######################################
#Let's remove inserttion and deletions#
#######################################

yes_vect <- c("A", "G", "T", "C")

hdl_eaf_alleles <- hdl_eaf_OK[which(hdl_eaf_OK$REF%in%yes_vect & hdl_eaf_OK$ALT%in%yes_vect),]

#############################
#Let's remove the MHC region#
#############################

summary(hdl_eaf_alleles$POS_b37) #all good, no issues here.

hdl_eaf_alleles_mhc <- hdl_eaf_alleles[which(as.numeric(hdl_eaf_alleles$CHROM) == 6),]
hdl_eaf_alleles_mhc <- hdl_eaf_alleles_mhc[which(as.numeric(hdl_eaf_alleles_mhc$POS_b37) >= 26000000),]
hdl_eaf_alleles_mhc <- hdl_eaf_alleles_mhc[which(as.numeric(hdl_eaf_alleles_mhc$POS_b37) <= 34000000),]

head((hdl_eaf_alleles_mhc$CHROM)) #perfect.
summary(as.numeric(hdl_eaf_alleles_mhc$POS_b37)) #perfect.

hdl_eaf_alleles_NO_mhc <- hdl_eaf_alleles[which(!(hdl_eaf_alleles$chr_pos%in%hdl_eaf_alleles_mhc$chr_pos)),] 

#Let's check if this is done properly...

length(hdl_eaf_alleles$POS_b37)-length(hdl_eaf_alleles_NO_mhc$POS_b37) #perfect matching.

###########################################################################################################
#Let's take the columns that we want and make them in the same format as fiadjbmi, which is honestly great#
###########################################################################################################

hdl_end <- hdl_eaf_alleles_NO_mhc %>%
  select(rsID, CHROM, POS_b37, ALT, REF, POOLED_ALT_AF, EFFECT_SIZE, SE, pvalue_GC, N,chr_pos)

colnames(hdl_end) <- c("variant", "chromosome","base_pair_location", "effect_allele", "other_allele", "effect_allele_frequency", "beta", "standard_error", "p_value", "sample_size", "chr_pos")

fwrite(hdl_end, "output/1_curated_gwas/hdl_curated_afr.txt")

#############
#WE ARE DONE#
#############