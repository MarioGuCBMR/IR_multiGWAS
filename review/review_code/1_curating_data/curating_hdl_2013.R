##############
#INTRODUCTION#
##############

#This code is to curate lipid data from Willer 2013.

###########
#Libraries#
###########

library(data.table)
library(tidyverse)

###################
#Loading functions#
###################

chr_parser <- function(chr_pos){
  
  #This function retrieves chromosomes in chr1 format.
  
  tmp <- strsplit(as.character(chr_pos), ":")[[1]][1]
  
  chr_ = strsplit(tmp, "chr")[[1]][2]
  
  return(as.numeric(chr_))
  
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

hdl <- fread("review/review_raw_data/jointGwasMc_HDL.txt.gz")

############################################################
#1. Let's do toupper to get the alleles in the right format#
############################################################

hdl$A1 = toupper(hdl$A1)
hdl$A2 = toupper(hdl$A2)

table(hdl$A1)

######################################################
#2. Remove those that have a EAF > 0.99 or EAF < 0.01#
######################################################

#Let's check first for any NAs:

summary(hdl$Freq.A1.1000G.EUR) #we have so many NAs..., I am going to load the original that has no NAs and try to recover the NAs from there:

hdl_ref=fread("output/1_curated_gwas/hdl_curated.txt")

#Let's perform some quick matches:

hdl_missing = hdl[which(is.na(hdl$Freq.A1.1000G.EUR)),]
hdl_missing_recoverable = hdl_missing[which(hdl_missing$SNP_hg19%in%hdl_ref$chr_pos),] #we only miss 5000, which might have been MHC or other, so happy we got this!
hdl_missing_recoverable = hdl_missing_recoverable[order(as.numeric(hdl_missing_recoverable$`P-value`)),]
hdl_missing_recoverable = hdl_missing_recoverable[which(duplicated(hdl_missing_recoverable$SNP_hg19) == FALSE),] #removing duplicates with higher p-values so that it is easier to do it downstream

#Let's get the reference

hdl_ref_match = hdl_ref[which(hdl_ref$chr_pos%in%hdl_missing_recoverable$SNP_hg19),] #we seem to have some duplicates, but we know how to solve that:

#Matching with unique IDs will remove triallelic SNPs or duplicates that we might find:

hdl_missing_recoverable$id = paste(hdl_missing_recoverable$SNP_hg19, hdl_missing_recoverable$A1, hdl_missing_recoverable$A2, sep = ":")
hdl_ref_match$id_1 = paste(hdl_ref_match$chr_pos, hdl_ref_match$effect_allele, hdl_ref_match$other_allele, sep = ":")
hdl_ref_match$id_2 = paste(hdl_ref_match$chr_pos, hdl_ref_match$other_allele, hdl_ref_match$effect_allele, sep = ":")

#Let's do the matching and make sure that the data got aligned:

hdl_ref_match = hdl_ref_match[which(hdl_ref_match$id_1%in%hdl_missing_recoverable$id | hdl_ref_match$id_2%in%hdl_missing_recoverable$id),] #yes!
hdl_missing_recoverable = hdl_missing_recoverable[which(hdl_missing_recoverable$id%in%hdl_ref_match$id_1 | hdl_missing_recoverable$id%in%hdl_ref_match$id_2),] #yes!!! We finally got the match. All alleles match

hdl_ref_match = hdl_ref_match[order(match(hdl_ref_match$chr_pos, hdl_missing_recoverable$SNP_hg19)),] #chr:pos, should be fine now!
print(length(which(hdl_ref_match$chr_pos == hdl_missing_recoverable$SNP_hg19))) #perfect

#Let's do a switch of alleles, which should go in line with the frequencies:

hdl_ref_match$a1_ordered = ifelse(hdl_ref_match$effect_allele != hdl_missing_recoverable$A1, hdl_ref_match$other_allele, hdl_ref_match$effect_allele)
hdl_ref_match$a2_ordered = ifelse(hdl_ref_match$effect_allele != hdl_missing_recoverable$A1, hdl_ref_match$effect_allele, hdl_ref_match$other_allele)
hdl_ref_match$aligned_eaf = ifelse(hdl_ref_match$effect_allele != hdl_missing_recoverable$A1, 1-as.numeric(hdl_ref_match$effect_allele_frequency), as.numeric(hdl_ref_match$effect_allele_frequency))

#Now the alleles should be matched:

print(length(which(hdl_ref_match$a1_ordered == hdl_missing_recoverable$A1))) #perfect!
print(length(which(hdl_ref_match$a2_ordered == hdl_missing_recoverable$A2))) #perfect!

#We can give the new effect allele frequency to hdl_missining_recoverable:

hdl_missing_recoverable$Freq.A1.1000G.EUR = as.numeric(hdl_ref_match$aligned_eaf)
hdl_missing_recoverable$rsid = hdl_ref_match$variant #I am changing some of the rsIDs, we have some issues with old ones having been updated and they are now indels

#Let's remove some columns to make thing OK

hdl_missing_recoverable = hdl_missing_recoverable %>%
  dplyr::select(-c("id"))

#Let's merge the data again and continue with our pipeline:

hdl_ok = hdl[which(is.na(hdl$Freq.A1.1000G.EUR) == FALSE),]

hdl = rbind(hdl_ok, hdl_missing_recoverable)

summary(as.numeric(hdl$Freq.A1.1000G.EUR)) #This looks great!!

#We are gonna remove those with MAF < 0.01.

hdl_eaf_OK <- hdl[which(as.numeric(hdl$Freq.A1.1000G.EUR) > 0.01),] #We go from 46M to...9.9M
hdl_eaf_OK <- hdl_eaf_OK[which(as.numeric(hdl_eaf_OK$Freq.A1.1000G.EUR) < 0.99),] #This remained the same.

summary(hdl_eaf_OK$Freq.A1.1000G.EUR) #perfect

#Let's get the opportunity to change some columns and get things right:

hdl_eaf_OK = hdl_eaf_OK %>%
  dplyr::select(-c("SNP_hg18"))

colnames(hdl_eaf_OK) = c("chr_pos", "variant", "effect_allele", "other_allele", "beta", "standard_error", "sample_size", "p_value", "effect_allele_frequency")
hdl_eaf_OK$chromosome=as.numeric(as.character(unlist(sapply(hdl_eaf_OK$chr_pos, chr_parser))))
hdl_eaf_OK$base_pair_location=as.numeric(as.character(unlist(sapply(hdl_eaf_OK$chr_pos, pos_parser))))

#############################
#Let's remove the MHC region#
#############################

summary(hdl_eaf_OK$base_pair_location) #all good, no issues here.

hdl_eaf_alleles_mhc <- hdl_eaf_OK[which(as.numeric(hdl_eaf_OK$chromosome) == 6),]
hdl_eaf_alleles_mhc <- hdl_eaf_alleles_mhc[which(as.numeric(hdl_eaf_alleles_mhc$base_pair_location) >= 26000000),]
hdl_eaf_alleles_mhc <- hdl_eaf_alleles_mhc[which(as.numeric(hdl_eaf_alleles_mhc$base_pair_location) <= 34000000),]

summary(as.numeric(hdl_eaf_alleles_mhc$chromosome)) #perfect.
summary(as.numeric(hdl_eaf_alleles_mhc$base_pair_location)) #perfect.

hdl_eaf_alleles_NO_mhc <- hdl_eaf_OK[which(!(hdl_eaf_OK$chr_pos%in%hdl_eaf_alleles_mhc$chr_pos)),] 

#Let's check if this is done properly...

length(hdl_eaf_OK$base_pair_location)-length(hdl_eaf_alleles_NO_mhc$base_pair_location) #perfect matching.

###########################################
#Let's save this data! We are finally done#
###########################################

dir.create("review/review_output/1_curated_data")

fwrite(hdl_eaf_alleles_NO_mhc, "review/review_output/1_curated_data/hdl_2013_curated.txt")

#############
#WE ARE DONE#
#############