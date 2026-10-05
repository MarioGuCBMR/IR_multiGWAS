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

tg <- fread("review/review_raw_data/jointGwasMc_TG.txt.gz")

############################################################
#1. Let's do toupper to get the alleles in the right format#
############################################################

tg$A1 = toupper(tg$A1)
tg$A2 = toupper(tg$A2)

table(tg$A1)
table(tg$A2) #no indels

######################################################
#2. Remove those that have a EAF > 0.99 or EAF < 0.01#
######################################################

#Let's check first for any NAs:

summary(tg$Freq.A1.1000G.EUR) #we have so many NAs..., I am going to load the original that has no NAs and try to recover the NAs from there:

tg_ref=fread("output/1_curated_gwas/tg_curated.txt")

#Let's perform some quick matches:

tg_missing = tg[which(is.na(tg$Freq.A1.1000G.EUR)),]
tg_missing_recoverable = tg_missing[which(tg_missing$SNP_hg19%in%tg_ref$chr_pos),] #we only miss 5000, which might have been MHC or other, so happy we got this!
tg_missing_recoverable = tg_missing_recoverable[order(as.numeric(tg_missing_recoverable$`P-value`)),]
tg_missing_recoverable = tg_missing_recoverable[which(duplicated(tg_missing_recoverable$SNP_hg19) == FALSE),] #removing duplicates with higher p-values so that it is easier to do it downstream

#Let's get the reference

tg_ref_match = tg_ref[which(tg_ref$chr_pos%in%tg_missing_recoverable$SNP_hg19),] #we seem to have some duplicates, but we know how to solve that:

#Matching with unique IDs will remove triallelic SNPs or duplicates that we might find:

tg_missing_recoverable$id = paste(tg_missing_recoverable$SNP_hg19, tg_missing_recoverable$A1, tg_missing_recoverable$A2, sep = ":")
tg_ref_match$id_1 = paste(tg_ref_match$chr_pos, tg_ref_match$effect_allele, tg_ref_match$other_allele, sep = ":")
tg_ref_match$id_2 = paste(tg_ref_match$chr_pos, tg_ref_match$other_allele, tg_ref_match$effect_allele, sep = ":")

#Let's do the matching and make sure that the data got aligned:

tg_ref_match = tg_ref_match[which(tg_ref_match$id_1%in%tg_missing_recoverable$id | tg_ref_match$id_2%in%tg_missing_recoverable$id),] #yes!
tg_missing_recoverable = tg_missing_recoverable[which(tg_missing_recoverable$id%in%tg_ref_match$id_1 | tg_missing_recoverable$id%in%tg_ref_match$id_2),] #yes!!! We finally got the match. All alleles match

tg_ref_match = tg_ref_match[order(match(tg_ref_match$chr_pos, tg_missing_recoverable$SNP_hg19)),] #chr:pos, should be fine now!
print(length(which(tg_ref_match$chr_pos == tg_missing_recoverable$SNP_hg19))) #perfect

#Let's do a switch of alleles, which should go in line with the frequencies:

tg_ref_match$a1_ordered = ifelse(tg_ref_match$effect_allele != tg_missing_recoverable$A1, tg_ref_match$other_allele, tg_ref_match$effect_allele)
tg_ref_match$a2_ordered = ifelse(tg_ref_match$effect_allele != tg_missing_recoverable$A1, tg_ref_match$effect_allele, tg_ref_match$other_allele)
tg_ref_match$aligned_eaf = ifelse(tg_ref_match$effect_allele != tg_missing_recoverable$A1, 1-as.numeric(tg_ref_match$effect_allele_frequency), as.numeric(tg_ref_match$effect_allele_frequency))

#Now the alleles should be matched:

print(length(which(tg_ref_match$a1_ordered == tg_missing_recoverable$A1))) #perfect!
print(length(which(tg_ref_match$a2_ordered == tg_missing_recoverable$A2))) #perfect!

#We can give the new effect allele frequency to tg_missining_recoverable:

tg_missing_recoverable$Freq.A1.1000G.EUR = as.numeric(tg_ref_match$aligned_eaf)
tg_missing_recoverable$rsid = tg_ref_match$variant #I am changing some of the rsIDs, we have some issues with old ones having been updated and they are now indels

#Let's remove some columns to make thing OK

tg_missing_recoverable = tg_missing_recoverable %>%
  dplyr::select(-c("id"))

#Let's merge the data again and continue with our pipeline:

tg_ok = tg[which(is.na(tg$Freq.A1.1000G.EUR) == FALSE),]

tg = rbind(tg_ok, tg_missing_recoverable)

summary(as.numeric(tg$Freq.A1.1000G.EUR)) #This looks great!!

#We are gonna remove those with MAF < 0.01.

tg_eaf_OK <- tg[which(as.numeric(tg$Freq.A1.1000G.EUR) > 0.01),] #We go from 46M to...9.9M
tg_eaf_OK <- tg_eaf_OK[which(as.numeric(tg_eaf_OK$Freq.A1.1000G.EUR) < 0.99),] #This remained the same.

summary(tg_eaf_OK$Freq.A1.1000G.EUR) #perfect

#Let's get the opportunity to change some columns and get things right:

tg_eaf_OK = tg_eaf_OK %>%
  dplyr::select(-c("SNP_hg18"))

colnames(tg_eaf_OK) = c("chr_pos", "variant", "effect_allele", "other_allele", "beta", "standard_error", "sample_size", "p_value", "effect_allele_frequency")
tg_eaf_OK$chromosome=as.numeric(as.character(unlist(sapply(tg_eaf_OK$chr_pos, chr_parser))))
tg_eaf_OK$base_pair_location=as.numeric(as.character(unlist(sapply(tg_eaf_OK$chr_pos, pos_parser))))

#############################
#Let's remove the MHC region#
#############################

summary(tg_eaf_OK$base_pair_location) #all good, no issues here.

tg_eaf_alleles_mhc <- tg_eaf_OK[which(as.numeric(tg_eaf_OK$chromosome) == 6),]
tg_eaf_alleles_mhc <- tg_eaf_alleles_mhc[which(as.numeric(tg_eaf_alleles_mhc$base_pair_location) >= 26000000),]
tg_eaf_alleles_mhc <- tg_eaf_alleles_mhc[which(as.numeric(tg_eaf_alleles_mhc$base_pair_location) <= 34000000),]

summary(as.numeric(tg_eaf_alleles_mhc$chromosome)) #perfect.
summary(as.numeric(tg_eaf_alleles_mhc$base_pair_location)) #perfect.

tg_eaf_alleles_NO_mhc <- tg_eaf_OK[which(!(tg_eaf_OK$chr_pos%in%tg_eaf_alleles_mhc$chr_pos)),] 

#Let's check if this is done properly...

length(tg_eaf_OK$base_pair_location)-length(tg_eaf_alleles_NO_mhc$base_pair_location) #perfect matching.

###########################################
#Let's save this data! We are finally done#
###########################################

dir.create("review/review_output/1_curated_data")

fwrite(tg_eaf_alleles_NO_mhc, "review/review_output/1_curated_data/tg_2013_curated.txt")

#############
#WE ARE DONE#
#############