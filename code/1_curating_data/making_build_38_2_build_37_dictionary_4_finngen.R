##############
#INTRODUCTION#
##############

#This is a code makes a dataframe with chr and positions in build 36 and 37 so that we can properly clean all the curated data. 

#HC sex-combined GWAS has all the postiions for the summary statistics, making it easier for us to retrieve the rest.

###########
#libraries#
###########

library(data.table)
library(tidyverse)

###################
#Loading functions#
###################

rsid_parser <- function(id){
  
  rsid <- str_split(id, ":")[[1]][1]
  
  return(rsid)
  
}

#################################################################################################
#STRATEGY: WE ARE GOING TO USE THE BUILD 37 DATA FROM THE META-ANALYSIS WITH UKBB AS FIRST CHECK#
#################################################################################################

project_path <- "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2023/" #change it with your own path.

setwd(project_path)

ref <- fread("../../mho_variants_common_info/raw_data/whradjbmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz")

ckd <- fread("raw_data/finngen_R9_N14_CHRONKIDNEYDIS.gz")
t2d <- fread("raw_data/finngen_R9_T2D.gz")
chd <- fread("raw_data/finngen_R9_I9_CHD.gz")
pcos <- fread("raw_data/finngen_R9_E4_PCOS_CONCORTIUM.gz")
nafld <- fread("raw_data/finngen_R9_NAFLD.gz")
hyptens <- fread("raw_data/finngen_R9_I9_HYPTENS.gz")

#############################################################################
#Let's only do the variants that are not found on CKD to create a dictionary#
#############################################################################

t2d_unmatch <- t2d[which(!(t2d$rsids%in%ckd$rsids)),] #only 60
chd_unmatch <- chd[which(!(chd$rsids%in%ckd$rsids)),] #only 101
pcos_unmatch <- pcos[which(!(pcos$rsids%in%ckd$rsids)),] #only 1!
nafld_unmatch <- nafld[which(!(nafld$rsids%in%ckd$rsids)),] #only 101
hyptens_unmatch <- hyptens[which(!(hyptens$rsids%in%ckd$rsids)),] #only 101

all_unmatch <- rbind(t2d_unmatch, chd_unmatch, pcos_unmatch, nafld_unmatch, hyptens_unmatch)

all_unmatch <- all_unmatch[which(duplicated(all_unmatch$rsids) == FALSE),] #101!! Always the same.

ckd_plus_extra <- rbind(ckd, all_unmatch)

#This time around we are going to do some filtering in terms of allele frequency to make things easier when matching.

ckd_plus_extra_maf <- ckd_plus_extra[which(as.numeric(ckd_plus_extra$af_alt) > 0.01 & (as.numeric(ckd_plus_extra$af_alt) < 0.99)),]

summary(ckd_plus_extra_maf$af_alt)

######################################
#Let's clean the reference one second#
######################################

ref <- ref[which(is.na(ref$INFO) == FALSE)] #these have the chromosome and base_pair_location wrong
ref$variant <- as.character(unlist(sapply(ref$SNP, rsid_parser)))

###################################
#Let's prepare the dictionary here#
###################################

dict <- ckd_plus_extra_maf %>%
  select("rsids", "#chrom", "pos")

colnames(dict) <- c("variant", "chromosome", "base_pair_location_38")

################################
#Let's perform a first matching#
################################

ref_2_dict <- ref[which(ref$variant%in%dict$variant),]
dict_2_ref <- dict[which(dict$variant%in%ref$variant),]
dict_2_ref <- dict_2_ref[which(duplicated(dict_2_ref$variant) == FALSE),]

#Reference has multi-allelic SNPs. 
#We can remove them, we just want the positions.

ref_2_dict <- ref_2_dict[which(duplicated(ref_2_dict$variant) == FALSE),]

dict_2_ref <- dict_2_ref[order(match(dict_2_ref$variant, ref_2_dict$variant)),]

#Let's check that they are correctly aligned:

length(which(dict_2_ref$variant == ref_2_dict$variant)) #perfect match

#Finally add:

dict_2_ref$base_pair_location_37 <- ref_2_dict$POS

#########################################################################
#Let's work with the rest of the variants - we will put them in liftover#
#########################################################################

dict_miss <- dict[which(!(dict$variant%in%ref$variant)),] #55K are missing. That is not so much!

data_4_liftover <- dict_miss

#########################################
#Let's try doing it the other way around#
#########################################

#LOAD "rtracklayer"
library("rtracklayer")

#Attribute the most recent build-appropriate human genome build to a name of your choise (snps in this case)
snps <- SNPlocs.Hsapiens.dbSNP155.GRCh37

#Create Dataframe for conversion

chromosomesVEC = as.character(seq(1,22, by = 1))
all_snps_info = snpsBySeqname(snps, chromosomesVEC) #76.4 freaking GB!!!
all_snps_info_DF = data.frame(all_snps_info)
head(all_snps_info_DF) #worked out!!

snps_df_match <- all_snps_info_DF[which(all_snps_info_DF$RefSNP_id%in%dict_miss$variant),]

dict_miss_match <- dict_miss[which(dict_miss$variant%in%snps_df_match$RefSNP_id),]
dict_miss_mismmatch <- dict_miss[which(!(dict_miss$variant%in%snps_df_match$RefSNP_id)),]

#Most of the mismatches are actually just that: some do not exist in build 37, which is our reference panel, so we won't find a thing there.
#Second, many of them are insertions or deletions, which makes it perfect since we are not including them.
#With that I am confident that we can proceed. 

dict_miss_match_non_dupl <- dict_miss_match[which(duplicated(dict_miss_match$variant) == FALSE),] #perfect match:

dict_miss_match_non_dupl <- dict_miss_match_non_dupl[order(match(dict_miss_match_non_dupl$variant, snps_df_match$RefSNP_id)),]

length(which(dict_miss_match_non_dupl$variant == snps_df_match$RefSNP_id)) #perfect match

dict_miss_match_non_dupl$base_pair_location_37 <- snps_df_match$pos

###################
#Saving dictionary#
###################

dict_end <- rbind(dict_2_ref, dict_miss_match_non_dupl)

fwrite(dict_end, "output/1_curated_data/build_38_2_37_for_finngen_only_ss.txt")
