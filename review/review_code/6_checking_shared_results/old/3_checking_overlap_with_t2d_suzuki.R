##############
#INTRODUCTION#
##############

#This code will read the 282 IR variants, get the even sample size G-SEM GWAS and do comparisons!!

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)
library(TwoSampleMR)
source("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/mrhorse-main/mr_horse.R")

##############
#Loading data#
##############

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/")

ir_variants = fread("manuscript/supplementary_material/supplementary_tables/drafts/supplementary_table_2.csv")

#Let's get the lipodystrophy like hits:

table(ir_variants$bmi_subgroup)

bmi_neg = ir_variants[which(ir_variants$bmi_subgroup == "BMI-decreasing"),]

#Let's recover them with proxies:

proxies <- fread("output/4_ir_loci_discovery/4_proxies/proxies_build_37_and_38_4_282_variants.txt")

yes_vect <- c("A", "G", "C", "T")

proxies <- proxies[which(proxies$query_snp_rsid%in%bmi_neg$variant),] #3188
proxies <- proxies[which(proxies$ref%in%yes_vect & proxies$alt%in%yes_vect),]

length(which(duplicated(proxies$query_snp_rsid) == FALSE)) #63

##########################################################
#Let's load the lipodystrophy-like hits from Suzuki et al#
##########################################################

suzuki = readxl::read_xlsx("raw_data/previous_loci/suzuki_et_al.xlsx")
table(suzuki$`Cluster assignment`) #45
suzuki_lipo = suzuki[which(suzuki$`Cluster assignment` == "Lipodystrophy"),]

#Great, let's perform the matching! 

proxies_match_lipo = proxies[which(proxies$rsID%in%suzuki_lipo$`Index SNV`),]

gene_info = ir_variants$annotation[which(ir_variants$variant%in%proxies_match_lipo$query_snp_rsid)]

length(unique(proxies_match_lipo$query_snp_rsid)) #11/45 as r2>0.8

#We already have the data with r2 <0.01

bmi_neg_loci = bmi_neg[which(bmi_neg$t2d_cluster == "Lipodystrophy"),] #19

#We have 11 as same signal and 19 defined as same locus. 

####################################################################################
#Let's try getting the proxies for the 45 SNPs in Suzuki with r2>0.01 and re-assess#
####################################################################################

for(index in seq(1, length(suzuki_lipo$`Index SNV`))){
  
  print(index)
  
  #STEP 1: get the rsID
  
  rsid = suzuki_lipo$`Index SNV`[index]
  
  #STEP 2: get the proxies
  
  proxies_tmp = LDlinkR::LDproxy(rsid, pop = "EUR", token="04cad4ca4374") #the rest of default settings are build 37 and r2, just what we need 
  
  proxies_tmp$lead = rsid
  
  if(!(exists("proxies_df"))){
    
    proxies_df = proxies_tmp
    
  } else {
    
    proxies_df = rbind(proxies_df, proxies_tmp)
    
  }

}

#Great, we have all the proxies here!
#With lead matching should be enough:

proxies_suzuki_match = proxies_df[which(proxies_df$RS_Number%in%bmi_neg$variant),] #25!! Well, this is the list that we will end up using!

########################################################################
#Let's load in T2D associations by Suzuki et al and explore the results#
########################################################################

suzuki=fread("output/1_curated_gwas/t2d_suzuki_curated.txt")

suzuki_match = suzuki[which(suzuki$chr_pos%in%bmi_neg$chr_pos),]

suzuki_match = suzuki_match[order(match(suzuki_match$chr_pos, bmi_neg$chr_pos)),]

print(length(which(suzuki_match$chr_pos == bmi_neg$chr_pos))) #perfect match

suzuki_match$beta_aligned = ifelse(suzuki_match$effect_allele != bmi_neg$effect_allele, as.numeric(suzuki_match$beta)*(-1), as.numeric(suzuki_match$beta))

summary(suzuki_match$beta_aligned)
summary(as.numeric(suzuki_match$p_value))

#Let's take a look at those that are significant:

length(which(as.numeric(suzuki_match$p_value) < 0.05)) #51/63 - quite good

#Let's see how many have T2D-increasing associations:

suzuki_nominal = suzuki_match[which(as.numeric(suzuki_match$p_value) <0.05),]

summary(suzuki_nominal$beta_aligned)
length(which(as.numeric(suzuki_nominal$beta_aligned) > 0)) #51/63 - quite good

