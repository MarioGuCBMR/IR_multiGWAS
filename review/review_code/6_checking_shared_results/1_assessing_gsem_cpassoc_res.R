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

##############
#Loading data#
##############

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/")

ir_variants = fread("manuscript/supplementary_material/supplementary_tables/drafts/supplementary_table_2.csv")
even_sample_gsem = fread("review/review_output/2_mv_gwas/fiadjbmi_hdl_tg_common_dwls_curated.txt")

#We won't be able to have a fully perfect match!
#Let's get the proxies and get those bad boys talking

even_match = even_sample_gsem[which(even_sample_gsem$chr_pos%in%ir_variants$chr_pos),] #243

#Let's get the mismatched:

ir_missing = ir_variants[which(!(ir_variants$chr_pos%in%even_match$chr_pos)),] #39

#Let's recover them with proxies:

proxies <- fread("output/4_ir_loci_discovery/4_proxies/proxies_build_37_and_38_4_282_variants.txt")

yes_vect <- c("A", "G", "C", "T")

proxies <- proxies[which(proxies$query_snp_rsid%in%ir_variants$variant),] #3188
proxies <- proxies[which(proxies$ref%in%yes_vect & proxies$alt%in%yes_vect),]

length(which(duplicated(proxies$query_snp_rsid) == FALSE)) #282

#Great, let's see how many proxies we got

proxies_missing = proxies[which(proxies$query_snp_rsid%in%ir_missing$variant),]

length(which(duplicated(proxies_missing$query_snp_rsid) == FALSE)) #39

even_recovered = even_sample_gsem[which(even_sample_gsem$variant%in%proxies_missing$rsID),] #103 proxies!!

#Let's get one for each, not matter which ones

proxies_recovered = proxies_missing[which(proxies_missing$rsID%in%even_recovered$variant),]

even_recovered = even_recovered[order(match(even_recovered$variant, proxies_recovered$rsID)),]

print(length(which(even_recovered$variant == proxies_recovered$rsID))) #perfect

even_recovered$lead = proxies_recovered$query_snp_rsid

#Let's get the best hit for all of them:

even_recovered = even_recovered[order(as.numeric(even_recovered$p_value)),]
even_recovered = even_recovered[which(duplicated(even_recovered$lead) == FALSE),] #21/39! That is great tbh

#Let's match the data:

even_match$lead = even_match$variant 

even_df = rbind(even_match, even_recovered) #264/282

#Let's get the real data boys

print(length(which(as.numeric(even_df$p_value) < 0.05))) #233/264 - 88%

#This is amazing ;)


##################################################################################
#Let's use this space to make some new columns that we will add to Suppl Table 2!#
##################################################################################

#Let's load the Supplementary Info to have the variants in the right order:

all_variants=readxl::read_xlsx("review/manuscript/03062026/Supplementary Tables_20March2026.xlsx", sheet = 3)
all_variants = as.data.frame(all_variants)
all_variants = all_variants[2:nrow(all_variants),] #remove extra first row that explains stuff
colnames(all_variants) = all_variants[1,]  #change column names
all_variants = all_variants[2:nrow(all_variants),] #remove columns names in the first line now that we have used them

#We do not need to many columns AND we have functions ready for a certain subset of variables.

all_variants = all_variants[,1:12]
colnames(all_variants) = c("variant", "chromosome", "base_pair_location", "minimum_allele_frequency", "effect_allele", "other_allele", "novel_or_reported", "bmi_subgroup", "VEP", "beta", "standard_error", "pvalue")

##############################
#Let's loop and find the data#
##############################

all_variants$pvalue.gsem_even_n =NA

for(index in seq(1, length(all_variants$variant))){
  
  #STEP 1: get the variant in that position:
  
  rsid = all_variants$variant[index]
  
  #STEP 2: match with the recovered dataframe:
  
  even_tmp = even_df[which(even_df$lead == rsid),]
  
  if(is_empty(even_tmp$variant)){
    
    next()
    
  }
  
  #STEP 3: get the p-value; above we already parsed for the best association
  
  all_variants$pvalue.gsem_even_n[index] = as.numeric(even_tmp$p_value)
  
}

#######################################
#Let's do the same with CPASSOC info!!#
#######################################

cpassoc_df = fread("review/review_output/2_mv_gwas/cpassoc/282_ir_fiadjbmi_hdl_tg_dwls_snps_w_cpassoc_even_sample_proxies.txt")

#Let's first filter for the best proxy for each variant:

cpassoc_df = cpassoc_df[order(as.numeric(cpassoc_df$pval.shet)),]
cpassoc_df = cpassoc_df[which(duplicated(cpassoc_df$query_snp_rsid) == FALSE),] #270 - we recovered a lot!!

all_variants$pvalue.cpassoc_even_n =NA

for(index in seq(1, length(all_variants$variant))){
  
  #STEP 1: get the variant in that position:
  
  rsid = all_variants$variant[index]
  
  #STEP 2: match with the recovered dataframe:
  
  even_tmp = cpassoc_df[which(cpassoc_df$query_snp_rsid == rsid),]
  
  if(is_empty(even_tmp$rsID)){
    
    next()
    
  }
  
  #STEP 3: get the p-value; above we already parsed for the best association
  
  all_variants$pvalue.cpassoc_even_n[index] = as.numeric(even_tmp$pval.shet)
  
}

##############################
#We are ready to save this!!!#
##############################

fwrite(all_variants, "review/manuscript/tables/282_ir_with_gsem_cpassoc_even_sample_pvals.txt")

######################
#Let's check the info#
######################

#How many variants are present in both analyses:

length(which(is.na(all_variants$pvalue.gsem_even_n) == FALSE &
               is.na(all_variants$pvalue.cpassoc_even_n) == FALSE)) #264 in total


length(which(as.numeric(all_variants$pvalue.gsem_even_n) < 0.05)) #233/264
length(which(as.numeric(all_variants$pvalue.cpassoc_even_n) < 0.05)) #237/264

length(which(as.numeric(all_variants$pvalue.gsem_even_n) < 0.05  |
               as.numeric(all_variants$pvalue.cpassoc_even_n) < 0.05)) #253


length(which(as.numeric(all_variants$pvalue.gsem_even_n) < 0.05  &
               as.numeric(all_variants$pvalue.cpassoc_even_n) < 0.05)) #217


#Let's see how many are novel: 

novel = all_variants[which(all_variants$novel_or_reported == "Novel" | is.na(all_variants$novel_or_reported)),] #70- amazing

#How many of them are testable:

novel = novel[which(is.na(novel$pvalue.cpassoc_even_n) == FALSE & is.na(novel$pvalue.gsem_even_n) == FALSE),] #66
              
novel = novel[which(as.numeric(novel$pvalue.cpassoc_even_n) < 0.05),] #54
novel = novel[which(as.numeric(novel$pvalue.gsem_even_n) < 0.05),] #52

################################################
#Let's do some additional counting just in case#
################################################

#Get novel first:

novel_check = ir_variants[which(ir_variants$ir_source == ""),] #Let's see if the operations above made sense
length(which(novel_check$variant%in%novel$variant)) #52/52 it is correct

#How many are in the even_df and pass that threshold?

novel_even = even_df[which(even_df$lead%in%novel_check$variant),] #66 - now it is correct

print(length(which(as.numeric(novel_even$p_value) < 0.05))) #56

#How many in the original are genome-wide associated?

print(length(which(as.numeric(ir_variants$p_value) < 5e-08))) #169/282 - 60%
print(length(which(as.numeric(novel_check$p_value) < 5e-08))) #20/70 - 29%

#How many do they pass the Lotta threshold?

print(length(which(as.numeric(ir_variants$p_value.fiadjbmi) < 5e-03 &
                     as.numeric(ir_variants$p_value.hdl) < 5e-03 &
                     as.numeric(ir_variants$p_value.tg) < 5e-03))) #88/282 - 31%

print(length(which(as.numeric(novel_check$p_value.fiadjbmi) < 5e-03 &
                     as.numeric(novel_check$p_value.hdl) < 5e-03 &
                     as.numeric(novel_check$p_value.tg) < 5e-03))) #22/70 - 31%

#And how many the P<0.05 for all of them?

print(length(which(as.numeric(ir_variants$p_value.fiadjbmi) < 5e-02 &
                     as.numeric(ir_variants$p_value.hdl) < 5e-02 &
                     as.numeric(ir_variants$p_value.tg) < 5e-02))) #143/282 - 50%

print(length(which(as.numeric(novel_check$p_value.fiadjbmi) < 5e-02 &
                     as.numeric(novel_check$p_value.hdl) < 5e-02 &
                     as.numeric(novel_check$p_value.tg) < 5e-02))) #32/70 - 45%

#And how many novel are genome-wide significant for at least one trait?

print(length(which(as.numeric(novel_check$p_value.fiadjbmi) < 5e-08 |
                     as.numeric(novel_check$p_value.hdl) < 5e-08 |
                     as.numeric(novel_check$p_value.tg) < 5e-08))) #54/70


