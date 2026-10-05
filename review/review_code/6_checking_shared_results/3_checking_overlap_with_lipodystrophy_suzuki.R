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

#Let's recover them with proxies:

proxies <- fread("output/4_ir_loci_discovery/4_proxies/proxies_build_37_and_38_4_282_variants.txt")

yes_vect <- c("A", "G", "C", "T")

proxies <- proxies[which(proxies$query_snp_rsid%in%ir_variants$variant),] #3188
proxies <- proxies[which(proxies$ref%in%yes_vect & proxies$alt%in%yes_vect),]

length(which(duplicated(proxies$query_snp_rsid) == FALSE)) #282

#Load Suzuki data:

suzuki = readxl::read_xlsx("raw_data/previous_loci/suzuki_et_al.xlsx")
suzuki = suzuki[-1,]
suzuki = as.data.frame(suzuki)

table(suzuki$`Cluster assignment`)

################################################################################################
#Let's first do it with the proxies for each subset. We will know if the proportions make sense#
################################################################################################

bmi_neutral_proxies = proxies[which(proxies$query_snp_rsid%in%ir_variants$variant[which(ir_variants$bmi_subgroup == "BMI-neutral")]),]
bmi_neg_proxies = proxies[which(proxies$query_snp_rsid%in%ir_variants$variant[which(ir_variants$bmi_subgroup == "BMI-decreasing")]),]
bmi_pos_proxies = proxies[which(proxies$query_snp_rsid%in%ir_variants$variant[which(ir_variants$bmi_subgroup == "BMI-increasing")]),]

print(length(unique(bmi_neutral_proxies$query_snp_rsid))) #141
print(length(unique(bmi_neg_proxies$query_snp_rsid))) #63
print(length(unique(bmi_pos_proxies$query_snp_rsid))) #78

suzuki_neutral = suzuki[which(suzuki$`Index SNV`%in%bmi_neutral_proxies$rsID),] #25
suzuki_neg = suzuki[which(suzuki$`Index SNV`%in%bmi_neg_proxies$rsID),] #15
suzuki_pos = suzuki[which(suzuki$`Index SNV`%in%bmi_pos_proxies$rsID),] #15

table(suzuki_neutral$`Cluster assignment`)
table(suzuki_neg$`Cluster assignment`)
table(suzuki_pos$`Cluster assignment`)

####################################################################################
#Let's try getting the proxies for the 45 SNPs in Suzuki with r2>0.01 and re-assess#
####################################################################################

ir_variants$t2d_signal=NA
ir_variants$r2_t2d_signal=NA

for(index in seq(1, length(ir_variants$variant))){
  
  print(index)
  
  #STEP 1: get the rsID
  
  rsid = ir_variants$variant[index]
  
  #STEP 2: get the proxies from the lead that are matching T2D
  
  proxies_tmp = LDlinkR::LDproxy(rsid, pop = "EUR", token="04cad4ca4374") #the rest of default settings are build 37 and r2, just what we need 
  proxies_tmp = proxies_tmp[which(as.numeric(proxies_tmp$R2) >= 0.1),]
  proxies_tmp = proxies_tmp[which(proxies_tmp$RS_Number%in%suzuki$`Index SNV`),]
  
  if(is_empty(proxies_tmp$RS_Number)){
    
    next()
    
  } else {
    
    suzuki_match = suzuki[which(suzuki$`Index SNV`%in%proxies_tmp$RS_Number),]
    
    ir_variants$t2d_signal[index]=paste(suzuki_match$`Cluster assignment`, collapse=";")
    
    signal = paste(proxies_tmp$RS_Number, " (r2=", proxies_tmp$R2, ")", sep = "")
    
    ir_variants$r2_t2d_signal[index]=paste(signal, collapse=";")
    
    
  }

}

fwrite(ir_variants, "review/manuscript/tables/282_ir_with_t2d_data.txt")

##############################################
#Let's try answering Tuomas on these analyses#
##############################################

same_signals = ir_variants[which(str_detect(ir_variants$r2_t2d_signal, "r2=0.8") |
                                 str_detect(ir_variants$r2_t2d_signal, "r2=0.9") | 
                                 str_detect(ir_variants$r2_t2d_signal, "r2=1")),]

table(same_signals$t2d_signal)

same_neutral = same_signals[which(same_signals$bmi_subgroup == "BMI-neutral"),]
table(same_neutral$t2d_signal)

same_decreasing = same_signals[which(same_signals$bmi_subgroup == "BMI-decreasing"),]
table(same_decreasing$t2d_signal)

same_increasing = same_signals[which(same_signals$bmi_subgroup == "BMI-increasing"),]
table(same_increasing$t2d_signal)

