##############
#INTRODUCTION#
##############

#This code will prepare the proxies input for Lotta et al 53 IR loci for Go-shifter.

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)

###################
#Loading functions#
###################

chr_parser = function(chr_pos){
  
  #STEP 1: strip away by :
  
  info_vect = unlist(str_split(chr_pos, ":"))
  
  chr = info_vect[1]

  #STEP 3: clean chr
  
  chr_clean = unlist(str_split(chr, "chr"))[2]
  
  return(as.numeric(chr_clean))
    
}

pos_parser = function(chr_pos){
  
  #STEP 1: strip away by :
  
  info_vect = unlist(str_split(chr_pos, ":"))
  
  pos_clean = info_vect[2]
  
  return(as.numeric(pos_clean))
  
}

#####################
#Let's read the data#
#####################

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/")

lotta <- readxl::read_xlsx("raw_data/previous_loci/lotta_et_al.xlsx")
lotta_clean = lotta[which(lotta$SNP != "rs6822892"),] #removing tri-alleic SNPs for which we cannot recover proxies

for(index in seq(1, length(lotta_clean$SNP))){
  
  print(index)
  
  #STEP 1: get the rsID
  
  rsid = lotta_clean$SNP[index]
  
  #STEP 2: get the proxies
  
  proxies_tmp = LDlinkR::LDproxy(rsid, pop = "EUR", token="04cad4ca4374") #the rest of default settings are build 37 and r2, just what we need 
  
  proxies_tmp$lead = rsid
  
  if(!(exists("proxies_df"))){
    
    proxies_df = proxies_tmp
    
  } else {
    
    proxies_df = rbind(proxies_df, proxies_tmp)
    
  }
  
}


#########################################
#Let's format and save the lead variants#
#########################################

leads = proxies_df[which(proxies_df$Distance == 0),] #52
leads$chromosome = as.numeric(as.character(unlist(sapply(leads$Coord, chr_parser))))
leads$base_pair_location = as.numeric(as.character(unlist(sapply(leads$Coord, pos_parser))))

leads <- leads %>%
  select(RS_Number, chromosome, base_pair_location)

colnames(leads) <- c("rsID","chr","pos")

dir.create("review/review_output/5_enrichment_analyses/")
dir.create("review/review_output/5_enrichment_analyses/2_binomial_lotta")
dir.create("review/review_output/5_enrichment_analyses/2_binomial_lotta/input")
dir.create("review/review_output/5_enrichment_analyses/2_binomial_lotta/output")
fwrite(leads, "review/review_output/5_enrichment_analyses/2_binomial_lotta/input/lotta_input.txt", sep = "\t")

###################################
#Let's get the input for our SNPs!#
###################################

ir_variants <- fread("output/4_ir_loci_discovery/3_novel_variants/3_novel_variants/282_ir_variants_with_ir_and_t2d_comparisons_w_closest_genes.txt")

ir_282_leads <- ir_variants %>%
  select(variant, chromosome, base_pair_location)

colnames(ir_282_leads) <- c("rsID","chr","pos")

fwrite(ir_282_leads, "review/review_output/5_enrichment_analyses/2_binomial_lotta/input/all_282_ir_input.txt", sep = "\t")

#Now let's split the according to the subsets - let's use the original data from the GoShifter analyses

og_bmi_ns = fread("output/5_enrichment_analyses/3_go_shifter/input_bmi_ns.txt")
og_bmi_neg = fread("output/5_enrichment_analyses/3_go_shifter/input_bmi_neg.txt")
og_bmi_pos = fread("output/5_enrichment_analyses/3_go_shifter/input_bmi_pos.txt")

#Let's split the data and save it:

ns_ir_282_leads = ir_282_leads[which(ir_282_leads$rsID%in%og_bmi_ns$SNP),] #141
neg_ir_282_leads = ir_282_leads[which(ir_282_leads$rsID%in%og_bmi_neg$SNP),] #63
pos_ir_282_leads = ir_282_leads[which(ir_282_leads$rsID%in%og_bmi_pos$SNP),] #141

fwrite(ns_ir_282_leads, "review/review_output/5_enrichment_analyses/2_binomial_lotta/input/ns_141_ir_input.txt", sep = "\t")
fwrite(neg_ir_282_leads, "review/review_output/5_enrichment_analyses/2_binomial_lotta/input/neg_63_ir_input.txt", sep = "\t")
fwrite(pos_ir_282_leads, "review/review_output/5_enrichment_analyses/2_binomial_lotta/input/pos_78_ir_input.txt", sep = "\t")