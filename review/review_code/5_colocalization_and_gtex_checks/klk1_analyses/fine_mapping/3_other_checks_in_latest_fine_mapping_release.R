##############
#INTRODUCTION#
##############

#This code reads the fine-mapped data for FIadjBMI and assesses whether our IR loci and or their proxies are there

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)

#######################
#Loading IR data first#
#######################

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/")

ir_variants = fread("manuscript/supplementary_material/supplementary_tables/drafts/supplementary_table_2.csv")

#Let's recover them with proxies:

proxies <- fread("output/4_ir_loci_discovery/4_proxies/proxies_build_37_and_38_4_282_variants.txt")

yes_vect <- c("A", "G", "C", "T")

proxies <- proxies[which(proxies$query_snp_rsid%in%ir_variants$variant),] #6182
proxies <- proxies[which(proxies$ref%in%yes_vect & proxies$alt%in%yes_vect),]

length(which(duplicated(proxies$query_snp_rsid) == FALSE)) #282

##############################
#Loading fine-mapped datasets# We are gonna focus on the KLK1 one though...
##############################

ir_variants$credible_set = NA

fi_st_agno = fread("review/review_raw_data/credible_sets_ST_agno_FI.tsv.gz")
fi_st_anno = fread("review/review_raw_data/credible_sets_ST_anno_FI.tsv.gz")
fi_mt_agno = fread("review/review_raw_data/credible_sets_MT_agno_FI.tsv.gz")
fi_mt_anno = fread("review/review_raw_data/credible_sets_MT_anno_FI.tsv.gz")
fi_ma_agno = fread("review/review_raw_data/credible_sets_MA_agno_FI.tsv.gz")

#Let's add chr_pos to all of them:

fi_st_agno$chr_pos = paste("chr", fi_st_agno$chr, ":", fi_st_agno$pos, sep = "")
fi_st_anno$chr_pos = paste("chr", fi_st_anno$chr, ":", fi_st_anno$pos, sep = "")
fi_mt_agno$chr_pos = paste("chr", fi_mt_agno$chr, ":", fi_mt_agno$pos, sep = "")
fi_mt_anno$chr_pos = paste("chr", fi_mt_anno$chr, ":", fi_mt_anno$pos, sep = "")
fi_ma_agno$chr_pos = paste("chr", fi_ma_agno$chr, ":", fi_ma_agno$pos, sep = "")

for(index in seq(1, length(ir_variants$variant))){
  
  #STEP 1: get data of IR SNP:
  
  ir_snp = ir_variants$variant[index]
  
  #STEP 2: get the damn proxies:
  
  proxy_tmp = proxies[which(proxies$query_snp_rsid == ir_snp),]
  
  #STEP 3: get the chr_pos in build 37 for the proxies:
  
  proxy_tmp$chr_pos = paste("chr", proxy_tmp$chr, ":", proxy_tmp$pos_hg19, sep = "")
  
  #STEP 4: let's do the overlaps real quick:
  
  ##########################
  #FI SINGLE-TRAIT AGNOSTIC#
  ##########################
  
  st_agno_match = fi_st_agno[which(fi_st_agno$chr_pos%in%proxy_tmp$chr_pos),]
  
  if(is_empty(st_agno_match$chr) == FALSE){
    
    id_ = unique(st_agno_match$credset_id)
    
    ir_variants$credible_set[index] = ifelse(is.na(ir_variants$credible_set[index]), id_, paste(ir_variants$credible_set[index], id_, ";"))
    
  }
  
  ############################
  #FI SINGLE-TRAIT ANNOTATION#
  ############################
  
  st_anno_match = fi_st_anno[which(fi_st_anno$chr_pos%in%proxy_tmp$chr_pos),]
  
  if(is_empty(st_anno_match$chr) == FALSE){
    
    id_ = unique(st_anno_match$credset_id)
    
    ir_variants$credible_set[index] = ifelse(is.na(ir_variants$credible_set[index]), id_, paste(ir_variants$credible_set[index], id_, ";"))
    
  }
  
  ##########################
  #FI SINGLE-TRAIT AGNOSTIC#
  ##########################
  
  mt_agno_match = fi_mt_agno[which(fi_mt_agno$chr_pos%in%proxy_tmp$chr_pos),]
  
  if(is_empty(mt_agno_match$chr) == FALSE){
    
    id_ = unique(mt_agno_match$credset_id)
    
    ir_variants$credible_set[index] = ifelse(is.na(ir_variants$credible_set[index]), id_, paste(ir_variants$credible_set[index], id_, ";"))
    
  }
  
  ###########################
  #FI MULTI-TRAIT ANNOTATION#
  ###########################
  
  mt_anno_match = fi_mt_anno[which(fi_mt_anno$chr_pos%in%proxy_tmp$chr_pos),]
  
  if(is_empty(mt_anno_match$chr) == FALSE){
    
    id_ = unique(mt_anno_match$credset_id)
    
    ir_variants$credible_set[index] = ifelse(is.na(ir_variants$credible_set[index]), id_, paste(ir_variants$credible_set[index], id_, ";"))
    
  }
  
  #############################
  #FI MULTI-ANCESTRYT AGNOSTIC#
  #############################
  
  ma_agno_match = fi_ma_agno[which(fi_ma_agno$chr_pos%in%proxy_tmp$chr_pos),]
  
  if(is_empty(ma_agno_match$chr) == FALSE){
    
    id_ = unique(ma_agno_match$credset_id)
    
    ir_variants$credible_set[index] = ifelse(is.na(ir_variants$credible_set[index]), id_, paste(ir_variants$credible_set[index], id_, ";"))
    
  }
  
}

#########################################################################################################################
#A decent, but lacking overlap - this is just finding the glycemic traits, not the ones that also affect dyslipidemia...#
#########################################################################################################################

#-> and ofc, KLK1 is not there. 
#-> we are just gonna have tone down the results

