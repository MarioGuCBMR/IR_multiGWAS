##############
#INTRODUCTION#
##############

#This code tries to handle the data from QTL analyses in the same narrative way as the paper section from QTLs!
#Many more things can be done, but this is the best that we can with 5000 words. 

###################
#Loading libraries#
###################

library(tidyverse)
library(data.table)

##############
#Loading data#
##############

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/")

ir_variants <- fread("output/4_ir_loci_discovery/3_novel_variants/4_colocalization/352_variants_with_all_info_plus_fiadjbmi_and_tg_hdl_coloc.txt")
cs2g_df <- fread("output/8_fine_mapping/all_loci/gene_prioritization/output/cS2G_prioritization_output_raw.csv") #148 variant-gene links
annotation_df <- fread("output/8_fine_mapping/all_loci/gene_prioritization/input/variants_4_cS2G_dictionary.txt")

#Let's get those with cS2G >0.1

cs2g_df <- cs2g_df[which(cs2g_df$Confidence_score >0.1),]

############################################
#STEP 1: filter for those that are in exons#
############################################

cs2g_df <- cs2g_df[which(cs2g_df$Exon == 1),] #50

#Let's filter for those that have no information:

cs2g_df <- cs2g_df[which(cs2g_df$GTEx == 0 & cs2g_df$eQTLGen == 0 & cs2g_df$EpiMap == 0 & cs2g_df$ABC == 0 & cs2g_df$Cicero == 0),]

#####################################################
#STEP 2: let's get more information on each data set#
#####################################################

cs2g_df$lead <- NA
cs2g_df$gene_id <- NA

for(index_proxy in seq(1, length(cs2g_df$Variants))){
  
  #STEP 1: get proxy:
  
  proxy_ <- cs2g_df$Variants[index_proxy]
  gene_ <- cs2g_df$Genes[index_proxy]
  
  #STEP 2: get the lead:
  
  lead_ <- annotation_df$lead_variant[which(annotation_df$chr_pos == proxy_)]
  
  #STEP 3: add it to the dataframe:
  
  cs2g_df$lead[index_proxy] <- lead_
  
  #STEP 4: let's updated the ENSG gene too:
  
  gene_info <-   tryCatch(
    otargen::geneInfo(gene_),
    error = function(e) NA,
    warning = function(w) NA
  )
  
  #STEP 5 add if possible:
  
  if(length(gene_info) == 1){ #this means it returned an NA
    
    next()
    
  } else {
    
    cs2g_df$gene_id[index_proxy] <- gene_info$id
    
  }
  
  
}

#############################################################################
#STEP 3: let's use this opportunity to add more information on these fuckers#
#############################################################################

cs2g_df$reported <- NA
cs2g_df$ir_source <- NA

for(index_proxy in seq(1, length(cs2g_df$Variants))){
  
  #STEP 1: get proxy:
  
  #print(index_proxy)
  
  lead <- cs2g_df$lead[index_proxy]
  
  #STEP 2: get the info:
  
  ir_variants_match <- ir_variants[which(ir_variants$variant == lead),]
  
  #STEP 3: get info:
  
  cs2g_df$reported[index_proxy] <- ir_variants_match$reported_ir
  cs2g_df$ir_source[index_proxy] <- ir_variants_match$ir_source
  
  
}

#######################################################################################################
#STEP 4: let's check the compiled data for enhancer-gene links and assess whether we have info on them#
#######################################################################################################

enhancer_df <- fread("output/7_functional_annotation/2_enhancer_gene_links/2_enhancer_gene_links_and_hic/proxies_with_stare_and_hic_data.txt")

enhancer_df <- enhancer_df[which(enhancer_df$query_snp_rsid%in%cs2g_df$lead),]

#Now let's loop, do we have any evidence that any of the proxies in the signal are in potential enhancer-promoter contact regions?

cs2g_df$enhancer_data <- NA

for(index_proxy  in seq(1, length(cs2g_df$Variants))){
  
  #STEP 1: get lead:
  
  lead_ <- cs2g_df$lead[index_proxy]
  
  #STEP 2: get enhancer data:
  
  enhancer_tmp <- enhancer_df[which(query_snp_rsid == lead_),]
  
  #STEP 3: check whether we have data on this:
  
  check_all <- ifelse(enhancer_tmp$stare_day_0 != "" | enhancer_tmp$stare_day_4 != "" | enhancer_tmp$stare_day_14 != "" | enhancer_tmp$adipocyte_bait_1 != "" | enhancer_tmp$adipocyte_bait_2 != "" | enhancer_tmp$adipocyte_bait_3 != "" | enhancer_tmp$adipocyte_bait_4 != "" | enhancer_tmp$differentiated_bait_1 != "" | enhancer_tmp$differentiated_bait_2 != "", "overlap", "none")
  
  #If we have at least one overlap, the results is overlap:
  
  if("overlap"%in%check_all){
    
    cs2g_df$enhancer_data[index_proxy] <- "overlap"
    
  } else {
    
    cs2g_df$enhancer_data[index_proxy] <- "none"
    
  }
  
}

#Just one variant in AKNA 3'UTR!!

cs2g_df <- cs2g_df[which(cs2g_df$enhancer_data == "none"),]

##################################################################
#STEP 3: let's check whether our variants present QTL signals too#
##################################################################

qtl_df <- readRDS("output/6_qtl_annotation/conditional_and_fine_mapped_asat_vat_liver_and_muscle_eqtls_sqtls.txt") 

qtl_df <- qtl_df[which(qtl_df$query_snp_rsid%in%cs2g_df$lead),]

#############################################
#Let's load the data for all tissues in GTEx#
#############################################

path_2_data <- "C:/Users/zlc436/Downloads/GTEx_Analysis_v10_apaQTL/GTEx_Analysis_v10_apaQTL_updated/"

list_of_files <- list.files(path_2_data)

list_of_files <- list_of_files[which(str_detect(list_of_files, "parquet"))]

for(file_ in list_of_files){
  
  df <- arrow::read_parquet(paste(path_2_data, file_, sep = ""))
  df <- df[which(str_detect(df$gene_id, "ENSG00000106948")),]
  
  df$tissue <- file_
  
  if(!(exists("final_df"))){
    
    final_df <- df
    
  } else {
    
    final_df <- rbind(final_df, df)
    
  }
  
}
