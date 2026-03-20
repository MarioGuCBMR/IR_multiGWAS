##############
#INTRODUCTION#
##############

#This code adds information about the variants that are in high LD with QTL in subcutaneous adipose tissue.

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)
library(otargen)

###################
#Loading functions#
###################

chr_parser <- function(snp){
  
  chr <- unlist(str_split(snp, "_")[[1]][1])
  
  chr <- unlist(str_split(chr, "chr")[[1]][2])
  
  return(chr)
  
}

pos_parser <- function(snp){
  
  pos <- unlist(str_split(snp, "_")[[1]][2])
  
  return(pos)
  
}

ref_parser <- function(snp){
  
  data <- str_split(snp, "_")[[1]]
  
  allele <- data[3]
  
  return(allele)
  
}

alt_parser <- function(snp){
  
  data <- str_split(snp, "_")[[1]]
  
  allele <- data[4]
  
  return(allele)
  
}

gene_parser <- function(gene){
  
  data <- str_split(gene, "[.]")[[1]]
  
  cleaned_ <- data[1]
  
  return(cleaned_)
  
}

pheno_parser <- function(pheno){
  
  data <- str_split(pheno, ":ENSG")[[1]]
  
  cleaned_ <- as.character(data[1])
  
  return(cleaned_)
  
}


get_afc_direction <- function(splicing_vector) {
  sapply(splicing_vector, function(splice) {
    if (is.na(splice)) return(NA)
    pheno <- substr(splice, nchar(splice), nchar(splice))
    if (pheno == "-") return(-1) else return(1)
  })
}

##############
#Loading data#
##############

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025//output")

ir_new <- fread("4_ir_loci_discovery/3_novel_variants/3_novel_variants/282_ir_variants_with_ir_and_t2d_comparisons_w_closest_genes.txt")

proxies <- readRDS("6_qtl_annotation/conditional_and_fine_mapped_asat_and_vat_eqtls.txt")

#Finally, let's load the QTL data:

qtl_df <- as.data.frame(arrow::read_parquet("../../gtex_v10/GTEx_v10_SuSiE_sQTL/Adipose_Visceral_Omentum.v10.sQTLs.SuSiE_summary.parquet"))

#Now let's add a chr_pos column to get info

proxies$chr_pos <- paste("chr", proxies$chr, ":", proxies$pos_hg38, sep ="")

#Now the same for the QTLs

qtl_df$chr <- unlist(sapply(qtl_df$variant_id, chr_parser))
qtl_df$pos <- unlist(sapply(qtl_df$variant_id, pos_parser))
qtl_df$ref <- unlist(sapply(qtl_df$variant_id, ref_parser))
qtl_df$alt <- unlist(sapply(qtl_df$variant_id, alt_parser))
qtl_df$gene_id_clean <- unlist(sapply(qtl_df$gene_id, gene_parser))
qtl_df$splicing <- unlist(sapply(qtl_df$phenotype_id, pheno_parser))
qtl_df$afc <- unlist(sapply(qtl_df$splicing, get_afc_direction))

qtl_df$chr_pos <- paste("chr", qtl_df$chr, ":", qtl_df$pos, sep = "")

#################################################################
#Now we can properly fill the data with the data that is missing#
#################################################################

proxies$gene_id.vat_sqtl <- NA
proxies$gene_name.vat_sqtl <- NA
proxies$splicing_region.vat <- NA

proxies$splicing_effect.vat <- NA
proxies$ir_allele.vat_sqtl <- NA
proxies$is_allele.vat_sqtl <- NA

proxies$pip.vat_sqtl <- NA

for(index in seq(1, length(proxies$rsID))){
  
  #STEP 1: get the info we want:
  
  variant= proxies$rsID[index]
  query = proxies$query_snp_rsid[index]
  
  chr_38 <- proxies$chr[index]
  bp_38 = proxies$pos_hg38[index]
  
  #Then we query the data in qtl_df
  
  qtl_tmp <- qtl_df[which(qtl_df$chr == chr_38 & qtl_df$pos == bp_38),]
  
  if(is_empty(qtl_tmp$gene_id_clean)){
    
    next()
    
  } else {
    
    #######################################################
    #Let's add the information from OTG to align the betas#
    #######################################################
    
    #STEP 1: retrieve info to do alignment (frequencies for lead and proxy)
    
    print(index)
    
    variant_otg <- tryCatch(otargen::variantInfo(variant), error=function(e) NA)
    lead_otg <- tryCatch(otargen::variantInfo(query), error=function(e) NA)
    
    ############################################################################################
    #Careful, if we have a mismatch for some reason, we will add a comment and jump on the next#
    ############################################################################################
    
    if(length(variant_otg) == 1| length(lead_otg) == 1){ #one of them have missing info:
      
      proxies$gene_id.vat_sqtl[index] <- paste(qtl_tmp$gene_id_clean, collapse = ";")
      proxies$gene_name.vat_sqtl[index] <- paste(qtl_tmp$gene_name, collapse = ";")
      proxies$splicing_region.vat[index] <- paste(qtl_tmp$splicing, collapse = ";")
      
      proxies$ir_allele.vat_sqtl[index] <- qtl_tmp$alt[1]
      proxies$is_allele.vat_sqtl[index] <- qtl_tmp$ref[1]
      
      proxies$splicing_effect.vat[index] <- paste(qtl_tmp$afc, collapse = ";")
      proxies$splicing_effect.vat[index] <- paste("non_aligned", proxies$splicing_effect.vat[index], sep= "_") #adding that the allele might be wrong! Needs to be checked
      
      proxies$pip.vat_sqtl[index] <- paste(qtl_tmp$pip, collapse = ";")
      #proxies$cs.vat_sqtl[index] <- paste(qtl_tmp$cs_id, collapse = ";")
      
      next()
      
    }
    
    ############################################
    #If all good, let's extract the AF for both#
    ############################################
    
    #STEP 1: first for the lead:
    
    ir_allele <- ir_new$effect_allele[which(ir_new$variant == query)]
    
    af_ref <- ifelse(ir_allele == lead_otg$altAllele, as.numeric(lead_otg$gnomadNFE), 1-as.numeric(lead_otg$gnomadNFE))
    
    af_ref_res <- ifelse(af_ref < 0.5, "min", "max")
    
    #STEP 2: let's check the frequency for the proxy:
    
    qtl_allele <- qtl_tmp$alt[1] # we just need one, all associations are aligned to the same allele. Hence, having only one of the alleles should be more than enough!
    
    af_qtl <- ifelse(qtl_allele == variant_otg$altAllele, as.numeric(variant_otg$gnomadNFE), 1-as.numeric(variant_otg$gnomadNFE))
    
    af_qtl_res <- ifelse(af_qtl < 0.5, "min", "max")
    
    #####################################################################################################
    #IF af_qtl_res and af_ref_res are the same, we have the beta aligned to the same correlated alleles!#
    #####################################################################################################
    
    if(af_qtl_res == af_ref_res){
      
      proxies$gene_id.vat_sqtl[index] <- paste(qtl_tmp$gene_id_clean, collapse = ";")
      proxies$gene_name.vat_sqtl[index] <- paste(qtl_tmp$gene_name, collapse = ";")
      proxies$splicing_region.vat[index] <- paste(qtl_tmp$splicing, collapse = ";")
      
      proxies$ir_allele.vat_sqtl[index] <- qtl_tmp$alt[1]
      proxies$is_allele.vat_sqtl[index] <- qtl_tmp$ref[1]
      
      proxies$splicing_effect.vat[index] <- paste(qtl_tmp$afc, collapse = ";")
      proxies$pip.vat_sqtl[index] <- paste(qtl_tmp$pip, collapse = ";")
      #proxies$cs.vat_sqtl[index] <- paste(qtl_tmp$cs_id, collapse = ";")
      
    } else {
      
      ###########################################
      #But if we don't then we need to change!!!#
      ###########################################
      
      print("exception")
      
      proxies$gene_id.vat_sqtl[index] <- paste(qtl_tmp$gene_id_clean, collapse = ";")
      proxies$gene_name.vat_sqtl[index] <- paste(qtl_tmp$gene_name, collapse = ";")
      proxies$splicing_region.vat[index] <- paste(qtl_tmp$splicing, collapse = ";")
      
      proxies$ir_allele.vat_sqtl[index] <- qtl_tmp$ref[1]
      proxies$is_allele.vat_sqtl[index] <- qtl_tmp$alt[1]
      
      proxies$splicing_effect.vat[index] <- paste(as.numeric(qtl_tmp$afc)*(-1), collapse = ";")
      proxies$pip.vat_sqtl[index] <- paste(qtl_tmp$pip, collapse = ";")
      #proxies$cs.vat_sqtl[index] <- paste(qtl_tmp$cs_id, collapse = ";")
      
    }
    
  }
  
}

saveRDS(proxies, "6_qtl_annotation/conditional_and_fine_mapped_asat_and_vat_eqtls_sqtls.txt")
