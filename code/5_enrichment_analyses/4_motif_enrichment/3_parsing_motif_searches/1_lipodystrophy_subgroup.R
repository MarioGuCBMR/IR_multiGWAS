##############
#INTRODUCTION#
##############

#This code reads the motif enrichment data and tries to identify the overlap of the IR peaks

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)

###################
#Loading functions#
###################

check_motif_overlap <- function(peak_start, peak_end, offset, motif_length, variant_position, strand) {
  #This first finds the region where the motif might be...
  
  # Check if strand is "+" or "-"
  
  if (strand == "+") {
    
    motif_start <- peak_start + offset
    
    motif_end <- motif_start + motif_length - 1
    
    } else if (strand == "-") {
      
    motif_end <- peak_end - offset
    
    motif_start <- motif_end - motif_length + 1
    
    } 
  
  # Check if variant overlaps with motif
  
  if (variant_position >= motif_start && variant_position <= motif_end) {
    
    return(TRUE)  # Variant overlaps motif
    
    } else {
      
      return(FALSE)  # No overlap
      
    }
  
}


##############
#Loading data#
##############

#First let's get the peaks at day 4, we need this data as reference:

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/")

peaks <- fread('raw_data/atac_seq_data/GSE178794_SGBS-day4_consensusPeaks-fromUnionRepPeaks_top100k-byMedianPval.bed.gz')

proxies <- fread("output/7_functional_annotation/1_atac_seq_overlaps/proxies_matched_w_roadmap_and_perrin_data_w_sgbs_diff.txt")

#Let's load the variants in this cluster, specifically,

bmi_neg <- fread("output/5_enrichment_analyses/3_go_shifter/input_bmi_neg.txt")

proxies <- proxies[which(proxies$query_snp_rsid%in%bmi_neg$SNP),]

#Let's get the motif matches:

matches <- fread("output/5_enrichment_analyses/4_homer/output/63_sgbs_day_4/cebp_motif_instances.txt")

#Let's get the ChIP data:

haploreg <- fread("output/4_ir_loci_discovery/4_proxies/282_proxies_core_15.txt", fill = TRUE)
haploreg_match <- haploreg[which(haploreg$rsID%in%proxies$rsID),]
haploreg_match <- haploreg_match[order(match(haploreg_match$rsID, proxies$rsID)),]

length(which(haploreg_match$rsID == proxies$rsID))

proxies$protein_binding <- haploreg_match$Proteins
proxies$motifs <- haploreg_match$Motifs

#####################################################
#Let's get the proxies that are dirsupting the motif#
#####################################################

proxies_cebpb <- proxies[which(str_detect(proxies$motifs, "CEBPB_")),] #33

proxies_match <- proxies[which(proxies$query_snp_rsid%in%proxies_cebpb$query_snp_rsid),]



proxies$source <- NA

for(index_proxy in seq(1, length(proxies$rsID))){
  
  print(index_proxy)
  
  #STEP 1: 
  
  chr_ <- paste("chr", proxies$chr[index_proxy], sep = "")
  pos_ <- proxies$pos_hg19[index_proxy]
  
  #STEP 2: get the matches:
  
  chip_match <- chip[which(chip$chr == chr_ & chip$start <= pos_ & chip$stop >= pos_),]
  chip_match$binding_region <- paste(chip_match$chr, "_", chip_match$start, "_", chip_match$stop, sep = "")
  
  #STEP 3: add the data
  
  if(is_empty(chip_match$chr)){
    
    next()
    
  } else {
    
    print("SUCCESS")
    
    proxies$cebpb_chip[index_proxy] <- paste(chip_match$binding_region, collapse = ";")
    proxies$source[index_proxy] <- paste(chip_match$source, collapse = ";")
    
  }
  
}

#########################################
#Let's assess how many of them are there#
#########################################

leads_cepbp <- unique(proxies$rsID[which(is.na(proxies$cebpb_chip) == FALSE)])

proxies_cebpb_match <- proxies[which(proxies$rsID%in%leads_cepbp),]

length(unique(proxies_cebpb_match$query_snp_rsid)) #23

#Let's take those that are open only after day 4:

leads_day_4 <- unique(proxies_cebpb_match$rsID[which(proxies_cebpb_match$mature_classification == "More accessible in adipocytes (day 4)")]) #10

proxies_day_4 <- proxies_cebpb_match[which(proxies_cebpb_match$rsID%in%leads_day_4),]

length(unique(proxies_day_4$query_snp_rsid)) #8

################################################
#Let's get information on these loci right away#
################################################

proxies_day_4$closest_gene <- NA

for(index_proxy in seq(1, length(proxies_day_4$rsID))){
  
  tmp <-  tryCatch(otargen::variantInfo(proxies_day_4$rsID[index_proxy]), error = function(e) NULL)
  
  proxies_day_4$closest_gene[index_proxy] <- paste(unique(tmp$nearestCodingGene.symbol), collapse=";")
  
  
}

######################################################################
#Let's add the QTL and the enhancer data, let's keep the ball rolling#
######################################################################

qtl_df <- readRDS("output/6_qtl_annotation/conditional_and_fine_mapped_asat_and_vat_eqtls_sqtls.txt")

qtl_df <- qtl_df[which(qtl_df$rsID%in%proxies_day_4$rsID),]
qtl_df <- qtl_df[which(is.na(qtl_df$beta.asat) == FALSE | is.na(qtl_df$afc.asat_gtex) == FALSE | is.na(qtl_df$pip.asat_sqtl) == FALSE),]

#Let's add the haploreg data:

haploreg <- fread("output/4_ir_loci_discovery/4_proxies/282_proxies_core_15.txt", fill = TRUE)
haploreg_match <- haploreg[which(haploreg$rsID%in%qtl_df$rsID),]
haploreg_match <- haploreg_match[order(match(haploreg_match$rsID, qtl_df$rsID)),]

###########################################################
#Let's do this again for all variants with proxies binding#
###########################################################

leads_cepbp <- unique(proxies$query_snp_rsid[which(is.na(proxies$cebpb_chip) == FALSE)])

proxies_cebpb_match <- proxies[which(proxies$query_snp_rsid%in%leads_cepbp),]

haploreg_match <- haploreg[which(haploreg$rsID%in%proxies_cebpb_match$rsID),]
haploreg_match <- haploreg_match[order(match(haploreg_match$rsID, proxies_cebpb_match$rsID)),]

###########################################################
#Let's do the same with the enhancer-promoter contact data#
###########################################################

enhancer_df <- fread("output/7_functional_annotation/2_enhancer_gene_links/2_enhancer_gene_links_and_hic/proxies_with_stare_and_hic_data.txt")

enhancer_df <- enhancer_df[which(enhancer_df$query_snp_rsid%in%proxies_day_4$query_snp_rsid),]

###################################################################################
#Let's take the original data from Haploreg and assess the motif data from there!!#
###################################################################################

haploreg <- fread("output/4_ir_loci_discovery/4_proxies/282_proxies_core_15.txt", fill = TRUE)
haploreg_match <- haploreg[which(haploreg$rsID%in%proxies_day_4$rsID),]
haploreg_match <- haploreg_match[order(match(haploreg_match$rsID, proxies_day_4$rsID)),]

length(which(proxies$rsID == haploreg_match$rsID)) #perfect

proxies$motif_data <- haploreg_match$Motifs

#########################################################
#Let's explore those that are matching with our criteria#
#########################################################

proxies_cepbp <- proxies[which(proxies$mature_classification == "More accessible in adipocytes (day 4)" & is.na(proxies$cebpb_chip) == FALSE),]
