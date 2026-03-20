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

my_otg_function <- function (variant_id, lead_snp) 
{
  cli::cli_progress_step("Connecting to database..", spinner = TRUE)
  otg_cli <- ghql::GraphqlClient$new(url = "https://api.genetics.opentargets.org/graphql")
  otg_qry <- ghql::Query$new()
  if (grepl(pattern = "rs\\d+", variant_id)) {
    query_searchid <- "query ConvertRSIDtoVID($queryString:String!) {\n    search(queryString:$queryString){\n      totalVariants\n      variants{\n        id\n        }\n      }\n    }"
    variables <- list(queryString = variant_id)
    otg_qry$query(name = "convertid", x = query_searchid)
    id_result <- jsonlite::fromJSON(otg_cli$exec(otg_qry$queries$convertid, 
                                                 variables), flatten = TRUE)$data
    input_variant_id <- id_result$search$variants$id
  }
  else if (grepl(pattern = "\\d+_\\d+_[a-zA-Z]+_[a-zA-Z]+", 
                 variant_id)) {
    input_variant_id <- variant_id
  }
  else {
    stop("\n Please provide a variant Id")
  }
  query <- "query v2gquery($variantId: String!){\n  genesForVariant(variantId: $variantId) {\n    gene{\n    id\n    symbol\n    }\n  variant\n  overallScore\n  qtls{\n    typeId\n    aggregatedScore\n    tissues{\n      tissue{\n        id\n        name\n      }\n      quantile\n      beta\n      pval\n    }\n  }\n  intervals{\n  typeId\n  sourceId\n  aggregatedScore\n  tissues{\n  tissue{\n  id\n  name\n  }\n  quantile\n  score\n  }\n  }\n  functionalPredictions{\n    typeId\n    sourceId\n    aggregatedScore\n    tissues{\n      tissue{\n        id\n        name\n      }\n      maxEffectLabel\n      maxEffectScore\n    }\n  }\n  distances{\n   typeId\n    sourceId\n    aggregatedScore\n    tissues{\n      tissue{\n        id\n        name\n      }\n      distance\n      score\n      quantile\n    }\n  }\n  }\n}"
  result_pkg = list()
  variables <- list(variantId = input_variant_id)
  otg_qry$query(name = "v2g_query", x = query)
  cli::cli_progress_step(paste0("Downloading data for ", variant_id, 
                                " ..."), spinner = TRUE)
  result <- jsonlite::fromJSON(otg_cli$exec(otg_qry$queries$v2g_query, 
                                            variables), flatten = TRUE)$data
  
  if(is.null(unlist(result))){
    
    return(NA)
    
  }
  
  result_df <- as.data.frame(result$genesForVariant) %>% dplyr::arrange(desc(overallScore))
  
  result_df$rsid <- variant_id
  result_df$lead_snp <- lead_snp
  
  
  if(is.null(unlist(result_df$qtls))){
    
    result_core <- result_df %>% dplyr::select(gene.symbol, 
                                               variant, overallScore, gene.id, rsid, lead_snp)
    
    result_pkg <- list(v2g = result_core)
    
    return(result_pkg)
    
  } 

  if (nrow(result_df) != 0) {
    result_core <- result_df %>% dplyr::select(gene.symbol, 
                                               variant, overallScore, gene.id)
    result_qtl <- result_df %>% dplyr::select(gene.symbol, 
                                              variant, qtls) %>% tidyr::unnest(qtls, names_sep = ".", 
                                                                               keep_empty = TRUE) %>% dplyr::rename(typeId = "qtls.typeId", 
                                                                                                                    aggregatedScore = "qtls.aggregatedScore")
    if ("qtls.tissues" %in% colnames(result_qtl)) {
      result_qtl <- result_qtl %>% tidyr::unnest(qtls.tissues, 
                                                 names_sep = "_", keep_empty = TRUE) %>% dplyr::rename(tissues_id = "qtls.tissues_tissue.id", 
                                                                                                       tissues_name = "qtls.tissues_tissue.name")
      base::colnames(result_qtl) <- stringr::str_replace_all(colnames(result_qtl), 
                                                             "qtls.", "")
      
      
    }
    
    result_qtl$rsid <- variant_id
    result_qtl$lead_snp <- lead_snp
    
    #result_intervals <- result_df %>% dplyr::select(gene.symbol, 
    #                                                variant, intervals) %>% tidyr::unnest(intervals, 
    #                                                                                      names_sep = ".", keep_empty = TRUE) %>% dplyr::rename(typeId = "intervals.typeId", 
    #                                                                                                                                            aggregatedScore = "intervals.aggregatedScore")
    #if ("intervals.tissues" %in% colnames(result_intervals)) {
    #  result_intervals <- result_intervals %>% tidyr::unnest(intervals.tissues, 
    #                                                         names_sep = "_", keep_empty = TRUE) %>% dplyr::rename(tissues_id = "intervals.tissues_tissue.id", 
    #                                                                                                               tissues_name = "intervals.tissues_tissue.name")
    #  base::colnames(result_intervals) <- stringr::str_replace_all(colnames(result_intervals), 
    #                                                               "intervals.", "")
    #}
    #result_distances <- result_df %>% dplyr::select(gene.symbol, 
    #                                                variant, distances) %>% tidyr::unnest(distances, 
    #                                                                                      names_sep = ".", keep_empty = TRUE) %>% dplyr::rename(typeId = "distances.typeId", 
    #                                                                                                                                            aggregatedScore = "distances.aggregatedScore")
    #if ("distances.tissues" %in% colnames(result_distances)) {
    #  result_distances <- result_distances %>% tidyr::unnest(distances.tissues, 
    #                                                         names_sep = "_", keep_empty = TRUE) %>% dplyr::rename(tissues_id = "distances.tissues_tissue.id", 
    #                                                                                                               tissues_name = "distances.tissues_tissue.name")
    #  base::colnames(result_distances) <- stringr::str_replace_all(colnames(result_distances), 
    #                                                               "distances.", "")
    #}
    #result_functionalPredictions <- result_df %>% dplyr::select(gene.symbol, 
    #                                                            variant, functionalPredictions) %>% tidyr::unnest(functionalPredictions, 
    #                                                                                                              names_sep = ".", keep_empty = TRUE) %>% dplyr::rename(typeId = "functionalPredictions.typeId", 
    #                                                                                                                                                                    aggregatedScore = "functionalPredictions.aggregatedScore")
    #if ("functionalPredictions.tissues" %in% colnames(result_functionalPredictions)) {
    #  result_functionalPredictions <- result_functionalPredictions %>% 
    #    tidyr::unnest(functionalPredictions.tissues, 
    #                  names_sep = "_", keep_empty = TRUE) %>% dplyr::rename(tissues_id = "functionalPredictions.tissues_tissue.id", 
    #                                                                        tissues_name = "functionalPredictions.tissues_tissue.name")
    #  base::colnames(result_functionalPredictions) <- stringr::str_replace_all(colnames(result_functionalPredictions), 
    #                                                                           "functionalPredictions.", "")
    #}
    
    result_core$rsid <- variant_id
    result_core$lead_snp <- lead_snp
    
    
    result_pkg <- list(v2g = result_core,  qtls = result_qtl)
  }
  cli::cli_progress_update()
  return(result_pkg)
}

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

##############
#Loading data#
##############

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025//output")

ir_new <- fread("4_ir_loci_discovery/3_novel_variants/3_novel_variants/282_ir_variants_with_ir_and_t2d_comparisons_w_closest_genes.txt")

proxies <- readRDS("6_qtl_annotation/conditional_and_fine_mapped_asat_vat_and_liver_eqtls_sqtls.txt")

#Finally, let's load the QTL data:

qtl_df <- arrow::read_parquet("../../gtex_v10/GTEx_v10_SuSiE_eQTL/Muscle_Skeletal.v10.eQTLs.SuSiE_summary.parquet")

#Now let's add a chr_pos column to get info

proxies$chr_pos <- paste("chr", proxies$chr, ":", proxies$pos_hg38, sep ="")

#Now the same for the QTLs

qtl_df$chr <- unlist(sapply(qtl_df$variant_id, chr_parser))
qtl_df$pos <- unlist(sapply(qtl_df$variant_id, pos_parser))
qtl_df$ref <- unlist(sapply(qtl_df$variant_id, ref_parser))
qtl_df$alt <- unlist(sapply(qtl_df$variant_id, alt_parser))
qtl_df$gene_id_clean <- unlist(sapply(qtl_df$phenotype_id, gene_parser))

qtl_df$chr_pos <- paste("chr", qtl_df$chr, ":", qtl_df$pos, sep = "")

#################################################################################################################################
#Not all variants have information, we only have info for one of the variants of each credible set..., we need to solve this out#
#################################################################################################################################

#We can loop over the phenotype_id and credible sets and make a quick fix:

qtl_df$phenotype_cs_id <- paste(qtl_df$phenotype_id, "_", qtl_df$cs_id, sep = "")

list_of_ids <- unique(qtl_df$phenotype_cs_id)

for(id_ in list_of_ids){
  
  qtl_tmp <- qtl_df[which(qtl_df$phenotype_cs_id == id_),]
  
  afc_ <- qtl_tmp$afc[which(is.na(qtl_tmp$afc) == FALSE)]
  afc_se_ <- qtl_tmp$afc_se[which(is.na(qtl_tmp$afc_se) == FALSE)]
  
  qtl_df$afc[which(qtl_df$phenotype_cs_id == id_)] <- afc_ #checked, this worked like a charm
  qtl_df$afc_se[which(qtl_df$phenotype_cs_id == id_)] <- afc_se_ #checked, this worked like a charm
  
}

#################################################################
#Now we can properly fill the data with the data that is missing#
#################################################################

proxies$gene_id.muscle_gtex <- NA
proxies$gene_name.muscle_gtex <- NA

proxies$ir_allele.muscle_gtex <- NA
proxies$is_allele.muscle_gtex <- NA

proxies$afc.muscle_gtex <- NA
proxies$standard_error.muscle_gtex <- NA
proxies$pip.muscle_gtex <- NA

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
      
      proxies$gene_id.muscle_gtex[index] <- paste(qtl_tmp$gene_id_clean, collapse = ";")
      proxies$gene_name.muscle_gtex[index] <- paste(qtl_tmp$gene_name, collapse = ";")
      
      proxies$ir_allele.muscle_gtex[index] <- qtl_tmp$alt[1]
      proxies$is_allele.muscle_gtex[index] <- qtl_tmp$ref[1]
      
      proxies$afc.muscle_gtex[index] <- paste(qtl_tmp$afc, collapse = ";")
      proxies$afc.muscle_gtex[index] <- paste("non_aligned", proxies$afc.muscle_gtex[index], sep= "_") #adding that the allele might be wrong! Needs to be checked
      
      proxies$standard_error.muscle_gtex[index] <- paste(qtl_tmp$afc_se, collapse = ";")
      proxies$pip.muscle_gtex[index] <- paste(qtl_tmp$pip, collapse = ";")
      #proxies$cs.muscle_gtex[index] <- paste(qtl_tmp$cs_id, collapse = ";")
      
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
      
      proxies$gene_id.muscle_gtex[index] <- paste(qtl_tmp$gene_id_clean, collapse = ";")
      proxies$gene_name.muscle_gtex[index] <- paste(qtl_tmp$gene_name, collapse = ";")
      
      proxies$ir_allele.muscle_gtex[index] <- qtl_tmp$alt[1]
      proxies$is_allele.muscle_gtex[index] <- qtl_tmp$ref[1]
      
      proxies$afc.muscle_gtex[index] <- paste(qtl_tmp$afc, collapse = ";")
      proxies$standard_error.muscle_gtex[index] <- paste(qtl_tmp$afc_se, collapse = ";")
      proxies$pip.muscle_gtex[index] <- paste(qtl_tmp$pip, collapse = ";")
      #proxies$cs.muscle_gtex[index] <- paste(qtl_tmp$cs_id, collapse = ";")
      
    } else {
      
      ###########################################
      #But if we don't then we need to change!!!#
      ###########################################
      
      print("exception")
      
      proxies$gene_id.muscle_gtex[index] <- paste(qtl_tmp$gene_id_clean, collapse = ";")
      proxies$gene_name.muscle_gtex[index] <- paste(qtl_tmp$gene_name, collapse = ";")
      
      proxies$ir_allele.muscle_gtex[index] <- qtl_tmp$ref[1]
      proxies$is_allele.muscle_gtex[index] <- qtl_tmp$alt[1]
      
      proxies$afc.muscle_gtex[index] <- paste(as.numeric(qtl_tmp$afc)*(-1), collapse = ";")
      proxies$standard_error.muscle_gtex[index] <- paste(qtl_tmp$afc_se, collapse = ";")
      proxies$pip.muscle_gtex[index] <- paste(qtl_tmp$pip, collapse = ";")
      #proxies$cs.muscle_gtex[index] <- paste(qtl_tmp$cs_id, collapse = ";")
      
    }
    
  }
  
}

saveRDS(proxies, "6_qtl_annotation/conditional_and_fine_mapped_asat_vat_liver_and_muscle_eqtls.txt")
