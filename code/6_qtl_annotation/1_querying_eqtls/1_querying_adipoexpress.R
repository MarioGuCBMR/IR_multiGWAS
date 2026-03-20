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
  
  return(chr)
  
}

pos_parser <- function(snp){
  
  pos <- unlist(str_split(snp, "_")[[1]][2])
  
  return(pos)
  
}

variant_parser <- function(snp){
  
  data <- str_split(snp, "_")[[1]]
  
  variant <- paste(data[1], "_", data[2], sep ="")
  
  return(variant)
  
}

other_allele_parser <- function(snp){
  
  data <- str_split(snp, "_")[[1]]
  
  allele <- data[3]
  
  return(allele)
  
}

effect_allele_parser <- function(snp){
  
  data <- str_split(snp, "_")[[1]]
  
  allele <- data[4]
  
  return(allele)
  
}

##############
#Loading data#
##############

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/output")

ir_new <- fread("4_ir_loci_discovery/3_novel_variants/3_novel_variants/282_ir_variants_with_ir_and_t2d_comparisons_w_closest_genes.txt")

proxies <- fread("4_ir_loci_discovery/4_proxies/proxies_build_37_and_38_4_282_variants.txt")

yes_vect <- c("A", "G", "C", "T")

proxies <- proxies[which(proxies$query_snp_rsid%in%ir_new$variant),] #3188
proxies <- proxies[which(proxies$ref%in%yes_vect & proxies$alt%in%yes_vect),]

length(which(duplicated(proxies$query_snp_rsid) == FALSE)) #282

#Finally, let's load the adipoexpress data:

qtl_df <- fread("../../IR_GSEM_2023/raw_data/SuppTable4_EURonly_eQTL_signals.txt")

###########################################################
#Let's first try to get all data associated to the proxies#
###########################################################

qtl_df$chr <- unlist(sapply(qtl_df$variant, chr_parser))
qtl_df$pos <- unlist(sapply(qtl_df$variant, pos_parser))
qtl_df$variant_clean <- unlist(sapply(qtl_df$variant, variant_parser))
qtl_df$other_allele <- unlist(sapply(qtl_df$variant, other_allele_parser))
qtl_df$effect_allele <- unlist(sapply(qtl_df$variant, effect_allele_parser))

################################################
#Let's try matching the variants to the proxies#
################################################

proxies$gene_id.asat <- NA
proxies$gene_name.asat <- NA
proxies$ir_allele.asat <- NA
proxies$is_allele.asat <- NA
proxies$beta.asat <- NA
proxies$standard_error.asat <- NA
proxies$p_value.asat <- NA

for(index in seq(1, length(proxies$rsID))){
    
    #STEP 1: get the info we want:
    
    variant= proxies$rsID[index]
    query = proxies$query_snp_rsid[index]
    
    chr_37 <- proxies$chr[index]
    bp_37 = proxies$pos_hg19[index]
    
    #Then we query the data in qtl_df
    
    qtl_tmp <- qtl_df[which(qtl_df$chr == chr_37 & qtl_df$pos == bp_37),]
    
    if(is_empty(qtl_tmp$ENSG)){
      
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
        
        proxies$gene_id.asat[index] <- paste(qtl_tmp$ENSG, collapse = ";")
        proxies$gene_name.asat[index] <- paste(qtl_tmp$gene, collapse = ";")
        
        proxies$ir_allele.asat[index] <- qtl_tmp$effect_allele[1]
        proxies$is_allele.asat[index] <- qtl_tmp$other_allele[1]
        
        proxies$beta.asat[index] <- paste(qtl_tmp$beta, collapse = ";")
        proxies$beta.asat[index] <- paste("non_aligned", proxies$beta.asat[index], sep= "_") #adding that the allele might be wrong! Needs to be checked
        
        proxies$standard_error.asat[index] <- paste(qtl_tmp$se, collapse = ";")
        proxies$p_value.asat[index] <- paste(qtl_tmp$pval_joint, collapse = ";")
        
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
      
      qtl_allele <- qtl_tmp$effect_allele[1] # we just need one, all associations are aligned to the same allele. Hence, having only one of the alleles should be more than enough!
      
      af_qtl <- ifelse(qtl_allele == variant_otg$altAllele, as.numeric(variant_otg$gnomadNFE), 1-as.numeric(variant_otg$gnomadNFE))
      
      af_qtl_res <- ifelse(af_qtl < 0.5, "min", "max")
      
      #####################################################################################################
      #IF af_qtl_res and af_ref_res are the same, we have the beta aligned to the same correlated alleles!#
      #####################################################################################################
      
      if(af_qtl_res == af_ref_res){
        
        proxies$gene_id.asat[index] <- paste(qtl_tmp$ENSG, collapse = ";")
        proxies$gene_name.asat[index] <- paste(qtl_tmp$gene, collapse = ";")
        
        proxies$ir_allele.asat[index] <- qtl_tmp$effect_allele[1]
        proxies$is_allele.asat[index] <- qtl_tmp$other_allele[1]
        
        proxies$beta.asat[index] <- paste(qtl_tmp$beta, collapse = ";")
        proxies$standard_error.asat[index] <- paste(qtl_tmp$se, collapse = ";")
        proxies$p_value.asat[index] <- paste(qtl_tmp$pval_joint, collapse = ";")
        
      } else {
        
        ###########################################
        #But if we don't then we need to change!!!#
        ###########################################
        
        print("exception")

        proxies$gene_id.asat[index] <- paste(qtl_tmp$ENSG, collapse = ";")
        proxies$gene_name.asat[index] <- paste(qtl_tmp$gene, collapse = ";")
        
        proxies$ir_allele.asat[index] <- qtl_tmp$other_allele[1]
        proxies$is_allele.asat[index] <- qtl_tmp$effect_allele[1]
        
        proxies$beta.asat[index] <- paste(qtl_tmp$beta*(-1), collapse = ";") #only thing we need to do is align the results here
        proxies$standard_error.asat[index] <- paste(qtl_tmp$se, collapse = ";")
        proxies$p_value.asat[index] <- paste(qtl_tmp$pval_joint, collapse = ";")
        
      }
    
    }
  
}

#####################
#Let's save the data#
#####################

dir.create("6_qtl_annotation")

saveRDS(proxies, "6_qtl_annotation/conditional_asat_eqtls.txt")
