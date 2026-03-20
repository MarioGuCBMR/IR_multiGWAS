##############
#INTRODUCTION#
##############

#This code performs a E2G pipeline that tries to identify enhancer-promoter interactions and links them to genes.

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)
library(MotifDb)

##############
#Loading data#
##############

#First let's load the data with enhancer-promoter evidence:

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025")

#STEP 1: get the BMI-decreasing hits:

enhancer_df <- fread("output/7_functional_annotation/2_enhancer_gene_links/2_enhancer_gene_links_and_hic/proxies_with_stare_and_hic_data.txt")

bmi_neg <- fread("output/5_enrichment_analyses/3_go_shifter/input_bmi_neg.txt")

enhancer_df <- enhancer_df[which(enhancer_df$query_snp_rsid%in%bmi_neg$SNP),]

#STEP 2: take the original IR variants, we will need this:

ir_variants <- fread("output/4_ir_loci_discovery/3_novel_variants/4_colocalization/352_variants_with_all_info_plus_fiadjbmi_and_tg_hdl_coloc.txt")

ir_variants <- ir_variants[which(ir_variants$variant%in%enhancer_df$query_snp_rsid),]

#STEP 3: get QTL data for later:

qtl_df <- readRDS("output/6_qtl_annotation/conditional_and_fine_mapped_asat_and_vat_eqtls_sqtls.txt")

#STEP 4: let's read the ChIP data:

#chip <- fread("raw_data/chip_atlas/SRX7807102.20.bed")

####################################################################################
#STEP 2: let's get all the proxies from the signals that meet the following criteria
####################################################################################

#They are in day 4 peaks:

enhancer_day_4 <- enhancer_df[which(enhancer_df$sgbs_day4 != ""),]

leads_day_4 <- unique(enhancer_day_4$query_snp_rsid) #42/63

#Let's load the variants in this cluster, 

enhancer_df$cebpb <- NA

for(index_proxy in seq(1, length(enhancer_df$rsID))){
  
  #STEP 1: get the proxy
  
  chr_ <- paste("chr", enhancer_df$chr[index_proxy],sep="")
  pos_ <- enhancer_df$pos_hg19[index_proxy]
  
  #STEP 2: find the match in chip data
  
  chip_match <- chip[which(chip$V1 == chr_ & as.numeric(chip$V2) <= pos_ & as.numeric(chip$V3) >= pos_),]
  chip_match$id <- paste(chip_match$V1, "_", chip_match$V2, "_", chip_match$V3, sep = "")
  
  if(is_empty(chip_match$V1)){
    
    next()
    
  } else {
    
    enhancer_df$cebpb[index_proxy] <- paste(chip_match$id, collapse=";")
    
  }
  
}

#Find the leads enriched for this:

enhancer_cebpb <- enhancer_df[which(is.na(enhancer_df$cebpb) == FALSE),]
leads_cebpb <- unique(enhancer_cebpb$query_snp_rsid) #8

#And now let's find the overlap:

lead_day_4_cebpb <- Reduce(intersect, list(leads_cebpb, leads_day_4)) #all of them!!

#################################################
#Let's do a quick look at these loci in QTL data#
#################################################

qtl_match <- qtl_df[which(qtl_df$query_snp_rsid%in%lead_day_4_cebpb),]

##############################################
#Let's run now motif breaker for these fellas#
##############################################

snps_4_break <- enhancer_df$rsID[which(enhancer_df$query_snp_rsid%in%leads_day_4)] #183

library(motifbreakR)
library(SNPlocs.Hsapiens.dbSNP155.GRCh37)
library(BSgenome.Hsapiens.UCSC.hg19)

#Let's extract the variants:

variants <- motifbreakR::snps.from.rsid(rsid = snps_4_break,
                                        dbSNP = SNPlocs.Hsapiens.dbSNP155.GRCh37,
                                        search.genome = BSgenome.Hsapiens.UCSC.hg19)

#Let's run motif breaker for only the human motifs

motif_db <- query(MotifDb, "Hsapiens")  # Filter for Hsapiens motifs - NOT EXCLUSIVE!

# Exclude "stamlabs" motifs
motif_db <- motif_db[!grepl("stamlab", names(motif_db))]

# Further filter for specific motif databases
pattern <- "(cisbp_1.02|HOCOMOCOv10|HOCOMOCOv11|hPDI|JASPAR_2014|JASPAR_CORE|jaspar2016|jaspar2018|jolma2013|SwissRegulon|UniPROBE)"

# Filter motifs based on the pattern matching the source names in the motif names
selected_motifs <- motif_db[grepl(pattern, names(motif_db))]

#We cannot do serial parametization here so let's change things a bit

library(BiocParallel)
register(SerialParam())

results <- motifbreakR(snpList = variants, filterp = TRUE,
                       pwmList = selected_motifs,
                       threshold = 1e-4,
                       method = "ic")

saveRDS(results, "output/7_functional_annotation/3_motif_analyses/ALL_day_ir_variants_disruption.RDS")

# #STEP 2: Let's explore the data:
# 
# motif_day_4_df <- as.data.frame(results, row.names = NULL)
# 
# #STEP 3: let's explore the strong data:
# 
# motif_day_4_df <- motif_day_4_df[which(motif_day_4_df$effect == "strong"),]
# 
# #Let's get the combination of IDs that truly match:
# 
# motif_day_4_df$id_ <- paste(motif_day_4_df$SNP_id, "_", motif_day_4_df$REF, "_", motif_day_4_df$ALT, sep  = "")
# enhancer_df$id_1 <- paste(enhancer_df$rsID, "_", enhancer_df$alt, "_", enhancer_df$ref, sep = "")
# enhancer_df$id_2 <- paste(enhancer_df$rsID, "_", enhancer_df$ref, "_", enhancer_df$alt, sep = "")
# 
# match_day_4_df <- motif_day_4_df[which(motif_day_4_df$id_%in%enhancer_df$id_1 | motif_day_4_df$id_%in%enhancer_df$id_2),]
# 
# #Let's remove duplicates now:
# 
# match_day_4_df$id_prot <- paste(match_day_4_df$SNP_id, "_", match_day_4_df$geneSymbol, sep = "")
# match_day_4_df <- match_day_4_df[order(abs(match_day_4_df$alleleEffectSize), decreasing = TRUE),]
# match_day_4_df <- match_day_4_df[which(duplicated(match_day_4_df$id_prot) == FALSE),] #782
# 
# check_day_4 <- as.data.frame(table(match_day_4_df$geneSymbol))
# check_day_4 <- check_day_4[order(check_day_4$Freq, decreasing = TRUE),]
# 
# #Let's add family information here, we can have sort the signal later:
# 
# check_day_4$family <- as.character(unlist(sapply(check_day_4$Var1, get_tf_family_annotation)))
# 
# ###################################
# #Let's try to add this information#
# ################################### IMPORTANT TO HAVE IT, BUT WE CAN DO IT LATER
# 
# #Let's get the results:
# 
# # proxies <- enhancer_df
# # 
# # proxies$protein_motif <- NA
# # proxies$strand_motif <- NA
# # proxies$ir_allele.motif <- NA
# # proxies$is_allele.motif <- NA
# # proxies$effect_binding <- NA
# # 
# # for(index in seq(1, length(proxies$rsID))){
# #   
# #   #STEP 1: get the info we want:
# #   
# #   variant= proxies$rsID[index]
# #   query = proxies$query_snp_rsid[index]
# #   
# #   #Let's find the matching:
# #   
# #   motif_tmp <- match_day_4_df[which(match_day_4_df$SNP_id == variant),]
# #   
# #   if(is_empty(motif_tmp$seqnames)){
# #     
# #     next()
# #     
# #   } else {
# #     
# #     #######################################################
# #     #Let's add the information from OTG to align the betas#
# #     #######################################################
# #     
# #     #STEP 1: retrieve info to do alignment (frequencies for lead and proxy)
# #     
# #     print(index)
# #     
# #     variant_otg <- tryCatch(otargen::variantInfo(variant), error=function(e) NA)
# #     lead_otg <- tryCatch(otargen::variantInfo(query), error=function(e) NA)
# #     
# #     ############################################################################################
# #     #Careful, if we have a mismatch for some reason, we will add a comment and jump on the next#
# #     ############################################################################################
# #     
# #     if(length(variant_otg) == 1| length(lead_otg) == 1){ #one of them have missing info:
# #       
# #       proxies$protein_motif[index] <- paste(motif_tmp$geneSymbol, collapse = ";")
# #       proxies$strand_motif[index] <- paste(motif_tmp$strand, collapse = ";")
# #       
# #       proxies$ir_allele.motif[index] <- motif_tmp$ALT[1]
# #       proxies$is_allele.motif[index] <- motif_tmp$REF[1]
# #       
# #       proxies$effect_binding[index] <- paste(motif_tmp$alleleEffectSize, collapse = ";")
# #       proxies$effect_binding[index] <- paste("non_aligned", proxies$effect_binding[index], sep= "_") #adding that the allele might be wrong! Needs to be checked
# # 
# #       next()
# #       
# #     }
# #     
# #     ############################################
# #     #If all good, let's extract the AF for both#
# #     ############################################
# #     
# #     #STEP 1: first for the lead:
# #     
# #     ir_allele <- ir_variants$effect_allele[which(ir_variants$variant == query)]
# #     
# #     af_ref <- ifelse(ir_allele == lead_otg$altAllele, as.numeric(lead_otg$gnomadNFE), 1-as.numeric(lead_otg$gnomadNFE))
# #     
# #     af_ref_res <- ifelse(af_ref < 0.5, "min", "max")
# #     
# #     #STEP 2: let's check the frequency for the proxy:
# #     
# #     qtl_allele <- motif_tmp$ALT[1] # we just need one, all associations are aligned to the same allele. Hence, having only one of the alleles should be more than enough!
# #     
# #     af_qtl <- ifelse(qtl_allele == variant_otg$altAllele, as.numeric(variant_otg$gnomadNFE), 1-as.numeric(variant_otg$gnomadNFE))
# #     
# #     af_qtl_res <- ifelse(af_qtl < 0.5, "min", "max")
# #     
# #     #####################################################################################################
# #     #IF af_qtl_res and af_ref_res are the same, we have the beta aligned to the same correlated alleles!#
# #     #####################################################################################################
# #     
# #     if(af_qtl_res == af_ref_res){
# #       
# #       proxies$protein_motif[index] <- paste(motif_tmp$geneSymbol, collapse = ";")
# #       proxies$strand_motif[index] <- paste(motif_tmp$strand, collapse = ";")
# #       
# #       proxies$ir_allele.motif[index] <- motif_tmp$ALT[1]
# #       proxies$is_allele.motif[index] <- motif_tmp$REF[1]
# #       
# #       proxies$effect_binding[index] <- paste(motif_tmp$alleleEffectSize, collapse = ";")
# # 
# #       
# #     } else {
# #       
# #       ###########################################
# #       #But if we don't then we need to change!!!#
# #       ###########################################
# #       
# #       print("exception")
# #       
# #       proxies$protein_motif[index] <- paste(motif_tmp$geneSymbol, collapse = ";")
# #       proxies$strand_motif[index] <- paste(motif_tmp$strand, collapse = ";")
# #       
# #       proxies$ir_allele.motif[index] <- motif_tmp$REF[1]
# #       proxies$is_allele.motif[index] <- motif_tmp$ALT[1]
# #       
# #       proxies$effect_binding[index] <- paste(as.numeric(motif_tmp$alleleEffectSize)*(-1), collapse = ";")
# #       
# #     }
# #     
# #   }
# #   
# # }
# 
# ###############################################################
# #Let's make a dataframe trying to compile all this information#
# ###############################################################
# 
# disrupted_tf <- as.data.frame(table(match_day_4_df$geneSymbol))
# disrupted_tf$family <- NA
# 
# #Great, now let's add the families here:
# 
# for(index in seq(1, length(disrupted_tf$Var1))){
#   
#   gene_ <- disrupted_tf$Var1[index]
#   
#   family_ <- get_tf_family_annotation(gene_)
#   
#   disrupted_tf$family[index] <- ifelse(is.null(family_), NA, family_)
#   
# }
# 
# #############################################################################
# #And now, with this info, let's go and try to find the best TFs for proxy!!!#
# #############################################################################
# 
# #The key is to loop by proxy and retrieve only those that are in the same family and have the best score:
# 
# for(index_query in seq(1, length(enhancer_df$query_snp_rsid))){
#   
#   #STEP 1: get the proxies:
#   
#   lead <- enhancer_df$query_snp_rsid[index_query]
#   proxies <-  enhancer_df$rsID[which(enhancer_df$query_snp_rsid == lead)]
#   
#   #STEP 2: get the motifs:
#   
#   motif_tmp <- match_day_4_df[which(match_day_4_df$SNP_id%in%proxies),]
#   
#   #STEP 3: get the families:
#   
#   motif_tmp$family <- as.character(unlist(sapply(motif_tmp$geneSymbol, get_tf_family_annotation)))
#   motif_tmp<- motif_tmp[which(is.na(motif_tmp$family) == FALSE),]
#   motif_tmp$family <- as.character(unlist(sapply(motif_tmp$family, family_cleaner)))
#   
#   
#   #We will remove the families that are missing to be sure we are not introducing biases:
#   #Great, now we are going to split the data according to the family and find the genes that match the familes.
# 
#   motif_family_reducted <- separate_rows(motif_tmp, family, sep = ";")
#   motif_family_reducted <- motif_family_reducted[which(is.na(motif_family_reducted$family) == FALSE & motif_family_reducted$family != "NA"),]
#   
#   motif_family_reducted$family <- as.character(unlist(sapply(motif_family_reducted$family, family_cleaner)))
#   
#   
# 
#   
# }
# 
# #######################################
# #Let's get the data for each dataframe#
# #######################################
# 
# for(index in seq(1, length(check_mature$Var1))){
# 
#   #STEP 1: get gene and rsID
#   
#   gene_ <- as.character(check_mature$Var1)[index]
#   
#   #STEP 2: get variants:
#   
#   rsid <- match_mature$SNP_id[which(match_mature$geneSymbol == gene_)]
#   
#   #STEP 3: get query SNPs:
#   
#   leads <- day_4$query_snp_rsid[which(day_4$rsID%in%rsid)]
#   
#   #STEP 4: add this info: 
#   
#   dummy_df <- rep(gene_, length(leads))
#   dummy_df <- cbind(dummy_df, leads)
#   
#   if(!(exists("final_df"))){
#     
#     final_df <- dummy_df
#     
#     } else {
#     
#     final_df <- rbind(final_df, dummy_df)
#     
#   }
#   
# }
# 
# 
# 
# #Let's go and classify each gene per family:
# 
# for(index in seq(1, length(final_df$dummy_df))){
#   
#   gene_ <- final_df$dummy_df[index]
#   
#   family_ <- get_tf_family_annotation(gene_)
#   
#   final_df$family[index] <- ifelse(is.null(family_), NA, family_)
#    
# }
# 
# ##############################
# #Now let's divide the data...#
# ##############################
# 
# day_4 <- proxies[which(proxies$mature_classification == "More accessible in adipocytes (day 4)")]
# all_days <- proxies[which(proxies$mature_classification == "Similarly accessible across days")]
# preadipocytes <- proxies[which(proxies$mature_classification == "More accessible in preadipocytes")]
# 
# match_mature <- match_day_4_df[which(match_day_4_df$SNP_id%in%day_4$rsID),] #782
# match_preadipocytes <- match_day_4_df[which(match_day_4_df$SNP_id%in%preadipocytes$rsID),] #782
# match_across <- match_day_4_df[which(match_day_4_df$SNP_id%in%all_days$rsID),] #782
# 
# check_mature <- as.data.frame(table(match_mature$geneSymbol))
# check_mature <- check_mature[order(check_mature$Freq, decreasing = TRUE),]
# 
# check_preadipocytes <- as.data.frame(table(match_preadipocytes$geneSymbol))
# check_preadipocytes <- check_preadipocytes[order(check_preadipocytes$Freq, decreasing = TRUE),]
# 
# check_across <- as.data.frame(table(match_across$geneSymbol))
# check_across <- check_across[order(check_across$Freq, decreasing = TRUE),]
# 
# 
# ####################################################
# #Let's perform enrichment analyses for each overlap#
# ####################################################
# 
# prot_day_4 <- as.character(check_mature$Var1)
# prot_pre <- as.character(check_preadipocytes$Var1)
# prot_across <- as.character(check_across$Var1)
# prot_all <- c(prot_day_4, prot_pre, prot_across)
# prot_all <- unique(prot_all)
# 
# #Let's do the enrichment:
# 
# library(STRINGdb)
# 
# prot_all <- data.frame(gene = prot_all, stringsAsFactors = FALSE)
# 
# string_db <- STRINGdb$new(version = "11.5", species = 9606, score_threshold = 400, input_directory = "")
# mapped_genes <- string_db$map(prot_all, "gene", removeUnmappedRows = TRUE)
# enrichment_all <- string_db$get_enrichment(mapped_genes$STRING_id)
# 
# #####################################
# #Let's do it only for those in day 4#
# #####################################
# 
# prot_day_4 <- data.frame(gene = prot_day_4, stringsAsFactors = FALSE)
# 
# mapped_genes <- string_db$map(prot_day_4, "gene", removeUnmappedRows = TRUE)
# enrichment_day_4 <- string_db$get_enrichment(mapped_genes$STRING_id)
# 
# ###################################################
# #Let's do it only for those in across adipogenesis#
# ###################################################
# 
# prot_across <- data.frame(gene = prot_across, stringsAsFactors = FALSE)
# 
# mapped_genes <- string_db$map(prot_across, "gene", removeUnmappedRows = TRUE)
# enrichment_across <- string_db$get_enrichment(mapped_genes$STRING_id)
# 
# #############################################
# #Let's do it only for those in preadipocytes#
# #############################################
# 
# prot_pre <- data.frame(gene = prot_pre, stringsAsFactors = FALSE)
# 
# mapped_genes <- string_db$map(prot_pre, "gene", removeUnmappedRows = TRUE)
# enrichment_pre <- string_db$get_enrichment(mapped_genes$STRING_id)
