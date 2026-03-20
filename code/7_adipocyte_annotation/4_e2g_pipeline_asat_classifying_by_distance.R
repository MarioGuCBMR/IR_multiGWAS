##############
#INTRODUCTION#
##############

#This code performs a E2G pipeline that tries to identify enhancer-promoter interactions and links them to genes.

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)

##############
#Loading data#
##############

#First let's load the data with enhancer-promoter evidence:

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025")

enhancer_df <- fread("output/7_functional_annotation/2_enhancer_gene_links/2_enhancer_gene_links_and_hic/proxies_with_stare_and_hic_data.txt")

ir_variants <- fread("output/4_ir_loci_discovery/3_novel_variants/4_colocalization/352_variants_with_all_info_plus_fiadjbmi_and_tg_hdl_coloc.txt")

ir_variants <- ir_variants[which(ir_variants$variant%in%enhancer_df$query_snp_rsid),]

qtl_df <- readRDS("output/6_qtl_annotation/conditional_and_fine_mapped_asat_and_vat_eqtls_sqtls.txt")

##############################################################################################
#STEP 0: let's check how many of the variants were connected to at least 1 gene through STARE#
##############################################################################################

check <- enhancer_df[which(enhancer_df$stare_adipose_nuclei != "" | enhancer_df$stare_day_0 != "" | enhancer_df$stare_day_4 != "" | enhancer_df$stare_day_14 != "" ),]

length(unique(check$query_snp_rsid)) #187

#################################################################################################################################################################
#STEP 0: let's not assume that does that do not have associations from conditional leads do not have a link. It is not fine-mapped, we do not have credible sets#
#################################################################################################################################################################

for(index_lead in seq(1, length(ir_variants$variant))){
  
  #STEP 1: get the lead_
  
  lead_ <- ir_variants$variant[index_lead]
  
  #STEP 2: get the qtl_df for this lead:
  
  qtl_tmp <- qtl_df[which(qtl_df$query_snp_rsid == lead_),]
  
  #STEP 3: let's merge the data for all of them:
  
  gene_id <- qtl_tmp$gene_id.asat
  gene_name <- qtl_tmp$gene_name.asat
  betas <- qtl_tmp$beta.asat
  ses <- qtl_tmp$standard_error.asat
  pvals <- qtl_tmp$p_value.asat
  
  #Let's clean the data:
  
  gene_id <- gene_id[which(is.na(gene_id) == FALSE)]
  gene_name <- gene_name[which(is.na(gene_name) == FALSE)]
  betas <- betas[which(is.na(betas) == FALSE)]
  ses <- ses[which(is.na(ses) == FALSE)]
  pvals <- pvals[which(is.na(pvals) == FALSE)]
  
  #Let's merge it:
  
  gene_id <- paste(gene_id, collapse = ";")
  gene_name <- paste(gene_name, collapse = ";")
  betas <- paste(betas, collapse = ";")
  ses <- paste(ses, collapse = ";")
  pvals <- paste(pvals, collapse = ";")
  
  #Let's change it:
  
  qtl_df$gene_id.asat[which(qtl_df$query_snp_rsid == lead_)] <- gene_id
  qtl_df$gene_name.asat[which(qtl_df$query_snp_rsid == lead_)] <- gene_name
  qtl_df$beta.asat[which(qtl_df$query_snp_rsid == lead_)] <- betas
  qtl_df$standard_error.asat[which(qtl_df$query_snp_rsid == lead_)] <- ses
  qtl_df$p_value.asat[which(qtl_df$query_snp_rsid == lead_)] <- pvals
  
}

############################################################################
#STEP 1: let's add info in a single dataframe so that it is easier to parse#
############################################################################

qtl_df_ <- qtl_df %>%
  select(c(1, 10:45))

#Let's merge both datasets:

merged_info <- merge(enhancer_df, qtl_df_, by.x = "rsID", by.y = "rsID", all.x = T) #this worked!!

#Let's filter for those that have data (careful these are truly saved like NAs!) No idea why these differences, they are all saved the same way.
#We are only working with those that have ASAT for now...

enhancer_linked_w_qtl <- merged_info[which(merged_info$gene_id.asat != "" |
                                             is.na(merged_info$gene_id.asat_gtex) ==FALSE |
                                             is.na(merged_info$gene_id.asat_sqtl) == FALSE | 
                                             is.na(merged_info$gene_id.vat_gtex) ==FALSE |
                                             is.na(merged_info$gene_id.vat_sqtl) == FALSE),] #3236

length(unique(enhancer_linked_w_qtl$query_snp_rsid)) #118

#####################################
#STEP 2: let's find the common genes#
#####################################

#Let's clean the data obtained from all adipose tissue data:

gene_ids_qtls <- c(enhancer_linked_w_qtl$gene_id.asat, enhancer_linked_w_qtl$gene_id.asat_gtex, enhancer_linked_w_qtl$gene_id.asat_sqtl, enhancer_linked_w_qtl$gene_id.vat_gtex, enhancer_linked_w_qtl$gene_id.vat_sqtl)
gene_ids_qtls <- gene_ids_qtls[which(is.na(gene_ids_qtls) == FALSE)]
gene_ids_qtls <- as.character(unlist(str_split(gene_ids_qtls, ";")))
gene_ids_qtls <- unique(gene_ids_qtls) #211
table(gene_ids_qtls) #perfect! No weird stuff

#Let's do the same for the links obtained with SGBS

gene_ids_stare <- c(enhancer_linked_w_qtl$gene_id_adipose_nuclei,enhancer_linked_w_qtl$gene_id_day_4, enhancer_linked_w_qtl$gene_id_day_14)
gene_ids_stare <- gene_ids_stare[which(gene_ids_stare != "")] #careful with differences in format of saving...
gene_ids_stare <- as.character(unlist(str_split(gene_ids_stare, ";")))
gene_ids_stare <- unique(gene_ids_stare) #664
table(gene_ids_stare) #perfect! No weird stuff

#OK, let's get the match:

common_gene_ids <- Reduce(intersect, list(gene_ids_qtls, gene_ids_stare)) #119

#######################################################
#STEP 2: let's work on finding the best genes per loci#
#######################################################

ir_variants <- ir_variants[which(ir_variants$variant%in%enhancer_linked_w_qtl$query_snp_rsid),] #118

ir_variants$gene_id_qlt_abc <- NA
ir_variants$gene_name_qlt_abc <- NA

for(index_lead in seq(1, length(ir_variants$variant))){
  
  #STEP 1: get the lead
  
  lead_ <- ir_variants$variant[index_lead]
  
  #STEP 2: get the proxies for the variant:
  
  proxy_tmp <- enhancer_linked_w_qtl[which(enhancer_linked_w_qtl$query_snp_rsid == lead_),]
  
  #STEP 3: let's combine all gene_ids found for this proxy and get only those that are in the common_gene_ids:
  
  gene_ids_proxy <- c(proxy_tmp$gene_id.asat, proxy_tmp$gene_id.asat_gtex, proxy_tmp$gene_id.asat_sqtl, proxy_tmp$gene_id.vat_gtex, proxy_tmp$gene_id.vat_sqtl) #it is only the shared ones, so this should work
  gene_ids_proxy <- unique(as.character(unlist(str_split(gene_ids_proxy, ";"))))
  gene_ids_proxy <- gene_ids_proxy[which(gene_ids_proxy%in%common_gene_ids)]
  
  #STEP 4: go through a loop and make sure we can have this info:
  
  common_gene_names <- c()
  
  #Let's get the names for each of them and make a freaking dictionary
  
  for(gene_ in gene_ids_proxy){ #just in case we have more than one
    
    try_otg <- tryCatch({
      otargen::geneticConstraintQuery(gene_)  # hypothetical otargen call
    }, error = function(e) {
      message(paste("Error with gene:", gene_, ":", e$message))
      return(NULL)
    })
    
    if(is.null(try_otg)){
      
      #print(index_gene)
      
      common_gene_names <- c(common_gene_names, gene_)
      
    } else {
      
      common_gene_names <- c(common_gene_names, try_otg$approvedSymbol[1]) #this provides data for 3 scores, we are using this as a trick to retrieve data
      
    }
    
  }
  
  #STEP 5: add the data:
  
  ir_variants$gene_id_qlt_abc[index_lead] <- paste(gene_ids_proxy, collapse = ";")
  ir_variants$gene_name_qlt_abc[index_lead] <- paste(common_gene_names, collapse = ";")
  
  
}

##########################################
#Let's clean the enhancer data for now...#
##########################################

#Let's remove those that have not worked out...

ir_variants <- ir_variants[which(ir_variants$gene_id_qlt_abc != ""),] #72! Great

#And let's split the dataframe in the way that things work out:

ir_variants <- separate_rows(ir_variants, gene_id_qlt_abc, gene_name_qlt_abc, sep = ";")

ir_variants$qtl_effect <- NA

#And now let's loop over the data to know whether the link is VAT or ASAT

for(index_gene in seq(1, length(ir_variants$gene_id_qlt_abc))){
  
  #STEP 1: get lead:
  
  lead_ <- ir_variants$variant[index_gene]
  gene_ <- ir_variants$gene_id_qlt_abc[index_gene]
  
  #STEP 2: get qtl data for the genes, asat and vat
  
  qtl_tmp <- qtl_df[which(qtl_df$query_snp_rsid == lead_),]
  
  #STEP 3: let's get the data for each analyses:
  
  asat_genes <- as.character(unlist(str_split(qtl_tmp$gene_id.asat, ";"))) 
  asat_gtex_genes <- as.character(unlist(str_split(qtl_tmp$gene_id.asat_gtex, ";")))
  asat_sqtl_genes <- as.character(unlist(str_split(qtl_tmp$gene_id.asat_sqtl, ";"))) 
  vat_gtex_genes <- as.character(unlist(str_split(qtl_tmp$gene_id.vat_gtex, ";"))) 
  vat_sqtl_genes <- as.character(unlist(str_split(qtl_tmp$gene_id.vat_sqtl, ";"))) 
  
  #Let's also get the betas for each of them: 
  
  asat_betas <- as.character(unlist(str_split(qtl_tmp$beta.asat, ";"))) 
  asat_gtex_betas <- as.character(unlist(str_split(qtl_tmp$afc.asat_gtex, ";")))
  asat_gtex_pips <- as.character(unlist(str_split(qtl_tmp$pip.asat_gtex, ";")))
  
  asat_sqtl_betas <- as.character(unlist(str_split(qtl_tmp$splicing_region, ";")))
  asat_sqtl_pips <- as.character(unlist(str_split(qtl_tmp$pip.asat_sqtl, ";")))
  
  vat_gtex_betas <- as.character(unlist(str_split(qtl_tmp$afc.vat_gtex, ";"))) 
  vat_gtex_pips <- as.character(unlist(str_split(qtl_tmp$pip.vat_gtex, ";"))) 
  
  vat_sqtl_betas <- as.character(unlist(str_split(qtl_tmp$splicing_region.vat, ";"))) 
  vat_sqtl_pips <- as.character(unlist(str_split(qtl_tmp$pip.vat_sqtl, ";"))) 
  
  #STEP 4: let's filter for those that actually match
  
  index_gene_asat <- which(asat_genes == gene_)
  index_gene_asat_gtex <- which(asat_gtex_genes == gene_)
  index_gene_sqtl <- which(asat_sqtl_genes == gene_)
  index_gene_vat <- which(vat_gtex_genes == gene_)
  index_gene_vat_sqtl <- which(vat_sqtl_genes == gene_)
  
  #STEP 5: let's filter for all of them
  
  asat_betas <- asat_betas[index_gene_asat][1] #just one effect per gene, no need for more
  
  asat_gtex_betas <- asat_gtex_betas[index_gene_asat_gtex][1] #just one effect per gene, no need for more
  #asat_sqtl_pips <- asat_sqtl_pips[index_gene_asat_gtex][1] #just one effect per gene, no need for more
  
  asat_sqtl_betas <- asat_sqtl_betas[index_gene_sqtl][1] #just one effect per gene, no need for more
  #asat_gtex_pips <- asat_gtex_pips[index_gene_asat_gtex][1] #just one effect per gene, no need for more
  
  vat_gtex_betas <- vat_gtex_betas[index_gene_vat][1] #just one effect per gene, no need for more
  #vat_gtex_pips <- vat_gtex_pips[index_gene_vat][1] #just one effect per gene, no need for more
  
  vat_sqtl_betas <- vat_sqtl_betas[index_gene_vat_sqtl][1] #just one effect per gene, no need for more
  #vat_sqtl_pips <- vat_sqtl_pips[index_gene_vat_sqtl][1] #just one effect per gene, no need for more
  
  #We will not add PIPs, but we will ensure the checks just in case...
  
  # if(is.na(asat_gtex_betas) == TRUE & is.na(asat_gtex_pips) == FALSE | asat_gtex_betas == "NA" & is.na(asat_gtex_pips) == FALSE){
  #   
  #   print("ASAT GTEx")
  #   print(gene_)
  #   
  # }
  # 
  # if(is.na(asat_sqtl_betas) == TRUE & is.na(asat_sqtl_pips) == FALSE | asat_sqtl_betas == "NA" & is.na(asat_sqtl_pips) == FALSE){
  #   
  #   print(index_gene)
  #   print("ASAT sQTL")
  #   print(gene_)
  #   
  # }
  # 
  # if(vat_gtex_betas == "NA" & is.na(vat_gtex_pips) == FALSE | is.na(vat_gtex_betas) == TRUE & is.na(vat_gtex_pips) == FALSE){
  #   
  #   print("VAT GTEx")
  #   print(gene_)
  #   
  # }
  # 
  # if(is.na(vat_sqtl_betas) == TRUE & is.na(vat_sqtl_pips) == FALSE | vat_sqtl_betas == "NA" & is.na(vat_sqtl_pips) == FALSE){
  #   
  #   print("VAT sQTL")
  #   print(gene_)
  #   
  # }
  # 

  final_score = paste("beta (ASAT)=", asat_betas, ";",
                      "afc (ASAT)=", asat_gtex_betas, ";",
                      "Splicing (ASAT)=", asat_sqtl_betas, ";",
                      "afc (VAT)=", vat_gtex_betas, ";",
                      "Splicing (VAT)=", vat_sqtl_betas)
  
  ir_variants$qtl_effect[index_gene] <- final_score
   
  
}

###################################################################################################################################################################
#Great, now we are going to get the proxies that overlap the most for each variant!! The proxies that have the most amount of annotations will win, for simplicity#
###################################################################################################################################################################

#Now we are going to go variant per variant 

ir_variants$proxy <- NA
ir_variants$hic_annotation <- NA
ir_variants$overlapping_annotations <- NA

for(index_gene in seq(1, length(ir_variants$gene_id_qlt_abc))){
  
  #STEP 1 get lead and gene:'
  
  lead_ <- ir_variants$variant[index_gene]
  gene_ <- ir_variants$gene_id_qlt_abc[index_gene]
  
  #STEP 2: let's the data from the matches:
  
  proxy_tmp <- enhancer_linked_w_qtl[which(query_snp_rsid == lead_),]
  proxy_tmp <- proxy_tmp[which(str_detect(proxy_tmp$gene_id.asat, gene_) |
                       str_detect(proxy_tmp$gene_id.asat_gtex, gene_) |
                       str_detect(proxy_tmp$gene_id.asat_sqtl, gene_) |
                       str_detect(proxy_tmp$gene_id.vat_gtex, gene_) |
                       str_detect(proxy_tmp$gene_id.vat_sqtl, gene_)),]
  
  #STEP 3: let's do a small check and add whether we have Hi-C data:
  
  check <- which(proxy_tmp$adipocyte_bait_1 != "" |
                             proxy_tmp$adipocyte_other_2 != "" |
                             proxy_tmp$adipocyte_bait_3 != "" |
                             proxy_tmp$adipocyte_other_4 != "" |
                             proxy_tmp$differentiated_bait_1 != "" |
                             proxy_tmp$differentiated_other_2 != "" |
                             proxy_tmp$adipocyte_bait_3 != "")
  
  if(!(is_empty(check))){
    
    ir_variants$hic_annotation[index_gene] <- "Yes"
    
    }
  
  #STEP 4: let's count the annotations for Hi-C data
  
  proxy_tmp$score <- 0
  proxy_tmp$score <- ifelse(proxy_tmp$tss_dist_adipose_nuclei != "", proxy_tmp$score+1, proxy_tmp$score)
  proxy_tmp$score <- ifelse(proxy_tmp$tss_dist_day_4 != "", proxy_tmp$score+1, proxy_tmp$score)
  proxy_tmp$score <- ifelse(proxy_tmp$tss_dist_day_14 != "", proxy_tmp$score+1, proxy_tmp$score)
  proxy_tmp$score <- ifelse(proxy_tmp$adipocyte_bait_1 != "", proxy_tmp$score+1, proxy_tmp$score) #match is done in bait_1, other_1 is the links
  proxy_tmp$score <- ifelse(proxy_tmp$adipocyte_other_2 != "", proxy_tmp$score+1, proxy_tmp$score) #match is done in other_2, bait_1 is the links
  proxy_tmp$score <- ifelse(proxy_tmp$adipocyte_bait_3 != "", proxy_tmp$score+1, proxy_tmp$score) #match is done in bait_3, other_3 is the links
  proxy_tmp$score <- ifelse(proxy_tmp$adipocyte_other_4 != "", proxy_tmp$score+1, proxy_tmp$score) #match is done in other_4, bait_5 is the links
  proxy_tmp$score <- ifelse(proxy_tmp$differentiated_bait_1 != "", proxy_tmp$score+1, proxy_tmp$score) #match is done in bait_1, other_1 is the links
  proxy_tmp$score <- ifelse(proxy_tmp$differentiated_other_2 != "", proxy_tmp$score+1, proxy_tmp$score) #match is done in other_2, bait_2 is the links
  
  max_score <- proxy_tmp$score[which.max(proxy_tmp$score)]
  
  proxy_tmp <- proxy_tmp[which(proxy_tmp$score == max_score),]
  proxy_tmp$chr_pos_b37 <- paste("chr", proxy_tmp$chr, ":", proxy_tmp$pos_hg19, sep = "")
  
  ir_variants$proxy[index_gene] <- paste(proxy_tmp$chr_pos_b37, collapse = ";")
  
  region_adipose_betas <- paste(unique(proxy_tmp$stare_adipose_nuclei), collapse="/") 
  region_day_4 <- paste(unique(proxy_tmp$stare_day_4), collapse="/") 
  region_day_14 <- paste(unique(proxy_tmp$stare_day_14), collapse="/") 
  
  region_hic_adipose_bait_1 <- paste(unique(proxy_tmp$adipocyte_bait_1), collapse="/")
  region_hic_adipose_other_1 <- paste(unique(proxy_tmp$adipocyte_other_2), collapse="/") 
  region_hic_adipose_1 <- paste(unique(c(region_hic_adipose_bait_1, region_hic_adipose_other_1)), collapse = "/")
  
  region_hic_adipose_bait_2 <- paste(unique(proxy_tmp$adipocyte_bait_3), collapse="/")
  region_hic_adipose_other_2 <- paste(unique(proxy_tmp$adipocyte_other_4), collapse="/") 
  region_hic_adipose_2 <- paste(unique(c(region_hic_adipose_bait_2, region_hic_adipose_other_2)), collapse = "/")
  
  region_hic_diff_bait <- paste(unique(proxy_tmp$differentiated_bait_1), collapse="/")
  region_hic_diff_other <- paste(unique(proxy_tmp$differentiated_other_2), collapse="/") 
  region_hic_diff <- paste(unique(c(region_hic_diff_bait, region_hic_diff_other)), collapse = "/")
  
  final_score = paste("Enhancer-gene link (adipose)=", region_adipose_betas, ";",
                      "Enhancer-gene link (Day 4)=", region_day_4, ";",
                      "Enhancer-gene link (Day 14)=", region_day_14, ";",
                      "Enhancer-promoter contact (adipocytes 1)=", region_hic_adipose_1, ";",
                      "Enhancer-promoter contact (adipocytes 2)=", region_hic_adipose_2, ";",
                      "Enhancer-promoter contact (differentiated SGBS)=", region_hic_diff)  
  
  ir_variants$overlapping_annotations[index_gene] <- final_score
  
  
}

##########################################################
#Let's assess in which subgroups these variants are found#
##########################################################

bmi_ns <- fread("output/5_enrichment_analyses/3_go_shifter/ld_bmi_ns.txt")
bmi_ns <- ir_variants[which(ir_variants$variant%in%bmi_ns$RsIdA),]

bmi_neg <- fread("output/5_enrichment_analyses/3_go_shifter/ld_bmi_neg.txt")
bmi_neg <- ir_variants[which(ir_variants$variant%in%bmi_neg$RsIdA),]

bmi_pos <- fread("output/5_enrichment_analyses/3_go_shifter/ld_bmi_pos.txt")
bmi_pos <- ir_variants[which(ir_variants$variant%in%bmi_pos$RsIdA),]

bmi_outlier <- fread("output/5_enrichment_analyses/3_go_shifter/ld_outliers.txt")
bmi_outlier <- ir_variants[which(ir_variants$variant%in%bmi_outlier$RsIdA),]

###########################
#Let's try adding the data#
###########################

ir_variants$cluster <- NA

for(index_lead in seq(1, length(ir_variants$variant))){
  
  #Variant:
  
  rsid <- ir_variants$variant[index_lead]

  if(rsid%in%bmi_ns$variant){
    
    ir_variants$cluster[index_lead] <- "BMI (P>0.05)"
    
  }
  
  if(rsid%in%bmi_neg$variant){
    
    ir_variants$cluster[index_lead] <- "BMI- (P<0.05)"
    
  }
  
  if(rsid%in%bmi_pos$variant){
    
    ir_variants$cluster[index_lead] <- "BMI+ (P<0.05)"
    
  }
  
  if(rsid%in%bmi_outlier$variant){
    
    ir_variants$cluster[index_lead] <- "BMI+;GSATadjBMI-"
    
  }
  
  
}


#########################################################
#Let's save the data in a format we are comfortable with#
#########################################################

working_df <- ir_variants %>%
  select(variant, reported_ir, ir_source, cluster, gene_name_qlt_abc, qtl_effect, proxy, hic_annotation, overlapping_annotations)

working_df <- working_df[order(working_df$cluster, working_df$gene_name_qlt_abc, decreasing = TRUE),]

fwrite(working_df, "output/7_functional_annotation/2_enhancer_gene_links/2_enhancer_gene_links_and_hic/enhancer_qtl_matches_for_publication.txt")

# novel_df <- working_df[which(working_df$reported_ir == "no" & working_df$ir_source == ""),]
# 
# fwrite(novel_df, "output/7_functional_annotation/2_enhancer_gene_links/2_enhancer_gene_links_and_hic/novel_enhancer_qtl_matches_for_table_1.txt")
# 
# 
# #################################################################################################################
# #Let's work with the novel ones just to show off, but the supplementary table needs to be checked and worked one#
# #################################################################################################################
# 
# #######################################################################################################################
# #Let's loop over the data and make a column with the common enhancer-gene and QTL genes so that we can work with those#
# #######################################################################################################################
# 
# #Let's loop over the variant:
# 
# enhancer_linked_w_qtl$enhancer_qtl_gene <- NA
# 
# for(query_rsid in enhancer_linked_w_qtl$query_snp_rsid){
#   
#   #STEP 1: get the proxy df:
#   
#   proxy_tmp <- enhancer_linked_w_qtl[which(enhancer_linked_w_qtl$query_snp_rsid == query_rsid),]
#   
#   #STEP 2: get the genes:
#   
#   genes <- c(proxy_tmp$gene_id.asat, proxy_tmp$gene_id.asat_gtex, proxy_tmp$gene_id.asat_sqtl, proxy_tmp$gene_id.vat_gtex, proxy_tmp$gene_id.vat_sqtl)
#   genes <- genes[which(is.na(genes) == FALSE)]
#   genes <- as.character(unlist(str_split(genes, ";")))
#   genes <- unique(genes)
#   
#   genes <- genes[which(genes%in%common_gene_ids)]
#   
#   #Let's get the gene names:
#   
#   positions_id_vect <- which(common_gene_ids%in%genes)
#   gene_names <- common_gene_names[positions_id_vect]
#   gene_names <- gene_names[order(gene_names)]
#   
#   enhancer_linked_w_qtl$enhancer_qtl_gene[which(enhancer_linked_w_qtl$query_snp_rsid == query_rsid)] <- paste(gene_names, collapse = ";")
#   
# }
# 
# ##################################################################################
# #Now let's filter for those variants that are in accessible peaks with STARE data# - an assess which ones are in promoters and such
# ##################################################################################
# 
# enhancer_linked_w_qtl_clean <- enhancer_linked_w_qtl[which(enhancer_linked_w_qtl$stare_adipose_nuclei != "" |
#                                                              enhancer_linked_w_qtl$stare_day_4 != "" |
#                                                              enhancer_linked_w_qtl$stare_day_14!= ""),]
# 
# enhancer_linked_w_qtl_clean <- enhancer_linked_w_qtl_clean[which(enhancer_linked_w_qtl_clean$enhancer_qtl_gene != ""),] #270
# 
# #Check how many independent associations do we have:
# 
# length(unique(enhancer_linked_w_qtl_clean$query_snp_rsid)) #61
# 
# #Let's get the data a cleaned:
# 
# working_df <- enhancer_linked_w_qtl_clean %>%
#   select(query_snp_rsid, rsID, enhancer_qtl_gene, gene_name_adipose_nuclei, tss_dist_adipose_nuclei, gene_name_day_4, tss_dist_day_4, gene_name_day_14, tss_dist_day_14, gene_name.asat, gene_name.asat_gtex, gene_name.asat_sqtl, gene_name.vat_gtex, gene_name.vat_sqtl)
# 
# #Now let's work only on the proxy with most annotations:
# 
# working_df$score <- 0
# working_df$score <- ifelse(working_df$tss_dist_adipose_nuclei != "", working_df$score+1, working_df$score)
# working_df$score <- ifelse(working_df$tss_dist_day_4 != "", working_df$score+1, working_df$score)
# working_df$score <- ifelse(working_df$tss_dist_day_14 != "", working_df$score+1, working_df$score)
# working_df$score <- ifelse(is.na(working_df$gene_name.asat) != TRUE, working_df$score+1, working_df$score)
# working_df$score <- ifelse(is.na(working_df$gene_name.asat_gtex) != TRUE, working_df$score+1, working_df$score)
# working_df$score <- ifelse(is.na(working_df$gene_name.asat_sqtl) != TRUE, working_df$score+1, working_df$score)
# working_df$score <- ifelse(is.na(working_df$gene_name.vat_gtex) != TRUE, working_df$score+1, working_df$score)
# working_df$score <- ifelse(is.na(working_df$gene_name.vat_sqtl) != TRUE, working_df$score+1, working_df$score)
# 
# working_df <- working_df[order(working_df$score, decreasing = TRUE),]
# working_df <- working_df[which(duplicated(working_df$enhancer_qtl_gene) == FALSE),]
# 
# ##################################################################################
# #Let's clean this data, it is going to make it so much easier to parse afterwards#
# ##################################################################################
# 
# for(index_proxy in seq(1, length(working_df$query_snp_rsid))){
#   
#   #STEP 1: let's get the genes per row and filter
#   
#   tmp_df <- working_df[index_proxy,]
#   
#   ref_genes <- as.character(unlist(str_split(tmp_df$enhancer_qtl_gene, ";")))
#   
#   genes_adipose <- as.character(unlist(str_split(tmp_df$gene_name_adipose_nuclei, ";")))
#   genes_day4 <- as.character(unlist(str_split(tmp_df$gene_name_day_4, ";")))
#   genes_day14 <- as.character(unlist(str_split(tmp_df$gene_name_day_14, ";")))
#   
#   #Let's get the indexes:
#   
#   index_adipose <- which(genes_adipose%in%ref_genes)
#   index_day4 <- which(genes_day4%in%ref_genes)
#   index_day14 <- which(genes_day14%in%ref_genes)
#   
#   genes_adipose <- genes_adipose[index_adipose]
#   genes_day4 <- genes_day4[index_day4]
#   genes_day14 <- genes_day14[index_day14]
#   
#   #Same with TSS:
#   
#   tss_adipose <- as.character(unlist(str_split(tmp_df$tss_dist_adipose_nuclei, ";")))
#   tss_day4 <- as.character(unlist(str_split(tmp_df$tss_dist_day_4, ";")))
#   tss_day14 <- as.character(unlist(str_split(tmp_df$tss_dist_day_14, ";")))
#   
#   tss_adipose <- tss_adipose[index_adipose]
#   tss_day4 <- tss_day4[index_day4]
#   tss_day14 <- tss_day14[index_day14]
#   
#   #And now let's subsitute:
#   
#   working_df$gene_name_adipose_nuclei[index_proxy] <- paste(genes_adipose, collapse = ";")
#   working_df$gene_name_day_4[index_proxy] <- paste(genes_day4, collapse = ";")
#   working_df$gene_name_day_14[index_proxy] <- paste(genes_day14, collapse = ";")
#   
#   working_df$tss_dist_adipose_nuclei[index_proxy] <- paste(tss_adipose, collapse = ";")
#   working_df$tss_dist_day_4[index_proxy] <- paste(tss_day4, collapse = ";")
#   working_df$tss_dist_day_14[index_proxy] <- paste(tss_day14, collapse = ";")
#   
# }
# 
# #################################################################################
# #Let's now separate the data across the subtypes and reveal the novel ones too!!#
# #################################################################################
# 
# ir_variants <- fread("output/4_ir_loci_discovery/3_novel_variants/3_novel_variants/282_ir_variants_with_ir_and_t2d_comparisons_w_closest_genes.txt")
# novel_variants <- ir_variants[which(ir_variants$reported_ir == "no"),]
# novel_variants <- novel_variants[which(novel_variants$ir_source != "/Oliveri_gw"),] #70, just in case
# 
# #Now let's add the data for each subset of variants:
# 
# #Let's load only the fat distribution variants:
# 
# bmi_ns <- fread("output/5_enrichment_analyses/3_go_shifter/ld_bmi_ns.txt")
# bmi_ns <- ir_variants[which(ir_variants$variant%in%bmi_ns$RsIdA),]
# 
# bmi_neg <- fread("output/5_enrichment_analyses/3_go_shifter/ld_bmi_neg.txt")
# bmi_neg <- ir_variants[which(ir_variants$variant%in%bmi_neg$RsIdA),]
# 
# bmi_pos <- fread("output/5_enrichment_analyses/3_go_shifter/ld_bmi_pos.txt")
# bmi_pos <- ir_variants[which(ir_variants$variant%in%bmi_pos$RsIdA),]
# 
# bmi_outlier <- fread("output/5_enrichment_analyses/3_go_shifter/ld_outliers.txt")
# bmi_outlier <- ir_variants[which(ir_variants$variant%in%bmi_outlier$RsIdA),]
# 
# ###########################
# #Let's try adding the data#
# ###########################
# 
# working_df$reported <- NA
# working_df$source <- NA
# working_df$cluster <- NA
# 
# for(index_proxy in seq(1, length(working_df$query_snp_rsid))){
#   
#   #Variant:
#   
#   rsid <- working_df$query_snp_rsid[index_proxy]
#   
#   ir_tmp <- ir_variants[which(ir_variants$variant == rsid),]
#   
#   working_df$reported[index_proxy] <- ir_tmp$reported_ir
#   working_df$source[index_proxy] <- ir_tmp$ir_source
#   
#   if(rsid%in%bmi_ns$variant){
#     
#     working_df$cluster[index_proxy] <- "BMI (P>0.05)"
#     
#   }
#   
#   if(rsid%in%bmi_neg$variant){
#     
#     working_df$cluster[index_proxy] <- "BMI- (P<0.05)"
#     
#   }
#   
#   if(rsid%in%bmi_pos$variant){
#     
#     working_df$cluster[index_proxy] <- "BMI+ (P<0.05)"
#     
#   }
#   
#   if(rsid%in%bmi_outlier$variant){
#     
#     working_df$cluster[index_proxy] <- "BMI+;GSATadjBMI-"
#     
#   }
#   
#   
# }
# 
# #############################################
# #Let's now include those that have Hi-C data#
# #############################################
# 
# working_bmi_neg <- working_df[which(working_df$cluster == "BMI- (P<0.05)"),]
# working_bmi_ns <- working_df[which(working_df$cluster == "BMI (P>0.05)"),]
