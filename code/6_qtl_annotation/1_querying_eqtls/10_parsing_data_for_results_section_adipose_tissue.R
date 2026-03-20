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
qtl_df <- readRDS("output/6_qtl_annotation/conditional_and_fine_mapped_asat_vat_liver_and_muscle_eqtls_sqtls.txt") 

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


######################################################################################
#STEP 1: let's see how many exclusive variant-gene associations we have for ASAT QTLs#
######################################################################################

asat_qtl_df <- qtl_df[which(qtl_df$gene_name.asat != "" | is.na(qtl_df$gene_name.asat_gtex) == FALSE | is.na(qtl_df$gene_name.asat_sqtl) == FALSE),]

#The best way to all variant-gene combos is doing a loop and making things sensibly like the following:

asat_qtl_df$variant_gene_links <- NA

for(index_proxy in seq(1, length(asat_qtl_df$rsID))){
  
  tmp_df <- asat_qtl_df[index_proxy,]
  
  #STEP 1: extract the genes for this particular case
  
  genes <- c(tmp_df$gene_id.asat, tmp_df$gene_id.asat_gtex, tmp_df$gene_id.asat_sqtl)
  genes <- as.character(unlist(str_split(genes, ";")))
  
  genes <- unique(genes)
  genes <- genes[which(is.na(genes) == FALSE)]
  genes <- genes[which(genes != "")]
  
  #STEP 2: make the rsID-gene link unique:
  
  variant_gene <- paste(tmp_df$query_snp_rsid, "_", genes, sep = "")
  
  asat_qtl_df$variant_gene_links[index_proxy] <- paste(variant_gene, collapse = ";") 
  
}

#Now let's get all unique variant-gene links

all_variant_gene_links_in_asat <- as.character(unlist(str_split(asat_qtl_df$variant_gene_links, ";")))
all_variant_gene_links_in_asat <- unique(all_variant_gene_links_in_asat) 
print(length(all_variant_gene_links_in_asat)) #187

#187 Among how many leads???

length(unique(asat_qtl_df$query_snp_rsid)) #109

###############################################################################
#Great, now we should say how many of these are from each BMI-specific cluster#
###############################################################################

novel_variants <- ir_variants[which(ir_variants$reported_ir == "no"),]
novel_variants <- novel_variants[which(novel_variants$ir_source != "/Oliveri_gw"),] #70, just in case

#Now let's add the data for each subset of variants:

#Let's load only the fat distribution variants:

bmi_ns <- fread("output/5_enrichment_analyses/3_go_shifter/ld_bmi_ns.txt")
bmi_ns <- ir_variants[which(ir_variants$variant%in%bmi_ns$RsIdA),]

bmi_neg <- fread("output/5_enrichment_analyses/3_go_shifter/ld_bmi_neg.txt")
bmi_neg <- ir_variants[which(ir_variants$variant%in%bmi_neg$RsIdA),]

bmi_pos <- fread("output/5_enrichment_analyses/3_go_shifter/ld_bmi_pos.txt")
bmi_pos <- ir_variants[which(ir_variants$variant%in%bmi_pos$RsIdA),]

bmi_outlier <- fread("output/5_enrichment_analyses/3_go_shifter/ld_outliers.txt")
bmi_outlier <- ir_variants[which(ir_variants$variant%in%bmi_outlier$RsIdA),]

#And now let's split the data:

bmi_ns_match <- asat_qtl_df[which(asat_qtl_df$query_snp_rsid%in%bmi_ns$variant),] 
bmi_neg_match <- asat_qtl_df[which(asat_qtl_df$query_snp_rsid%in%bmi_neg$variant),] 
bmi_pos_match <- asat_qtl_df[which(asat_qtl_df$query_snp_rsid%in%bmi_pos$variant),] 
bmi_outlier_match <- asat_qtl_df[which(asat_qtl_df$query_snp_rsid%in%bmi_outlier$variant),]

#And finally, let's get the specific variants:

#BMI-NS or GSAT-specific lipodystrophy-like effects:

bmi_ns_gene_links_in_asat <- as.character(unlist(str_split(bmi_ns_match$variant_gene_links, ";")))
bmi_ns_gene_links_in_asat <- unique(bmi_ns_gene_links_in_asat) #111 variant-gene links
length(unique(bmi_ns_match$query_snp_rsid)) #56 IR leads
length(which(unique(bmi_ns_match$query_snp_rsid)%in%novel_variants$variant)) #15 novel

#BMI-NEG or generalized lipodystrophy-like effects:

bmi_neg_gene_links_in_asat <- as.character(unlist(str_split(bmi_neg_match$variant_gene_links, ";")))
bmi_neg_gene_links_in_asat <- unique(bmi_neg_gene_links_in_asat) #68 variant-gene links
length(unique(bmi_neg_match$query_snp_rsid)) #28 IR leads
length(which(unique(bmi_neg_match$query_snp_rsid)%in%novel_variants$variant)) #4 novel

#BMI-POS or generalized lipodystrophy-like effects:

bmi_pos_gene_links_in_asat <- as.character(unlist(str_split(bmi_pos_match$variant_gene_links, ";")))
bmi_pos_gene_links_in_asat <- unique(bmi_pos_gene_links_in_asat) #40 variant-gene links
length(unique(bmi_pos_match$query_snp_rsid)) #17 IR leads
length(which(unique(bmi_pos_match$query_snp_rsid)%in%novel_variants$variant)) #6 novel

#Outliers:

bmi_outlier_gene_links_in_asat <- as.character(unlist(str_split(bmi_outlier_match$variant_gene_links, ";")))
bmi_outlier_gene_links_in_asat <- unique(bmi_outlier_gene_links_in_asat) #40 variant-gene links
length(unique(bmi_outlier_match$query_snp_rsid)) #8 IR leads
length(which(unique(bmi_outlier_match$query_snp_rsid)%in%novel_variants$variant)) #2 novel

##################################################################################################
#Alright!! Next section is all about VAT, let's go and assess the data compared to adipose tissue#
##################################################################################################

vat_qtl_df <- qtl_df[which(is.na(qtl_df$gene_name.vat_gtex) == FALSE | is.na(qtl_df$gene_name.vat_sqtl) == FALSE),]

#The best way to all variant-gene combos is doing a loop and making things sensibly like the following:

vat_qtl_df$variant_gene_links <- NA

for(index_proxy in seq(1, length(vat_qtl_df$rsID))){
  
  tmp_df <- vat_qtl_df[index_proxy,]
  
  #STEP 1: extract the genes for this particular case
  
  genes <- c(tmp_df$gene_id.vat_gtex, tmp_df$gene_id.vat_sqtl)
  genes <- as.character(unlist(str_split(genes, ";")))
  genes <- unique(genes)
  genes <- genes[which(is.na(genes) == FALSE)]
  
  #STEP 2: make the rsID-gene link unique:
  
  variant_gene <- paste(tmp_df$query_snp_rsid, "_", genes, sep = "")
  
  vat_qtl_df$variant_gene_links[index_proxy] <- paste(variant_gene, collapse = ";") 
  
}

#And finally, let's extract them all and assess:

all_variant_gene_links_in_vat <- as.character(unlist(str_split(vat_qtl_df$variant_gene_links, ";")))
all_variant_gene_links_in_vat <- unique(all_variant_gene_links_in_vat) #91
print(length(all_variant_gene_links_in_vat)) #91

#122 Among how many leads???

length(unique(vat_qtl_df$query_snp_rsid)) #62

#Now!!! How many of these have already been reported in ASAT!!

shared_vat_variant_gene_links <- all_variant_gene_links_in_vat[which(all_variant_gene_links_in_vat%in%all_variant_gene_links_in_asat)] #64
print(length(unique(shared_vat_variant_gene_links)))

#Let's get noe those that are exclusive...

exclusive_vat_variant_gene_links <- all_variant_gene_links_in_vat[which(!(all_variant_gene_links_in_vat%in%all_variant_gene_links_in_asat))] #27
print(length(unique(exclusive_vat_variant_gene_links)))

#Be really careful with how we handle the dataframe:

vat_qtl_df <- separate_rows(vat_qtl_df, variant_gene_links, sep=";")

exclusive_vat_qtl_df <- vat_qtl_df[which(vat_qtl_df$variant_gene_links%in%exclusive_vat_variant_gene_links),]

length(unique(exclusive_vat_qtl_df$query_snp_rsid)) #25

###########################################################################################################
#Let's define then, properly the SAT-exclusive data and check the novel hits and the groups they belong to#
###########################################################################################################

#How many exclusive links?

exclusive_asat_variant_gene_links <- all_variant_gene_links_in_asat[which(!(all_variant_gene_links_in_asat%in%all_variant_gene_links_in_vat))]

print(length(unique(exclusive_asat_variant_gene_links))) #123

#Let's process the dataframe so that we can work better...

asat_qtl_df <- separate_rows(asat_qtl_df, variant_gene_links, sep=";")

exclusive_asat_qtl_df <- asat_qtl_df[which(asat_qtl_df$variant_gene_links%in%exclusive_asat_variant_gene_links),]

length(unique(exclusive_asat_qtl_df$query_snp_rsid)) #86

#Let's focus on the novel ones:

novel_exclusive_asat <- exclusive_asat_qtl_df[which(exclusive_asat_qtl_df$query_snp_rsid%in%novel_variants$variant),]

#Now let's split this 

novel_exclusive_asat_bmi_ns <- novel_exclusive_asat[which(novel_exclusive_asat$query_snp_rsid%in%bmi_ns$variant),]
novel_exclusive_asat_bmi_neg <- novel_exclusive_asat[which(novel_exclusive_asat$query_snp_rsid%in%bmi_neg$variant),]
novel_exclusive_asat_bmi_pos <- novel_exclusive_asat[which(novel_exclusive_asat$query_snp_rsid%in%bmi_pos$variant),]
novel_exclusive_asat_bmi_outlier <- novel_exclusive_asat[which(novel_exclusive_asat$query_snp_rsid%in%bmi_outlier$variant),]

#Let's get the data for BMI-independent

novel_bmi_ns <- unique(novel_exclusive_asat_bmi_ns$variant_gene_links)
novel_bmi_ns_genes <- unlist(str_split(novel_bmi_ns, "_"))
novel_bmi_ns_genes <- novel_bmi_ns_genes[which(str_detect(novel_bmi_ns_genes, "ENSG"))]

for(gene in novel_bmi_ns_genes){
  
  print(gene)
  
  tmp_df <- tryCatch(
    otargen::geneInfo(gene),
    error = function(e) {
      "NA"
    })
    
  if(length(tmp_df) == 1){
    
    next()
    
  }
  
  print(tmp_df$symbol)
  
} #9  in total, one could not be found, but in the df we show it is SACS-AS1

#Now, let's do it for BMI-neg

novel_bmi_neg <- unique(novel_exclusive_asat_bmi_neg$variant_gene_links)
novel_bmi_neg_genes <- unlist(str_split(novel_bmi_neg, "_"))
novel_bmi_neg_genes <- novel_bmi_neg_genes[which(str_detect(novel_bmi_neg_genes, "ENSG"))]

for(gene in novel_bmi_neg_genes){
  
  print(gene)
  
  tmp_df <- tryCatch(
    otargen::geneInfo(gene),
    error = function(e) {
      "NA"
    })
  
  if(length(tmp_df) == 1){
    
    next()
    
  }
  
  print(tmp_df$symbol)
  
} #3 in total...

#Now for BMI-increasing

novel_bmi_pos <- unique(novel_exclusive_asat_bmi_pos$variant_gene_links)
novel_bmi_pos_genes <- unlist(str_split(novel_bmi_pos, "_"))
novel_bmi_pos_genes <- novel_bmi_pos_genes[which(str_detect(novel_bmi_pos_genes, "ENSG"))]

for(gene in novel_bmi_pos_genes){
  
  print(gene)
  
  tmp_df <- tryCatch(
    otargen::geneInfo(gene),
    error = function(e) {
      "NA"
    })
  
  if(length(tmp_df) == 1){
    
    next()
    
  }
  
  print(tmp_df$symbol)
  
} #3 in total...

#Finally, outliers:

novel_bmi_outlier <- unique(novel_exclusive_asat_bmi_outlier$variant_gene_links)
novel_bmi_outlier <- unlist(str_split(novel_bmi_outlier, "_"))
novel_bmi_outlier <- novel_bmi_outlier[which(str_detect(novel_bmi_outlier, "ENSG"))]

for(gene in novel_bmi_outlier){
  
  print(gene)
  
  tmp_df <- tryCatch(
    otargen::geneInfo(gene),
    error = function(e) {
      "NA"
    })
  
  if(length(tmp_df) == 1){
    
    next()
    
  }
  
  print(tmp_df$symbol)
  
} #3 in total...

