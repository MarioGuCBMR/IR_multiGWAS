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
ir_variants = ir_variants[which(ir_variants$ir_direction == "yes"),]
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

#STEP 1: remove those that are not present...

asat_qtl_df <- qtl_df[which(qtl_df$gene_name.asat != "" | is.na(qtl_df$gene_name.asat_gtex) == FALSE | is.na(qtl_df$gene_name.asat_sqtl) == FALSE),]

#The best way to find all variant-gene combos is doing a loop and combining lead SNP per gene

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

#Let's take this change to split rows for one second..., which will come in handy afterwards

asat_qtl_df <- separate_rows(asat_qtl_df, variant_gene_links, sep=";")

########################################################################################################################
#STEP 2: let's read the data as we have it now and try to understand why we have a missing link and where it comes from#
########################################################################################################################

#Let's get the lead SNPs info just in case:

all_variants=readxl::read_xlsx("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/review/manuscript/03062026/Supplementary Tables_20March2026.xlsx", sheet = 3)
all_variants = as.data.frame(all_variants)
all_variants = all_variants[2:nrow(all_variants),] #remove extra first row that explains stuff
colnames(all_variants) = all_variants[1,]  #change column names
all_variants = all_variants[2:nrow(all_variants),] #remove columns names in the first line now that we have used them

#We do not need to many columns AND we have functions ready for a certain subset of variables.

all_variants = all_variants[,1:12]
colnames(all_variants) = c("variant", "chromosome", "base_pair_location", "minimum_allele_frequency", "effect_allele", "other_allele", "novel_or_reported", "bmi_subgroup", "VEP", "beta", "standard_error", "pvalue")
all_variants$sample_size = 2e06 #an approximation, not needed for the actual analyses, just to make the functions work

#PERFECT!!
#Now let's move to get the proxies for each SNP:

proxies <- fread("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/output/4_ir_loci_discovery/4_proxies/proxies_build_37_and_38_4_282_variants.txt")

yes_vect <- c("A", "G", "C", "T")

proxies <- proxies[which(proxies$query_snp_rsid%in%all_variants$variant),] #3188
proxies <- proxies[which(proxies$ref%in%yes_vect & proxies$alt%in%yes_vect),]

length(which(duplicated(proxies$query_snp_rsid) == FALSE)) #282

#PERFECT! Now we need the full sumstats, of course:

iradjbmi = fread("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/output/2_mv_gwas/fiadjbmi_hdl_tg_common_dwls_curated.txt")

#Last thing: the list of genes:

qtls=readxl::read_xlsx("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/review/manuscript/03062026/Supplementary Tables_20March2026.xlsx", sheet = 12)

#Let's pass this through the same filters:

qtls = as.data.frame(qtls)
qtls = qtls[2:nrow(qtls),] #remove extra first row that explains stuff
colnames(qtls) = qtls[1,]  #change column names
qtls = qtls[2:nrow(qtls),] #remove columns names in the first line now that we have used them

#The columns here did not work as good as I thought due to the table's structure. 
#No matter, we are gonna try changing some of the columns:

colnames(qtls)[1:11] = c("Lead IR SNP", "Variant", "R2", "Chromosome", "Position (b37)", "Position (b38)", "Effect Allele Frequency", "IR allele", "IS allele", "Gene ID", "Gene Symbol")

#Last thing - let's get the co-localizing high-confidence data:

coloc=readxl::read_xlsx("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/review/manuscript/10082026/Extended_Data_11August2026.xlsx", sheet = 1)

#Let's clean this info a bit to...

coloc=as.data.frame(coloc)
coloc = coloc[2:nrow(coloc),] #remove extra first row that explains stuff
colnames(coloc) = coloc[1,]  #change column names
coloc = coloc[2:nrow(coloc),] #remove columns names in the first line now that we have used them

#The names got a bit shambled, but we can work with this just fine.

colnames(coloc)[1:12] #these are normal cols
colnames(coloc)[13:20] = paste(colnames(coloc)[13:20], ".sat_marginal", sep = "")
colnames(coloc)[21:28] = paste(colnames(coloc)[21:28], ".sat_conditional", sep = "")
colnames(coloc)[29:36] = paste(colnames(coloc)[29:36], ".sat_gtex", sep = "")
colnames(coloc)[37:44] = paste(colnames(coloc)[37:44], ".sat_sqtl", sep = "")
colnames(coloc)[45:52] = paste(colnames(coloc)[45:52], ".vat_gtex", sep = "")
colnames(coloc)[53:60] = paste(colnames(coloc)[53:60], ".vat_sqtl", sep = "")
colnames(coloc)[61:68] = paste(colnames(coloc)[61:68], ".liver_gtex", sep = "")
colnames(coloc)[69:76] = paste(colnames(coloc)[69:76], ".liver_sqtl", sep = "")
colnames(coloc)[77:84] = paste(colnames(coloc)[77:84], ".muscle_gtex", sep = "")
colnames(coloc)[85:92] = paste(colnames(coloc)[85:92], ".muscle_sqtl", sep = "")
colnames(coloc)[93:99] = paste(colnames(coloc)[93:99], ".pqtl", sep = "")

#############################################################
#Let's get the genes for all analyses in AdipoXpress tissues#
#############################################################

#Let's get the high confidence data:

asat_genes_all = qtls %>%
  dplyr::select("Variant", "Gene ID", "P-value (ASAT cis-eQTL)", "PIP (ASAT cis-eQTL)", "PIP (ASAT cis-sQTL)")

asat_genes_all = asat_genes_all[which(is.na(asat_genes_all$`PIP (ASAT cis-eQTL)`) == FALSE | 
                                        is.na(asat_genes_all$`PIP (ASAT cis-sQTL)`) == FALSE |
                                        is.na(asat_genes_all$`P-value (ASAT cis-eQTL)`) == FALSE),] #233

#Great, let's get the high-confidence genes list too: 

sat_coloc_genes = coloc[which(is.na(coloc$`High-confidence genes.sat_marginal`) ==FALSE | 
                                is.na(coloc$`High-confidence genes.sat_conditional`) ==FALSE | 
                                is.na(coloc$`High-confidence genes.sat_gtex`) ==FALSE |
                                is.na(coloc$`High-confidence genes.sat_sqtl`) ==FALSE),]

#We are gonna use the dataframe above to match stuff, but for the time being let's make a vector with all different genes. 

coloc_genes_vect = c(unlist(str_split(sat_coloc_genes$`High-confidence genes.sat_marginal`, ";")),
                     unlist(str_split(sat_coloc_genes$`High-confidence genes.sat_conditional`, ";")),
                     unlist(str_split(sat_coloc_genes$`High-confidence genes.sat_gtex`, ";")),
                     unlist(str_split(sat_coloc_genes$`High-confidence genes.sat_sqtl`, ";")))

coloc_genes_vect = unique(coloc_genes_vect[which(is.na(coloc_genes_vect) == FALSE)]) #175

#OK, let's get the loci linked to these by compacting all the possible genes. 

ir_leads =as.data.frame(all_variants[,1])
colnames(ir_leads) = "leads"

ir_leads$genes=NA
ir_leads$high_confidence_genes=NA
ir_leads$gene_counts=0
ir_leads$high_confidence_gene_counts=0

#And now do the classic proxy search to append the data

for(index in seq(1, length(ir_leads$leads))){
  
  #STEP 1: get the proxies for this leads
  
  proxies_tmp = proxies[which(proxies$query_snp_rsid == ir_leads$leads[index]),]
  
  #STEP 2: lets' get the genes linked to this gene:
  
  asat_genes_all_tmp = asat_genes_all[which(asat_genes_all$Variant%in%proxies_tmp$rsID),]
  
  if(is_empty(asat_genes_all_tmp$Variant)){
    
    next()
    
  }
  
  print(index)
  
  #Let's get the data.
  
  ir_leads$genes[index] = paste(asat_genes_all_tmp$`Gene ID`, collapse=";")
  ir_leads$gene_counts[index] = nrow(as.data.frame(asat_genes_all_tmp))
  
  #STEP 3: we need symbols for our genes, since that is the info we have for the high-confidence hits:
  
  ensdb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
  
  symbol_map <- AnnotationDbi::mapIds(
    ensdb,
    keys      = asat_genes_all_tmp$`Gene ID`,
    column    = "SYMBOL",
    keytype   = "GENEID",
    multiVals = "first"
  )
  
  #STEP 4: let's see how many of these genes are in the coloc_list:
  
  coloc_genes_vect_match = coloc_genes_vect[which(coloc_genes_vect%in%unlist(symbol_map) | coloc_genes_vect%in%names(symbol_map))] #this allows to capture the ones that did not go through!
  
  if(is_empty(coloc_genes_vect_match)){
    
    next()
    
  }  
  
  ir_leads$high_confidence_genes[index] = paste(coloc_genes_vect_match, collapse=";")
  ir_leads$high_confidence_gene_counts[index] = length(coloc_genes_vect_match)
  
}

#With this info we can get the number of loci that were linked to SAT Vs the co-localized.
#The number of genes is a bit more complex, since some of them are repeated due to be secondary signals. 
#Let's get those numbers:

ir_sat=ir_leads[which(is.na(ir_leads$genes) == FALSE),] #107 - we found before 109. Maybe because of the way we treated ENSGs? We need to check

###################################################################################################
#We really are finding different genes in both versions - let's figure out which are the real ones#
###################################################################################################

length(which(unique(asat_qtl_df$query_snp_rsid)%in%ir_sat$leads)) #107/109 - correct. 

ir_missing = asat_qtl_df[which(!(asat_qtl_df$query_snp_rsid%in%ir_sat$leads)),]


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

#This was great to separate stuff, but let's focus... on adding info to the main dataframe:

ir_variants$bmi_cluster <- NA
ir_variants$bmi_cluster <- ifelse(ir_variants$variant%in%bmi_ns$variant, "BMI-independent", ir_variants$bmi_cluster)
ir_variants$bmi_cluster <- ifelse(ir_variants$variant%in%bmi_neg$variant, "BMI-decreasing", ir_variants$bmi_cluster)
ir_variants$bmi_cluster <- ifelse(ir_variants$variant%in%bmi_pos$variant, "BMI-increasing", ir_variants$bmi_cluster)
ir_variants$bmi_cluster <- ifelse(ir_variants$variant%in%bmi_outlier$variant, "BMI-outlier", ir_variants$bmi_cluster)

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

#Q1: how many variant-gene links?

print(length(all_variant_gene_links_in_vat)) #91

#Q2: Among how many leads???

length(unique(vat_qtl_df$query_snp_rsid)) #62

#Now!!! How many of these have already been reported in ASAT!!

shared_vat_variant_gene_links <- all_variant_gene_links_in_vat[which(all_variant_gene_links_in_vat%in%all_variant_gene_links_in_asat)] #64

#How many are shared?

print(length(unique(shared_vat_variant_gene_links))) #64

#Let's get now those that are exclusive...

exclusive_vat_variant_gene_links <- all_variant_gene_links_in_vat[which(!(all_variant_gene_links_in_vat%in%all_variant_gene_links_in_asat))] #27
print(length(unique(exclusive_vat_variant_gene_links))) #27

#Be really careful with how we handle the dataframe:

vat_qtl_df <- separate_rows(vat_qtl_df, variant_gene_links, sep=";")

###############################
#Let's do the same for liver!!#
###############################

liver_qtl_df <- qtl_df[which(is.na(qtl_df$gene_name.liver_gtex) == FALSE | is.na(qtl_df$gene_name.liver_sqtl) == FALSE),]

#The best way to all variant-gene combos is doing a loop and making things sensibly like the following:

liver_qtl_df$variant_gene_links <- NA

for(index_proxy in seq(1, length(liver_qtl_df$rsID))){
  
  tmp_df <- liver_qtl_df[index_proxy,]
  
  #STEP 1: extract the genes for this particular case
  
  genes <- c(tmp_df$gene_id.liver_gtex, tmp_df$gene_id.liver_sqtl)
  genes <- as.character(unlist(str_split(genes, ";")))
  genes <- unique(genes)
  genes <- genes[which(is.na(genes) == FALSE)]
  
  #STEP 2: make the rsID-gene link unique:
  
  variant_gene <- paste(tmp_df$query_snp_rsid, "_", genes, sep = "")
  
  liver_qtl_df$variant_gene_links[index_proxy] <- paste(variant_gene, collapse = ";") 
  
}

#And finally, let's extract them all and assess:

all_variant_gene_links_in_liver <- as.character(unlist(str_split(liver_qtl_df$variant_gene_links, ";")))
all_variant_gene_links_in_liver <- unique(all_variant_gene_links_in_liver) #27

#How many variant-gene links?

print(length(all_variant_gene_links_in_liver)) #27
 
#In how many many leads???

length(unique(liver_qtl_df$query_snp_rsid)) #22

#Now!!! How many of these have already been reported in ASAT/VAT?

shared_liver_variant_gene_links_asat_only <- all_variant_gene_links_in_liver[which(all_variant_gene_links_in_liver%in%all_variant_gene_links_in_asat)] #82

length(unique(shared_liver_variant_gene_links_asat_only)) #16

#Let's separate the hits...

liver_qtl_df <- separate_rows(liver_qtl_df, variant_gene_links, sep=";")

################################
#Let's do the same for muscle!!#
################################

muscle_qtl_df <- qtl_df[which(is.na(qtl_df$gene_name.muscle_gtex) == FALSE | is.na(qtl_df$gene_name.muscle_sqtl) == FALSE),]

#The best way to all variant-gene combos is doing a loop and making things sensibly like the following:

muscle_qtl_df$variant_gene_links <- NA

for(index_proxy in seq(1, length(muscle_qtl_df$rsID))){
  
  tmp_df <- muscle_qtl_df[index_proxy,]
  
  #STEP 1: extract the genes for this particular case
  
  genes <- c(tmp_df$gene_id.muscle_gtex, tmp_df$gene_id.muscle_sqtl)
  genes <- as.character(unlist(str_split(genes, ";")))
  genes <- unique(genes)
  genes <- genes[which(is.na(genes) == FALSE)]
  
  #STEP 2: make the rsID-gene link unique:
  
  variant_gene <- paste(tmp_df$query_snp_rsid, "_", genes, sep = "")
  
  muscle_qtl_df$variant_gene_links[index_proxy] <- paste(variant_gene, collapse = ";") 
  
}

#And finally, let's extract them all and assess:

all_variant_gene_links_in_muscle <- as.character(unlist(str_split(muscle_qtl_df$variant_gene_links, ";")))
all_variant_gene_links_in_muscle <- unique(all_variant_gene_links_in_muscle) #27
length(all_variant_gene_links_in_muscle) #69

#122 Among how many leads???

length(unique(muscle_qtl_df$query_snp_rsid)) #53

#Now!!! How many of these have already been reported in ASAT!!

shared_muscle_variant_gene_links_asat <- all_variant_gene_links_in_muscle[which(all_variant_gene_links_in_muscle%in%all_variant_gene_links_in_asat)] #82

length(unique(shared_muscle_variant_gene_links_asat)) #38

#Let's separate the hits...
muscle_qtl_df <- separate_rows(muscle_qtl_df, variant_gene_links, sep=";")

###################################
#Let's obtain tissue specific hits#
###################################

#First, let's get all exclusive hits for each tissue:

exclusive_asat_variant_gene_links <- all_variant_gene_links_in_asat[which(!(all_variant_gene_links_in_asat%in%all_variant_gene_links_in_muscle) &
                                                                                !(all_variant_gene_links_in_asat%in%all_variant_gene_links_in_vat) &
                                                                                !(all_variant_gene_links_in_asat%in%all_variant_gene_links_in_liver))] 

print(length(unique(exclusive_asat_variant_gene_links))) #108

exclusive_vat_variant_gene_links <- all_variant_gene_links_in_vat[which(!(all_variant_gene_links_in_vat%in%all_variant_gene_links_in_asat) &
                                                                                !(all_variant_gene_links_in_vat%in%all_variant_gene_links_in_muscle) &
                                                                                !(all_variant_gene_links_in_vat%in%all_variant_gene_links_in_liver))] 

print(length(unique(exclusive_vat_variant_gene_links))) #22

exclusive_liver_variant_gene_links <- all_variant_gene_links_in_liver[which(!(all_variant_gene_links_in_liver%in%all_variant_gene_links_in_asat) &
                                                                                !(all_variant_gene_links_in_liver%in%all_variant_gene_links_in_vat) &
                                                                                !(all_variant_gene_links_in_liver%in%all_variant_gene_links_in_muscle))] 


print(length(unique(exclusive_liver_variant_gene_links))) #7


exclusive_muscle_variant_gene_links <- all_variant_gene_links_in_muscle[which(!(all_variant_gene_links_in_muscle%in%all_variant_gene_links_in_asat) &
                                                                                !(all_variant_gene_links_in_muscle%in%all_variant_gene_links_in_vat) &
                                                                                !(all_variant_gene_links_in_muscle%in%all_variant_gene_links_in_liver))] 

print(length(unique(exclusive_muscle_variant_gene_links))) #26

#############################################################
#Let's make dataframes for these - and obtain the novel ones#
#############################################################

exclusive_asat_qtl_df <- asat_qtl_df[which(asat_qtl_df$variant_gene_links%in%exclusive_asat_variant_gene_links),]
exclusive_vat_qtl_df <- vat_qtl_df[which(vat_qtl_df$variant_gene_links%in%exclusive_vat_variant_gene_links),]
exclusive_liver_qtl_df <- liver_qtl_df[which(liver_qtl_df$variant_gene_links%in%exclusive_liver_variant_gene_links),]
exclusive_muscle_qtl_df <- muscle_qtl_df[which(muscle_qtl_df$variant_gene_links%in%exclusive_muscle_variant_gene_links),]

#Let's get the novel ones:

novel_excl_asat <- exclusive_asat_qtl_df[which(exclusive_asat_qtl_df$query_snp_rsid%in%novel_variants$variant),]
print(length(unique(novel_excl_asat$query_snp_rsid))) #15
print(length(unique(novel_excl_asat$variant_gene_links))) #15

novel_excl_vat <- exclusive_vat_qtl_df[which(exclusive_vat_qtl_df$query_snp_rsid%in%novel_variants$variant),]
print(length(unique(novel_excl_vat$query_snp_rsid))) #4
print(length(unique(novel_excl_vat$variant_gene_links))) #4

novel_excl_liver <- exclusive_liver_qtl_df[which(exclusive_liver_qtl_df$query_snp_rsid%in%novel_variants$variant),]
print(length(unique(novel_excl_liver$query_snp_rsid))) #2
print(length(unique(novel_excl_liver$variant_gene_links))) #2

novel_excl_muscle <- exclusive_muscle_qtl_df[which(exclusive_muscle_qtl_df$query_snp_rsid%in%novel_variants$variant),]
print(length(unique(novel_excl_muscle$query_snp_rsid))) #7
print(length(unique(novel_excl_muscle$variant_gene_links))) #8

################################################################
#Let's format the data so that we can have a sensible dataframe#
################################################################

novel_variant_gene_links_asat <- unique(novel_excl_asat$variant_gene_links)

#Let's split asat_qtl_df first...

asat_qtl_adipoxpress <- asat_qtl_df[which(asat_qtl_df$gene_id.asat != ""),]
asat_qtl_adipoxpress <- separate_rows(asat_qtl_adipoxpress, gene_id.asat, gene_name.asat, beta.asat, p_value.asat, sep=";")

#Now for gtex QTL:

asat_qtl_gtex_eqtl <- asat_qtl_df[which(asat_qtl_df$gene_id.asat_gtex != ""),]
asat_qtl_gtex_eqtl <- separate_rows(asat_qtl_gtex_eqtl, gene_id.asat_gtex, gene_name.asat_gtex, afc.asat_gtex, pip.asat_gtex, sep=";")

asat_qtl_gtex_sqtl <- asat_qtl_df[which(asat_qtl_df$gene_id.asat_sqtl != ""),]
asat_qtl_gtex_sqtl <- separate_rows(asat_qtl_gtex_sqtl, gene_id.asat_sqtl, gene_name.asat_sqtl, splicing_effect, pip.asat_sqtl, sep=";")

#Let' add an analysis column so that it is easier:

asat_qtl_adipoxpress$analysis <- "AdipoXpress cis-eQTL"
asat_qtl_gtex_eqtl$analysis <- "GTEx V10 cis-eQTL"
asat_qtl_gtex_sqtl$analysis <- "GTEx V10 cis-sQTL"

for(index in seq(1, length(novel_variant_gene_links_asat))){
  
  #STEP 0: get the variant-gene link
  
  v2g_tmp <- novel_variant_gene_links_asat[index]
  
  #STEP 1: get the rsID and gene:
  
  rsid <- as.character(unlist(str_split(v2g_tmp, "_"))[1])
  gene <- as.character(unlist(str_split(v2g_tmp, "_"))[2])
  
  #STEP 2: get info on the variant 
  
  variant_info <- ir_variants[which(ir_variants$variant == rsid),]
  
  #STEP 3: get QTL data
  
  adipoxpress <- asat_qtl_adipoxpress[which(asat_qtl_adipoxpress$query_snp_rsid == rsid & asat_qtl_adipoxpress$gene_id.asat == gene),]
  gtex_eqtl <- asat_qtl_gtex_eqtl[which(asat_qtl_gtex_eqtl$query_snp_rsid == rsid & asat_qtl_gtex_eqtl$gene_id.asat_gtex == gene),]
  gtex_sqtl <- asat_qtl_gtex_sqtl[which(asat_qtl_gtex_sqtl$query_snp_rsid == rsid & asat_qtl_gtex_sqtl$gene_id.asat_sqtl == gene),]
  
  #Let's get the best info of all dataframes:
  
  yes_vect <- c("A", "G", "C", "T")
  
  adipoxpress <- adipoxpress[which(adipoxpress$ir_allele.asat%in%yes_vect),]
  adipoxpress <- adipoxpress[which(duplicated(adipoxpress$rsID) == FALSE),] #to get just one
  adipoxpress <- adipoxpress[which.min(adipoxpress$p_value.asat),] #to get just one
  
  gtex_eqtl <- gtex_eqtl[which(gtex_eqtl$ir_allele.asat_gtex%in%yes_vect),]
  gtex_eqtl <- gtex_eqtl[which(duplicated(gtex_eqtl$rsID) == FALSE),] 
  gtex_eqtl <- gtex_eqtl[which.max(gtex_eqtl$pip.asat_gtex),] #to get just one
  
  gtex_sqtl <- gtex_sqtl[which(gtex_sqtl$ir_allele.asat_sqtl%in%yes_vect),]
  gtex_sqtl <- gtex_sqtl[which(duplicated(gtex_sqtl$rsID) == FALSE),] #to get just one
  gtex_sqtl <- gtex_sqtl[which.max(gtex_sqtl$pip.asat_sqtl),]
  
  #STEP 4: let's combine all information in one dataframe
  
  gene_name <- paste(adipoxpress$gene_name.asat, gtex_eqtl$gene_name.asat_gtex, gtex_sqtl$gene_name.asat_sqtl, sep = ";")
  ir_allele <- paste(adipoxpress$ir_allele.asat, gtex_eqtl$ir_allele.asat_gtex, gtex_sqtl$ir_allele.asat_sqtl, sep = ";")
  qtl_effect <- paste(adipoxpress$beta.asat, gtex_eqtl$afc.asat_gtex, gtex_sqtl$splicing_region, sep = ";")
  association <- paste(adipoxpress$p_value.asat, gtex_eqtl$pip.asat_gtex, gtex_sqtl$pip.asat_sqtl, sep = ";")
  analysis <- paste(adipoxpress$analysis, gtex_eqtl$analysis, gtex_sqtl$analysis, sep = ";")
  
  #STEP 5: let's make our dataframe:
  
  if(!(exists("asat_excl_novel_summary"))){
    
    info_df <- as.data.frame(t(c(variant_info$variant, variant_info$chr_pos, variant_info$bmi_cluster, gene_name, qtl_effect, association, analysis)))
    colnames(info_df) <- c("variant", "chr_pos", "cluster", "gene", "effect", "pval", "study")
    asat_excl_novel_summary <- info_df
    
  } else {
    
    info_df <- as.data.frame(t(c(variant_info$variant, variant_info$chr_pos, variant_info$bmi_cluster, gene_name, qtl_effect, association, analysis)))
    colnames(info_df) <- c("variant", "chr_pos", "cluster", "gene", "effect", "pval", "study")
    asat_excl_novel_summary <- rbind(asat_excl_novel_summary, info_df)
    
  }
  
}

###########################
#Let's do the same for VAT#
###########################

novel_variant_gene_links_vat <- unique(novel_excl_vat$variant_gene_links)

#Now for gtex QTL:

vat_qtl_gtex_eqtl <- vat_qtl_df[which(vat_qtl_df$gene_id.vat_gtex != ""),]
vat_qtl_gtex_eqtl <- separate_rows(vat_qtl_gtex_eqtl, gene_id.vat_gtex, gene_name.vat_gtex, afc.vat_gtex, pip.vat_gtex, sep=";")

vat_qtl_gtex_sqtl <- vat_qtl_df[which(vat_qtl_df$gene_id.vat_sqtl != ""),]
vat_qtl_gtex_sqtl <- separate_rows(vat_qtl_gtex_sqtl, gene_id.vat_sqtl, gene_name.vat_sqtl, splicing_effect.vat, pip.vat_sqtl, sep=";")

#Let' add an analysis column so that it is easier:

vat_qtl_gtex_eqtl$analysis <- "GTEx V10 cis-eQTL"
vat_qtl_gtex_sqtl$analysis <- "GTEx V10 cis-sQTL"

for(index in seq(1, length(novel_variant_gene_links_vat))){
  
  #STEP 0: get the variant-gene link
  
  v2g_tmp <- novel_variant_gene_links_vat[index]
  
  #STEP 1: get the rsID and gene:
  
  rsid <- as.character(unlist(str_split(v2g_tmp, "_"))[1])
  gene <- as.character(unlist(str_split(v2g_tmp, "_"))[2])
  
  #STEP 2: get info on the variant 
  
  variant_info <- ir_variants[which(ir_variants$variant == rsid),]
  
  #STEP 3: get QTL data
  
  #adipoxpress <- vat_qtl_adipoxpress[which(vat_qtl_adipoxpress$query_snp_rsid == rsid & vat_qtl_adipoxpress$gene_id.vat == gene),]
  gtex_eqtl <- vat_qtl_gtex_eqtl[which(vat_qtl_gtex_eqtl$query_snp_rsid == rsid & vat_qtl_gtex_eqtl$gene_id.vat_gtex == gene),]
  gtex_sqtl <- vat_qtl_gtex_sqtl[which(vat_qtl_gtex_sqtl$query_snp_rsid == rsid & vat_qtl_gtex_sqtl$gene_id.vat_sqtl == gene),]
  
  #Let's get the best info of all dataframes:
  
  yes_vect <- c("A", "G", "C", "T")
  
  # adipoxpress <- adipoxpress[which(adipoxpress$ir_allele.vat%in%yes_vect),]
  # adipoxpress <- adipoxpress[which(duplicated(adipoxpress$rsID) == FALSE),] #to get just one
  # adipoxpress <- adipoxpress[which.min(adipoxpress$p_value.vat),] #to get just one
  
  gtex_eqtl <- gtex_eqtl[which(gtex_eqtl$ir_allele.vat_gtex%in%yes_vect),]
  gtex_eqtl <- gtex_eqtl[which(duplicated(gtex_eqtl$rsID) == FALSE),] 
  gtex_eqtl <- gtex_eqtl[which.max(gtex_eqtl$pip.vat_gtex),] #to get just one
  
  gtex_sqtl <- gtex_sqtl[which(gtex_sqtl$ir_allele.vat_sqtl%in%yes_vect),]
  gtex_sqtl <- gtex_sqtl[which(duplicated(gtex_sqtl$rsID) == FALSE),] #to get just one
  gtex_sqtl <- gtex_sqtl[which.max(gtex_sqtl$pip.vat_sqtl),]
  
  #STEP 4: let's combine all information in one dataframe
  
  gene_name <- paste(gtex_eqtl$gene_name.vat_gtex, gtex_sqtl$gene_name.vat_sqtl, sep = ";")
  ir_allele <- paste(gtex_eqtl$ir_allele.vat_gtex, gtex_sqtl$ir_allele.vat_sqtl, sep = ";")
  qtl_effect <- paste(gtex_eqtl$afc.vat_gtex, gtex_sqtl$splicing_region.vat, sep = ";")
  association <- paste(gtex_eqtl$pip.vat_gtex, gtex_sqtl$pip.vat_sqtl, sep = ";")
  analysis <- paste(gtex_eqtl$analysis, gtex_sqtl$analysis, sep = ";")
  
  #STEP 5: let's make our dataframe:
  
  if(!(exists("vat_excl_novel_summary"))){
    
    info_df <- as.data.frame(t(c(variant_info$variant, variant_info$chr_pos, variant_info$bmi_cluster, gene_name, qtl_effect, association, analysis)))
    colnames(info_df) <- c("variant", "chr_pos", "cluster", "gene", "effect", "pval", "study")
    vat_excl_novel_summary <- info_df
    
  } else {
    
    info_df <- as.data.frame(t(c(variant_info$variant, variant_info$chr_pos, variant_info$bmi_cluster, gene_name, qtl_effect, association, analysis)))
    colnames(info_df) <- c("variant", "chr_pos", "cluster", "gene", "effect", "pval", "study")
    vat_excl_novel_summary <- rbind(vat_excl_novel_summary, info_df)
    
  }
  
}

#######################################
#Let's do the same for skeletal muscle#
#######################################

novel_variant_gene_links_muscle <- unique(novel_excl_muscle$variant_gene_links)

#Now for gtex QTL:

muscle_qtl_gtex_eqtl <- muscle_qtl_df[which(muscle_qtl_df$gene_id.muscle_gtex != ""),]
muscle_qtl_gtex_eqtl <- separate_rows(muscle_qtl_gtex_eqtl, gene_id.muscle_gtex, gene_name.muscle_gtex, afc.muscle_gtex, pip.muscle_gtex, sep=";")

muscle_qtl_gtex_sqtl <- muscle_qtl_df[which(muscle_qtl_df$gene_id.muscle_sqtl != ""),]
muscle_qtl_gtex_sqtl <- separate_rows(muscle_qtl_gtex_sqtl, gene_id.muscle_sqtl, gene_name.muscle_sqtl, splicing_effect.muscle, pip.muscle_sqtl, sep=";")

#Let' add an analysis column so that it is easier:

muscle_qtl_gtex_eqtl$analysis <- "GTEx V10 cis-eQTL"
muscle_qtl_gtex_sqtl$analysis <- "GTEx V10 cis-sQTL"

for(index in seq(1, length(novel_variant_gene_links_muscle))){
  
  #STEP 0: get the variant-gene link
  
  v2g_tmp <- novel_variant_gene_links_muscle[index]
  
  #STEP 1: get the rsID and gene:
  
  rsid <- as.character(unlist(str_split(v2g_tmp, "_"))[1])
  gene <- as.character(unlist(str_split(v2g_tmp, "_"))[2])
  
  #STEP 2: get info on the variant 
  
  variant_info <- ir_variants[which(ir_variants$variant == rsid),]
  
  #STEP 3: get QTL data
  
  #adipoxpress <- muscle_qtl_adipoxpress[which(muscle_qtl_adipoxpress$query_snp_rsid == rsid & muscle_qtl_adipoxpress$gene_id.muscle == gene),]
  gtex_eqtl <- muscle_qtl_gtex_eqtl[which(muscle_qtl_gtex_eqtl$query_snp_rsid == rsid & muscle_qtl_gtex_eqtl$gene_id.muscle_gtex == gene),]
  gtex_sqtl <- muscle_qtl_gtex_sqtl[which(muscle_qtl_gtex_sqtl$query_snp_rsid == rsid & muscle_qtl_gtex_sqtl$gene_id.muscle_sqtl == gene),]
  
  #Let's get the best info of all dataframes:
  
  yes_vect <- c("A", "G", "C", "T")
  
  # adipoxpress <- adipoxpress[which(adipoxpress$ir_allele.muscle%in%yes_vect),]
  # adipoxpress <- adipoxpress[which(duplicated(adipoxpress$rsID) == FALSE),] #to get just one
  # adipoxpress <- adipoxpress[which.min(adipoxpress$p_value.muscle),] #to get just one
  
  gtex_eqtl <- gtex_eqtl[which(gtex_eqtl$ir_allele.muscle_gtex%in%yes_vect),]
  gtex_eqtl <- gtex_eqtl[which(duplicated(gtex_eqtl$rsID) == FALSE),] 
  gtex_eqtl <- gtex_eqtl[which.max(gtex_eqtl$pip.muscle_gtex),] #to get just one
  
  gtex_sqtl <- gtex_sqtl[which(gtex_sqtl$ir_allele.muscle_sqtl%in%yes_vect),]
  gtex_sqtl <- gtex_sqtl[which(duplicated(gtex_sqtl$rsID) == FALSE),] #to get just one
  gtex_sqtl <- gtex_sqtl[which.max(gtex_sqtl$pip.muscle_sqtl),]
  
  #STEP 4: let's combine all information in one dataframe
  
  gene_name <- paste(gtex_eqtl$gene_name.muscle_gtex, gtex_sqtl$gene_name.muscle_sqtl, sep = ";")
  ir_allele <- paste(gtex_eqtl$ir_allele.muscle_gtex, gtex_sqtl$ir_allele.muscle_sqtl, sep = ";")
  qtl_effect <- paste(gtex_eqtl$afc.muscle_gtex, gtex_sqtl$splicing_region.muscle, sep = ";")
  association <- paste(gtex_eqtl$pip.muscle_gtex, gtex_sqtl$pip.muscle_sqtl, sep = ";")
  analysis <- paste(gtex_eqtl$analysis, gtex_sqtl$analysis, sep = ";")
  
  #STEP 5: let's make our dataframe:
  
  if(!(exists("muscle_excl_novel_summary"))){
    
    info_df <- as.data.frame(t(c(variant_info$variant, variant_info$chr_pos, variant_info$bmi_cluster, gene_name, qtl_effect, association, analysis)))
    colnames(info_df) <- c("variant", "chr_pos", "cluster", "gene", "effect", "pval", "study")
    muscle_excl_novel_summary <- info_df
    
  } else {
    
    info_df <- as.data.frame(t(c(variant_info$variant, variant_info$chr_pos, variant_info$bmi_cluster, gene_name, qtl_effect, association, analysis)))
    colnames(info_df) <- c("variant", "chr_pos", "cluster", "gene", "effect", "pval", "study")
    muscle_excl_novel_summary <- rbind(muscle_excl_novel_summary, info_df)
    
  }
  
}

#########################
#Finally, let's do liver#
#########################

novel_variant_gene_links_liver <- unique(novel_excl_liver$variant_gene_links)

#Now for gtex QTL:

liver_qtl_gtex_eqtl <- liver_qtl_df[which(liver_qtl_df$gene_id.liver_gtex != ""),]
liver_qtl_gtex_eqtl <- separate_rows(liver_qtl_gtex_eqtl, gene_id.liver_gtex, gene_name.liver_gtex, afc.liver_gtex, pip.liver_gtex, sep=";")

liver_qtl_gtex_sqtl <- liver_qtl_df[which(liver_qtl_df$gene_id.liver_sqtl != ""),]
liver_qtl_gtex_sqtl <- separate_rows(liver_qtl_gtex_sqtl, gene_id.liver_sqtl, gene_name.liver_sqtl, splicing_effect.liver, pip.liver_sqtl, sep=";")

#Let' add an analysis column so that it is easier:

liver_qtl_gtex_eqtl$analysis <- "GTEx V10 cis-eQTL"
liver_qtl_gtex_sqtl$analysis <- "GTEx V10 cis-sQTL"

for(index in seq(1, length(novel_variant_gene_links_liver))){
  
  #STEP 0: get the variant-gene link
  
  v2g_tmp <- novel_variant_gene_links_liver[index]
  
  #STEP 1: get the rsID and gene:
  
  rsid <- as.character(unlist(str_split(v2g_tmp, "_"))[1])
  gene <- as.character(unlist(str_split(v2g_tmp, "_"))[2])
  
  #STEP 2: get info on the variant 
  
  variant_info <- ir_variants[which(ir_variants$variant == rsid),]
  
  #STEP 3: get QTL data
  
  #adipoxpress <- liver_qtl_adipoxpress[which(liver_qtl_adipoxpress$query_snp_rsid == rsid & liver_qtl_adipoxpress$gene_id.liver == gene),]
  gtex_eqtl <- liver_qtl_gtex_eqtl[which(liver_qtl_gtex_eqtl$query_snp_rsid == rsid & liver_qtl_gtex_eqtl$gene_id.liver_gtex == gene),]
  gtex_sqtl <- liver_qtl_gtex_sqtl[which(liver_qtl_gtex_sqtl$query_snp_rsid == rsid & liver_qtl_gtex_sqtl$gene_id.liver_sqtl == gene),]
  
  #Let's get the best info of all dataframes:
  
  yes_vect <- c("A", "G", "C", "T")
  
  # adipoxpress <- adipoxpress[which(adipoxpress$ir_allele.liver%in%yes_vect),]
  # adipoxpress <- adipoxpress[which(duplicated(adipoxpress$rsID) == FALSE),] #to get just one
  # adipoxpress <- adipoxpress[which.min(adipoxpress$p_value.liver),] #to get just one
  
  gtex_eqtl <- gtex_eqtl[which(gtex_eqtl$ir_allele.liver_gtex%in%yes_vect),]
  gtex_eqtl <- gtex_eqtl[which(duplicated(gtex_eqtl$rsID) == FALSE),] 
  gtex_eqtl <- gtex_eqtl[which.max(gtex_eqtl$pip.liver_gtex),] #to get just one
  
  gtex_sqtl <- gtex_sqtl[which(gtex_sqtl$ir_allele.liver_sqtl%in%yes_vect),]
  gtex_sqtl <- gtex_sqtl[which(duplicated(gtex_sqtl$rsID) == FALSE),] #to get just one
  gtex_sqtl <- gtex_sqtl[which.max(gtex_sqtl$pip.liver_sqtl),]
  
  #STEP 4: let's combine all information in one dataframe
  
  gene_name <- paste(gtex_eqtl$gene_name.liver_gtex, gtex_sqtl$gene_name.liver_sqtl, sep = ";")
  ir_allele <- paste(gtex_eqtl$ir_allele.liver_gtex, gtex_sqtl$ir_allele.liver_sqtl, sep = ";")
  qtl_effect <- paste(gtex_eqtl$afc.liver_gtex, gtex_sqtl$splicing_region.liver, sep = ";")
  association <- paste(gtex_eqtl$pip.liver_gtex, gtex_sqtl$pip.liver_sqtl, sep = ";")
  analysis <- paste(gtex_eqtl$analysis, gtex_sqtl$analysis, sep = ";")
  
  #STEP 5: let's make our dataframe:
  
  if(!(exists("liver_excl_novel_summary"))){
    
    info_df <- as.data.frame(t(c(variant_info$variant, variant_info$chr_pos, variant_info$bmi_cluster, gene_name, qtl_effect, association, analysis)))
    colnames(info_df) <- c("variant", "chr_pos", "cluster", "gene", "effect", "pval", "study")
    liver_excl_novel_summary <- info_df
    
  } else {
    
    info_df <- as.data.frame(t(c(variant_info$variant, variant_info$chr_pos, variant_info$bmi_cluster, gene_name, qtl_effect, association, analysis)))
    colnames(info_df) <- c("variant", "chr_pos", "cluster", "gene", "effect", "pval", "study")
    liver_excl_novel_summary <- rbind(liver_excl_novel_summary, info_df)
    
  }
  
}

##########################################
#Let's merge the novel hits and go for it#
##########################################

asat_excl_novel_summary$tissue <- "ASAT"
vat_excl_novel_summary$tissue <- "VAT"
muscle_excl_novel_summary$tissue <- "Muscle"
liver_excl_novel_summary$tissue <- "Liver"

novel_hits_df <- rbind(asat_excl_novel_summary, vat_excl_novel_summary, muscle_excl_novel_summary, liver_excl_novel_summary)

fwrite(novel_hits_df, "output/6_qtl_annotation/summary_novel_hits.txt", sep = "\t")
