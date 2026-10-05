##############
#INTRODUCTION#
##############

#This code reads the genes reported in the Supplementary Tables and compares the hits across the 4 IR tissues

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)
library(TwoSampleMR)
library(hyprcoloc)
library(GenomicRanges)
library(rtracklayer)

###################
#Loading functions#
###################

gene_parser = function(gene_data){
  
  #STEP 1: take the string and split it:
  
  vect_=unlist(str_split(gene_data, "[.]"))
  
  id_=vect_[1]
  
  return(id_)
  
}


chr_parser = function(variant_id){
  
  #STEP 1: take the string and split it:
  
  vect_=unlist(str_split(variant_id, "_"))
  
  chr_=vect_[1]
  
  return(chr_)
  
}

pos_parser = function(variant_id){
  
  #STEP 1: take the string and split it:
  
  vect_=unlist(str_split(variant_id, "_"))
  
  pos_=vect_[2]
  
  return(pos_)
  
}

ref_parser = function(variant_id){
  
  #STEP 1: take the string and split it:
  
  vect_=unlist(str_split(variant_id, "_"))
  
  ref_=vect_[3]
  
  return(ref_)
  
}

alt_parser = function(variant_id){
  
  #STEP 1: take the string and split it:
  
  vect_=unlist(str_split(variant_id, "_"))
  
  alt_=vect_[4]
  
  return(alt_)
  
}


trait_aligner <- function(query_ss, other_ss){
  
  #We are gonna put the example so that we can know what happens:
  
  #whradjbmi_ss <- exp_df_found_pos
  #other_ss <- hdl_005
  #query_ss <- triangulated_fi_all
  #other_ss <- ebmd_ss
  
  #STEP 0: let's run the matching with TwoSampleMR so we need the data:
  
  #Let's first check if we have the effect_allele_frequency column:
  
  check <- which(colnames(query_ss) == "effect_allele_frequency")
  
  if(is_empty(check)){
    
    query_ss$effect_allele_frequency <- NA
    
  }
  
  #Now we can proceed
  
  exposure <- query_ss %>%
    select(chr_pos, chromosome, base_pair_location, effect_allele, other_allele,
           effect_allele_frequency, beta, standard_error, p_value, sample_size, variant)
  
  colnames(exposure) <- c("SNP", "chr.exposure", "pos.exposure", "effect_allele.exposure",
                          "other_allele.exposure", "eaf.exposure", "beta.exposure", "se.exposure",
                          "pval.exposure", "samplesize.exposure", "rsid.exposure")
  
  exposure$exposure <- "fat_distr"
  exposure$id.exposure <- "fat_distr"
  
  #Now with the outcome:
  
  check <- which(colnames(other_ss) == "effect_allele_frequency")
  
  if(is_empty(check)){
    
    other_ss$effect_allele_frequency <- NA
    
  }
  
  #Now we can proceed
  
  outcome <- other_ss %>%
    select(chr_pos, chromosome, base_pair_location, effect_allele, other_allele, effect_allele_frequency, beta, standard_error, p_value, sample_size, variant)
  
  colnames(outcome) <- c("SNP", "chr.outcome", "pos.outcome", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome", "pval.outcome", "samplesize.outcome", "rsid.outcome")
  
  outcome$outcome <- "Outcome"
  outcome$id.outcome <- "Outcome"
  
  #STEP 1: match the data. This will probably fail with tri-allelic SNPs. Here I think we are gonna be OK...
  
  merged_df <- harmonise_data(exposure, outcome, action=3)
  
  merged_df <- merged_df[which(merged_df$remove == FALSE),] #removing incompatible alleles
  #merged_df <- merged_df[which(merged_df$mr_keep),] #removing incompatible alleles
  
  #I checked that all is working great. whradjbmi betas is positive. The rest is great.
  
  #STEP 3: reorganize the dataframe, cuz we need something clean:
  
  other_ss_aligned <- merged_df %>%
    select("rsid.outcome", "chr.outcome", "pos.outcome", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome", "pval.outcome", "samplesize.outcome", "SNP")
  
  colnames(other_ss_aligned) <- c("variant", "chromosome", "base_pair_location",  "effect_allele", "other_allele", "effect_allele_frequency", "beta", "standard_error", "p_value", "sample_size", "chr_pos")
  
  return(other_ss_aligned)
  
}

organizing_common_data <- function(ss_df, common_snps){
  
  ss_end <- ss_df[which(ss_df$chr_pos%in%common_snps),]
  ss_ordered <- ss_end[order(as.numeric(ss_end$chromosome), as.numeric(ss_end$base_pair_location)),]
  
  return(ss_ordered)
  
}

duplicate_remover <- function(aligned_df){
  
  #Let's remove the duplicates, because with the way that the data is built there is no way to separate those...
  
  duplicated_snps=aligned_df$chr_pos[which(duplicated(aligned_df$chr_pos))]
  
  if(is_empty(duplicated_snps)){
    
    return(aligned_df)
    
  }
  
  aligned_df=aligned_df[which(!(aligned_df$chr_pos%in%duplicated_snps)),]
  
  return(aligned_df)
  
}

#############################
#STEP 0: Let's load the data#
#############################

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

qtls=data.table::fread("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/manuscript/supplementary_material/supplementary_tables/drafts/qtl_matches_03092026.csv")

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
ir_sat_coloc = ir_leads[which(is.na(ir_leads$high_confidence_genes) == FALSE),]

###################################################
#Let's do the comparisons for all the GTEx tissues#
###################################################

asat_genes = qtls %>%
  dplyr::select("Gene ID", "PIP (ASAT cis-eQTL)", "PIP (ASAT cis-sQTL)")

asat_genes = asat_genes[which(is.na(asat_genes$`PIP (ASAT cis-eQTL)`) == FALSE | 
                                is.na(asat_genes$`PIP (ASAT cis-sQTL)`) == FALSE),] #114

#Let's do the same with VAT.

vat_genes = qtls %>%
  dplyr::select("Gene ID", "PIP (VAT cis-eQTL)", "PIP (VAT cis-sQTL)")

vat_genes = vat_genes[which(is.na(vat_genes$`PIP (VAT cis-eQTL)`) == FALSE | 
                                is.na(vat_genes$`PIP (VAT cis-sQTL)`) == FALSE),] #94

#Now the muscle:

muscle_genes = qtls %>%
  dplyr::select("Gene ID", "PIP (Skeletal muscle cis-eQTL)", "PIP (Skeletal muscle cis-sQTL)")

muscle_genes = muscle_genes[which(is.na(muscle_genes$"PIP (Skeletal muscle cis-eQTL)") == FALSE | 
                              is.na(muscle_genes$"PIP (Skeletal muscle cis-sQTL)") == FALSE),] #73
#Now the liver:

liver_genes = qtls %>%
  dplyr::select("Gene ID", "PIP (Liver cis-eQTL)", "PIP (Liver cis-sQTL)")

liver_genes = liver_genes[which(is.na(liver_genes$"PIP (Liver cis-eQTL)") == FALSE | 
                                    is.na(liver_genes$"PIP (Liver cis-sQTL)") == FALSE),] #28

#Careful, due to how we edited the data, we have duplicates due to proxies!!
#Now that we have the hits, we can just take the data as a vect

asat_genes_vect = unique(asat_genes$`Gene ID`) #107
vat_genes_vect = unique(vat_genes$`Gene ID`) #91
muscle_genes_vect = unique(muscle_genes$`Gene ID`) #69
liver_genes_vect = unique(liver_genes$`Gene ID`) #27

################################################################
#Let's get only those that do not replicate in any other tissue#
################################################################

#First we are gonna do a SAT vs VAT comparison:

sat_vs_vat = asat_genes_vect[which(!(asat_genes_vect%in%vat_genes_vect))] #37

vat_vs_sat = vat_genes_vect[which(!(vat_genes_vect%in%asat_genes_vect))] #26

muscle_only = muscle_genes_vect[which(!(muscle_genes_vect%in%vat_genes_vect) &
                                    !(muscle_genes_vect%in%asat_genes_vect) &
                                    !(muscle_genes_vect%in%liver_genes_vect))] #27

liver_only = liver_genes_vect[which(!(liver_genes_vect%in%vat_genes_vect) &
                                        !(liver_genes_vect%in%asat_genes_vect) &
                                        !(liver_genes_vect%in%muscle_genes_vect))] #10

#Let's check those that are found in all:

all_hits = Reduce(intersect, list(asat_genes_vect, vat_genes_vect, muscle_genes_vect, liver_genes_vect)) #8

###############################################################################
#We might want to do this for the co-localization, but maybe in another moment#
###############################################################################