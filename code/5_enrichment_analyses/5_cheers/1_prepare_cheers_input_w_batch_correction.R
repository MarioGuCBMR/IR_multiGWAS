##############
#INTRODUCTION#
##############

#This code will prepare the input data for CHEERS.

###################
#Load the packages#
###################

library(tidyverse)
library(data.table)
library(DESeq2)
library(limma)
library(data.table)

###################
#Loading functions#
###################

adding_b37 <- function(proxies_){
  #This code takes the haploreg output and just changes the column names so that they fit the CHEERs format.
  
  colnames(proxies_) <- c("rsID", "chr_37", "pos_hg38", "query_snp_rsid", "alt", "ref", "EUR", "r2", "bp_37")
  
  return(proxies_)
  
}

formatting_data_4_cheers <- function(data_df){
  
  #STEP 0: be careful we might have weird rsID;
  
  data_df <- data_df[which(data_df$rsID != ""),]
  
  #STEP 1: let's add the chromosome and base_pair_location to all the proxies:
  
  data_df_b37 <- adding_b37_otg(data_df)

  #STEP 2: let's clean this data from NAs
  
  data_df_b37_clean <- data_df_b37[which(is.na(data_df_b37$bp_37) == FALSE),]
  
  #sTEP 3: reformat the data:
  
  data_df_b37_clean <- data_df_b37_clean %>% dplyr::select(rsID, chr_37, bp_37)
  data_df_b37_clean$chr_37=paste0("chr",data_df_b37_clean$chr_37)
  
  #STEP 4: return the cleaned dataframe:
  
  return(data_df_b37_clean)
  
}

######################
#STEP 1: Loading data#
######################

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/")

# Feature counts input for each replicate (sample):
peaks <- read_tsv('raw_data/atac_seq_data/GSE178794_mergedConsensusPeaks_forDifferentialTesting_SGBS-days0-2-4-14.bed')
load('raw_data/atac_seq_data/GSE178794_DESeq2Datasets_forGEO_ATAC_SGBS-days0-2-4-14.RData')

##############################################################
#STEP 2: Controlling for BATCH effect + reformatting the data#
##############################################################

#https://pmc.ncbi.nlm.nih.gov/articles/PMC7518324/
#https://github.com/zhangyuqing/ComBat-seq 
#devtools::install_github("zhangyuqing/sva-devel")

library(sva)
#In ComBat-seq, user may specify biological covariates, whose signals will be preserved in the adjusted data. If the user would like to specify one biological variable, they may use the group parameter:
adjusted_counts <- ComBat_seq(assay(peaks.gc.dds), batch=peaks.gc.dds$Batch, group=peaks.gc.dds$condition)
assay(peaks.gc.dds) <- adjusted_counts

# Get the counts for each sample (omit Day2) in the required format and write a file

dir.create("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/output/5_enrichment_analyses/5_cheers")
dir.create("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/output/5_enrichment_analyses/5_cheers/1_cheers_input")
dir.create("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/output/5_enrichment_analyses/5_cheers/1_cheers_input/peak_data")
dir.create("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/output/5_enrichment_analyses/5_cheers/1_cheers_input/peak_data/combat_formatted_normalized_data")


setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/output/5_enrichment_analyses/5_cheers/1_cheers_input/peak_data/")


SGBS_day4_batch2_rep1 = as.data.frame(assay(peaks.gc.dds)[,"SGBS_day4_batch2_rep1"]) %>% dplyr::rename("counts"="assay(peaks.gc.dds)[, \"SGBS_day4_batch2_rep1\"]") %>% rownames_to_column("peak_ID")  %>% merge(x = ., y = peaks, by = "peak_ID") %>% dplyr::select("#chr", "start", "stop", "counts") %>% write.table(., file = "SGBS_day4_batch2_rep1_combat_corrected.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

SGBS_day4_batch2_rep2 = as.data.frame(assay(peaks.gc.dds)[,"SGBS_day4_batch2_rep2"]) %>% dplyr::rename("counts"="assay(peaks.gc.dds)[, \"SGBS_day4_batch2_rep2\"]") %>% rownames_to_column("peak_ID")  %>% merge(x = ., y = peaks, by = "peak_ID") %>% dplyr::select("#chr", "start", "stop", "counts") %>% write.table(., file = "SGBS_day4_batch2_rep2_combat_corrected.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

SGBS_day4_batch2_rep3 = as.data.frame(assay(peaks.gc.dds)[,"SGBS_day4_batch2_rep3"]) %>% dplyr::rename("counts"="assay(peaks.gc.dds)[, \"SGBS_day4_batch2_rep3\"]") %>% rownames_to_column("peak_ID")  %>% merge(x = ., y = peaks, by = "peak_ID") %>% dplyr::select("#chr", "start", "stop", "counts") %>% write.table(., file = "SGBS_day4_batch2_rep3_combat_corrected.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

SGBS_day4_batch2_rep4 = as.data.frame(assay(peaks.gc.dds)[,"SGBS_day4_batch2_rep4"]) %>% dplyr::rename("counts"="assay(peaks.gc.dds)[, \"SGBS_day4_batch2_rep4\"]") %>% rownames_to_column("peak_ID")  %>% merge(x = ., y = peaks, by = "peak_ID") %>% dplyr::select("#chr", "start", "stop", "counts") %>% write.table(., file = "SGBS_day4_batch2_rep4_combat_corrected.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

SGBS_day4_batch2_rep5 =  as.data.frame(assay(peaks.gc.dds)[,"SGBS_day4_batch2_rep5"]) %>% dplyr::rename("counts"="assay(peaks.gc.dds)[, \"SGBS_day4_batch2_rep5\"]") %>% rownames_to_column("peak_ID")  %>% merge(x = ., y = peaks, by = "peak_ID") %>% dplyr::select("#chr", "start", "stop", "counts") %>% write.table(., file = "SGBS_day4_batch2_rep5_combat_corrected.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

SGBS_day4_batch2_rep6 = as.data.frame(assay(peaks.gc.dds)[,"SGBS_day4_batch2_rep6"]) %>% dplyr::rename("counts"="assay(peaks.gc.dds)[, \"SGBS_day4_batch2_rep6\"]") %>% rownames_to_column("peak_ID")  %>% merge(x = ., y = peaks, by = "peak_ID") %>% dplyr::select("#chr", "start", "stop", "counts") %>% write.table(., file = "SGBS_day4_batch2_rep6_combat_corrected.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

SGBS_day0_batch2_rep1 = as.data.frame(assay(peaks.gc.dds)[,"SGBS_day0_batch2_rep1"]) %>% dplyr::rename("counts"="assay(peaks.gc.dds)[, \"SGBS_day0_batch2_rep1\"]") %>% rownames_to_column("peak_ID")  %>% merge(x = ., y = peaks, by = "peak_ID") %>% dplyr::select("#chr", "start", "stop", "counts") %>% write.table(., file = "SGBS_day0_batch2_rep1_combat_corrected.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

SGBS_day0_batch2_rep2 = as.data.frame(assay(peaks.gc.dds)[,"SGBS_day0_batch2_rep2"]) %>% dplyr::rename("counts"="assay(peaks.gc.dds)[, \"SGBS_day0_batch2_rep2\"]") %>% rownames_to_column("peak_ID")  %>% merge(x = ., y = peaks, by = "peak_ID") %>% dplyr::select("#chr", "start", "stop", "counts") %>% write.table(., file = "SGBS_day0_batch2_rep2_combat_corrected.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

SGBS_day0_batch2_rep3 = as.data.frame(assay(peaks.gc.dds)[,"SGBS_day0_batch2_rep3"]) %>% dplyr::rename("counts"="assay(peaks.gc.dds)[, \"SGBS_day0_batch2_rep3\"]") %>% rownames_to_column("peak_ID")  %>% merge(x = ., y = peaks, by = "peak_ID") %>% dplyr::select("#chr", "start", "stop", "counts") %>% write.table(., file = "SGBS_day0_batch2_rep3_combat_corrected.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

SGBS_day0_batch2_rep4 = as.data.frame(assay(peaks.gc.dds)[,"SGBS_day0_batch2_rep4"]) %>% dplyr::rename("counts"="assay(peaks.gc.dds)[, \"SGBS_day0_batch2_rep4\"]") %>% rownames_to_column("peak_ID")  %>% merge(x = ., y = peaks, by = "peak_ID") %>% dplyr::select("#chr", "start", "stop", "counts") %>% write.table(., file = "SGBS_day0_batch2_rep4_combat_corrected.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

SGBS_day0_batch2_rep5 = as.data.frame(assay(peaks.gc.dds)[,"SGBS_day0_batch2_rep5"]) %>% dplyr::rename("counts"="assay(peaks.gc.dds)[, \"SGBS_day0_batch2_rep5\"]") %>% rownames_to_column("peak_ID")  %>% merge(x = ., y = peaks, by = "peak_ID") %>% dplyr::select("#chr", "start", "stop", "counts") %>% write.table(., file = "SGBS_day0_batch2_rep5_combat_corrected.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

SGBS_day0_batch2_rep6 = as.data.frame(assay(peaks.gc.dds)[,"SGBS_day0_batch2_rep6"]) %>% dplyr::rename("counts"="assay(peaks.gc.dds)[, \"SGBS_day0_batch2_rep6\"]") %>% rownames_to_column("peak_ID")  %>% merge(x = ., y = peaks, by = "peak_ID") %>% dplyr::select("#chr", "start", "stop", "counts") %>% write.table(., file = "SGBS_day0_batch2_rep6_combat_corrected.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

SGBS_day14_batch1_rep1 = as.data.frame(assay(peaks.gc.dds)[,"SGBS_day14_batch1_rep1"]) %>% dplyr::rename("counts"="assay(peaks.gc.dds)[, \"SGBS_day14_batch1_rep1\"]") %>% rownames_to_column("peak_ID")  %>% merge(x = ., y = peaks, by = "peak_ID") %>% dplyr::select("#chr", "start", "stop", "counts") %>% write.table(., file = "SGBS_day14_batch1_rep1_combat_corrected.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

SGBS_day14_batch1_rep2 = as.data.frame(assay(peaks.gc.dds)[,"SGBS_day14_batch1_rep2"]) %>% dplyr::rename("counts"="assay(peaks.gc.dds)[, \"SGBS_day14_batch1_rep2\"]") %>% rownames_to_column("peak_ID")  %>% merge(x = ., y = peaks, by = "peak_ID") %>% dplyr::select("#chr", "start", "stop", "counts") %>% write.table(., file = "SGBS_day14_batch1_rep2_combat_corrected.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

SGBS_day14_batch1_rep3 = as.data.frame(assay(peaks.gc.dds)[,"SGBS_day14_batch1_rep3"]) %>% dplyr::rename("counts"="assay(peaks.gc.dds)[, \"SGBS_day14_batch1_rep3\"]") %>% rownames_to_column("peak_ID")  %>% merge(x = ., y = peaks, by = "peak_ID") %>% dplyr::select("#chr", "start", "stop", "counts") %>% write.table(., file = "SGBS_day14_batch1_rep3_combat_corrected.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

SGBS_day0_batch1_rep1 = as.data.frame(assay(peaks.gc.dds)[,"SGBS_day0_batch1_rep1"]) %>% dplyr::rename("counts"="assay(peaks.gc.dds)[, \"SGBS_day0_batch1_rep1\"]") %>% rownames_to_column("peak_ID")  %>% merge(x = ., y = peaks, by = "peak_ID") %>% dplyr::select("#chr", "start", "stop", "counts") %>% write.table(., file = "SGBS_day0_batch1_rep1_combat_corrected.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

SGBS_day0_batch1_rep2 = as.data.frame(assay(peaks.gc.dds)[,"SGBS_day0_batch1_rep2"]) %>% dplyr::rename("counts"="assay(peaks.gc.dds)[, \"SGBS_day0_batch1_rep2\"]") %>% rownames_to_column("peak_ID")  %>% merge(x = ., y = peaks, by = "peak_ID") %>% dplyr::select("#chr", "start", "stop", "counts") %>% write.table(., file = "SGBS_day0_batch1_rep2_combat_corrected.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

SGBS_day4_batch3_rep1 = as.data.frame(assay(peaks.gc.dds)[,"SGBS_day4_batch3_rep1"]) %>% dplyr::rename("counts"="assay(peaks.gc.dds)[, \"SGBS_day4_batch3_rep1\"]") %>% rownames_to_column("peak_ID")  %>% merge(x = ., y = peaks, by = "peak_ID") %>% dplyr::select("#chr", "start", "stop", "counts") %>% write.table(., file = "SGBS_day4_batch3_rep1_combat_corrected.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

SGBS_day4_batch3_rep2 = as.data.frame(assay(peaks.gc.dds)[,"SGBS_day4_batch3_rep2"]) %>% dplyr::rename("counts"="assay(peaks.gc.dds)[, \"SGBS_day4_batch3_rep2\"]") %>% rownames_to_column("peak_ID")  %>% merge(x = ., y = peaks, by = "peak_ID") %>% dplyr::select("#chr", "start", "stop", "counts") %>% write.table(., file = "SGBS_day4_batch3_rep2_combat_corrected.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

SGBS_day14_batch3_rep1 = as.data.frame(assay(peaks.gc.dds)[,"SGBS_day14_batch3_rep1"]) %>% dplyr::rename("counts"="assay(peaks.gc.dds)[, \"SGBS_day14_batch3_rep1\"]") %>% rownames_to_column("peak_ID")  %>% merge(x = ., y = peaks, by = "peak_ID") %>% dplyr::select("#chr", "start", "stop", "counts") %>% write.table(., file = "SGBS_day14_batch3_rep1_combat_corrected.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

SGBS_day14_batch3_rep2 = as.data.frame(assay(peaks.gc.dds)[,"SGBS_day14_batch3_rep2"]) %>% dplyr::rename("counts"="assay(peaks.gc.dds)[, \"SGBS_day14_batch3_rep2\"]") %>% rownames_to_column("peak_ID")  %>% merge(x = ., y = peaks, by = "peak_ID") %>% dplyr::select("#chr", "start", "stop", "counts") %>% write.table(., file = "SGBS_day14_batch3_rep2_combat_corrected.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

SGBS_day0_batch4_rep1 = as.data.frame(assay(peaks.gc.dds)[,"SGBS_day0_batch4_rep1"]) %>% dplyr::rename("counts"="assay(peaks.gc.dds)[, \"SGBS_day0_batch4_rep1\"]") %>% rownames_to_column("peak_ID")  %>% merge(x = ., y = peaks, by = "peak_ID") %>% dplyr::select("#chr", "start", "stop", "counts") %>% write.table(., file = "SGBS_day0_batch4_rep1_combat_corrected.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

SGBS_day0_batch4_rep2 = as.data.frame(assay(peaks.gc.dds)[,"SGBS_day0_batch4_rep2"]) %>% dplyr::rename("counts"="assay(peaks.gc.dds)[, \"SGBS_day0_batch4_rep2\"]") %>% rownames_to_column("peak_ID")  %>% merge(x = ., y = peaks, by = "peak_ID") %>% dplyr::select("#chr", "start", "stop", "counts") %>% write.table(., file = "SGBS_day0_batch4_rep2_combat_corrected.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

SGBS_day4_batch4_rep1 = as.data.frame(assay(peaks.gc.dds)[,"SGBS_day4_batch4_rep1"]) %>% dplyr::rename("counts"="assay(peaks.gc.dds)[, \"SGBS_day4_batch4_rep1\"]") %>% rownames_to_column("peak_ID")  %>% merge(x = ., y = peaks, by = "peak_ID") %>% dplyr::select("#chr", "start", "stop", "counts") %>% write.table(., file = "SGBS_day4_batch4_rep1_combat_corrected.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

SGBS_day4_batch4_rep2 =  as.data.frame(assay(peaks.gc.dds)[,"SGBS_day4_batch4_rep2"]) %>% dplyr::rename("counts"="assay(peaks.gc.dds)[, \"SGBS_day4_batch4_rep2\"]") %>% rownames_to_column("peak_ID")  %>% merge(x = ., y = peaks, by = "peak_ID") %>% dplyr::select("#chr", "start", "stop", "counts") %>% write.table(., file = "SGBS_day4_batch4_rep2_combat_corrected.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

###########################################################
#STEP 3: let's load all proxies for all data and save them#
###########################################################

# Proxies and FM structure: ID chr position
proxies <- fread("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/output/4_ir_loci_discovery/4_proxies/proxies_build_37_and_38_4_282_variants.txt")

##############################################################################################
#Let's get the lead SNPs for each of the clusters, so that we can run CHEERS for each of them#
##############################################################################################

#STEP 1: we can load the data for each cluster utilized in go-shifter_

bmi_ns <- fread("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/output/5_enrichment_analyses/3_go_shifter/ld_bmi_ns.txt")
bmi_neg <- fread("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/output/5_enrichment_analyses/3_go_shifter/ld_bmi_neg.txt")
bmi_pos <-  fread("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/output/5_enrichment_analyses/3_go_shifter/ld_bmi_pos.txt")
all_variants <- fread("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/output/5_enrichment_analyses/3_go_shifter/ld.txt")

#STEP 2: let's get proxies for each of variants in the clusters:

proxies_bmi_ns <- proxies[which(proxies$query_snp_rsid%in%bmi_ns$RsIdA),]
proxies_bmi_neg <- proxies[which(proxies$query_snp_rsid%in%bmi_neg$RsIdA),]
proxies_bmi_pos <- proxies[which(!(proxies$query_snp_rsid%in%bmi_neg$RsIdA | proxies$query_snp_rsid%in%bmi_ns$RsIdA)),]
proxies_all <- proxies[which(proxies$query_snp_rsid%in%all_variants$RsIdA),] #kind of unnecessary, but great sensitivity test

#Let's do a quick check:

length(unique(proxies_bmi_ns$query_snp_rsid)) #141
length(unique(proxies_bmi_neg$query_snp_rsid)) #63
length(unique(proxies_bmi_pos$query_snp_rsid)) #78
length(unique(proxies_all$query_snp_rsid)) #82/82 

#STEP 3: we used to add chr_37 and bp_37 here, but with our new approach..., we already have them! Just need to have the right column names

proxies_bmi_ns_b37 <- adding_b37(proxies_bmi_ns)
proxies_bmi_neg_b37 <- adding_b37(proxies_bmi_neg)
proxies_bmi_pos_b37 <- adding_b37(proxies_bmi_pos)
proxies_all_b37 <- adding_b37(proxies_all)

#STEP 4: let's clean this data:

proxies_bmi_ns_b37_clean <- proxies_bmi_ns_b37[which(is.na(proxies_bmi_ns_b37$bp_37) == FALSE),]
proxies_bmi_neg_b37_clean <- proxies_bmi_neg_b37[which(is.na(proxies_bmi_neg_b37$bp_37) == FALSE),]
proxies_bmi_pos_b37_clean <- proxies_bmi_pos_b37[which(is.na(proxies_bmi_pos_b37$bp_37) == FALSE),]
proxies_all_b37_clean <- proxies_all_b37[which(is.na(proxies_all_b37$bp_37) == FALSE),]

#Now let's see if we have all leads, if we miss any we put it manually, only for the leads though

length(unique(proxies_bmi_ns_b37_clean$query_snp_rsid)) #141
length(unique(proxies_bmi_neg_b37_clean$query_snp_rsid)) #63
length(unique(proxies_bmi_pos_b37_clean$query_snp_rsid)) #78
length(unique(proxies_all_b37_clean$query_snp_rsid)) #282

###########################################
#Let's format the data for all the proxies#
###########################################

#CLUSTER 1:

proxies_bmi_ns_b37_clean <- proxies_bmi_ns_b37_clean %>% dplyr::select(rsID, chr_37, bp_37)
proxies_bmi_ns_b37_clean$chr_37=paste0("chr",proxies_bmi_ns_b37_clean$chr_37)

dir.create("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/output/5_enrichment_analyses/5_cheers/1_cheers_input/variant_data")

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/output/5_enrichment_analyses/5_cheers/1_cheers_input/variant_data/")

write.table(proxies_bmi_ns_b37_clean, file = "proxies_bmi_ns.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

#CLUSTER 2:

proxies_bmi_neg_b37_clean <- proxies_bmi_neg_b37_clean %>% dplyr::select(rsID, chr_37, bp_37)
proxies_bmi_neg_b37_clean$chr_37=paste0("chr",proxies_bmi_neg_b37_clean$chr_37)

write.table(proxies_bmi_neg_b37_clean, file = "proxies_bmi_neg.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

#CLUSTER 3:

proxies_bmi_pos_b37_clean <- proxies_bmi_pos_b37_clean %>% dplyr::select(rsID, chr_37, bp_37)
proxies_bmi_pos_b37_clean$chr_37=paste0("chr",proxies_bmi_pos_b37_clean$chr_37)

write.table(proxies_bmi_pos_b37_clean, file = "proxies_bmi_pos.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

#Now all results:

proxies_all_b37_clean <- proxies_all_b37_clean %>% dplyr::select(rsID, chr_37, bp_37)
proxies_all_b37_clean$chr_37=paste0("chr",proxies_all_b37_clean$chr_37)

write.table(proxies_all_b37_clean, file = "proxies_full.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")