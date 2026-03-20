##############
#INTRODUCTION#
##############

#This code parses enhancer-promoter connections from STARE.

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)
library(otargen)

###################
#Loading functions#
###################

parsing_gene_id <- function(gene){
  
  gene_ <- as.character(unlist(str_split(gene, "[.]"))[1])
  
  return(gene_)
  
}

#################################################################################
#Let's load the IR variants and the dataframe we have with all SGBS data already#
#################################################################################

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025")

#Let's load the most updated datasets here:

ir_variants <- fread("output/4_ir_loci_discovery/3_novel_variants/4_colocalization/352_variants_with_all_info_plus_fiadjbmi_and_tg_hdl_coloc.txt")
proxies <- fread("output/7_functional_annotation/1_atac_seq_overlaps/proxies_matched_w_roadmap_and_perrin_data_w_sgbs_diff.txt")

#Now let's focus only on the 282 IR variants:

ir_variants <- ir_variants[which(ir_variants$variant%in%proxies$rsID),] #282! Perfect

#####################
#Let's load the data#
#####################

stare_e063 <- fread("output/7_functional_annotation/2_enhancer_gene_links/1_stare/output/roadmap_e063_ABCpp_scoredInteractions_value.txt.gz")
stare_day0 <- fread("output/7_functional_annotation/2_enhancer_gene_links/1_stare/output/sgbs_day_0_ABCpp_scoredInteractions_value.txt.gz")
stare_day4 <- fread("output/7_functional_annotation/2_enhancer_gene_links/1_stare/output/sgbs_day_4_ABCpp_scoredInteractions_value.txt.gz")
stare_day14 <- fread("output/7_functional_annotation/2_enhancer_gene_links/1_stare/output/sgbs_day_14_ABCpp_scoredInteractions_value.txt.gz")

stare_e063$gene_id <- as.character(unlist(sapply(stare_e063$`Ensembl ID`, parsing_gene_id)))
stare_day0$gene_id <- as.character(unlist(sapply(stare_day0$`Ensembl ID`, parsing_gene_id)))
stare_day4$gene_id <- as.character(unlist(sapply(stare_day4$`Ensembl ID`, parsing_gene_id)))
stare_day14$gene_id <- as.character(unlist(sapply(stare_day14$`Ensembl ID`, parsing_gene_id)))

#Let's add a quick annotation section to help us parse the data afterwards cuz it gets quite crazy...

stare_e063$annotation <- NA
stare_e063$annotation <- ifelse(as.numeric(stare_e063$`TSS-dist`) <= 250, "Promoter", stare_e063$annotation)
stare_e063$annotation <- ifelse(as.numeric(stare_e063$`TSS-dist`) <= 250, "Promoter", stare_e063$annotation)
stare_e063$annotation <- ifelse(as.numeric(stare_e063$`TSS-dist`) <= 1000 & as.numeric(stare_e063$`TSS-dist`) > 250, "Flanking promoter", stare_e063$annotation)
stare_e063$annotation <- ifelse(as.numeric(stare_e063$`TSS-dist`) <= 2000 & as.numeric(stare_e063$`TSS-dist`) > 1000, "Upstream gene", stare_e063$annotation)
stare_e063$annotation <- ifelse(as.numeric(stare_e063$`TSS-dist`) > 2000, "Distal", stare_e063$annotation)

stare_day0$annotation <- NA
stare_day0$annotation <- ifelse(as.numeric(stare_day0$`TSS-dist`) <= 250, "Promoter", stare_day0$annotation)
stare_day0$annotation <- ifelse(as.numeric(stare_day0$`TSS-dist`) <= 250, "Promoter", stare_day0$annotation)
stare_day0$annotation <- ifelse(as.numeric(stare_day0$`TSS-dist`) <= 1000 & as.numeric(stare_day0$`TSS-dist`) > 250, "Flanking promoter", stare_day0$annotation)
stare_day0$annotation <- ifelse(as.numeric(stare_day0$`TSS-dist`) <= 2000 & as.numeric(stare_day0$`TSS-dist`) > 1000, "Upstream gene", stare_day0$annotation)
stare_day0$annotation <- ifelse(as.numeric(stare_day0$`TSS-dist`) > 2000, "Distal", stare_day0$annotation)


stare_day4$annotation <- NA
stare_day4$annotation <- ifelse(as.numeric(stare_day4$`TSS-dist`) <= 250, "Promoter", stare_day4$annotation)
stare_day4$annotation <- ifelse(as.numeric(stare_day4$`TSS-dist`) <= 250, "Promoter", stare_day4$annotation)
stare_day4$annotation <- ifelse(as.numeric(stare_day4$`TSS-dist`) <= 1000 & as.numeric(stare_day4$`TSS-dist`) > 250, "Flanking promoter", stare_day4$annotation)
stare_day4$annotation <- ifelse(as.numeric(stare_day4$`TSS-dist`) <= 2000 & as.numeric(stare_day4$`TSS-dist`) > 1000, "Upstream gene", stare_day4$annotation)
stare_day4$annotation <- ifelse(as.numeric(stare_day4$`TSS-dist`) > 2000, "Distal", stare_day4$annotation)


stare_day14$annotation <- NA
stare_day14$annotation <- ifelse(as.numeric(stare_day14$`TSS-dist`) <= 250, "Promoter", stare_day14$annotation)
stare_day14$annotation <- ifelse(as.numeric(stare_day14$`TSS-dist`) <= 250, "Promoter", stare_day14$annotation)
stare_day14$annotation <- ifelse(as.numeric(stare_day14$`TSS-dist`) <= 1000 & as.numeric(stare_day14$`TSS-dist`) > 250, "Flanking promoter", stare_day14$annotation)
stare_day14$annotation <- ifelse(as.numeric(stare_day14$`TSS-dist`) <= 2000 & as.numeric(stare_day14$`TSS-dist`) > 1000, "Upstream gene", stare_day14$annotation)
stare_day14$annotation <- ifelse(as.numeric(stare_day14$`TSS-dist`) > 2000, "Distal", stare_day14$annotation)

#########################################
#Let's match the data for adipose tissue#
#########################################

proxies$stare_adipose_nuclei <- NA
proxies$tss_dist_adipose_nuclei<- NA
proxies$annotation_adipose_nuclei<- NA
proxies$gene_id_adipose_nuclei <- NA
proxies$gene_name_adipose_nuclei <- NA
proxies$score_adipose_nuclei <- NA

for(index_proxy in seq(1, length(proxies$rsID))){
  
  chr_ <- proxies$chr[index_proxy]
  pos_ <- proxies$pos_hg19[index_proxy]
  
  peak_tmp <- stare_e063[which(as.numeric(stare_e063$`#chr`) == as.numeric(chr_) & as.numeric(stare_e063$start) <= as.numeric(pos_) & as.numeric(stare_e063$end) >= as.numeric(pos_))]
  peak_tmp$peak_region <- paste("chr", peak_tmp$'#chr', ":", peak_tmp$start, "_", peak_tmp$end, sep ="")
  
  if(is_empty(peak_tmp$start)){
    
    next()
    
  } else {
    
    print(index_proxy)
    
    proxies$stare_adipose_nuclei[index_proxy] <- paste(unique(peak_tmp$peak_region), collapse =";")
    proxies$tss_dist_adipose_nuclei[index_proxy] <- paste(peak_tmp$'TSS-dist', collapse =";")
    proxies$annotation_adipose_nuclei[index_proxy] <- paste(peak_tmp$annotation, collapse =";")
    proxies$gene_id_adipose_nuclei[index_proxy] <- paste(peak_tmp$gene_id, collapse =";")
    proxies$gene_name_adipose_nuclei[index_proxy] <- paste(peak_tmp$'Gene Name', collapse =";")
    proxies$score_adipose_nuclei[index_proxy] <- paste(peak_tmp$`ABC-Score`, collapse =";")
    
  }
  
}

######################
#Let's match the data#
######################

proxies$stare_day_0 <- NA
proxies$tss_dist_day_0 <- NA
proxies$annotation_day_0<- NA
proxies$gene_id_day_0 <- NA
proxies$gene_name_day_0 <- NA
proxies$score_day_0 <- NA

for(index_proxy in seq(1, length(proxies$rsID))){
  
  chr_ <- proxies$chr[index_proxy]
  pos_ <- proxies$pos_hg19[index_proxy]
  
  peak_tmp <- stare_day0[which(as.numeric(stare_day0$`#chr`) == as.numeric(chr_) & as.numeric(stare_day0$start) <= as.numeric(pos_) & as.numeric(stare_day0$end) >= as.numeric(pos_))]
  peak_tmp$peak_region <- paste("chr", peak_tmp$'#chr', ":", peak_tmp$start, "_", peak_tmp$end, sep ="")
  
  if(is_empty(peak_tmp$start)){
    
    next()
    
  } else {
    
    print(index_proxy)
    
    proxies$stare_day_0[index_proxy] <- paste(unique(peak_tmp$peak_region), collapse =";")
    proxies$tss_dist_day_0[index_proxy] <- paste(peak_tmp$'TSS-dist', collapse =";")
    proxies$annotation_day_0[index_proxy] <- paste(peak_tmp$annotation, collapse =";")
    proxies$gene_id_day_0[index_proxy] <- paste(peak_tmp$gene_id, collapse =";")
    proxies$gene_name_day_0[index_proxy] <- paste(peak_tmp$'Gene Name', collapse =";")
    proxies$score_day_0[index_proxy] <- paste(peak_tmp$`ABC-Score`, collapse =";")

  }
  
}

################################
#Let's match the data for Day 4#
################################

proxies$stare_day_4 <- NA
proxies$tss_dist_day_4 <- NA
proxies$annotation_day_4<- NA
proxies$gene_id_day_4 <- NA
proxies$gene_name_day_4 <- NA
proxies$score_day_4 <- NA

for(index_proxy in seq(1, length(proxies$rsID))){
  
  chr_ <- proxies$chr[index_proxy]
  pos_ <- proxies$pos_hg19[index_proxy]
  
  peak_tmp <- stare_day4[which(as.numeric(stare_day4$`#chr`) == as.numeric(chr_) & as.numeric(stare_day4$start) <= as.numeric(pos_) & as.numeric(stare_day4$end) >= as.numeric(pos_))]
  peak_tmp$peak_region <- paste("chr", peak_tmp$'#chr', ":", peak_tmp$start, "_", peak_tmp$end, sep ="")
  
  if(is_empty(peak_tmp$start)){
    
    next()
    
  } else {
    
    print(index_proxy)
    
    proxies$stare_day_4[index_proxy] <- paste(unique(peak_tmp$peak_region), collapse =";")
    proxies$tss_dist_day_4[index_proxy] <- paste(peak_tmp$'TSS-dist', collapse =";")
    proxies$annotation_day_4[index_proxy] <- paste(peak_tmp$annotation, collapse =";")
    proxies$gene_id_day_4[index_proxy] <- paste(peak_tmp$gene_id, collapse =";")
    proxies$gene_name_day_4[index_proxy] <- paste(peak_tmp$'Gene Name', collapse =";")
    proxies$score_day_4[index_proxy] <- paste(peak_tmp$`ABC-Score`, collapse =";")
    
  }
  
}

################
#Also at day 14#
################

proxies$stare_day_14 <- NA
proxies$tss_dist_day_14 <- NA
proxies$annotation_day_14 <- NA
proxies$gene_id_day_14 <- NA
proxies$gene_name_day_14 <- NA
proxies$score_day_14 <- NA

for(index_proxy in seq(1, length(proxies$rsID))){
  
  chr_ <- proxies$chr[index_proxy]
  pos_ <- proxies$pos_hg19[index_proxy]
  
  peak_tmp <- stare_day14[which(as.numeric(stare_day14$`#chr`) == as.numeric(chr_) & as.numeric(stare_day14$start) <= as.numeric(pos_) & as.numeric(stare_day14$end) >= as.numeric(pos_))]
  peak_tmp$peak_region <- paste("chr", peak_tmp$'#chr', ":", peak_tmp$start, "_", peak_tmp$end, sep ="")
  
  if(is_empty(peak_tmp$start)){
    
    next()
    
  } else {
    
    print(index_proxy)
    
    proxies$stare_day_14[index_proxy] <- paste(unique(peak_tmp$peak_region), collapse =";")
    proxies$tss_dist_day_14[index_proxy] <- paste(peak_tmp$'TSS-dist', collapse =";")
    proxies$annotation_day_14[index_proxy] <- paste(peak_tmp$annotation, collapse =";")
    proxies$gene_id_day_14[index_proxy] <- paste(peak_tmp$gene_id, collapse =";")
    proxies$gene_name_day_14[index_proxy] <- paste(peak_tmp$'Gene Name', collapse =";")
    proxies$score_day_14[index_proxy] <- paste(peak_tmp$'ABC-Score', collapse =";")
    
  }
  
}

##############################################################################################
#Let's take a break here and divide the data according to their proximity to promoter effects#
##############################################################################################

fwrite(proxies, "output/7_functional_annotation/2_enhancer_gene_links/2_enhancer_gene_links_and_hic/proxies_with_stare_data.txt")

#########################################################################################################
#Let's use this opportunity to add the data from the Hi-C data that already exists for mature adipocytes#
#########################################################################################################

hic1 <- read.table(file = 'raw_data/hic_data/GSM3004355_BSA1_Down_hg19_cisInt.txt',header = T) %>% as.tibble(.) %>% dplyr::filter(CHiCAGOscore >= 5) %>% dplyr::filter(baitChr != "chrX")
hic2 <- read.table(file = 'raw_data/hic_data/GSE129574_Ad_pCHiC_washU.txt') %>%
  separate(., col = V1, into = c('baitChr', 'baitStart','baitEnd'), sep = ',') %>%
  separate(., col = V2, into = c('otherEndChr', 'otherEndStart','otherEndEnd'), sep = ',') %>%
  as_tibble(.) %>%
  mutate(baitChr = paste0('chr',baitChr),
         otherEndChr = paste0('chr',otherEndChr),
         baitStart = as.integer(baitStart),
         baitEnd = as.integer(baitEnd),
         otherEndStart = as.integer(otherEndStart),
         otherEndEnd = as.integer(otherEndEnd)) %>%
  dplyr::rename('CHiCAGOscore'= V3) %>%
  dplyr::filter(CHiCAGOscore >= 5) %>%
  dplyr::filter(baitChr !="chrX")

#Let's match the data for hic1, first for the variants that are overlapping one end:

proxies$adipocyte_bait_1 <- NA
proxies$adipocyte_other_1 <- NA

for(index_proxy in seq(1, length(proxies$rsID))){
  
  chr_ <- proxies$chr[index_proxy]
  pos_ <- proxies$pos_hg19[index_proxy]
  
  peak_tmp <- hic1[which(hic1$baitChr == paste("chr", chr_, sep = "") & as.numeric(hic1$baitStart) <= as.numeric(pos_) & as.numeric(hic1$baitEnd) >= as.numeric(pos_)),]
  peak_tmp$bait_region <- paste(peak_tmp$baitChr, ":", peak_tmp$baitStart, "_", peak_tmp$baitEnd, sep ="")
  peak_tmp$other_region <- paste(peak_tmp$otherEndChr , ":", peak_tmp$otherEndStart, "_", peak_tmp$otherEndEnd, sep ="")
  
  if(is_empty(peak_tmp$baitStart)){
    
    next()
    
  } else {
    
    print(index_proxy)
    
    proxies$adipocyte_bait_1[index_proxy] <- paste(unique(peak_tmp$bait_region), collapse =";")
    proxies$adipocyte_other_1[index_proxy] <- paste(unique(peak_tmp$other_region), collapse =";")
    
  }
  
}

proxies$adipocyte_bait_2 <- NA
proxies$adipocyte_other_2 <- NA

for(index_proxy in seq(1, length(proxies$rsID))){
  
  chr_ <- proxies$chr[index_proxy]
  pos_ <- proxies$pos_hg19[index_proxy]
  
  peak_tmp <- hic1[which(hic1$baitChr == paste("chr", chr_, sep = "") & as.numeric(hic1$otherEndStart) <= as.numeric(pos_) & as.numeric(hic1$otherEndEnd) >= as.numeric(pos_)),]
  peak_tmp$bait_region <- paste(peak_tmp$baitChr, ":", peak_tmp$baitStart, "_", peak_tmp$baitEnd, sep ="")
  peak_tmp$other_region <- paste(peak_tmp$otherEndChr , ":", peak_tmp$otherEndStart, "_", peak_tmp$otherEndEnd, sep ="")
  
  if(is_empty(peak_tmp$baitStart)){
    
    next()
    
  } else {
    
    print(index_proxy)
    
    proxies$adipocyte_bait_2[index_proxy] <- paste(unique(peak_tmp$bait_region), collapse =";")
    proxies$adipocyte_other_2[index_proxy] <- paste(unique(peak_tmp$other_region), collapse =";")
    
  }
  
}

#Let's do the same for hic_2!!

proxies$adipocyte_bait_3 <- NA
proxies$adipocyte_other_3 <- NA

for(index_proxy in seq(1, length(proxies$rsID))){
  
  chr_ <- proxies$chr[index_proxy]
  pos_ <- proxies$pos_hg19[index_proxy]
  
  peak_tmp <- hic2[which(hic2$baitChr == paste("chr", chr_, sep = "") & as.numeric(hic2$baitStart) <= as.numeric(pos_) & as.numeric(hic2$baitEnd) >= as.numeric(pos_)),]
  peak_tmp$bait_region <- paste(peak_tmp$baitChr, ":", peak_tmp$baitStart, "_", peak_tmp$baitEnd, sep ="")
  peak_tmp$other_region <- paste(peak_tmp$otherEndChr , ":", peak_tmp$otherEndStart, "_", peak_tmp$otherEndEnd, sep ="")
  
  if(is_empty(peak_tmp$baitStart)){
    
    next()
    
  } else {
    
    print(index_proxy)
    
    proxies$adipocyte_bait_3[index_proxy] <- paste(unique(peak_tmp$bait_region), collapse =";")
    proxies$adipocyte_other_3[index_proxy] <- paste(unique(peak_tmp$other_region), collapse =";")
    
  }
  
}

proxies$adipocyte_bait_4 <- NA
proxies$adipocyte_other_4 <- NA

for(index_proxy in seq(1, length(proxies$rsID))){
  
  chr_ <- proxies$chr[index_proxy]
  pos_ <- proxies$pos_hg19[index_proxy]
  
  peak_tmp <- hic2[which(hic2$baitChr == paste("chr", chr_, sep = "") & as.numeric(hic2$otherEndStart) <= as.numeric(pos_) & as.numeric(hic2$otherEndEnd) >= as.numeric(pos_)),]
  peak_tmp$bait_region <- paste(peak_tmp$baitChr, ":", peak_tmp$baitStart, "_", peak_tmp$baitEnd, sep ="")
  peak_tmp$other_region <- paste(peak_tmp$otherEndChr , ":", peak_tmp$otherEndStart, "_", peak_tmp$otherEndEnd, sep ="")
  
  if(is_empty(peak_tmp$baitStart)){
    
    next()
    
  } else {
    
    print(index_proxy)
    
    proxies$adipocyte_bait_4[index_proxy] <- paste(unique(peak_tmp$bait_region), collapse =";")
    proxies$adipocyte_other_4[index_proxy] <- paste(unique(peak_tmp$other_region), collapse =";")
    
  }
  
}

##########################################################
#Let's also take the chance to add the data in SGBS cells#
##########################################################

preadipocytes_hic <- fread('raw_data/hic_data/GSE262496_SGBS_Undiff_4frag.ibed.gz') 
adipocytes_hic <- fread('raw_data/hic_data/GSE262496_SGBS_Diff_4frag.ibed.gz') 

#This is build 38!! Careful

#First let's filter for scores >= 5

preadipocytes_hic <- preadipocytes_hic[which(as.numeric(preadipocytes_hic$score) >= 5),]
adipocytes_hic <- adipocytes_hic[which(as.numeric(adipocytes_hic$score) >= 5),]

#All of them are already scored like that. Great, let's change the names so that we can recycle the loops better...

colnames(preadipocytes_hic)[which(colnames(preadipocytes_hic) == "bait_chr")] <- "baitChr"
colnames(preadipocytes_hic)[which(colnames(preadipocytes_hic) == "bait_start")] <- "baitStart"
colnames(preadipocytes_hic)[which(colnames(preadipocytes_hic) == "bait_end")] <- "baitEnd"

colnames(preadipocytes_hic)[which(colnames(preadipocytes_hic) == "otherEnd_chr")] <- "otherEndChr"
colnames(preadipocytes_hic)[which(colnames(preadipocytes_hic) == "otherEnd_start")] <- "otherEndStart"
colnames(preadipocytes_hic)[which(colnames(preadipocytes_hic) == "otherEnd_end")] <- "otherEndEnd"

#And now the same for mature adipocytes

colnames(adipocytes_hic)[which(colnames(adipocytes_hic) == "bait_chr")] <- "baitChr"
colnames(adipocytes_hic)[which(colnames(adipocytes_hic) == "bait_start")] <- "baitStart"
colnames(adipocytes_hic)[which(colnames(adipocytes_hic) == "bait_end")] <- "baitEnd"

colnames(adipocytes_hic)[which(colnames(adipocytes_hic) == "otherEnd_chr")] <- "otherEndChr"
colnames(adipocytes_hic)[which(colnames(adipocytes_hic) == "otherEnd_start")] <- "otherEndStart"
colnames(adipocytes_hic)[which(colnames(adipocytes_hic) == "otherEnd_end")] <- "otherEndEnd"

####################################################
#Let's run the loop for undifferentiated SGBS cells#
####################################################

proxies$undifferentiated_bait_1 <- NA
proxies$undifferentiated_other_1 <- NA

for(index_proxy in seq(1, length(proxies$rsID))){
  
  chr_ <- proxies$chr[index_proxy]
  pos_ <- proxies$pos_hg38[index_proxy]
  
  peak_tmp <- preadipocytes_hic[which(preadipocytes_hic$baitChr == paste("chr", chr_, sep = "") & as.numeric(preadipocytes_hic$baitStart) <= as.numeric(pos_) & as.numeric(preadipocytes_hic$baitEnd) >= as.numeric(pos_)),]
  peak_tmp$bait_region <- paste(peak_tmp$baitChr, ":", peak_tmp$baitStart, "_", peak_tmp$baitEnd, sep ="")
  peak_tmp$other_region <- paste(peak_tmp$otherEndChr , ":", peak_tmp$otherEndStart, "_", peak_tmp$otherEndEnd, sep ="")
  
  if(is_empty(peak_tmp$baitStart)){
    
    next()
    
  } else {
    
    print(index_proxy)
    
    proxies$undifferentiated_bait_1[index_proxy] <- paste(unique(peak_tmp$bait_region), collapse =";")
    proxies$undifferentiated_other_1[index_proxy] <- paste(unique(peak_tmp$other_region), collapse =";")
    
  }
  
}

proxies$undifferentiated_bait_2 <- NA
proxies$undifferentiated_other_2 <- NA

for(index_proxy in seq(1, length(proxies$rsID))){
  
  chr_ <- proxies$chr[index_proxy]
  pos_ <- proxies$pos_hg38[index_proxy]
  
  peak_tmp <- preadipocytes_hic[which(preadipocytes_hic$baitChr == paste("chr", chr_, sep = "") & as.numeric(preadipocytes_hic$otherEndStart) <= as.numeric(pos_) & as.numeric(preadipocytes_hic$otherEndEnd) >= as.numeric(pos_)),]
  peak_tmp$bait_region <- paste(peak_tmp$baitChr, ":", peak_tmp$baitStart, "_", peak_tmp$baitEnd, sep ="")
  peak_tmp$other_region <- paste(peak_tmp$otherEndChr , ":", peak_tmp$otherEndStart, "_", peak_tmp$otherEndEnd, sep ="")
  
  if(is_empty(peak_tmp$baitStart)){
    
    next()
    
  } else {
    
    print(index_proxy)
    
    proxies$undifferentiated_bait_2[index_proxy] <- paste(unique(peak_tmp$bait_region), collapse =";")
    proxies$undifferentiated_other_2[index_proxy] <- paste(unique(peak_tmp$other_region), collapse =";")
    
  }
  
}

##################################################
#Let's run the loop for differentiated SGBS cells#
##################################################

proxies$differentiated_bait_1 <- NA
proxies$differentiated_other_1 <- NA

for(index_proxy in seq(1, length(proxies$rsID))){
  
  chr_ <- proxies$chr[index_proxy]
  pos_ <- proxies$pos_hg38[index_proxy]
  
  peak_tmp <- adipocytes_hic[which(adipocytes_hic$baitChr == paste("chr", chr_, sep = "") & as.numeric(adipocytes_hic$baitStart) <= as.numeric(pos_) & as.numeric(adipocytes_hic$baitEnd) >= as.numeric(pos_)),]
  peak_tmp$bait_region <- paste(peak_tmp$baitChr, ":", peak_tmp$baitStart, "_", peak_tmp$baitEnd, sep ="")
  peak_tmp$other_region <- paste(peak_tmp$otherEndChr , ":", peak_tmp$otherEndStart, "_", peak_tmp$otherEndEnd, sep ="")
  
  if(is_empty(peak_tmp$baitStart)){
    
    next()
    
  } else {
    
    print(index_proxy)
    
    proxies$differentiated_bait_1[index_proxy] <- paste(unique(peak_tmp$bait_region), collapse =";")
    proxies$differentiated_other_1[index_proxy] <- paste(unique(peak_tmp$other_region), collapse =";")
    
  }
  
}

proxies$differentiated_bait_2 <- NA
proxies$differentiated_other_2 <- NA

for(index_proxy in seq(1, length(proxies$rsID))){
  
  chr_ <- proxies$chr[index_proxy]
  pos_ <- proxies$pos_hg38[index_proxy]
  
  peak_tmp <- adipocytes_hic[which(adipocytes_hic$baitChr == paste("chr", chr_, sep = "") & as.numeric(adipocytes_hic$otherEndStart) <= as.numeric(pos_) & as.numeric(adipocytes_hic$otherEndEnd) >= as.numeric(pos_)),]
  peak_tmp$bait_region <- paste(peak_tmp$baitChr, ":", peak_tmp$baitStart, "_", peak_tmp$baitEnd, sep ="")
  peak_tmp$other_region <- paste(peak_tmp$otherEndChr , ":", peak_tmp$otherEndStart, "_", peak_tmp$otherEndEnd, sep ="")
  
  if(is_empty(peak_tmp$baitStart)){
    
    next()
    
  } else {
    
    print(index_proxy)
    
    proxies$differentiated_bait_2[index_proxy] <- paste(unique(peak_tmp$bait_region), collapse =";")
    proxies$differentiated_other_2[index_proxy] <- paste(unique(peak_tmp$other_region), collapse =";")
    
  }
  
}

######################
#Let's save this data#
######################

dir.create("output/7_functional_annotation/2_enhancer_gene_links/2_enhancer_gene_links_and_hic/")

fwrite(proxies, "output/7_functional_annotation/2_enhancer_gene_links/2_enhancer_gene_links_and_hic/proxies_with_stare_and_hic_data.txt")

