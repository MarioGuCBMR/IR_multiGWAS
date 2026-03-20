##############
#INTRODUCTION#
##############

#This code adds information of proxies that overlap ATAC-seq peaks in adipose consensus data.

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)
library(otargen)

###################
#Loading functions#
###################

parse_roadmap_data <- function(file_){
  #This function returns only peaks in promoters and enhancers from a certain roadmap data
  
  print(file_)
  
  #STEP 1: read the file:
  
  e_tmp <- data.table::fread(file_)
  
  #STEP 2: parse for the annotation we want:
  
  promoter_and_enhancer_15 <- c("1_TssA", "2_TssAFlnk", "6_EnhG", "7_Enh", "10_TssBiv", "11_BivFlnk", "12_EnhBiv") #Currin et al (2019)
  
  e_tmp <- e_tmp[which(e_tmp$V4%in%promoter_and_enhancer_15),]
  
  colnames(e_tmp) <- c("#chr", "start", "end", "state", "0", ".","start_2", "end_2", "other")
  
  #STEP 3: save the data: 
  
  e_clean <- e_tmp %>%
    dplyr::select("#chr","start","end")
  
  return(e_clean)  
}

##############
#Loading data#
##############

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025")

#Let's load the most updated datasets here:

ir_variants <- fread("output/4_ir_loci_discovery/3_novel_variants/4_colocalization/352_variants_with_all_info_plus_fiadjbmi_and_tg_hdl_coloc.txt")
proxies <- fread("output/4_ir_loci_discovery/4_proxies/proxies_build_37_and_38_4_282_variants.txt")

#Now let's focus only on the 282 IR variants:

ir_variants <- ir_variants[which(ir_variants$variant%in%proxies$rsID),] #282! Perfect

#And also let's filter the proxies just in case:

yes_vect <- c("A", "G", "C", "T")

proxies <- proxies[which(proxies$query_snp_rsid%in%ir_variants$variant),] #3015
proxies_curated <- proxies[which(proxies$ref%in%yes_vect & proxies$alt%in%yes_vect),] #perfect

############################
#Loading the ATAC-seq peaks#
############################

adipose_tissue_consensus <- fread("raw_data/atac_seq_data/GSE178794_adiposeTissue_consensusPeaks-from11samples_unionPeaksInAtLeast3samples.bed.gz")
roadmap_adipose_bulk <- parse_roadmap_data("raw_data/atac_seq_data/E063_15_coreMarks_dense.bed.gz")
roadmap_adipocyte <- parse_roadmap_data("raw_data/atac_seq_data/E023_15_coreMarks_dense.bed.gz")
sgbs_day0 <- fread("raw_data/atac_seq_data/GSE178794_SGBS-day0_consensusPeaks-fromUnionRepPeaks_top100k-byMedianPval.bed.gz")
sgbs_day4 <- fread("raw_data/atac_seq_data/GSE178794_SGBS-day4_consensusPeaks-fromUnionRepPeaks_top100k-byMedianPval.bed.gz")
sgbs_day14 <- fread("raw_data/atac_seq_data/GSE178794_SGBS-day14_consensusPeaks-fromUnionRepPeaks_top100k-byMedianPval.bed.gz")


##########################################################
#Let's match the proxy data with adipose tissue consensus#
##########################################################

proxies_curated$adipose_tissue_peak <- NA

for(index_proxy in seq(1, length(proxies_curated$rsID))){
  
  chr_ <- proxies_curated$chr[index_proxy]
  pos_ <- proxies_curated$pos_hg19[index_proxy]
  
  peak_tmp <- adipose_tissue_consensus[which(adipose_tissue_consensus$`#chr` == paste("chr", chr_, sep = "") & as.numeric(adipose_tissue_consensus$start) <= as.numeric(pos_) & as.numeric(adipose_tissue_consensus$stop) >= as.numeric(pos_))]
  peak_tmp$peak_region <- paste(peak_tmp$'#chr', ":", peak_tmp$start, "_", peak_tmp$stop, sep ="")
  
  if(is_empty(peak_tmp$start)){
    
    next()
    
  } else {
    
    print(index_proxy)
    
    proxies_curated$adipose_tissue_peak[index_proxy] <- paste(peak_tmp$peak_region, collapse =";")
  
  }
  
}

##################################################################
#Let's match the proxy data with adipose tissue bulk from roadmap#
##################################################################

proxies_curated$adipose_bulk_roadmap <- NA

for(index_proxy in seq(1, length(proxies_curated$rsID))){
  
  chr_ <- proxies_curated$chr[index_proxy]
  pos_ <- proxies_curated$pos_hg19[index_proxy]
  
  peak_tmp <- roadmap_adipocyte[which(roadmap_adipocyte$`#chr` == paste("chr", chr_, sep = "") & as.numeric(roadmap_adipocyte$start) <= as.numeric(pos_) & as.numeric(roadmap_adipocyte$end) >= as.numeric(pos_))]
  peak_tmp$peak_region <- paste(peak_tmp$'#chr', ":", peak_tmp$start, "_", peak_tmp$end, sep ="")
  
  if(is_empty(peak_tmp$start)){
    
    next()
    
  } else {
    
    print(index_proxy)
    
    proxies_curated$adipose_bulk_roadmap[index_proxy] <- paste(peak_tmp$peak_region, collapse =";")
    
  }
  
}


#####################################################################
#Let's match the proxy data with in-vivo adipocyte data from ROADMAP#
#####################################################################

proxies_curated$adipocyte_roadmap <- NA

for(index_proxy in seq(1, length(proxies_curated$rsID))){
  
  chr_ <- proxies_curated$chr[index_proxy]
  pos_ <- proxies_curated$pos_hg19[index_proxy]
  
  peak_tmp <- roadmap_adipose_bulk[which(roadmap_adipose_bulk$`#chr` == paste("chr", chr_, sep = "") & as.numeric(roadmap_adipose_bulk$start) <= as.numeric(pos_) & as.numeric(roadmap_adipose_bulk$end) >= as.numeric(pos_))]
  peak_tmp$peak_region <- paste(peak_tmp$'#chr', ":", peak_tmp$start, "_", peak_tmp$end, sep ="")
  
  if(is_empty(peak_tmp$start)){
    
    next()
    
  } else {
    
    print(index_proxy)
    
    proxies_curated$adipocyte_roadmap[index_proxy] <- paste(peak_tmp$peak_region, collapse =";")
    
  }
  
}

###############################################
#Let's match the proxy data with SGBS at day 0#
###############################################

proxies_curated$sgbs_day0 <- NA

for(index_proxy in seq(1, length(proxies_curated$rsID))){
  
  chr_ <- proxies_curated$chr[index_proxy]
  pos_ <- proxies_curated$pos_hg19[index_proxy]
  
  peak_tmp <- sgbs_day0[which(sgbs_day0$`#chr` == paste("chr", chr_, sep = "") & as.numeric(sgbs_day0$start) <= as.numeric(pos_) & as.numeric(sgbs_day0$stop) >= as.numeric(pos_))]
  peak_tmp$peak_region <- paste(peak_tmp$'#chr', ":", peak_tmp$start, "_", peak_tmp$stop, sep ="")
  
  if(is_empty(peak_tmp$start)){
    
    next()
    
  } else {
    
    print(index_proxy)
    
    proxies_curated$sgbs_day0[index_proxy] <- paste(peak_tmp$peak_region, collapse =";")
    
  }
  
}

###############################################
#Let's match the proxy data with SGBS at day 4#
###############################################

proxies_curated$sgbs_day4 <- NA

for(index_proxy in seq(1, length(proxies_curated$rsID))){
  
  chr_ <- proxies_curated$chr[index_proxy]
  pos_ <- proxies_curated$pos_hg19[index_proxy]
  
  peak_tmp <- sgbs_day4[which(sgbs_day4$`#chr` == paste("chr", chr_, sep = "") & as.numeric(sgbs_day4$start) <= as.numeric(pos_) & as.numeric(sgbs_day4$stop) >= as.numeric(pos_))]
  peak_tmp$peak_region <- paste(peak_tmp$'#chr', ":", peak_tmp$start, "_", peak_tmp$stop, sep ="")
  
  if(is_empty(peak_tmp$start)){
    
    next()
    
  } else {
    
    print(index_proxy)
    
    proxies_curated$sgbs_day4[index_proxy] <- paste(peak_tmp$peak_region, collapse =";")
    
  }
  
}

###############################################
#Let's match the proxy data with SGBS at day 14#
###############################################

proxies_curated$sgbs_day14 <- NA

for(index_proxy in seq(1, length(proxies_curated$rsID))){
  
  chr_ <- proxies_curated$chr[index_proxy]
  pos_ <- proxies_curated$pos_hg19[index_proxy]
  
  peak_tmp <- sgbs_day14[which(sgbs_day14$`#chr` == paste("chr", chr_, sep = "") & as.numeric(sgbs_day14$start) <= as.numeric(pos_) & as.numeric(sgbs_day14$stop) >= as.numeric(pos_))]
  peak_tmp$peak_region <- paste(peak_tmp$'#chr', ":", peak_tmp$start, "_", peak_tmp$stop, sep ="")
  
  if(is_empty(peak_tmp$start)){
    
    next()
    
  } else {
    
    print(index_proxy)
    
    proxies_curated$sgbs_day14[index_proxy] <- paste(peak_tmp$peak_region, collapse =";")
    
  }
  
}

###############################################################################################################################
#Finally, let's assess whether we can assess which peaks are distinctly mature adipocytes and which might be exclusively Day 4#
###############################################################################################################################

atac_seq <- fread("raw_data/sgbs_differential_data/normalized_counts_chromatin_accessibility_ATAC.csv")
peak_info_og <- fread("raw_data/sgbs_differential_data/peak_info.csv")
peak_info_og$chromosome <- peak_info_og$`#chr`
peak_info_og$peak_region <- paste(peak_info_og$'#chr', ":", peak_info_og$start, "_", peak_info_og$stop, sep ="")

#First, let's just arrange the original peak data so that we can have, more or less, 

D0D4 <- fread("raw_data/sgbs_differential_data/atac_0_4.csv")
D0D14 <- fread("raw_data/sgbs_differential_data/atac_0_14.csv")
D4D14 <- fread("raw_data/sgbs_differential_data/atac_4_14.csv")

#We have duplicates for each variant so let's overlap the data:

peak_info_og$day_0_vs_4 <- NA
peak_info_og$day_0_vs_14 <- NA
peak_info_og$day_4_vs_14 <- NA

################################################################################################################################
#Let's match the consensus peaks that match in all conditions first, then we can reduce the matches and add the interpretations#
################################################################################################################################

proxies_curated$sgbs_consensus_all_days <- NA

for(index_proxy in seq(1, length(proxies_curated$rsID))){
  
  chr_ <- proxies_curated$chr[index_proxy]
  pos_ <- proxies_curated$pos_hg19[index_proxy]
  
  peak_tmp <- peak_info_og[which(peak_info_og$`#chr` == paste("chr", chr_, sep = "") & as.numeric(peak_info_og$start) <= as.numeric(pos_) & as.numeric(peak_info_og$stop) >= as.numeric(pos_))]
  peak_tmp$peak_region <- paste(peak_tmp$'#chr', ":", peak_tmp$start, "_", peak_tmp$stop, sep ="")
  
  if(is_empty(peak_tmp$start)){
    
    next()
    
  } else {
    
    print(index_proxy)
    
    proxies_curated$sgbs_consensus_all_days[index_proxy] <- paste(peak_tmp$peak_region, collapse =";")
    
  }
  
}

################################################################
#Now we can add the info we want on the differences across days#
################################################################

peak_info_match <- peak_info_og[which(peak_info_og$peak_region%in%proxies_curated$sgbs_consensus_all_days),]

for(peak_index in seq(1, length(peak_info_match$peak_ID))){
  
  print(peak_index)
  
  ################
  #DAY 0 Vs DAY 4#
  ################
  
  #STEP 1: get the data for that peak
  
  d0d4_match <- D0D4[which(D0D4$peakID == peak_info_match$peak_ID[peak_index]),]
  
  #STEP 2: interpret:
  
  d0d4_interpretation_a = ifelse(as.numeric(d0d4_match$padj) <0.05, "change", "No")
  d0d4_interpretation_b = ifelse(d0d4_interpretation_a == "change" & as.numeric(d0d4_match$log2FoldChange) > 0 & as.numeric(d0d4_match$padj) < 0.05, "Up", d0d4_interpretation_a)
  d0d4_interpretation_b = ifelse(d0d4_interpretation_a == "change" & as.numeric(d0d4_match$log2FoldChange) < 0 & as.numeric(d0d4_match$padj) < 0.05, "Down", d0d4_interpretation_b)
  
  #STEP 3: add info:
  
  peak_info_match$day_0_vs_4[peak_index] = d0d4_interpretation_b
  
  #################
  #DAY 4 Vs DAY 14#
  #################
  
  #STEP 1: matching
  
  d4d14_match <- D4D14[which(D4D14$peakID == peak_info_match$peak_ID[peak_index]),]
  
  #STEP 2: interpreting:
  
  d4d14_interpretation_a = ifelse(d4d14_match$padj <0.05, "change", "No")
  d4d14_interpretation_b = ifelse(d4d14_interpretation_a == "change" & as.numeric(d4d14_match$log2FoldChange) > 0 & as.numeric(d4d14_match$padj) < 0.05, "Up", d4d14_interpretation_a)
  d4d14_interpretation_b = ifelse(d4d14_interpretation_a == "change" & as.numeric(d4d14_match$log2FoldChange) < 0 & as.numeric(d4d14_match$padj) < 0.05, "Down", d4d14_interpretation_b)
  
  #STEP 3: adding:
  
  peak_info_match$day_4_vs_14[peak_index] = d4d14_interpretation_b
  
  #################
  #DAY 0 Vs DAY 14#
  #################
  
  #STEP 1: matching
  
  d0d14_match <- D0D14[which(D0D14$peakID == peak_info_match$peak_ID[peak_index]),]
  
  #STEP 2: interpreting:
  
  d0d14_interpretation_a = ifelse(d0d14_match$padj <0.05, "change", "No")
  d0d14_interpretation_b = ifelse(d0d14_interpretation_a == "change" & as.numeric(d0d14_match$log2FoldChange) > 0 & as.numeric(d0d14_match$padj) < 0.05, "Up", d0d14_interpretation_a)
  d0d14_interpretation_b = ifelse(d0d14_interpretation_a == "change" & as.numeric(d0d14_match$log2FoldChange) < 0 & as.numeric(d0d14_match$padj) < 0.05, "Down", d0d14_interpretation_b)
  
  #STEP 3: adding:
  
  peak_info_match$day_0_vs_14[peak_index] = d0d14_interpretation_b
  
}

peak_info_match$interpretation <- paste(peak_info_match$day_0_vs_4, "_", peak_info_match$day_0_vs_14, "_", peak_info_match$day_4_vs_14, sep = "")

table(peak_info_match$interpretation)

#Down_Down_No: immature Decreases after day 0
#Down_No_No = immature  Decreases after day 0 ()
#No_Down_No = immature

consistent_vect <- c("No_No_No")
preadipocyte_vect <- c("Down_Down_No", "Down_No_No", "No_Down_No")
adipocyte_vect <- c("Up_No_Down", "Up_Up_No", "No_Up_No", "Up_No_No")
day4_vect <- c("Up_No_Down", "Up_Up_No", "Up_No_No")

##########################################################################################
#We have all the data to add to our dear peaks and proceed to do some cool analyses later#
##########################################################################################

proxies_curated$pairwise_diff_sgbs <- NA
proxies_curated$mature_classification <- NA

for(index_proxy in seq(1, length(proxies_curated$rsID))){
  
  region <- proxies_curated$sgbs_consensus_all_days[index_proxy]

  peak_tmp <- peak_info_match[which(peak_info_match$peak_region == region),]

  if(is_empty(peak_tmp$start)){
    
    next()
    
  } else {
    
    #STEP 1 assess which interpretation:
    
    classification <- "Similarly accessible across days"
    classification <- ifelse(peak_tmp$interpretation%in%adipocyte_vect, "More accessible in adipocytes", classification)
    classification <- ifelse(peak_tmp$interpretation%in%preadipocyte_vect, "More accessible in preadipocytes", classification)
    
    day_4 <- ifelse(classification == "More accessible in adipocytes" &  peak_tmp$interpretation%in%day4_vect, "More accessible in adipocytes (day 4)", NA)
    
    if(is.na(day_4)){
      
      proxies_curated$pairwise_diff_sgbs[index_proxy] <- peak_tmp$interpretation
      proxies_curated$mature_classification[index_proxy] <- classification
      
    } else {
      
      proxies_curated$pairwise_diff_sgbs[index_proxy] <- peak_tmp$interpretation
      proxies_curated$mature_classification[index_proxy] <- day_4
      
    }

  }
  
}

###################
#Let's save this!!#
###################

dir.create("output/7_functional_annotation/")
dir.create("output/7_functional_annotation/1_atac_seq_overlaps")

fwrite(proxies_curated, "output/7_functional_annotation/1_atac_seq_overlaps/proxies_matched_w_roadmap_and_perrin_data_w_sgbs_diff.txt")
