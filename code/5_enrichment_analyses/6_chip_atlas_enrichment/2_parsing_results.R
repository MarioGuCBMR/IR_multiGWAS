##############
#INTRODUCTION#
##############

#This code prepares data for foot-printing enrichment analyses with ChIP-atlas.

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)
library(MotifDb)
library(igraph)

###################
#Loading functions#
###################

cleaning_peaks <- function(enrichment_check_1){
  
  chr <- unlist(str_split(enrichment_check_1, ":"))
  chr_ <- chr[which(str_detect(chr, "chr"))]
  pos_ <- chr[which(str_detect(chr, "_"))]
  
  pos_ <- str_split(pos_, "_")
  
  start_ <- c()
  end_ <- c()
  
  for(i in seq(1, length(pos_))){
    
    start_tmp <- pos_[[i]][1]
    end_tmp <- pos_[[i]][2]
    
    start_ <- c(start_, start_tmp)
    end_ <- c(end_, end_tmp)
    
  }
  
  data_4_chip <- cbind(chr_, start_, end_)
  
  return(data_4_chip)
  
  
}



##############
#Loading data#
##############

#STEP 1: get the data:

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025")

peaks_neg_clean <- fread("output/5_enrichment_analyses/6_chip_atlas/output/day_4_IR_vs_non_IR_loci.tsv")
peaks_ns_clean <- fread("output/5_enrichment_analyses/6_chip_atlas/output/day_14_IR_vs_non_IR_loci.tsv")

####################################
#Let's compare the results, then...#
####################################

colnames(peaks_neg_clean) <- c("ID", "antigen_class", "antigen", "cell_class", "cell", "number_of_peaks", "overlaps_a", "overlaps_b", "logP", "logQ", "fold_enrichment")

cells_ok <- c("Adipocytes", "SGBS")

peaks_neg_clean <- peaks_neg_clean[which(peaks_neg_clean$cell%in%cells_ok),]
peaks_neg_clean <- peaks_neg_clean[which(peaks_neg_clean$overlaps_a != "0/72"),]
peaks_neg_clean$q_value <- 10^(peaks_neg_clean$logQ)

peaks_neg_clean <- peaks_neg_clean[which(peaks_neg_clean$q_value < 0.05),] 

#We are gonna use SRX7807102!!!

####################################
#Let's compare the results, then...#
####################################

colnames(peaks_ns_clean) <- c("ID", "antigen_class", "antigen", "cell_class", "cell", "number_of_peaks", "overlaps_a", "overlaps_b", "logP", "logQ", "fold_enrichment")

cells_ok <- c("Adipocytes", "SGBS")

peaks_ns_clean <- peaks_ns_clean[which(peaks_ns_clean$cell%in%cells_ok),]
peaks_ns_clean <- peaks_ns_clean[which(peaks_ns_clean$overlaps_a != "0/126"),]
peaks_ns_clean$q_value <- 10^(peaks_ns_clean$logQ)

peaks_ns_clean <- peaks_ns_clean[which(peaks_ns_clean$q_value < 0.05),] #super interesting, but we will pass since we did not replicate them with HOMER.
