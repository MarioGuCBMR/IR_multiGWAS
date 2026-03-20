##############
#INTRODUCTION#
##############

#This code reads ATAC-seq data and prepares it to perform adapted ABC model STARE.

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)

#####################
#Let's read the data#
#####################

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/")

e063 <- fread("raw_data/other_chip_seq_data/E063-H3K27ac.narrowPeak.gz")

sgbs_day0 <- fread("raw_data/atac_seq_data/GSE178794_SGBS-day0_consensusPeaks-fromUnionRepPeaks_top100k-byMedianPval.bed.gz")
sgbs_day4 <- fread("raw_data/atac_seq_data/GSE178794_SGBS-day4_consensusPeaks-fromUnionRepPeaks_top100k-byMedianPval.bed.gz")
sgbs_day14 <- fread("raw_data/atac_seq_data/GSE178794_SGBS-day14_consensusPeaks-fromUnionRepPeaks_top100k-byMedianPval.bed.gz")

#############################################
#Let's filter the data with different values#
#############################################

e063 <- e063 %>%
  dplyr::select("V1", "V2", "V3", "V9")

sgbs_day0 <- sgbs_day0 %>%
  dplyr::select("#chr", "start", "stop", "median_qValue")

sgbs_day4 <- sgbs_day4 %>%
  dplyr::select("#chr", "start", "stop", "median_qValue")

sgbs_day14 <- sgbs_day14 %>%
  dplyr::select("#chr", "start", "stop", "median_qValue")

#Let's name the columns:

colnames(e063) <- c("#CHR", "start", "end", "value")
colnames(sgbs_day0) <- c("#CHR", "start", "end", "value")
colnames(sgbs_day4) <- c("#CHR", "start", "end", "value")
colnames(sgbs_day14) <- c("#CHR", "start", "end", "value")

#####################
#Let's save the data#
#####################

dir.create("output/7_functional_annotation/2_enhancer_gene_links/")
dir.create("output/7_functional_annotation/2_enhancer_gene_links/1_stare")
dir.create("output/7_functional_annotation/2_enhancer_gene_links/1_stare/input")
dir.create("output/7_functional_annotation/2_enhancer_gene_links/1_stare/output")

data.table::fwrite(e063, "output/7_functional_annotation/2_enhancer_gene_links/1_stare/input/roadmap_e063_4_stare.bed", sep="\t")
data.table::fwrite(sgbs_day0, "output/7_functional_annotation/2_enhancer_gene_links/1_stare/input/sgbs_day_0_4_stare.bed", sep="\t")
data.table::fwrite(sgbs_day4, "output/7_functional_annotation/2_enhancer_gene_links/1_stare/input/sgbs_day_4_4_stare.bed", sep="\t")
data.table::fwrite(sgbs_day14, "output/7_functional_annotation/2_enhancer_gene_links/1_stare/input/sgbs_day_14_4_stare.bed", sep="\t")
