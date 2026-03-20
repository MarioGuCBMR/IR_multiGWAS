##############
#INTRODUCTION#
##############

#This is from consensus adipose tissue ATAC-seq data from Perrin et al.

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)

#####################
#Let's read the data#
#####################

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/")

consensus_adipose <- data.table::fread("raw_data/atac_seq_data/GSE178794_adiposeTissue_consensusPeaks-from11samples_unionPeaksInAtLeast3samples.bed.gz")
consensus_day0 <- data.table::fread("raw_data/atac_seq_data/GSE178794_SGBS-day0_consensusPeaks-fromUnionRepPeaks_top100k-byMedianPval.bed.gz")
consensus_day4 <- data.table::fread("raw_data/atac_seq_data/GSE178794_SGBS-day4_consensusPeaks-fromUnionRepPeaks_top100k-byMedianPval.bed.gz")
consensus_day14 <- data.table::fread("raw_data/atac_seq_data/GSE178794_SGBS-day14_consensusPeaks-fromUnionRepPeaks_top100k-byMedianPval.bed.gz")

#############################################
#Let's filter the data with different values#
#############################################

consensus_adipose <- consensus_adipose[,c(1,2,3)]
consensus_day0 <- consensus_day0[,c(1,2,3)]
consensus_day4 <- consensus_day4[,c(1,2,3)]
consensus_day14 <- consensus_day14[,c(1,2,3)]

#####################
#Let's save the data# #this format will help to run STARE
#####################

fwrite(unique(consensus_adipose), "output/5_enrichment_analyses/3_go_shifter/curated_atac_seq/peaks_adipose_tissue.bed", sep = "\t", col.names = F)
fwrite(unique(consensus_day0), "output/5_enrichment_analyses/3_go_shifter/curated_atac_seq/peaks_sgbs_day0.bed", sep = "\t", col.names = F)
fwrite(unique(consensus_day4), "output/5_enrichment_analyses/3_go_shifter/curated_atac_seq/peaks_sgbs_day4.bed", sep = "\t", col.names = F)
fwrite(unique(consensus_day14), "output/5_enrichment_analyses/3_go_shifter/curated_atac_seq/peaks_sgbs_day14.bed", sep = "\t", col.names = F)


