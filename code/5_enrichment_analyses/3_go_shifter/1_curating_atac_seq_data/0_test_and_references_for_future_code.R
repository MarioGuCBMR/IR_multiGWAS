##############
#INTRODUCTION#
##############

#This code reads ROADMAP Epigen data for adipose nuclei H3K27ac ATAC-seq data downloaded the 6/11/2024 from: https://egg2.wustl.edu/roadmap/data/byFileType/peaks/consolidated/narrowPeak/
#These are consolidated narrow peaks.

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)

#####################
#Let's read the data#
#####################

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/")

e023 <- data.table::fread("raw_data/atac_seq_data/E023_15_coreMarks_dense.bed")
e063 <- data.table::fread("raw_data/atac_seq_data/E063_15_coreMarks_dense.bed")
e063_18_marks <- data.table::fread("raw_data/atac_seq_data/E063_18_core_K27ac_dense.bed.gz")
# e063_25_marks <- data.table::fread("raw_data/atac_seq_data/E063_25_imputed12marks_dense.bed.gz")
# e063_25_marks <- e063_25_marks[which(str_detect(e063_25_marks$V4, "Enh") | str_detect(e063_25_marks$V4, "Tss") | str_detect(e063_25_marks$V4, "Prom") ),]

#Let's filter for promoter and enhancer:

promoter_and_enhancer_15 <- c("1_TssA", "2_TssAFlnk", "6_EnhG", "7_Enh", "10_TssBiv", "11_BivFlnk", "12_EnhBiv") #Currin et al (2019)
promoter_and_enhancer_18 <- c("1_TssA", "2_TssFlnk", "3_TssFlnkU", "4_TssFlnkD", "14_TssBiv", "7_EnhG1", "8_EnhG2", "9_EnhA1", "10_EnhA2", "11_EnhWk", "15_EnhBiv") #Currin et al (2019)

e023 <- e023[which(e023$V4%in%promoter_and_enhancer_15),]
e063 <- e063[which(e063$V4%in%promoter_and_enhancer_15),]
e063_18_marks <- e063_18_marks[which(e063_18_marks$V4%in%promoter_and_enhancer_18),]

#Let's change the names:

colnames(e023) <- c("#chr", "start", "end", "state", "0", ".","start_2", "end_2", "other")
colnames(e063) <- c("#chr", "start", "end", "state", "0", ".","start_2", "end_2", "other")
colnames(e063_18_marks) <- c("#chr", "start", "end", "state", "0", ".","start_2", "end_2", "other")
#colnames(e063_25_marks) <- c("#chr", "start", "end", "state", "0", ".","start_2", "end_2", "other")

#Past trials with H3K27ac which does not really capture every change...

#narrow_peak <- data.table::fread("raw_data/atac_seq_data/E063-H3K27ac.narrowPeak")
#broad_peak <- data.table::fread("raw_data/atac_seq_data/E063-H3K27ac.broadPeak")
#gapped_peak <- data.table::fread("raw_data/atac_seq_data/E063-H3K27ac.gappedPeak")

#############################################
#Let's filter the data with different values#
#############################################

#narrow_peak <- narrow_peak[,c(1,2,3,9)]
#broad_peak <- broad_peak[,c(1,2,3,9)]
#gapped_peak <- gapped_peak[,c(1,2,3,9)]

#Let's name the columns:

#colnames(narrow_peak) <- c("#CHR", "start", "end", "value")
#colnames(broad_peak) <- c("#CHR", "start", "end", "value")
#colnames(gapped_peak) <- c("#CHR", "start", "end", "value")

#####################
#Let's save the data# #this format will help to run STARE
#####################

dir.create("output/5_enrichment_analyses/")
#dir.create("output/5_enrichment_analyses/1_grs")
#dir.create("output/5_enrichment_analyses/2_depict")
dir.create("output/5_enrichment_analyses/3_go_shifter")
dir.create("output/5_enrichment_analyses/3_go_shifter/curated_atac_seq")

#data.table::fwrite(narrow_peak, "output/5_enrichment_analyses/3_go_shifter/curated_atac_seq/narrow_peak_4_stare.bed", sep="\t")
#data.table::fwrite(broad_peak, "output/5_enrichment_analyses/3_go_shifter/curated_atac_seq/broad_peak_4_stare.bed", sep="\t")
#data.table::fwrite(gapped_peak, "output/5_enrichment_analyses/3_go_shifter/curated_atac_seq/gapped_peak_4_stare.bed", sep="\t")

##############################################
#Let's change the format a bit for go-shifter#
##############################################

e023 <- e023 %>%
  dplyr::select("#chr","start","end")

fwrite(unique(e023), "output/5_enrichment_analyses/3_go_shifter/curated_atac_seq/peaks_roadmap_e023.bed", sep = "\t", col.names = F)

#same for e063:

e063 <- e063 %>%
  dplyr::select("#chr","start","end")

fwrite(unique(e063), "output/5_enrichment_analyses/3_go_shifter/curated_atac_seq/peaks_roadmap_e063.bed", sep = "\t", col.names = F)

#And one more:

e063_18_marks <- e063_18_marks %>%
  dplyr::select("#chr","start","end")

fwrite(unique(e063_18_marks), "output/5_enrichment_analyses/3_go_shifter/curated_atac_seq/peaks_roadmap_e063_18_marks.bed", sep = "\t", col.names = F)

#Finally:

#e063_25_marks <- e063_25_marks %>%
#  dplyr::select("#chr","start","end")

#fwrite(unique(e063_25_marks), "output/5_enrichment_analyses/3_go_shifter/curated_atac_seq/peaks_roadmap_e063_25_marks.bed", sep = "\t", col.names = F)
