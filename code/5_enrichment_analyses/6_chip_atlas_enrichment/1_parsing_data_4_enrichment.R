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

#STEP 1: get the data as in MOTIF DF with the same settings as the ones used in motifbreakr:

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025")

enhancer_df <- fread("output/7_functional_annotation/2_enhancer_gene_links/2_enhancer_gene_links_and_hic/proxies_with_stare_and_hic_data.txt")

peaks_4 <- fread('raw_data/atac_seq_data/GSE178794_SGBS-day4_consensusPeaks-fromUnionRepPeaks_top100k-byMedianPval.bed.gz')

peaks_4$id <- paste(peaks_4$`#chr`, ":", peaks_4$start, "_", peaks_4$stop, sep = "")   

peaks_14 <- fread('raw_data/atac_seq_data/GSE178794_SGBS-day14_consensusPeaks-fromUnionRepPeaks_top100k-byMedianPval.bed.gz')

peaks_14$id <- paste(peaks_14$`#chr`, ":", peaks_14$start, "_", peaks_14$stop, sep = "") 

###################################################################################################################
#STEP 1: before running the analyses for motifs, we are going to do an enrichment of foot-printing with ChIP-atlas#
###################################################################################################################

bmi_neg <- fread("output/5_enrichment_analyses/3_go_shifter/input_bmi_neg.txt")

enhancer_neg <- enhancer_df[which(enhancer_df$query_snp_rsid%in%bmi_neg$SNP),]
enhancer_neg <- enhancer_neg[which(enhancer_neg$sgbs_day4 != ""),]

#Let's take the peaks:

peaks_neg <- unique(enhancer_neg$sgbs_day4) #72

#Now let's take controls:

rest_of_peaks <- unique(peaks_4$id)
rest_of_peaks <- rest_of_peaks[which(!(rest_of_peaks%in%peaks_neg))]
rest_of_peaks <- rest_of_peaks[which(rest_of_peaks != "")]

#Finally, let's get individuals peaks...

peaks_neg_clean <- as.data.frame(cleaning_peaks(peaks_neg))

rest_of_peaks_day_4 <- peaks_4[which(peaks_4$id%in%rest_of_peaks),]

rest_of_peaks_day_4 <- rest_of_peaks_day_4 %>%
  select("#chr", "start", "stop")

dir.create("output/5_enrichment_analyses/6_chip_atlas")
dir.create("output/5_enrichment_analyses/6_chip_atlas/input")
dir.create("output/5_enrichment_analyses/6_chip_atlas/output")

fwrite(peaks_neg_clean, "output/5_enrichment_analyses/6_chip_atlas/input/day_4_overlapping_bmi_neg_4_enrichment.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote=FALSE)
fwrite(rest_of_peaks_day_4, "output/5_enrichment_analyses/6_chip_atlas/input/day_4_NOT_overlapping_bmi_neg_4_enrichment.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote=FALSE)

###################################################################################################################
#STEP 1: before running the analyses for motifs, we are going to do an enrichment of foot-printing with ChIP-atlas#
###################################################################################################################

bmi_ns <- fread("output/5_enrichment_analyses/3_go_shifter/input_bmi_ns.txt")

enhancer_ns <- enhancer_df[which(enhancer_df$query_snp_rsid%in%bmi_ns$SNP),]
enhancer_ns <- enhancer_ns[which(enhancer_ns$sgbs_day14 != ""),]

#Let's take the peaks:

peaks_ns <- unique(enhancer_ns$sgbs_day14) #126

#Now let's take controls:

rest_of_peaks <- unique(peaks_14$id)
rest_of_peaks <- rest_of_peaks[which(!(rest_of_peaks%in%peaks_ns))]
rest_of_peaks <- rest_of_peaks[which(rest_of_peaks != "")]

#Finally, let's get individuals peaks...

peaks_ns_clean <- as.data.frame(cleaning_peaks(peaks_ns))

rest_of_peaks_day_14 <- peaks_14[which(peaks_14$id%in%rest_of_peaks),]

rest_of_peaks_day_14 <- rest_of_peaks_day_14 %>%
  select("#chr", "start", "stop")

dir.create("output/5_enrichment_analyses/6_chip_atlas")
dir.create("output/5_enrichment_analyses/6_chip_atlas/input")
dir.create("output/5_enrichment_analyses/6_chip_atlas/output")

fwrite(peaks_ns_clean, "output/5_enrichment_analyses/6_chip_atlas/input/day_14_overlapping_bmi_ns_4_enrichment.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote=FALSE)
fwrite(rest_of_peaks_day_14, "output/5_enrichment_analyses/6_chip_atlas/input/day_14_NOT_overlapping_bmi_ns_4_enrichment.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote=FALSE)

###########################################################################################################################
#We are going to have to do these analyses for all peaks for the 282 IR variants for the final version of the manuscript!!#
###########################################################################################################################

#Let's go with adipose tissue - let's load the raw data first:

peaks_bulk <- fread('raw_data/atac_seq_data/GSE178794_adiposeTissue_consensusPeaks-from11samples_unionPeaksInAtLeast3samples.bed.gz')
peaks_bulk$id <- paste(peaks_bulk$`#chr`, ":", peaks_bulk$start, "_", peaks_bulk$stop, sep = "") 

#Now let's filter for the proxies overlapping adipose consensus:

enhancer_adipose <- enhancer_df[which(enhancer_df$adipose_tissue_peak != ""),] #218

#Let's get the peaks:

peaks_adipose <- unique(enhancer_adipose$adipose_tissue_peak) #157 - same as HOMER. Makes all the sense in the world

#Now let's take controls:

rest_of_peaks <- unique(peaks_bulk$id)
rest_of_peaks <- rest_of_peaks[which(!(rest_of_peaks%in%peaks_adipose))]
rest_of_peaks <- rest_of_peaks[which(rest_of_peaks != "")]

#Finally, let's get dataframes...

#First the dataframe for the peaks that overlaps

peaks_adipose_clean <- as.data.frame(cleaning_peaks(peaks_adipose))

rest_of_peaks_bulk <- peaks_bulk[which(peaks_bulk$id%in%rest_of_peaks),] #doing it like this we filter the original dataframe with the non-overlapping

rest_of_peaks_bulk <- rest_of_peaks_bulk %>%
  select("#chr", "start", "stop")

dir.create("output/5_enrichment_analyses/6_chip_atlas")
dir.create("output/5_enrichment_analyses/6_chip_atlas/input")
dir.create("output/5_enrichment_analyses/6_chip_atlas/output")

fwrite(peaks_adipose_clean, "output/5_enrichment_analyses/6_chip_atlas/input/adipose_bulk_overlapping_282_ir_4_enrichment.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote=FALSE)
fwrite(rest_of_peaks_bulk, "output/5_enrichment_analyses/6_chip_atlas/input/adipose_bulk_NOT_overlapping_282_ir_4_enrichment.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote=FALSE)

###########################################################################################################################
#We are going to have to do these analyses for all peaks for the 282 IR variants for the final version of the manuscript!!# DAY 0
###########################################################################################################################

#Let's go with adipose tissue - let's load the raw data first:

peaks_day0 <- fread('raw_data/atac_seq_data/GSE178794_SGBS-day0_consensusPeaks-fromUnionRepPeaks_top100k-byMedianPval.bed.gz')
peaks_day0$id <- paste(peaks_day0$`#chr`, ":", peaks_day0$start, "_", peaks_day0$stop, sep = "") 

#Now let's filter for the proxies overlapping adipose consensus:

enhancer_day0 <- enhancer_df[which(enhancer_df$sgbs_day0 != ""),] #354

#Let's get the peaks:

peaks_day0_overlap <- unique(enhancer_day0$sgbs_day0) #248 - same as HOMER.

#Now let's take controls:

rest_of_peaks <- unique(peaks_day0$id)
rest_of_peaks <- rest_of_peaks[which(!(rest_of_peaks%in%peaks_day0_overlap))]
rest_of_peaks <- rest_of_peaks[which(rest_of_peaks != "")]

#Finally, let's get dataframes...

#First the dataframe for the peaks that overlaps

peaks_day0_clean <- as.data.frame(cleaning_peaks(peaks_day0_overlap))

rest_of_peaks_day0 <- peaks_day0[which(peaks_day0$id%in%rest_of_peaks),] #doing it like this we filter the original dataframe with the non-overlapping

rest_of_peaks_day0 <- rest_of_peaks_day0 %>%
  select("#chr", "start", "stop")

dir.create("output/5_enrichment_analyses/6_chip_atlas")
dir.create("output/5_enrichment_analyses/6_chip_atlas/input")
dir.create("output/5_enrichment_analyses/6_chip_atlas/output")

fwrite(peaks_day0_clean, "output/5_enrichment_analyses/6_chip_atlas/input/day0_overlapping_282_ir_4_enrichment.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote=FALSE)
fwrite(rest_of_peaks_day0, "output/5_enrichment_analyses/6_chip_atlas/input/day0_NOT_overlapping_282_ir_4_enrichment.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote=FALSE)

###########################################################################################################################
#We are going to have to do these analyses for all peaks for the 282 IR variants for the final version of the manuscript!!# DAY 4
###########################################################################################################################

#Let's go with adipose tissue - let's load the raw data first:

peaks_day4 <- fread('raw_data/atac_seq_data/GSE178794_SGBS-day4_consensusPeaks-fromUnionRepPeaks_top100k-byMedianPval.bed.gz')
peaks_day4$id <- paste(peaks_day4$`#chr`, ":", peaks_day4$start, "_", peaks_day4$stop, sep = "") 

#Now let's filter for the proxies overlapping adipose consensus:

enhancer_day4 <- enhancer_df[which(enhancer_df$sgbs_day4 != ""),] #354

#Let's get the peaks:

peaks_day4_overlap <- unique(enhancer_day4$sgbs_day4) #277 - same as HOMER.

#Now let's take controls:

rest_of_peaks <- unique(peaks_day4$id)
rest_of_peaks <- rest_of_peaks[which(!(rest_of_peaks%in%peaks_day4_overlap))]
rest_of_peaks <- rest_of_peaks[which(rest_of_peaks != "")]

#Finally, let's get dataframes...

#First the dataframe for the peaks that overlaps

peaks_day4_clean <- as.data.frame(cleaning_peaks(peaks_day4_overlap))

rest_of_peaks_day4 <- peaks_day4[which(peaks_day4$id%in%rest_of_peaks),] #doing it like this we filter the original dataframe with the non-overlapping

rest_of_peaks_day4 <- rest_of_peaks_day4 %>%
  select("#chr", "start", "stop")

dir.create("output/5_enrichment_analyses/6_chip_atlas")
dir.create("output/5_enrichment_analyses/6_chip_atlas/input")
dir.create("output/5_enrichment_analyses/6_chip_atlas/output")

fwrite(peaks_day4_clean, "output/5_enrichment_analyses/6_chip_atlas/input/day4_overlapping_282_ir_4_enrichment.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote=FALSE)
fwrite(rest_of_peaks_day4, "output/5_enrichment_analyses/6_chip_atlas/input/day4_NOT_overlapping_282_ir_4_enrichment.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote=FALSE)

###########################################################################################################################
#We are going to have to do these analyses for all peaks for the 282 IR variants for the final version of the manuscript!!# DAY 14
###########################################################################################################################

#Let's go with adipose tissue - let's load the raw data first:

peaks_day14 <- fread('raw_data/atac_seq_data/GSE178794_SGBS-day14_consensusPeaks-fromUnionRepPeaks_top100k-byMedianPval.bed.gz')
peaks_day14$id <- paste(peaks_day14$`#chr`, ":", peaks_day14$start, "_", peaks_day14$stop, sep = "") 

#Now let's filter for the proxies overlapping adipose consensus:

enhancer_day14 <- enhancer_df[which(enhancer_df$sgbs_day14 != ""),] #3514

#Let's get the peaks:

peaks_day14_overlap <- unique(enhancer_day14$sgbs_day14) #272 - same as HOMER.

#Now let's take controls:

rest_of_peaks <- unique(peaks_day14$id)
rest_of_peaks <- rest_of_peaks[which(!(rest_of_peaks%in%peaks_day14_overlap))]
rest_of_peaks <- rest_of_peaks[which(rest_of_peaks != "")]

#Finally, let's get dataframes...

#First the dataframe for the peaks that overlaps

peaks_day14_clean <- as.data.frame(cleaning_peaks(peaks_day14_overlap))

rest_of_peaks_day14 <- peaks_day14[which(peaks_day14$id%in%rest_of_peaks),] #doing it like this we filter the original dataframe with the non-overlapping

rest_of_peaks_day14 <- rest_of_peaks_day14 %>%
  select("#chr", "start", "stop")

dir.create("output/5_enrichment_analyses/6_chip_atlas")
dir.create("output/5_enrichment_analyses/6_chip_atlas/input")
dir.create("output/5_enrichment_analyses/6_chip_atlas/output")

fwrite(peaks_day14_clean, "output/5_enrichment_analyses/6_chip_atlas/input/day14_overlapping_282_ir_4_enrichment.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote=FALSE)
fwrite(rest_of_peaks_day14, "output/5_enrichment_analyses/6_chip_atlas/input/day14_NOT_overlapping_282_ir_4_enrichment.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote=FALSE)

###########################################################################################################################
#STEP sensitivy: before running the analyses for motifs, we are going to do an enrichment of foot-printing with ChIP-atlas#
###########################################################################################################################

bmi_neg <- fread("output/5_enrichment_analyses/3_go_shifter/input_bmi_neg.txt")

enhancer_neg <- enhancer_df[which(enhancer_df$query_snp_rsid%in%bmi_neg$SNP),]
enhancer_neg <- enhancer_neg[which(enhancer_neg$sgbs_day4 != ""),]
enhancer_neg_unique <- enhancer_neg[order(enhancer_neg$rsID, enhancer_neg$r2),]
enhancer_neg_unique <- enhancer_neg_unique[which(duplicated(enhancer_neg_unique$query_snp_rsid) == FALSE),]

#Let's take the peaks:

peaks_neg_unique <- unique(enhancer_neg_unique$sgbs_day4) #42
peaks_neg_all <- unique(enhancer_neg$sgbs_day4) #72

#Now let's take controls:

rest_of_peaks <- unique(peaks_4$id)
rest_of_peaks <- rest_of_peaks[which(!(rest_of_peaks%in%peaks_neg))]
rest_of_peaks <- rest_of_peaks[which(rest_of_peaks != "")]

#Finally, let's get individuals peaks...

peaks_neg_clean <- as.data.frame(cleaning_peaks(peaks_neg_unique))

rest_of_peaks_day_4 <- peaks_4[which(peaks_4$id%in%rest_of_peaks),]

rest_of_peaks_day_4 <- rest_of_peaks_day_4 %>%
  select("#chr", "start", "stop")

dir.create("output/5_enrichment_analyses/6_chip_atlas")
dir.create("output/5_enrichment_analyses/6_chip_atlas/input")
dir.create("output/5_enrichment_analyses/6_chip_atlas/output")

fwrite(peaks_neg_clean, "output/5_enrichment_analyses/6_chip_atlas/input/day_4_overlapping_bmi_neg_4_enrichment_sensitivity.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote=FALSE)
fwrite(rest_of_peaks_day_4, "output/5_enrichment_analyses/6_chip_atlas/input/day_4_NOT_overlapping_bmi_neg_4_enrichment_sensitivity.bed", sep = "\t", row.names = FALSE, col.names = FALSE, quote=FALSE)
