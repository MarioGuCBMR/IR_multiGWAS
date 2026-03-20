##############
#INTRODUCTION#
##############

#This code prepares the data for running HOMER TF enrichment.
#Here we will use the proxies for all 282 IR variants.
#Additionally, we will use the peaks for bulk adipose tissue data.

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)
library(GenomicRanges)

###################
#Loading libraries#
###################

overlap_GR_SNPs <- function(IR_proxies, GR_peaks){
  #This function take chr_pos from Coord column in the proxies DF and tries overlapping them with the Genomic Ranges data.
  
  #STEP 1: arrange the data from the coordinates into Genomic Ranges format:
  
  variants_position <- IR_proxies %>%
    separate(.,Coord,into = c("chr","start"),sep = ":") %>%
    makeGRangesFromDataFrame(.,
                             seqnames.field = "chr",
                             start.field = "start",
                             end.field = "start")
  
  #STEP 2: find the overlap:
  
  out <- findOverlaps(query = variants_position,
                      subject = GR_peaks,
                      ignore.strand = T)
  return(out)
}

############################################
#Let' load the peak data and treat it first#
############################################

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/")

peaks <- fread('raw_data/atac_seq_data/GSE178794_SGBS-day14_consensusPeaks-fromUnionRepPeaks_top100k-byMedianPval.bed.gz')

# Make the GRanges object for ATACseq peaks
GR_peaks_ <- makeGRangesFromDataFrame(df = peaks,
                                     keep.extra.columns = T,
                                     seqnames.field = "#chr",
                                     start.field = "start",
                                     end.field = "stop")


###############################
#Now let's load the proxy data#
###############################

proxies <- fread("output/4_ir_loci_discovery/4_proxies/proxies_build_37_and_38_4_282_variants.txt")

#Let's load the variants in this cluster, especifically,

bmi_ns <- fread("output/5_enrichment_analyses/3_go_shifter/input_bmi_ns.txt")

proxies <- proxies[which(proxies$query_snp_rsid%in%bmi_ns$SNP),]

length(unique(proxies$query_snp_rsid)) #141 - perfect

#Let's take forward only the coordinates.

proxies <- proxies %>% dplyr::select(chr, pos_hg19)
colnames(proxies) <- c("chr_37", "bp_37")

#Let's change the column names so that we can do an overlap:

proxies$chr_37 <- paste0("chr",proxies$chr_37)
proxies$Coord <- paste0(proxies$chr_37,":", proxies$bp_37)

proxies <- proxies %>% dplyr::select(Coord)

#########################
#Let's find the overlap!#
#########################

overlapped_data <- overlap_GR_SNPs(proxies, GR_peaks_)
overlapped_data_info <- data.frame(proxies[queryHits(overlapped_data),], GR_peaks_[subjectHits(overlapped_data),])

#Let's extract the overlapping information:

IR_peaks <- overlapped_data_info %>% dplyr::select( "seqnames","start","end","peak_ID")
IR_peaks$col5 <- c(".")
IR_peaks$col6 <- c(".")

#And save this data:

dir.create("output/5_enrichment_analyses/4_homer")
dir.create("output/5_enrichment_analyses/4_homer/input/")
dir.create("output/5_enrichment_analyses/4_homer/input/141_sgbs_day_14")
dir.create("output/5_enrichment_analyses/4_homer/output/141_sgbs_day_14")

fwrite(unique(IR_peaks), "output/5_enrichment_analyses/4_homer/input/141_sgbs_day_14/IR_peaks.txt", sep = "\t", col.names = F)

#################################
#Now let's find the non-overlaps#
#################################

# control regions: non-overlapping IR variants

#Let's first format the input data:

non_IR_peaks <-  peaks
colnames(non_IR_peaks) <- c("chr_peak","start_peak","end_peak","peakID", "median_foldEnrichment", "median_pValue", "median_qValue", "num_samples")
non_IR_peaks <- non_IR_peaks %>% dplyr::select("chr_peak","start_peak","end_peak","peakID")
non_IR_peaks$col5 <- c(".")
non_IR_peaks$col6 <- c(".")

#Remove the ones that are overlapping:

non_IR_peaks <- non_IR_peaks[!non_IR_peaks$peakID %in% IR_peaks$peak_ID,]

#And finally save them:

fwrite(unique(non_IR_peaks), "output/5_enrichment_analyses/4_homer/input/141_sgbs_day_14/non_IR_peaks.txt", sep = "\t", col.names = F)
