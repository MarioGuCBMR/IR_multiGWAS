##############
#INTRODUCTION#
##############

#This code performs CPASSOC with our three traits!!

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)

path_2_cpassoc <- "N:/SUN-CBMR-Kilpelainen-Group/insulin_resistance_variants_common_info/programs/"

setwd(path_2_cpassoc)

source("CPASSOC/FunctionSet.R")

##############
#Loading data#
##############

#Let's change the working directory:

path_2_files <-  "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/" #change it with your own path.

setwd(path_2_files)

final_sumstats <- fread("output/4_ir_loci_discovery/1_cpassoc/merged_fiadjbmi_hdl_tg_sumstats_4_cpassoc.txt")
corMatrix <- readRDS("output/4_ir_loci_discovery/1_cpassoc/correlation_matrix_between_traits.RDS")
#triangulated_snps <- fread("output/4_ir_loci_discovery/2_ir_variants/282_ir_fiadjbmi_hdl_tg_dwls_snps.txt")

triangulated_snps <- final_sumstats[which(as.numeric(final_sumstats$beta.fiadjbmi) > 0 & as.numeric(final_sumstats$beta.hdl) < 0 & as.numeric(final_sumstats$beta.tg) > 0),] #aligned already to FIadjbmi+ - 3M variants approximately! This will be running for quite a while, but that is OK

#########################################################
#Now let's perform some pruning to get this data rolling#
#########################################################

#Let's get the z-score and the sample size to run this properly.

final_sumstats$z_score.fiadjbmi <- as.numeric(final_sumstats$beta.fiadjbmi)/as.numeric(final_sumstats$standard_error.fiadjbmi)
final_sumstats$z_score.hdl <- as.numeric(final_sumstats$beta.hdl)/as.numeric(final_sumstats$standard_error.hdl)
final_sumstats$z_score.tg <- as.numeric(final_sumstats$beta.tg)/as.numeric(final_sumstats$standard_error.tg)

Samplesize <-c(median(final_sumstats$sample_size.fiadjbmi), median(final_sumstats$sample_size.hdl), median(final_sumstats$sample_size.tg))

#And now let's get this done for the SNPs that we need the triangulation for:

Sumstats <- final_sumstats[which(final_sumstats$chr_pos%in%triangulated_snps$chr_pos),]

Sumstats_ordered <- Sumstats[order(match(Sumstats$chr_pos, triangulated_snps$chr_pos)),]

length(which(Sumstats_ordered$chr_pos == triangulated_snps$chr_pos))

############################
#Preparing data for CPASSOC#
############################

Sumstat <- Sumstats_ordered %>%
  dplyr::select(z_score.fiadjbmi, z_score.hdl, z_score.tg)

#estimate parameters of gamma distribution.
para = EstimateGamma(N = 1E6, Samplesize, corMatrix) #with this amount of N it loops enough times to converge

#1.2795933 1.7898744 0.1075788 -> these are the parameters approximately. They converge quite nicely.

Test.shet<-SHet(Sumstat,Samplesize,corMatrix)
p.shet = pgamma(q = Test.shet-para[3], shape = para[1], scale = para[2], lower.tail = F)
  
#And now let's load the data in our SNPs:

triangulated_snps$shet <- Test.shet
triangulated_snps$pval.shet <- p.shet

##########################################################################
#And now let's filter the data for those that are significant for CPASSOC#
##########################################################################

fwrite(triangulated_snps, "output/4_ir_loci_discovery/2_ir_variants/all_fiadjbmi_hdl_tg_dwls_snps_w_cpassoc.txt")

gw_variants <- triangulated_snps[which(as.numeric(triangulated_snps$pval.shet) < 5e-08),] #80,285

#Let's load the common variant gwas:

raw_gsem <- fread("output/2_mv_gwas/fiadjbmi_hdl_tg_dwls.txt")
11gsem <- fread("output/2_mv_gwas/fiadjbmi_hdl_tg_common_dwls_curated.txt")

#Let's take the matching variants:

gsem_match <- gsem[which(gsem$variant%in%gw_variants$variant),]

length(which(gsem_match$p_value < 5e-08)) #3054
length(which(gsem_match$p_value < 1e-06)) #5273

