##############
#INTRODUCTION#
##############

#This code performs DWLS common factor

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)
library(GenomicSEM)

##############################
#Let's obtain the LDSC output#
##############################

path_2_input <- "output/2_mv_gwas"

setwd(path_2_input)

LDSCoutput <- readRDS("fiadjbmi_hdl_tg.rds")

all_sumstats <- fread("fiadjbmi_hdl_tg_4_mvgwas.txt")

########################
#Let's run the analysis#
########################

common_factor <- commonfactorGWAS(covstruc = LDSCoutput, SNPs = all_sumstats, estimation = "DWLS", cores = NULL, toler = 1e-100, SNPSE = 0.0005, parallel = FALSE, GC="conserv")

fwrite(common_factor, "2_mv_gwas/fiadjbmi_hdl_tg_dwls.txt")
