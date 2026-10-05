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

path_2_input <- "output/2_models/fiadjbmi_hdl_tg"

setwd(path_2_input)

LDSCoutput <- readRDS("3_gc/fiadjbmi_hdl_tg/fiadjbmi_hdl_tg.rds")

all_sumstats <- fread("../../../3_mv_gwas/fiadjbmi_hdl_tg/fiadjbmi_hdl_tg_4_mvgwas.txt")

########################
#Let's run the analysis#
########################

common_factor <- commonfactorGWAS(covstruc = LDSCoutput, SNPs = all_sumstats, estimation = "DWLS", cores = NULL, toler = 1e-100, SNPSE = 0.0005, parallel = FALSE, GC="conserv")

fwrite(common_factor, "../../../3_mv_gwas/fiadjbmi_hdl_tg/fiadjbmi_hdl_tg_dwls.txt")
