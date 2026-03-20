##############
#INTRODUCTION#
##############

#This code checks the actual effective sample size for DWLS and MLR methods and checks them out.

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)

###################################
#Loading curated datasets for both#
###################################

path_2_input <- "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/"

setwd(path_2_input)

ir_dwls <- fread("output/2_mv_gwas/fiadjbmi_hdl_tg_common_dwls_curated.txt")
ir_mlr <- fread("output/2_mv_gwas/fiadjbmi_hdl_tg_ml.txt")

#################################
#Let's for a moment clean the ML#
#################################

head(ir_mlr)

ir_mlr <- ir_mlr[which(as.numeric(ir_mlr$Q_pval) > 5e-08),]
ir_mlr <- ir_mlr[which(as.numeric(ir_mlr$fail) == 0 & as.numeric(ir_mlr$warning) == 0),]

#Importantly, just take into account the SNPs that pass the filers:

ir_mlr <- ir_mlr[which(as.numeric(ir_mlr$MAF) >= 0.1 & as.numeric(ir_mlr$MAF) <= 0.4),]

individual_neff_num <- (ir_mlr$Z_Estimate/ir_mlr$est)^2
individual_neff_den <- 2*(ir_mlr$MAF)*(1-ir_mlr$MAF)
neff <- as.numeric(individual_neff_num)/as.numeric(individual_neff_den)
sum(neff)/length(neff) #520.9780!! What??

######################
#Let's do it for DWLS#
######################

ir_dwls <- ir_dwls[which(as.numeric(ir_dwls$minimum_allele_frequency) >= 0.1 & as.numeric(ir_dwls$minimum_allele_frequency) <= 0.4),]

dwls_z <- as.numeric(ir_dwls$beta)/as.numeric(ir_dwls$standard_error)

individual_neff_num <- (dwls_z/ir_dwls$beta)^2
individual_neff_den <- 2*(ir_dwls$minimum_allele_frequency)*(1-ir_dwls$minimum_allele_frequency)
neff <- as.numeric(individual_neff_num)/as.numeric(individual_neff_den)
sum(neff)/length(neff) #520.9780!! What??
