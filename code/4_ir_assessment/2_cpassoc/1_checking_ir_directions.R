##############
#INTRODUCTION#
##############

#This code compares the results for IR with and without FIadjBMI.

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)
library(GenomicSEM)
library(vautils)
library(vroom)
library(tidyverse)
library(ggplot2)
library(ggrepel)

##############################
#Let's obtain the LDSC output#
##############################

path_2_input <- "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/"

setwd(path_2_input)

#Let's load the FIadjBMI-HDL-TG data:

fiadjbmi_hdl_tg_mv_ss <- fread("output/2_mv_gwas/fiadjbmi_hdl_tg_common_dwls_curated.txt")
fiadjbmi_hdl_tg_ind <- fread("output/3_ld_clumping/2_clumped_data/clumped_data/fiadjbmi_hdl_tg_common_dwls_clumped.txt") #352

#Let's load the FIadjBMI, FI, HDL and TG original data aligned from CPASSOC and do proper comparisons...

cpassoc_fiadjbmi_hdl_tg <- fread("output/4_ir_loci_discovery/1_cpassoc/merged_fiadjbmi_hdl_tg_sumstats_4_cpassoc.txt")

#Let's use TG/HDL and fill all the info here:

tg_hdl_ratio <- fread("output/1_curated_gwas/tg_hdl_ratio_curated.txt")

######################################################################################################
#STEP 1: let's align the 352 independent loci to the IR-increasing allele, it will make things easier#
######################################################################################################

ir_variants <- fiadjbmi_hdl_tg_mv_ss[which(fiadjbmi_hdl_tg_mv_ss$variant%in%fiadjbmi_hdl_tg_ind$rsid),]

ir_pos <- ir_variants 

new_a1 <- ifelse(as.numeric(ir_pos$beta) < 0, ir_pos$other_allele, ir_pos$effect_allele)
new_a2 <- ifelse(as.numeric(ir_pos$beta) < 0, ir_pos$effect_allele, ir_pos$other_allele)
new_beta <- ifelse(as.numeric(ir_pos$beta) < 0, ir_pos$beta*(-1), as.numeric(ir_pos$beta))

ir_pos$effect_allele <- new_a1
ir_pos$other_allele <- new_a2
ir_pos$beta <- new_beta

###############################################################################################
#Let's merge the data so that we can fill in with the sumstats of FIadjBMI, TG, HDL and TG/HDL#
###############################################################################################

#Let's first match the cpassoc data:

cpassoc_match <- cpassoc_fiadjbmi_hdl_tg[which(cpassoc_fiadjbmi_hdl_tg$variant%in%ir_pos$variant),] #all of them

cpassoc_ordered <- cpassoc_match[order(match(cpassoc_match$variant, ir_pos$variant)),]
 
length(which(cpassoc_ordered$variant == ir_pos$variant)) #352
length(which(cpassoc_ordered$effect_allele == ir_pos$effect_allele)) #only 282!! Let's chance the data for those that do not match:

new_a1 <- ifelse(cpassoc_ordered$effect_allele != ir_pos$effect_allele, cpassoc_ordered$other_allele, cpassoc_ordered$effect_allele)
new_a2 <- ifelse(cpassoc_ordered$effect_allele != ir_pos$effect_allele, cpassoc_ordered$effect_allele, cpassoc_ordered$other_allele)

new_fiadjbmi <- ifelse(cpassoc_ordered$effect_allele != ir_pos$effect_allele, as.numeric(cpassoc_ordered$beta.fiadjbmi)*(-1), as.numeric(cpassoc_ordered$beta.fiadjbmi))
new_hdl <- ifelse(cpassoc_ordered$effect_allele != ir_pos$effect_allele, as.numeric(cpassoc_ordered$beta.hdl)*(-1), as.numeric(cpassoc_ordered$beta.hdl))
new_tg <- ifelse(cpassoc_ordered$effect_allele != ir_pos$effect_allele, as.numeric(cpassoc_ordered$beta.tg)*(-1), as.numeric(cpassoc_ordered$beta.tg))

cpassoc_aligned <- cpassoc_ordered

cpassoc_aligned$effect_allele <- new_a1
cpassoc_aligned$other_allele <- new_a2

cpassoc_aligned$beta.fiadjbmi <- new_fiadjbmi
cpassoc_aligned$beta.hdl <- new_hdl
cpassoc_aligned$beta.tg <- new_tg

length(which(cpassoc_aligned$variant == ir_pos$variant)) #352
length(which(cpassoc_aligned$effect_allele == ir_pos$effect_allele)) #352

#######################################################
#Let's do the same with TG/HDL and add the information#
#######################################################

tg_hdl_match <- tg_hdl_ratio[which(tg_hdl_ratio$rs_id%in%ir_pos$variant),] #all of them

tg_hdl_ordered <- tg_hdl_match[order(match(tg_hdl_match$rs_id, ir_pos$variant)),]

length(which(tg_hdl_ordered$rs_id == ir_pos$variant)) #352
length(which(tg_hdl_ordered$effect_allele == ir_pos$effect_allele)) #352!!! All of them already ordered!!!

################################
#Let's add all this information#
################################

ir_pos$beta.fiadjbmi <- cpassoc_aligned$beta.fiadjbmi
ir_pos$standard_error.fiadjbmi <- cpassoc_aligned$standard_error.fiadjbmi
ir_pos$p_value.fiadjbmi <- cpassoc_aligned$p_value.fiadjbmi
ir_pos$sample_size.fiadjbmi <- cpassoc_aligned$sample_size.fiadjbmi

ir_pos$beta.hdl<- cpassoc_aligned$beta.hdl
ir_pos$standard_error.hdl <- cpassoc_aligned$standard_error.hdl
ir_pos$p_value.hdl <- cpassoc_aligned$p_value.hdl
ir_pos$sample_size.hdl <- cpassoc_aligned$sample_size.hdl

ir_pos$beta.tg <- cpassoc_aligned$beta.tg
ir_pos$standard_error.tg <- cpassoc_aligned$standard_error.tg
ir_pos$p_value.tg <- cpassoc_aligned$p_value.tg
ir_pos$sample_size.tg <- cpassoc_aligned$sample_size.tg

ir_pos$beta.tg_hdl <- tg_hdl_ordered$beta
ir_pos$standard_error.tg_hdl <- tg_hdl_ordered$standard_error
ir_pos$p_value.tg_hdl <- tg_hdl_ordered$p_value
ir_pos$sample_size.tg_hdl <- tg_hdl_ordered$n

#Let's add just one column adding the info whether the data follows the expected pattern of associations:

ir_pos$ir_direction <- ifelse(as.numeric(ir_pos$beta.fiadjbmi) > 0 & as.numeric(ir_pos$beta.hdl) < 0 & as.numeric(ir_pos$beta.tg) > 0 | 
                                as.numeric(ir_pos$beta.fiadjbmi) < 0 & as.numeric(ir_pos$beta.hdl) > 0 & as.numeric(ir_pos$beta.tg) < 0, "yes", "no")

#####################
#Let's save the data#
#####################

dir.create("output/4_ir_loci_discovery/2_ir_variants")

fwrite(ir_pos, "output/4_ir_loci_discovery/2_ir_variants/282_ir_fiadjbmi_hdl_tg_dwls_snps.txt")
