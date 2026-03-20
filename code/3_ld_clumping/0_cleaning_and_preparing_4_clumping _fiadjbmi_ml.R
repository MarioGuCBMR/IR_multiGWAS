##############
#INTRODUCTION#
##############

#This code gets the significant variants for each trait.

#################
#Loading library#
#################

library(data.table)
library(tidyverse)

###################
#Loading functions#
###################

formatting_data_4_clumping <- function(ss_gw){
  #This function formats the data for the clumping.
  
  #STEP 2: now make the SNP column which is SNP:A1:A2
  
  rsid <- ifelse(str_detect(ss_gw$variant, "rs") == TRUE, ss_gw$variant, ss_gw$chr_pos)
  SNP <- paste(rsid, ":", ss_gw$other_allele, ":", ss_gw$effect_allele, sep = "")
  
  #Now we select the columns that we want:
  
  ss_gw$SNP <- SNP
  ss_gw$rsid <- rsid
  
  final_df <- ss_gw %>%
    select(SNP, p_value, chr_pos, effect_allele, other_allele, rsid)
  
  colnames(final_df) <- c("SNP", "pval", "chr_pos", "effect_allele", "other_allele", "rsid")
  
  return(final_df)
  
}

##############
#Loading data#
##############

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/")

fiadjbmi_hdl_tg <- fread("output/2_mv_gwas/fiadjbmi_hdl_tg_ml.txt")

#######################################
#First let's filter for weird variants#
#######################################

summary(fiadjbmi_hdl_tg$MAF) #perfect MAF filter

which(as.numeric(fiadjbmi_hdl_tg$CHR) == 6 & as.numeric(fiadjbmi_hdl_tg$BP) >= 26000000 & as.numeric(fiadjbmi_hdl_tg$BP) <= 34000000) #none!!

#Perfect, then let's filter for those that did not work!!

fiadjbmi_hdl_tg_clean <- fiadjbmi_hdl_tg[which(fiadjbmi_hdl_tg$fail == 0 & fiadjbmi_hdl_tg$warning == 0),] #perfect. 

#And now let's separate those that are not from the common factor!!

fiadjbmi_hdl_tg_common <- fiadjbmi_hdl_tg_clean[which(as.numeric(fiadjbmi_hdl_tg_clean$Q_pval) > 5e-08)]
fiadjbmi_hdl_tg_het <- fiadjbmi_hdl_tg_clean[which(as.numeric(fiadjbmi_hdl_tg_clean$Q_pval) < 5e-08)] #147,879

#Then, let's change the columns...

colnames(fiadjbmi_hdl_tg_common) <- c("variant", "chromosome", "base_pair_location", "minimum_allele_frequency", "effect_allele", "other_allele", "i", "lhs", "op", "rhs", "beta", "standard_error", "z_score", "p_value", "Q", "Q_df", "Q_p_value", "fail", "warning")
colnames(fiadjbmi_hdl_tg_het) <- c("variant", "chromosome", "base_pair_location", "minimum_allele_frequency", "effect_allele", "other_allele", "i", "lhs", "op", "rhs", "beta", "standard_error", "z_score", "p_value", "Q", "Q_df", "Q_p_value", "fail", "warning")

#Finally, let's rearrange the data:

fiadjbmi_hdl_tg_common_clean <- fiadjbmi_hdl_tg_common %>%
  select("variant", "chromosome", "base_pair_location", "minimum_allele_frequency", "effect_allele", "other_allele", "beta", "standard_error", "p_value")

fiadjbmi_hdl_tg_het_clean <- fiadjbmi_hdl_tg_het %>%
  select("variant", "chromosome", "base_pair_location", "minimum_allele_frequency", "effect_allele", "other_allele", "beta", "standard_error", "p_value")

#First for the common variants:

fiadjbmi_hdl_tg_common_clean$chr_pos <- paste("chr", fiadjbmi_hdl_tg_common_clean$chromosome, ":", fiadjbmi_hdl_tg_common_clean$base_pair_location, sep = "")
fiadjbmi_hdl_tg_common_clean$sample_size <- mean(1/((2*fiadjbmi_hdl_tg_clean$MAF*(1-fiadjbmi_hdl_tg_clean$MAF))*fiadjbmi_hdl_tg_clean$se_c^2))

#Now for the heterogeneic variants:

fiadjbmi_hdl_tg_het_clean$chr_pos <- paste("chr", fiadjbmi_hdl_tg_het_clean$chromosome, ":", fiadjbmi_hdl_tg_het_clean$base_pair_location, sep = "")
fiadjbmi_hdl_tg_het_clean$sample_size <- mean(1/((2*fiadjbmi_hdl_tg_clean$MAF*(1-fiadjbmi_hdl_tg_clean$MAF))*fiadjbmi_hdl_tg_clean$se_c^2))

#Now let's save the data...

fwrite(fiadjbmi_hdl_tg_common_clean, "output/2_mv_gwas/fiadjbmi_hdl_tg_common_ml_curated.txt")
fwrite(fiadjbmi_hdl_tg_het_clean, "output/2_mv_gwas/fiadjbmi_hdl_tg_het_ml_curated.txt")

########################################
#Now we set up the data for ld-clumping#
########################################

fiadjbmi_hdl_tg_common_4_clumping <- fiadjbmi_hdl_tg_common_clean[which(as.numeric(fiadjbmi_hdl_tg_common_clean$p_value) < 1e-06),] #WE DO NOT HAVE ANY