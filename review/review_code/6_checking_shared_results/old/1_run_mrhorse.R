##############
#INTRODUCTION#
##############

#This code runs MVMR with FIadjBMI for WHRadjBMI and our PCOS and other mediator traits.
#First it runs MR-HORSE and, then, it adjusted for FIadjBMI.

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)
library(TwoSampleMR)
source("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/mrhorse-main/mr_horse.R")

###################
#Loading functions#
###################

trait_aligner <- function(query_ss, other_ss){
  
  #We are gonna put the example so that we can know what happens:
  
  #whradjbmi_ss <- exp_df_found_pos
  #other_ss <- hdl_005
  #query_ss <- triangulated_fi_all
  #other_ss <- ebmd_ss
  
  #STEP 0: let's run the matching with TwoSampleMR so we need the data:
  
  #Let's first check if we have the effect_allele_frequency column:
  
  check <- which(colnames(query_ss) == "effect_allele_frequency")
  
  if(is_empty(check)){
    
    query_ss$effect_allele_frequency <- NA
    
  }
  
  #Now we can proceed
  
  exposure <- query_ss %>%
    select(chr_pos, chromosome, base_pair_location, effect_allele, other_allele,
           effect_allele_frequency, beta, standard_error, p_value, sample_size, variant)
  
  colnames(exposure) <- c("SNP", "chr.exposure", "pos.exposure", "effect_allele.exposure",
                          "other_allele.exposure", "eaf.exposure", "beta.exposure", "se.exposure",
                          "pval.exposure", "samplesize.exposure", "rsid.exposure")
  
  exposure$exposure <- "fat_distr"
  exposure$id.exposure <- "fat_distr"
  
  #Now with the outcome:
  
  check <- which(colnames(other_ss) == "effect_allele_frequency")
  
  if(is_empty(check)){
    
    other_ss$effect_allele_frequency <- NA
    
  }
  
  #Now we can proceed
  
  outcome <- other_ss %>%
    select(chr_pos, chromosome, base_pair_location, effect_allele, other_allele, effect_allele_frequency, beta, standard_error, p_value, sample_size, variant)
  
  colnames(outcome) <- c("SNP", "chr.outcome", "pos.outcome", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome", "pval.outcome", "samplesize.outcome", "rsid.outcome")
  
  outcome$outcome <- "Outcome"
  outcome$id.outcome <- "Outcome"
  
  #STEP 1: match the data. This will probably fail with tri-allelic SNPs. Here I think we are gonna be OK...
  
  merged_df <- harmonise_data(exposure, outcome, action=3)
  
  merged_df <- merged_df[which(merged_df$remove == FALSE),] #removing incompatible alleles
  #merged_df <- merged_df[which(merged_df$mr_keep),] #removing incompatible alleles
  
  #I checked that all is working great. whradjbmi betas is positive. The rest is great.
  
  #STEP 3: reorganize the dataframe, cuz we need something clean:
  
  other_ss_aligned <- merged_df %>%
    select("rsid.outcome", "chr.outcome", "pos.outcome", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome", "pval.outcome", "samplesize.outcome", "SNP")
  
  colnames(other_ss_aligned) <- c("variant", "chromosome", "base_pair_location",  "effect_allele", "other_allele", "effect_allele_frequency", "beta", "standard_error", "p_value", "sample_size", "chr_pos")
  
  return(other_ss_aligned)
  
}

organizing_common_data <- function(ss_df, common_snps){
  
  ss_end <- ss_df[which(ss_df$chr_pos%in%common_snps),]
  ss_ordered <- ss_end[order(as.numeric(ss_end$chromosome), as.numeric(ss_end$base_pair_location)),]
  
  return(ss_ordered)
  
}

compare_mr_horse_mvmr_no_fi_or_tg_hdl <- function(fat_pcos_df, fat_fiadjbmi_df, exposure_, mediator_, outcome_, ref){
  
  #STEP 0: let's see what the hell is going on:
  
  #fat_pcos_df <- bmi_pcos_broad
  #fat_fi_df <- fi
  
  #STEP 1: compute mr_horse analyses for fat_pcos_df:
  
  fat_pcos_4_mr_horse <- fat_pcos_df %>%
    dplyr::select(beta.exposure, se.exposure, beta.outcome, se.outcome)
  
  colnames(fat_pcos_4_mr_horse) <- c("betaX", "betaXse", "betaY", "betaYse")
  
  res_mr_horse <- mr_horse(fat_pcos_4_mr_horse)
  
  uni_trait_res <- res_mr_horse$MR_Estimate
  uni_trait_res$exposure <- exposure_
  uni_trait_res$mediator <- "None"
  uni_trait_res$outcome <- outcome_
  uni_trait_res$analysis <- "MR replication"
  uni_trait_res$Parameter <- NA
  
  #STEP 2: rerun analyses BUT with the variants that match with FI, FIadjBMI and TG/HDL
  
  fat_match=ref[which(ref$chr_pos%in%fat_pcos_df$chr_pos.exposure),]
  
  #Let's use these to align the FI, FIadjBMI and TG/HDL data!
  
  fiadjbmi_match=trait_aligner(fat_match, fat_fiadjbmi_df)

  common_snps = Reduce(intersect, list(fat_pcos_df$chr_pos.exposure, fiadjbmi_match$chr_pos))
  
  fat_match=fat_match[which(fat_match$chr_pos%in%common_snps),]
  pcos_match = fat_pcos_df[which(fat_pcos_df$chr_pos.exposure%in%common_snps),]
  #fi_match=fi_match[which(fi_match$chr_pos%in%common_snps),]
  fiadjbmi_match=fiadjbmi_match[which(fiadjbmi_match$chr_pos%in%common_snps),]
  #tg_hdl_match=tg_hdl_match[which(tg_hdl_match$chr_pos%in%common_snps),]
  
  #Let's organize the data:
  
  #fi_match = fi_match[order(match(fi_match$chr_pos, fat_match$chr_pos)),]
  pcos_match = pcos_match[order(match(pcos_match$chr_pos.exposure, fat_match$chr_pos)),]
  fiadjbmi_match = fiadjbmi_match[order(match(fiadjbmi_match$chr_pos, fat_match$chr_pos)),]
  #tg_hdl_match = tg_hdl_match[order(match(tg_hdl_match$chr_pos, fat_match$chr_pos)),]
  
  #print(length(which(fat_match$chr_pos == fi_match$chr_pos)))
  print(length(which(fat_match$chr_pos == pcos_match$chr_pos.exposure)))
  print(length(which(fat_match$chr_pos == fiadjbmi_match$chr_pos)))
  #print(length(which(fat_match$chr_pos == tg_hdl_match$chr_pos)))
  
  fat_match$beta.pcos=pcos_match$beta.outcome
  fat_match$standard_error.pcos=pcos_match$se.outcome
  
  #fat_match$beta.fi = fi_match$beta
  #fat_match$standard_error.fi = fi_match$standard_error
  
  fat_match$beta.fiadjbmi = fiadjbmi_match$beta
  fat_match$standard_error.fiadjbmi = fiadjbmi_match$standard_error
  
  #fat_match$beta.tg_hdl = tg_hdl_match$beta
  #fat_match$standard_error.tg_hdl = tg_hdl_match$standard_error
  
  fat_pcos_4_mvmr_horse <- fat_match %>%
    dplyr::select(beta, standard_error, beta.pcos, standard_error.pcos, beta.fiadjbmi, standard_error.fiadjbmi)
  
  colnames(fat_pcos_4_mvmr_horse) <- c("betaX1", "betaX1se", "betaY", "betaYse", "betaX2", "betaX2se")
  
  res_mvmr_horse <- mvmr_horse(as.data.frame(fat_pcos_4_mvmr_horse))
  
  multi_trait_res <- res_mvmr_horse$MR_Estimate[1,] #we are taking the association for WHRadjBMI only
  multi_trait_res$exposure <- exposure_
  multi_trait_res$mediator <- "FI"
  multi_trait_res$outcome <- outcome_
  multi_trait_res$analysis <- "MVMR"
  
  #Let's combine the results now:
  
  final_df <- rbind(uni_trait_res, multi_trait_res)
  
  return(list(final_df, res_mvmr_horse))
  
}

#######################
#Loading exposure data#
#######################

path_2_input <- "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/PCOS_2026/"

setwd(path_2_input)

#Third for BMI+WHR+:

whradjbmi_pcos_doctor <- readRDS("output/2_replicating_liu_et_al/2_mr/2_mr/whradjbmi/doctor/instruments_after_outlier_removal_df.RDS")
whradjbmi_pcos_broad <- readRDS("output/2_replicating_liu_et_al/2_mr/2_mr/whradjbmi/broad/instruments_after_outlier_removal_df.RDS")
whradjbmi_pcos_consortium <- readRDS("output/2_replicating_liu_et_al/2_mr/2_mr/whradjbmi/consortium/instruments_after_outlier_removal_df.RDS")
whradjbmi_pcos_meta_analysis <- readRDS("output/2_replicating_liu_et_al/2_mr/2_mr/whradjbmi/meta_analysis/instruments_after_steiger_df.RDS")
whradjbmi_pcos_venkatesh <- readRDS("output/2_replicating_liu_et_al/2_mr/2_mr/whradjbmi/venkatesh/instruments_after_steiger_df.RDS")
whradjbmi_pcos_adj_age <- readRDS("output/2_replicating_liu_et_al/2_mr/2_mr/whradjbmi/age_adjusted/instruments_after_steiger_df.RDS")
whradjbmi_pcos_adj_age_bmi <- readRDS("output/2_replicating_liu_et_al/2_mr/2_mr/whradjbmi/age_and_bmi_adjusted/instruments_after_outlier_removal_df.RDS")

#Let's now take the FI, FIadjBMI and TG/HDL data:

fiadjbmi <- fread("../IR_GSEM_2023//output/1_curated_data/fiadjbmi_curated.txt")

###########################################################
#Let's align the data first, if not we are gonna go crazy!#
###########################################################

#Let's get the BMI, WHR and WHRadjBMI original data, cuz it has the right columns to align the data in a relatively easy way, compared to the MR results!

whradjbmi=fread("output/1_curated_data/whradjbmi_curated_female.txt")

#Let's do it for the analyses that we are the most interested in: WHRadjBMI first:

whradjbmi_pos=whradjbmi

new_a1=ifelse(as.numeric(whradjbmi$beta) < 0, whradjbmi$other_allele, whradjbmi$effect_allele)
new_a2=ifelse(as.numeric(whradjbmi$beta) < 0, whradjbmi$effect_allele, whradjbmi$other_allele)
new_beta=ifelse(as.numeric(whradjbmi$beta) < 0, as.numeric(whradjbmi$beta)*(-1), as.numeric(whradjbmi$beta))
new_eaf=ifelse(as.numeric(whradjbmi$beta) < 0, 1-as.numeric(whradjbmi$effect_allele_frequency), as.numeric(whradjbmi$effect_allele_frequency))

whradjbmi_pos$effect_allele = new_a1
whradjbmi_pos$other_allele = new_a2
whradjbmi_pos$beta = new_beta
whradjbmi_pos$effect_allele_frequency = new_eaf

#Worked amazingly!

#######################################
#Let's first get the data for BMI+WHR-#
#######################################

whradjbmi_fiadjbmi_doctor <- compare_mr_horse_mvmr_no_fi_or_tg_hdl(fat_pcos_df = whradjbmi_pcos_doctor,
                                             fat_fiadjbmi_df=fiadjbmi,
                                             exposure_ = "WHRadjBMI", 
                                             mediator_ = "IR", 
                                             outcome_ = "PCOS (FinnGen)", 
                                             ref=whradjbmi_pos)

whradjbmi_fiadjbmi_broad <- compare_mr_horse_mvmr_no_fi_or_tg_hdl(fat_pcos_df = whradjbmi_pcos_broad,
                                            fat_fiadjbmi_df=fiadjbmi,
                                            exposure_ = "WHRadjBMI", 
                                            mediator_ = "IR", 
                                             outcome_ = "PCOS (Broad)", 
                                            ref=whradjbmi_pos)

whradjbmi_fiadjbmi_consortium <- compare_mr_horse_mvmr_no_fi_or_tg_hdl(fat_pcos_df = whradjbmi_pcos_consortium,
                                                 fat_fiadjbmi_df=fiadjbmi,
                                                 exposure_ = "WHRadjBMI", 
                                                 mediator_ = "IR", 
                                             outcome_ = "PCOS (Consortium)", 
                                             ref=whradjbmi_pos)

whradjbmi_fiadjbmi_meta_analysis <- compare_mr_horse_mvmr_no_fi_or_tg_hdl(fat_pcos_df = whradjbmi_pcos_meta_analysis,
                                                    fat_fiadjbmi_df=fiadjbmi,
                                                    exposure_ = "WHRadjBMI", 
                                                    mediator_ = "IR", 
                                             outcome_ = "PCOS (Day et al)", 
                                             ref=whradjbmi_pos)

whradjbmi_fiadjbmi_venkatesh <- compare_mr_horse_mvmr_no_fi_or_tg_hdl(fat_pcos_df = whradjbmi_pcos_venkatesh,
                                                fat_fiadjbmi_df=fiadjbmi,
                                                exposure_ = "WHRadjBMI", 
                                                mediator_ = "IR", 
                                                outcome_ = "PCOS (Venkatesh et al)", 
                                                ref=whradjbmi_pos)

whradjbmi_fiadjbmi_age <- compare_mr_horse_mvmr_no_fi_or_tg_hdl(fat_pcos_df = whradjbmi_pcos_adj_age,
                                          fat_fiadjbmi_df=fiadjbmi,
                                          exposure_ = "WHRadjBMI", 
                                          mediator_ = "IR", 
                                          outcome_ = "PCOS (Tyrmi et al)", 
                                          ref=whradjbmi_pos)

whradjbmi_fiadjbmi_age_bmi <- compare_mr_horse_mvmr_no_fi_or_tg_hdl(fat_pcos_df = whradjbmi_pcos_adj_age_bmi,
                                              fat_fiadjbmi_df=fiadjbmi,
                                              exposure_ = "WHRadjBMI", 
                                              mediator_ = "IR", 
                                              outcome_ = "PCOSadjBMI", 
                                              ref=whradjbmi_pos)

#Let's save this data just in case:

dir.create("output/3_prs_and_mr/3_mvmr")

saveRDS(whradjbmi_fiadjbmi_doctor, "output/3_prs_and_mr/3_mvmr/whradjbmi_fiadjbmi_pcos.RDS") 
saveRDS(whradjbmi_fiadjbmi_broad, "output/3_prs_and_mr/3_mvmr/whradjbmi_fiadjbmi_pcos_broad.RDS") 
saveRDS(whradjbmi_fiadjbmi_consortium, "output/3_prs_and_mr/3_mvmr/whradjbmi_fiadjbmi_pcos_consortium.RDS") 
saveRDS(whradjbmi_fiadjbmi_meta_analysis, "output/3_prs_and_mr/3_mvmr/whradjbmi_fiadjbmi_pcos_meta_analysis.RDS") 
saveRDS(whradjbmi_fiadjbmi_venkatesh, "output/3_prs_and_mr/3_mvmr/whradjbmi_fiadjbmi_venkatesh_meta_analysis.RDS") 
saveRDS(whradjbmi_fiadjbmi_age, "output/3_prs_and_mr/3_mvmr/whradjbmi_fiadjbmi_pcos_tyrmi.RDS") 
saveRDS(whradjbmi_fiadjbmi_age_bmi, "output/3_prs_and_mr/3_mvmr/whradjbmi_fiadjbmi_pcosadjbmi.RDS") 

#Let's check the results:

whradjbmi_fiadjbmi_doctor[[1]]
# Estimate    SD 2.5% quantile 97.5% quantile  Rhat  exposure mediator        outcome       analysis Parameter
# 1    0.339 0.093          0.16          0.525 1.000 WHRadjBMI     None PCOS (FinnGen) MR replication      <NA>
#   2    0.259 0.144         -0.02          0.550 1.005 WHRadjBMI       FI PCOS (FinnGen)           MVMR  theta[1]

whradjbmi_fiadjbmi_broad[[1]]

#Estimate    SD 2.5% quantile 97.5% quantile  Rhat  exposure mediator      outcome       analysis Parameter
#1    0.138 0.069         0.002          0.273 1.000 WHRadjBMI     None PCOS (Broad) MR replication      <NA>
#  2    0.102 0.106        -0.106          0.311 1.001 WHRadjBMI       FI PCOS (Broad)           MVMR  theta[1]

whradjbmi_fiadjbmi_consortium[[1]]

#Estimate    SD 2.5% quantile 97.5% quantile  Rhat  exposure mediator      outcome       analysis Parameter
#1    0.138 0.069         0.002          0.273 1.000 WHRadjBMI     None PCOS (Broad) MR replication      <NA>
#  2    0.102 0.106        -0.106          0.311 1.001 WHRadjBMI       FI PCOS (Broad)           MVMR  theta[1]

whradjbmi_fiadjbmi_meta_analysis[[1]]

#Estimate    SD 2.5% quantile 97.5% quantile  Rhat  exposure mediator          outcome       analysis Parameter
#1    0.002 0.094        -0.189          0.186 1.000 WHRadjBMI     None PCOS (Day et al) MR replication      <NA>
#  2    0.038 0.151        -0.260          0.339 1.001 WHRadjBMI       FI PCOS (Day et al)           MVMR  theta[1]

whradjbmi_fiadjbmi_venkatesh[[1]]

#Estimate    SD 2.5% quantile 97.5% quantile  Rhat  exposure mediator                outcome       analysis Parameter
#1    0.236 0.085         0.065          0.397 1.000 WHRadjBMI     None PCOS (Venkatesh et al) MR replication      <NA>
#2   -0.055 0.146        -0.340          0.231 1.001 WHRadjBMI       FI PCOS (Venkatesh et al)           MVMR  theta[1]

whradjbmi_fiadjbmi_age[[1]]

# Estimate    SD 2.5% quantile 97.5% quantile  Rhat  exposure mediator            outcome       analysis Parameter
# 1    0.228 0.076         0.080          0.376 1.003 WHRadjBMI     None PCOS (Tyrmi et al) MR replication      <NA>
#   2    0.146 0.122        -0.093          0.386 1.001 WHRadjBMI       FI PCOS (Tyrmi et al)           MVMR  theta[1]

whradjbmi_fiadjbmi_age_bmi[[1]]

# Estimate    SD 2.5% quantile 97.5% quantile  Rhat  exposure mediator    outcome       analysis Parameter
# 1    0.327 0.093         0.146          0.508 1.001 WHRadjBMI     None PCOSadjBMI MR replication      <NA>
#   2    0.312 0.142         0.042          0.587 1.006 WHRadjBMI       FI PCOSadjBMI           MVMR  theta[1]