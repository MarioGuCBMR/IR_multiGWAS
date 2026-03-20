##############
#INTRODUCTION#
##############

#This code performs the clustering of all IR variants with TG/HDL, but we will just only focus on the 83 novel in the downstream analysis.

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)
library(navmix)
library(TwoSampleMR)
library(hyprcoloc)

###################
#Loading functions#
###################

trait_aligner <- function(query_ss, other_ss){
  
  #We are gonna put the example so that we can know what happens:
  
  #fiadjbmi_ss <- exp_df_found_pos
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
  
  #I checked that all is working great. FIadjBMI betas is positive. The rest is great.
  
  #STEP 3: reorganize the dataframe, cuz we need something clean:
  
  other_ss_aligned <- merged_df %>%
    select("rsid.outcome", "chr.outcome", "pos.outcome", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome", "pval.outcome", "samplesize.outcome", "SNP")
  
  colnames(other_ss_aligned) <- c("variant", "chromosome", "base_pair_location",  "effect_allele", "other_allele", "effect_allele_frequency", "beta", "standard_error", "p_value", "sample_size", "chr_pos")
  
  return(other_ss_aligned)
  
}

organizing_common_data <- function(ss_df, common_df){
  
  ss_end <- ss_df[which(ss_df$chr_pos%in%common_snps),]
  ss_ordered <- ss_end[order(as.numeric(ss_end$chromosome), as.numeric(ss_end$base_pair_location)),]
  
  return(ss_ordered)
  
}


###############
#Load the data#
###############

path_2_input <- "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/"

setwd(path_2_input)

#First let's load the variants:

ir_variants <- fread("output/4_ir_loci_discovery/3_novel_variants/3_novel_variants/282_ir_variants_with_ir_and_t2d_comparisons_w_closest_genes.txt")

#Now let's load the data:

iradjbmi <- data.table::fread("output/2_mv_gwas/fiadjbmi_hdl_tg_common_dwls_curated.txt")
tg_hdl_ratio <- data.table::fread("output/1_curated_gwas/tg_hdl_ratio_curated.txt")

tg_hdl_ratio$variant <- tg_hdl_ratio$rs_id
tg_hdl_ratio$sample_size <- tg_hdl_ratio$n

#Finally, let's get the proxies:

#Finally, let's get the proxies #please go to next section regarding finding proxies to continue with this: N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/code/4_identify_proxies

proxies <- fread("output/4_ir_loci_discovery/4_proxies/proxies_build_37_and_38_4_282_variants.txt")

length(which(proxies$rsID%in%ir_variants$variant)) #282

proxies <- proxies[which(proxies$query_snp_rsid%in%ir_variants$variant),] #all of them!!

######################################
#STEP 1: let's load over all the SNPs#
######################################

for(i in seq(1, length(ir_variants$chr_pos))){
  
  print(i)
  
  #STEP 1.1: set up the loci that we are going to work with:
  
  chr <- ir_variants$chromosome[i]
  start <- ifelse(as.numeric(ir_variants$base_pair_location[i]-500000)<0, 0, as.numeric(ir_variants$base_pair_location[i]-500000))
  end <- as.numeric(ir_variants$base_pair_location[i]+500000)
  
  #STEP 1.2 find all SNPs in that region in BMI and align to the increasing allele:
  
  iradjbmi_tmp <- iradjbmi[which(iradjbmi$chromosome == chr & iradjbmi$base_pair_location >= start & iradjbmi$base_pair_location <= end),]
  
  new_beta <- ifelse(as.numeric(iradjbmi_tmp$beta) < 0, as.numeric(iradjbmi_tmp$beta)*(-1), as.numeric(iradjbmi_tmp$beta))
  new_a1 <- ifelse(as.numeric(iradjbmi_tmp$beta) < 0, iradjbmi_tmp$other_allele, iradjbmi_tmp$effect_allele)
  new_a2 <- ifelse(as.numeric(iradjbmi_tmp$beta) < 0, iradjbmi_tmp$effect_allele, iradjbmi_tmp$other_allele)
  new_eaf <- ifelse(as.numeric(iradjbmi_tmp$beta) < 0, 1-as.numeric(iradjbmi_tmp$effect_allele_frequency), as.numeric(iradjbmi_tmp$effect_allele_frequency))
  
  iradjbmi_tmp$effect_allele <- new_a1
  iradjbmi_tmp$other_allele <- new_a2
  iradjbmi_tmp$beta <- new_beta
  iradjbmi_tmp$effect_allele_frequency <- new_eaf
  
  #STEP 1.3: align the data for all the flavours:
  
  tg_hdl_ratio_2_iradjbmi <- trait_aligner(iradjbmi_tmp, tg_hdl_ratio) 
  #fiadjbmi_2_iradjbmi <- trait_aligner(iradjbmi_tmp, fiadjbmi) 
  
  #############################################################
  #STEP 3: let's just do it for the variants that are matching#
  #############################################################
  
  common_snps <- Reduce(intersect, list(iradjbmi_tmp$chr_pos, tg_hdl_ratio_2_iradjbmi$chr_pos))
  
  iradjbmi_matched <- iradjbmi_tmp[which(iradjbmi_tmp$chr_pos%in%common_snps),]
  iradjbmi_matched_ordered <- organizing_common_data(iradjbmi_matched, common_snps)
  
  tg_hdl_ratio_2_iradjbmi <- tg_hdl_ratio_2_iradjbmi[which(tg_hdl_ratio_2_iradjbmi$chr_pos%in%common_snps),]
  tg_hdl_ratio_2_iradjbmi <- organizing_common_data(tg_hdl_ratio_2_iradjbmi, common_snps)
  
  print(length(which(tg_hdl_ratio_2_iradjbmi$chr_pos == iradjbmi_matched_ordered$chr_pos)))
  print(length(which(tg_hdl_ratio_2_iradjbmi$effect_allele == iradjbmi_matched_ordered$effect_allele)))

  #####################################################
  #STEP 4: we are ready to set it all up for hyprcoloc#
  #####################################################
  
  betas <- as.matrix(as.data.frame(cbind(iradjbmi_matched_ordered$beta, tg_hdl_ratio_2_iradjbmi$beta)))
  ses <- as.matrix(as.data.frame(cbind(iradjbmi_matched_ordered$standard_error, tg_hdl_ratio_2_iradjbmi$standard_error)))
  traits <- c("IRadjBMI", "TG/HDL")
  colnames(betas) <- traits
  colnames(ses) <- traits
  rownames(betas) <- iradjbmi_matched_ordered$variant
  rownames(ses) <- iradjbmi_matched_ordered$variant
  
  res <- hyprcoloc::hyprcoloc(effect.est = betas, effect.se = ses, snp.id = iradjbmi_matched_ordered$variant, trait.names = traits)
  res_df <- res$results
  res_df$lead_snp <- ir_variants$variant[i]
  
  print(res_df)
  
  if(!(exists("coloc_df"))){
    
    coloc_df <- res_df
    
  } else {
    
    coloc_df <- rbind(coloc_df, res_df)
    
  }
  
}

#################################
#STEP 2: lets's save raw results#
#################################

dir.create("output/4_ir_loci_discovery/3_novel_variants/4_colocalization/")

fwrite(coloc_df, "output/4_ir_loci_discovery/3_novel_variants/4_colocalization/282_coloc_tg_hdl.txt")

#######################################################################################################
#Let's update the data with a merge so that we can have a supplementary table with all possible info!!#
#######################################################################################################

ir_variants_ <- fread("output/4_ir_loci_discovery/3_novel_variants/4_colocalization/352_variants_with_all_info_plus_fiadjbmi_coloc.txt")

ir_variants_$post_prob_coloc.tg_hdl <- NA
ir_variants_$driver_coloc.tg_hdl <- NA

for(i in seq(1, length(ir_variants_$variant))){
  
  #STEP 1: get the variant:
  
  rsid <- ir_variants_$variant[i]
  
  #STEP 2: get the info:
  
  coloc_tmp <- coloc_df[which(coloc_df$lead_snp == rsid),]
  
  if(is_empty(coloc_tmp$iteration)){
    
    next()
    
  } else {
    
    ir_variants_$post_prob_coloc.tg_hdl[i]  <- coloc_tmp$posterior_prob
    
    #Let's check whether the candidate proxy is a proxy of our signal and add that info:
    
    proxies_tmp <- proxies[which(proxies$query_snp_rsid == rsid),]
    
    candidate_snp <- coloc_tmp$candidate_snp
    
    check <- ifelse(candidate_snp%in%proxies_tmp$rsID, "(r2 >= 0.8)", "(r2 < 0.8)")
    
    candidate_snp <- paste(candidate_snp, check, sep = " ")
    
    ir_variants_$driver_coloc.tg_hdl[i]  <- candidate_snp
    
  }
  
}

#################
#This is done!!!#
#################

fwrite(ir_variants_, "output/4_ir_loci_discovery/3_novel_variants/4_colocalization/352_variants_with_all_info_plus_fiadjbmi_and_tg_hdl_coloc.txt")

#Let's do a small check:

novel <- ir_variants_[which(ir_variants_$reported_ir == "no"),]
table(novel$driver_coloc.tg_hdl)
