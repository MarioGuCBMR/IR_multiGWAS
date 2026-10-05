##############
#INTRODUCTION#
##############

#Let's run hyprcoloc for fiadjbmi, hdl and tg for all our IR loci

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

chr_parser = function(variant_id){
  
  #STEP 1: take the string and split it:
  
  vect_=unlist(str_split(variant_id, "_"))
  
  chr_=vect_[1]
  
  return(chr_)
  
}

pos_parser = function(variant_id){
  
  #STEP 1: take the string and split it:
  
  vect_=unlist(str_split(variant_id, "_"))
  
  pos_=vect_[2]
  
  return(pos_)
  
}

ref_parser = function(variant_id){
  
  #STEP 1: take the string and split it:
  
  vect_=unlist(str_split(variant_id, "_"))
  
  ref_=vect_[3]
  
  return(ref_)
  
}

alt_parser = function(variant_id){
  
  #STEP 1: take the string and split it:
  
  vect_=unlist(str_split(variant_id, "_"))
  
  alt_=vect_[4]
  
  return(alt_)
  
}


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
  
  #I checked that all is working great. fiadjbmi betas is positive. The rest is great.
  
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

#####################
#Let's load the data#
#####################

all_variants=fread("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/output/4_ir_loci_discovery/2_ir_variants/282_ir_fiadjbmi_hdl_tg_dwls_snps_w_cpassoc.txt") #353 hits!
all_variants = all_variants[which(all_variants$ir_direction == "yes"),]

#Let's load the data for all traits and assess overall co-localization

fiadjbmi=data.table::fread("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/output/1_curated_gwas/fiadjbmi_curated.txt")
hdl=data.table::fread("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/review/review_output/1_curated_data/hdl_2013_curated.txt")
tg=data.table::fread("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/review/review_output/1_curated_data/tg_2013_curated.txt")

#Let's just in case filter the fiadjbmi data to remove rare variants and MHC region

summary(as.numeric(fiadjbmi$effect_allele_frequency))
fiadjbmi = fiadjbmi[which(as.numeric(fiadjbmi$effect_allele_frequency) > 0.01 & as.numeric(fiadjbmi$effect_allele_frequency) < 0.99),]
summary(as.numeric(fiadjbmi$effect_allele_frequency))

mhc = fiadjbmi[which(fiadjbmi$chromosome == 6 & fiadjbmi$base_pair_location >= 26000000 & fiadjbmi$base_pair_location <= 34000000),]
fiadjbmi = fiadjbmi[which(!(fiadjbmi$chr_pos%in%mhc$chr_pos)),]

#########################################################
#STEP 1: let's go and loop over the leads and the traits#
#########################################################

for(index in seq(all_variants$variant)){ #we are looping through each index SNP
  
  #And retrieving fiadjbmi data:
  
  print(index)
  
  #STEP 1: let's get the data for each:
  
  lead_ = all_variants$chr_pos[index]
  chr_=all_variants$chromosome[index]
  start_ = as.numeric(all_variants$base_pair_location[index])-500000
  end_ = as.numeric(all_variants$base_pair_location[index])+500000
  
  #STEP 2: let's filter the data:
  
  fiadjbmi_tmp = fiadjbmi[which(as.numeric(fiadjbmi$chromosome) == chr_ & as.numeric(fiadjbmi$base_pair_location) > start_ & as.numeric(fiadjbmi$base_pair_location) < end_),]

  #Let's align it to fiadjbmi+ allele:
    
  fiadjbmi_pos <- fiadjbmi_tmp
    
  #Let's align to positive:
    
  new_a1 <- ifelse(fiadjbmi_pos$beta < 0, fiadjbmi_pos$other_allele, fiadjbmi_pos$effect_allele)
  new_a2 <- ifelse(fiadjbmi_pos$beta < 0, fiadjbmi_pos$effect_allele, fiadjbmi_pos$other_allele)
  new_beta <- ifelse(fiadjbmi_pos$beta < 0, as.numeric(fiadjbmi_pos$beta)*(-1), as.numeric(fiadjbmi_pos$beta))
  new_eaf <- ifelse(fiadjbmi_pos$beta < 0, 1-as.numeric(fiadjbmi_pos$effect_allele_frequency), as.numeric(fiadjbmi_pos$effect_allele_frequency))
    
  fiadjbmi_pos$effect_allele = new_a1
  fiadjbmi_pos$other_allele = new_a2
  fiadjbmi_pos$beta = new_beta
  fiadjbmi_pos$effect_allele_frequency = new_eaf #just tested! Switch work great!!
  
  fiadjbmi_pos = fiadjbmi_pos[order(fiadjbmi_pos$p_value),]
  fiadjbmi_pos=fiadjbmi_pos[which(duplicated(fiadjbmi_pos$variant)==FALSE),]
  fiadjbmi_pos=fiadjbmi_pos[order(fiadjbmi_pos$chromosome, fiadjbmi_pos$base_pair_location),]
  
  #Let's match our data:
  
  hdl_tmp = hdl[which(hdl$chr_pos%in%fiadjbmi_pos$chr_pos),]
  tg_tmp = tg[which(tg$chr_pos%in%fiadjbmi_pos$chr_pos),]

  #Let's align the data:
  
  hdl_aligned = trait_aligner(fiadjbmi_pos, hdl_tmp)
  tg_aligned = trait_aligner(fiadjbmi_pos, tg_tmp)
  
  #Get the shared SNPs:
  
  #these are very old dataframes so be careful:
  
  fiadjbmi_pos = fiadjbmi_pos[which(duplicated(fiadjbmi_pos$chr_pos)==FALSE),]
  hdl_aligned = hdl_aligned[which(duplicated(hdl_aligned$chr_pos)==FALSE),]
  tg_aligned = tg_aligned[which(duplicated(tg_aligned$chr_pos)==FALSE),]
  
  
  shared_snps = Reduce(intersect, list(hdl_aligned$chr_pos,
                                       tg_aligned$chr_pos))
  
  fiadjbmi_pos = fiadjbmi_pos[which(fiadjbmi_pos$chr_pos%in%shared_snps),]
  hdl_aligned = hdl_aligned[which(hdl_aligned$chr_pos%in%shared_snps),]
  tg_aligned = tg_aligned[which(tg_aligned$chr_pos%in%shared_snps),]
  
  #this is very old data so some of them might actually have 
  
  
  #Let's organize the SNPs:
  
  fiadjbmi_pos = fiadjbmi_pos[order(as.numeric(fiadjbmi_pos$chromosome), as.numeric(fiadjbmi_pos$base_pair_location)),] 
  hdl_aligned = hdl_aligned[order(as.numeric(hdl_aligned$chromosome), as.numeric(hdl_aligned$base_pair_location)),] 
  tg_aligned = tg_aligned[order(as.numeric(tg_aligned$chromosome), as.numeric(tg_aligned$base_pair_location)),] 
 
  print(length(which(fiadjbmi_pos$effect_allele == hdl_aligned$effect_allele)))
  print(length(which(fiadjbmi_pos$effect_allele == tg_aligned$effect_allele)))
  
  #We are all set!
  #Let's perform hyprocoloc
    
  betas <- as.matrix(as.data.frame(cbind(as.numeric(fiadjbmi_pos$beta), 
                                         as.numeric(hdl_aligned$beta),
                                         as.numeric(tg_aligned$beta)
                                         )))
  
  ses <- as.matrix(as.data.frame(cbind(as.numeric(fiadjbmi_pos$standard_error), 
                                       as.numeric(hdl_aligned$standard_error),
                                       as.numeric(tg_aligned$standard_error)
  )))
  
  traits <- c("fiadjbmi", "hdl", "tg")
  colnames(betas) <- traits
  colnames(ses) <- traits
  rownames(betas) <- fiadjbmi_pos$variant
  rownames(ses) <- fiadjbmi_pos$variant
    
  res <- hyprcoloc::hyprcoloc(effect.est = betas, effect.se = ses, snp.id = row.names(betas), trait.names = traits)
  res_df <- res$results
  res_df$lead_snp <- lead_
    
  print(res_df)
    
  if(!(exists("coloc_df"))){
      
    coloc_df <- res_df
      
  } else {
      
    coloc_df <- rbind(coloc_df, res_df)
      
  }
    
}

dir.create("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/review/review_output/4_colocalization")

saveRDS(coloc_df, "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/review/review_output/4_colocalization/raw_hyprcoloc_res.RDS")

###################################
#Let's filter now for TG/HDL hits!#
###################################

coloc_high_confidence = coloc_df[which(coloc_df$posterior_prob > 0.5),]
coloc_high_confidence = coloc_high_confidence[which(coloc_high_confidence$regional_prob > 0.5),]

table(coloc_high_confidence$traits)

#Let's take the LD between lead and candiate SNP and whether they are novel or not, which we can easily add later

coloc_fiadjbmi_high_confidence$r2 = NA

for(index in seq(1, length(coloc_fiadjbmi_high_confidence$lead_snp))){
  
  #STEP 1: get the lead and the potential SNP:
  
  lead_ = coloc_fiadjbmi_high_confidence$lead_snp[index]
  candidate_ = coloc_fiadjbmi_high_confidence$candidate_snp[index]
  
  if(lead_ == candidate_){
    
    coloc_fiadjbmi_high_confidence$r2[index] = 1
    
  } else {
    
    ld_res=LDlinkR::LDpair(var1 = lead_, var2 = candidate_,pop = "EUR", token = "04cad4ca4374")
    
    coloc_fiadjbmi_high_confidence$r2[index] = ld_res$r2
    
  }
  
}

#Let's filter for those that are driving our associations:

coloc_fiadjbmi_high_confidence = coloc_fiadjbmi_high_confidence[which(coloc_fiadjbmi_high_confidence$r2 > 0.5),]

##################################################################################
#Let's use this opportunity to add as high-confidence hits and their associations#
##################################################################################

all_variants$fiadjbmi_gw = ifelse(all_variants$p_value.fiadjbmi < 5e-08, "yes", "no")
all_variants$outcome_gw = ifelse(all_variants$p_value.tgadjbmi < 5e-08 |
                                all_variants$p_value.gsatadjbmi < 5e-08 |
                                all_variants$p_value.vatadjbmi < 5e-08 |
                                all_variants$p_value.liver_fat < 5e-08, "yes", "no")
all_variants$novel = ifelse(all_variants$fiadjbmi_gw == "no" & all_variants$outcome_gw == "no", "yes", "no")

#Add data:

all_variants$placo_coloc = "fiadjbmi"
all_variants$placo_coloc = ifelse(all_variants$p_placo.hdl < 5e-08, paste(all_variants$placo_coloc, ", hdl", sep = ""), all_variants$placo_coloc)
all_variants$placo_coloc = ifelse(all_variants$p_placo.tgadjbmi < 5e-08, paste(all_variants$placo_coloc, ", tg", sep = ""), all_variants$placo_coloc)
all_variants$placo_coloc = ifelse(all_variants$p_placo.gsatadjbmi < 5e-08, paste(all_variants$placo_coloc, ", gsat", sep = ""), all_variants$placo_coloc)
all_variants$placo_coloc = ifelse(all_variants$p_placo.vatadjbmi < 5e-08, paste(all_variants$placo_coloc, ", vat", sep = ""), all_variants$placo_coloc)
all_variants$placo_coloc = ifelse(all_variants$p_placo.liver_fat < 5e-08, paste(all_variants$placo_coloc, ", liver_fat", sep = ""), all_variants$placo_coloc)

#Let's add the co-localization data:

all_variants$hyprcoloc = NA

for(index in seq(1, length(all_variants$variant))){
  
  #STEP 1: get the lead:  
  
  lead_ = all_variants$variant[index]
  
  #STEP 2: find the matching:
  
  coloc_match = coloc_fiadjbmi_high_confidence[which(coloc_fiadjbmi_high_confidence$lead_snp == lead_),]
  
  if(is_empty(coloc_match$iteration)){
    
    next()
    
  } else {
    
    all_variants$hyprcoloc[index] = coloc_match$traits
  }
  
}

#############################################
#All right this data frame is SO MUCH BETTER#
#############################################

high_confidence = all_variants[which(all_variants$novel == "no" & is.na(all_variants$hyprcoloc) == FALSE),] #28
novel = all_variants[which(all_variants$novel == "yes"),] #126

#Perfect, this can boost our possibilities of finding loci that are similar in their biology

fwrite(all_variants, "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_PLACO_2026/output/4_colocalization/fiadjbmi/353_variants_classified_with_placo_and_hyprcoloc.txt")
