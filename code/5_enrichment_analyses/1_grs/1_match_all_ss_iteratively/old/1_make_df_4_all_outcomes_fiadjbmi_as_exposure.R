##############
#INTRODUCTION#
##############

#This code will run unweighted GRS for all idps traits using initiation smoking and smoking initiation exposure IVs.

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)
library(TwoSampleMR)

###################
#Loading functions#
###################

source("N:/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/Mediation_analysis/code/0_functions/functions_4_mediation.R")

chr_cleaner <- function(chr_pos){
  
  chr_ <- str_split(chr_pos, ":")[[1]][1]
  chr_ <- str_split(chr_, "chr")[[1]][2]
  
  return(as.numeric(chr_))
}

pos_cleaner <- function(chr_pos){
  
  pos_ <- str_split(chr_pos, ":")[[1]][2]
  
  return(as.numeric(pos_))

}

data_aligner <- function(query_ss, other_ss){
  
  ################################################################################
  #This code uses: exposure and proxy dataframes that should be loaded beforehand#
  ################################################################################
  
  #########################################################################################
  #STEP 0: let's run the matching with TwoSampleMR so we need the data in a certain format#
  #########################################################################################
  
  #Let's first check if we have the effect_allele_frequency column:
  
  check <- which(colnames(query_ss) == "effect_allele_frequency")
  
  if(is_empty(check)){
    
    query_ss$effect_allele_frequency <- NA
    
  }
  
  #Now we can proceed
  
  exposure <- query_ss %>%
    select(variant, chromosome, base_pair_location, effect_allele, other_allele, effect_allele_frequency, beta, standard_error, p_value, sample_size, chr_pos)
  
  colnames(exposure) <- c("SNP", "chr.exposure", "pos.exposure", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure", "beta.exposure", "se.exposure", "pval.exposure", "samplesize.exposure", "chr_pos.exposure")
  
  exposure$exposure <- "exposure"
  exposure$id.exposure <- "exposure"
  
  #Now with the outcome:
  
  check <- which(colnames(other_ss) == "effect_allele_frequency")
  
  if(is_empty(check)){
    
    other_ss$effect_allele_frequency <- NA
    
  }
  
  #Now we can proceed
  
  if("sample_size_cases"%in%colnames(other_ss)){
    
    outcome <- other_ss %>%
      select(variant, chromosome, base_pair_location, effect_allele, other_allele, effect_allele_frequency, beta, standard_error, p_value, sample_size, sample_size_cases, sample_size_controls, prevalence, chr_pos)
    
    colnames(outcome) <- c("SNP", "chr.outcome", "pos.outcome", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome", "pval.outcome", "samplesize.outcome", "ncase.outcome", "ncontrol.outcome", "prevalence.outcome", "chr_pos.outcome")
    
  } else {
    
    outcome <- other_ss %>%
      select(variant, chromosome, base_pair_location, effect_allele, other_allele, effect_allele_frequency, beta, standard_error, p_value, sample_size, chr_pos)
    
    colnames(outcome) <- c("SNP", "chr.outcome", "pos.outcome", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome", "pval.outcome", "samplesize.outcome", "chr_pos.outcome")
    
  }
  
  outcome$outcome <- "outcome"
  outcome$id.outcome <- "outcome"
  
  ############################################################################################################
  #STEP 1: match the data. This will probably fail with tri-allelic SNPs. Here I think we are gonna be OK... #
  ############################################################################################################
  
  merged_df <- harmonise_data(exposure, outcome, action=3)
  print(dim(merged_df))
  
  merged_df <- merged_df[which(merged_df$remove == FALSE),] #removing incompatible alleles
  
  ######################################
  #STEP 2: do we need to query proxies?#
  ######################################
  
  missing_snps <- query_ss[which(!(query_ss$chr_pos%in%merged_df$chr_pos.exposure)),] 
  
  if(length(missing_snps$chr_pos) == 0){
    
    return(merged_df)
    
  } else {
    
    #We have to deal with proxies!!
    
    #1: check proxies for missing variants:
    
    proxies_4_missing <- proxies[which(proxies$query_snp_rsid%in%missing_snps$variant),]
    
    #2: check which are available in outcome:
    
    proxies_in_outcome <- proxies_4_missing[which(proxies_4_missing$rsID%in%other_ss$variant),]
    
    #3. Let's order the proxies by exposure data:
    
    exposure_match <- exposure_df[which(exposure_df$chr_pos%in%proxies_in_outcome$chr_pos),]
    
    #4. Let's add lead SNP so that we can properly add data:
    
    proxies_match <- proxies_in_outcome[which(proxies_in_outcome$chr_pos%in%exposure_match$chr_pos),]
    proxies_ordered <- proxies_match[order(match(proxies_match$chr_pos, exposure_match$chr_pos)),]
    
    length(which(proxies_ordered$chr_pos == exposure_match$chr_pos))
    
    exposure_match$lead_snp <- proxies_ordered$query_snp_rsid
    exposure_match$variant <- proxies_ordered$rsID #adding RSID to have ALL INFO. Maybe the summary statistics only has CHR:POS.
    exposure_match$rsq <- proxies_ordered$r2 #adding RSID to have ALL INFO. Maybe the summary statistics only has CHR:POS.
    
    #5. Let's get only one hit per signal, retaining only the best:
    
    exposure_best <- exposure_match[order(as.numeric(exposure_match$rsq), decreasing = TRUE),]
    exposure_best <- exposure_best[which(duplicated(exposure_best$lead_snp) == FALSE),] #7!
    
    #6. Let's align to the positive allele:
    
    tmp_df <- exposure_best

    #And align them to the positive allele:
    
    new_a1 <- ifelse(as.numeric(tmp_df$beta < 0), tmp_df$other_allele, tmp_df$effect_allele)
    new_a2 <- ifelse(as.numeric(tmp_df$beta < 0), tmp_df$effect_allele, tmp_df$other_allele)
    new_beta <- ifelse(as.numeric(tmp_df$beta < 0), as.numeric(tmp_df$beta)*(-1), as.numeric(tmp_df$beta))
    new_eaf <- ifelse(as.numeric(tmp_df$beta < 0), 1-as.numeric(tmp_df$effect_allele_frequency), as.numeric(tmp_df$effect_allele_frequency))
    
    tmp_df$effect_allele <- new_a1
    tmp_df$other_allele <- new_a2
    tmp_df$beta <- new_beta
    tmp_df$effect_allele_frequency <- new_eaf
    
    #7. Let's get the data ready for merging...
    
    exposure_proxies <- tmp_df %>%
      select(variant, chromosome, base_pair_location, effect_allele, other_allele, effect_allele_frequency, beta, standard_error, p_value, sample_size, chr_pos)
    
    colnames(exposure_proxies) <- c("SNP", "chr.exposure", "pos.exposure", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure", "beta.exposure", "se.exposure", "pval.exposure", "samplesize.exposure", "chr_pos.exposure")
    
    exposure_proxies$exposure <- "exposure"
    exposure_proxies$id.exposure <- "exposure"
    
    #And let's merge with the outcome, that should have ALL the data:
    
    merged_proxies_df <- harmonise_data(exposure_proxies, outcome, action=3)
    print(dim(merged_proxies_df))
    
    merged_proxies_df <- merged_proxies_df[which(merged_proxies_df$remove == FALSE),] #removing incompatible alleles
    
    #8. Let's add this info to the previous version:
    
    merged_df <- rbind(merged_df, merged_proxies_df)
    
    #9. Let's double-check for independence:
    
    #merged_df <- filt_ld(merged_df)
    
    return(merged_df)
    
  }
  
}

recursively_matching_data <- function(exp_df, list_of_phenos_, trait_names_, output_path, exposure){
  
  #STEP 0: let's make a list of traits that we need to be careful with:
  
  #STEP 1: let's loop and perform the matching:
  
  for(index_trait in seq(1, length(list_of_phenos_))){
    
    print(index_trait)
    
    trait_df <- list_of_phenos_[index_trait][[1]]
    
    #Let's see if we use data_aligner_37 or other:
    
    trait_name <- trait_names_[index_trait]
    
    print(trait_name)
    
    matched_data <- data_aligner(exp_df, trait_df)
    
    fwrite(matched_data, paste(output_path, exposure, trait_name, ".txt", sep =""))    
    
  }
  
}


#######################
#Loading exposure data#
#######################

#Let's get a path where the clusters are so that we can loop throught them:

path_2_input <- "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2023/"

setwd(path_2_input)

iradjbmi <- fread("output/3_mv_gwas/fiadjbmi_hdl_tg/fiadjbmi_hdl_tg_commmon_dwls_curated.txt")
all_variants <- fread("output/4_ld_clumping/2_clumped_data/fiadjbmi_hdl_tg/common/clumped_data/fiadjbmi_hdl_tg_common_dwls_clumped.txt")
ir_variants <- fread("output/5_comparing_fiadjbmi_fi/2_sensitivity_analysis/282_old_and_new_ir_compared_w_tg_hdl.txt")
non_ir_variants <- all_variants[which(!(all_variants$rsid%in%ir_variants$variant)),]
fiadjbmi <- fread("output/1_curated_data/fiadjbmi_curated.txt")

#Also the whole dataframe to get proxies if necessary:

exposure_df <- fiadjbmi

proxies <- fread("output/5_comparing_fiadjbmi_fi/4_mr/all_variants/0_proxies/haploreg_output.txt", fill = TRUE)

proxies <- proxies[which(proxies$query_snp_rsid%in%all_variants$rsid),]

length(which(duplicated(proxies$query_snp_rsid) == FALSE)) #352/352

#####################################################################################
#Let's get the data in build37 for the proxies cuz it makes our lives easier, really#
#####################################################################################

# proxies$chr_pos <- NA
# 
# for(index_rsid in seq(1, length(proxies$rsID))){
#   
#   snp_df <- tryCatch(otargen::variantInfo(proxies$rsID[index_rsid]), error = function(e) { skip_to_next <<- TRUE})
#   
#   if(length(snp_df) == 1){
#     
#     next()
#   }
#   
#   proxies$chr_pos[index_rsid] <- paste("chr", snp_df$chromosomeB37, ":", snp_df$positionB37, sep="")
#   
#   
# }
# 
# #Finally, let's clean the proxies with NA in their chr_pos
# 
# proxies <- proxies[which(is.na(proxies$chr_pos) == FALSE),] #Worked!
# 
# fwrite(proxies, "output/5_comparing_fiadjbmi_fi/4_mr/all_variants/0_proxies/proxies_df.txt")

proxies <- fread("output/5_comparing_fiadjbmi_fi/4_mr/all_variants/0_proxies/proxies_df.txt")

###############################
#Let's prepare the output data#
###############################

output_path <- "output/5_comparing_fiadjbmi_fi/4_mr/all_variants/1_expo_outcome_df_fiadjbmi/"

dir.create(output_path)

######################
#Loading outcome data#
######################

fiadjbmi <- fread("output/1_curated_data/fiadjbmi_curated.txt")
hdl <- fread("output/1_curated_data/hdl_curated.txt")
tg <- fread("output/1_curated_data/tg_curated.txt")
isiadjbmi <- fread("../../Team projects/Hermina&Mario&MariaJose/output/1_curated_data/isiadjbmi_curated.txt")
tg_hdl_ratio <- fread("output/1_curated_data/tg_hdl_ratio_curated.txt")

tg_hdl_ratio$variant <- tg_hdl_ratio$rs_id
tg_hdl_ratio$sample_size <- tg_hdl_ratio$n

##################################################
#For the data I will need some help for ISIadjBMI#
##################################################

#Let's get RSIDs for ISIadjBMI

fiadjbmi_match <- fiadjbmi[which(fiadjbmi$chr_pos%in%isiadjbmi$chr_pos),]
fiadjbmi_match <- fiadjbmi_match[which(duplicated(fiadjbmi_match$chr_pos) == FALSE),]

isiadjbmi_match <- isiadjbmi[which(isiadjbmi$chr_pos%in%fiadjbmi_match$chr_pos),]
isiadjbmi_match <- isiadjbmi_match[which(duplicated(isiadjbmi_match$chr_pos) == FALSE),]

isiadjbmi_match_ordered <- isiadjbmi_match[order(match(isiadjbmi_match$chr_pos, fiadjbmi_match$chr_pos)),]

length(which(isiadjbmi_match_ordered$chr_pos == fiadjbmi_match$chr_pos)) #perfect

isiadjbmi_match_ordered$variant <- fiadjbmi_match$variant

###############################################
#We are ready to load all the pertinent data!!#
###############################################

bmi_ss <- fread("output/1_curated_data/bmi_giant_curated.txt")
whr_ss <- fread("output/1_curated_data/whr_giant_curated.txt")
whradjbmi_ss <- fread("output/1_curated_data/whradjbmi_giant_curated.txt")
wc_ss <- fread("output/1_curated_data/wc_giant_curated.txt")
wcadjbmi_ss <- fread("output/1_curated_data/wcadjbmi_giant_curated.txt")
hc_ss <- fread("output/1_curated_data/hc_giant_curated.txt")
hcadjbmi_ss <- fread("output/1_curated_data/hcadjbmi_giant_curated.txt")

#And now the diseases:

t2d_ss <- fread("output/1_curated_data/t2d_curated.txt")
chd_ss <- fread("output/1_curated_data/chd_curated.txt")
nafld_ss <- fread("output/1_curated_data/nafld_curated.txt")
pcos_ss <- fread("output/1_curated_data/pcos_curated.txt")
ckd_ss <- fread("output/1_curated_data/ckd_curated.txt")
hypertension_ss <- fread("output/1_curated_data/hypertension_curated.txt")

######################################################################
#STEP 1: Let's clean the data a little bit to get it through the SNPs#
######################################################################

t2d_ss$sample_size <- t2d_ss$sample_cases+t2d_ss$sample_controls
chd_ss$sample_size <- chd_ss$sample_cases+chd_ss$sample_controls
nafld_ss$sample_size <- nafld_ss$sample_cases+nafld_ss$sample_controls
pcos_ss$sample_size <- pcos_ss$sample_cases+pcos_ss$sample_controls
ckd_ss$sample_size <- ckd_ss$sample_cases+ckd_ss$sample_controls
hypertension_ss$sample_size <- hypertension_ss$sample_cases+hypertension_ss$sample_controls

#Let's get the percentage:

t2d_ss$prevalence<- as.numeric(t2d_ss$prevalence)/100
chd_ss$prevalence<- as.numeric(chd_ss$prevalence)/100
nafld_ss$prevalence<- as.numeric(nafld_ss$prevalence)/100
pcos_ss$prevalence<- as.numeric(pcos_ss$prevalence)/100
ckd_ss$prevalence<- as.numeric(ckd_ss$prevalence)/100
hypertension_ss$prevalence<- as.numeric(hypertension_ss$prevalence)/100

######################################
#STEP 2: Let's load the data properly#
######################################

list_of_phenos <- list(fiadjbmi,
                       hdl,
                       tg,
                       isiadjbmi_match_ordered,
                       tg_hdl_ratio,
                       bmi_ss,
                       hc_ss,
                       hcadjbmi_ss,
                       wc_ss,
                       wcadjbmi_ss,
                       whr_ss,
                       whradjbmi_ss,
                       t2d_ss,
                       chd_ss,
                       nafld_ss,
                       pcos_ss,
                       ckd_ss,
                       hypertension_ss)

trait_names <- c("FIadjBMI",
                 "HDL",
                 "TG",
                 "ISIadjBMI",
                 "tg_hdl_ratio",
                 "BMI",
                 "HC",
                 "HCadjBMI",
                 "WC",
                 "WCadjBMI",
                 "WHR",
                 "WHRadjBMI",
                 "T2D",
                 "CHD",
                 "NAFLD", 
                 "PCOS",
                 "CKD",
                 "Hypertension")

type_of_trait <- c(rep("Insulin Resistance", 5), rep("Anthropometric", 7), rep("Cardiometabolic complication", 6))

########################################################
#STEP 1: let's align the cluster to the positive allele#
########################################################

data_df <- fiadjbmi[which(fiadjbmi$chr_pos%in%all_variants$chr_pos),]

#And align them to the positive allele:
  
new_a1 <- ifelse(as.numeric(data_df$beta) < 0, data_df$other_allele, data_df$effect_allele)
new_a2 <- ifelse(as.numeric(data_df$beta) < 0, data_df$effect_allele, data_df$other_allele)
new_beta <- ifelse(as.numeric(data_df$beta) < 0, as.numeric(data_df$beta)*(-1), as.numeric(data_df$beta))
new_eaf <- ifelse(as.numeric(data_df$beta) < 0, 1-as.numeric(data_df$effect_allele_frequency), as.numeric(data_df$effect_allele_frequency))

data_df$effect_allele <- new_a1
data_df$other_allele <- new_a2
data_df$beta <- new_beta
data_df$effect_allele_frequency <- new_eaf

########################################
#STEP 1: recursively match all variants#
########################################

recursively_matching_data(exp_df = data_df, list_of_phenos_ =list_of_phenos, trait_names_ = trait_names, output_path <- output_path, exposure = "fiadjbmi_all_variants_")

#######################################
#STEP 2: recursively match ir variants#
#######################################

ir_df <- data_df[which(data_df$chr_pos%in%ir_variants$chr_pos),]

recursively_matching_data(exp_df = ir_df, list_of_phenos_ =list_of_phenos, trait_names_ = trait_names, output_path <- output_path, exposure = "fiadjbmi_ir_variants_")

###########################################
#STEP 3: recursively match non ir variants#
###########################################

non_ir_df <- data_df[which(!(data_df$chr_pos%in%ir_variants$chr_pos)),]

recursively_matching_data(exp_df = non_ir_df, list_of_phenos_ =list_of_phenos, trait_names_ = trait_names, output_path <- output_path, exposure = "fiadjbmi_non_ir_variants_")

#############################################
#STEP 4: recursively match novel ir variants#
#############################################

novel_ir <- ir_variants[which(ir_variants$new == "yes" & ir_variants$tg_hdl_ratio == "no"),]

novel_ir_df <- data_df[which(data_df$chr_pos%in%novel_ir$chr_pos),]

recursively_matching_data(exp_df = novel_ir_df, list_of_phenos_ =list_of_phenos, trait_names_ = trait_names, output_path <- output_path, exposure = "fiadjbmi_novel_ir_variants_")

################################################
#STEP 5: recursively match reported ir variants#
################################################

reported_ir <- ir_variants[which(!(ir_variants$new == "yes" & ir_variants$tg_hdl_ratio == "no")),]

reported_ir_df <- data_df[which(data_df$chr_pos%in%reported_ir$chr_pos),]

recursively_matching_data(exp_df = reported_ir_df, list_of_phenos_ =list_of_phenos, trait_names_ = trait_names, output_path <- output_path, exposure = "fiadjbmi_reported_ir_variants_")
