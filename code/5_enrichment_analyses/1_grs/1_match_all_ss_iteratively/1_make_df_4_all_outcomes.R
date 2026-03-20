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
  
  if("sample_cases"%in%colnames(other_ss)){
    
    outcome <- other_ss %>%
      select(variant, chromosome, base_pair_location, effect_allele, other_allele, effect_allele_frequency, beta, standard_error, p_value, sample_size, sample_cases, sample_controls, prevalence, chr_pos)
    
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
    
    merged_proxies_df <- harmonise_data(exposure_proxies, outcome, action=2)
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

path_2_input <- "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/"

setwd(path_2_input)

iradjbmi <- fread("output/2_mv_gwas/fiadjbmi_hdl_tg_common_dwls_curated.txt")
all_variants <- fread("output/4_ir_loci_discovery/2_ir_variants/282_ir_fiadjbmi_hdl_tg_dwls_snps_w_cpassoc.txt")
ir_variants <- fread("output/4_ir_loci_discovery/3_novel_variants/3_novel_variants/282_ir_variants_with_ir_and_t2d_comparisons_w_closest_genes.txt")
#non_ir_variants <- all_variants[which(!(all_variants$variant%in%ir_variants$variant)),]

#Also the whole dataframe to get proxies if necessary:

exposure_df <- iradjbmi

#Let's add the proxies for all the variants, making it a bit easier for us:

proxies <- fread("output/4_ir_loci_discovery/4_proxies/proxies_build_37_and_38_4_282_variants.txt")
proxies <- proxies[which(proxies$query_snp_rsid%in%ir_variants$variant),]

length(which(duplicated(proxies$query_snp_rsid) == FALSE)) #282/282

#####################################################################################
#Let's get the data in build37 for the proxies cuz it makes our lives easier, really#
#####################################################################################

proxies$chr_pos <- paste("chr", proxies$chr, ":", proxies$pos_hg19, sep= "") #let's add the chr_pos 

###############################
#Let's prepare the output data#
###############################

dir.create("output/5_enrichment_analyses/1_prs")
dir.create("output/5_enrichment_analyses/1_prs/1_expo_outcome_df")

output_path <- "output/5_enrichment_analyses/1_prs/1_expo_outcome_df/"

######################
#Loading outcome data#
######################

fiadjbmi <- fread("output/1_curated_gwas/fiadjbmi_curated.txt")
hdl <- fread("output/1_curated_gwas/hdl_curated.txt")
tg <- fread("output/1_curated_gwas/tg_curated.txt")
tg_hdl_ratio <- fread("output/1_curated_gwas/tg_hdl_ratio_curated.txt")
isiadjbmi <- fread("output/1_curated_gwas/isiadjbmi_curated.txt")
ifcadjbmi <- fread("output/1_curated_gwas/ifcadjbmi_curated.txt")
fgadjbmi <- fread("output/1_curated_gwas/fgadjbmi_curated.txt")
thgadjbmi <- fread("output/1_curated_gwas/thgadjbmi_curated.txt") 

##################################################
#For the data I will need some help for ISIadjBMI#
##################################################

yes_vect <- c("A", "G", "T", "C")
fiadjbmi <- fiadjbmi[which(fiadjbmi$effect_allele%in%yes_vect & fiadjbmi$other_allele%in%yes_vect),]

fiadjbmi_match <- fiadjbmi[which(fiadjbmi$chr_pos%in%isiadjbmi$chr_pos),]
fiadjbmi_match <- fiadjbmi_match[which(duplicated(fiadjbmi_match$chr_pos) == FALSE),]

isiadjbmi_match <- isiadjbmi[which(isiadjbmi$chr_pos%in%fiadjbmi_match$chr_pos),]
isiadjbmi_match <- isiadjbmi_match[which(duplicated(isiadjbmi_match$chr_pos) == FALSE),]

isiadjbmi_match_ordered <- isiadjbmi_match[order(match(isiadjbmi_match$chr_pos, fiadjbmi_match$chr_pos)),]

length(which(isiadjbmi_match_ordered$chr_pos == fiadjbmi_match$chr_pos)) #perfect

isiadjbmi_match_ordered$variant <- fiadjbmi_match$variant

#Now the same for IFCadjBMI

fiadjbmi_match <- fiadjbmi[which(fiadjbmi$chr_pos%in%ifcadjbmi$chr_pos),]
fiadjbmi_match <- fiadjbmi_match[which(duplicated(fiadjbmi_match$chr_pos) == FALSE),]

ifcadjbmi_match <- ifcadjbmi[which(ifcadjbmi$chr_pos%in%fiadjbmi_match$chr_pos),]
ifcadjbmi_match <- ifcadjbmi_match[which(duplicated(ifcadjbmi_match$chr_pos) == FALSE),]

ifcadjbmi_match_ordered <- ifcadjbmi_match[order(match(ifcadjbmi_match$chr_pos, fiadjbmi_match$chr_pos)),]

length(which(ifcadjbmi_match_ordered$chr_pos == fiadjbmi_match$chr_pos)) #perfect

ifcadjbmi_match_ordered$variant <- fiadjbmi_match$variant

#########################################################################
#For the munging I am going to need some help for FGadjBMI and 2hGadjBMI#
#########################################################################

#Let's get RSIDs for ISIadjBMI

yes_vect <- c("A", "G", "T", "C")
fiadjbmi <- fiadjbmi[which(fiadjbmi$effect_allele%in%yes_vect & fiadjbmi$other_allele%in%yes_vect),]

fiadjbmi_match <- fiadjbmi[which(fiadjbmi$chr_pos%in%fgadjbmi$chr_pos),]
fiadjbmi_match <- fiadjbmi_match[which(duplicated(fiadjbmi_match$chr_pos) == FALSE),]

fgadjbmi_match <- fgadjbmi[which(fgadjbmi$chr_pos%in%fiadjbmi_match$chr_pos),]
fgadjbmi_match <- fgadjbmi_match[which(duplicated(fgadjbmi_match$chr_pos) == FALSE),]

fgadjbmi_match_ordered <- fgadjbmi_match[order(match(fgadjbmi_match$chr_pos, fiadjbmi_match$chr_pos)),]

length(which(fgadjbmi_match_ordered$chr_pos == fiadjbmi_match$chr_pos)) #perfect

fgadjbmi_match_ordered$variant <- fiadjbmi_match$variant

#Now the same for thgadjbmi

fiadjbmi_match <- fiadjbmi[which(fiadjbmi$chr_pos%in%thgadjbmi$chr_pos),]
fiadjbmi_match <- fiadjbmi_match[which(duplicated(fiadjbmi_match$chr_pos) == FALSE),]

thgadjbmi_match <- thgadjbmi[which(thgadjbmi$chr_pos%in%fiadjbmi_match$chr_pos),]
thgadjbmi_match <- thgadjbmi_match[which(duplicated(thgadjbmi_match$chr_pos) == FALSE),]

thgadjbmi_match_ordered <- thgadjbmi_match[order(match(thgadjbmi_match$chr_pos, fiadjbmi_match$chr_pos)),]

length(which(thgadjbmi_match_ordered$chr_pos == fiadjbmi_match$chr_pos)) #perfect

thgadjbmi_match_ordered$variant <- fiadjbmi_match$variant

#######################################
#Also let's do small changes to TG/HDL#
#######################################

tg_hdl_ratio$sample_size <- tg_hdl_ratio$n
tg_hdl_ratio$variant <- tg_hdl_ratio$rs_id

###############################################
#We are ready to load all the pertinent data!!#
###############################################

bmi_ss <- fread("output/1_curated_gwas/bmi_curated.txt")
whr_ss <- fread("output/1_curated_gwas/whr_curated.txt")
whradjbmi_ss <- fread("output/1_curated_gwas/whradjbmi_curated.txt")
wc_ss <- fread("output/1_curated_gwas/wc_curated.txt")
wcadjbmi_ss <- fread("output/1_curated_gwas/wcadjbmi_giant_curated.txt")
hc_ss <- fread("output/1_curated_gwas/hc_curated.txt")
hcadjbmi_ss <- fread("output/1_curated_gwas/hcadjbmi_giant_curated.txt")

#And now the diseases:

t2d_ss <- fread("output/1_curated_gwas/t2d_curated.txt")
chd_ss <- fread("output/1_curated_gwas/chd_curated.txt")
nafld_ss <- fread("output/1_curated_gwas/nafld_curated.txt")
pcos_ss <- fread("output/1_curated_gwas/pcos_curated.txt")
ckd_ss <- fread("output/1_curated_gwas/ckd_curated.txt")
hypertension_ss <- fread("output/1_curated_gwas/hypertension_curated.txt")

######################################################################
#STEP 1: Let's clean the data a little bit to get it through the SNPs#
######################################################################

# #Let's do a small change for T2D Suzuki::
# colnames(t2d_suzuki) <- c("variant", "chromosome", "base_pair_location", "effect_allele", "other_allele", "effect_allele_frequency", "beta", "standard_error", "p_value", "sample_size", "sample_cases", "sample_controls","chr_pos")
# 
# #Let's do the same as with ISIadjBMI so that it runs properly:
# 
# fiadjbmi_match <- fiadjbmi[which(fiadjbmi$chr_pos%in%t2d_suzuki$chr_pos),]
# fiadjbmi_match <- fiadjbmi_match[which(duplicated(fiadjbmi_match$chr_pos) == FALSE),]
# 
# t2d_suzuki <- t2d_suzuki[which(t2d_suzuki$chr_pos%in%fiadjbmi_match$chr_pos),]
# t2d_suzuki <- t2d_suzuki[which(duplicated(t2d_suzuki$chr_pos) == FALSE),]
# 
# t2d_suzuki <- t2d_suzuki[order(match(t2d_suzuki$chr_pos, fiadjbmi_match$chr_pos)),]
# 
# length(which(t2d_suzuki$chr_pos == fiadjbmi_match$chr_pos)) #perfect
# 
# t2d_suzuki$variant <- fiadjbmi_match$variant

#And now add the data

t2d_ss$sample_size <- t2d_ss$sample_cases+t2d_ss$sample_controls
chd_ss$sample_size <- chd_ss$sample_cases+chd_ss$sample_controls
nafld_ss$sample_size <- nafld_ss$sample_cases+nafld_ss$sample_controls
pcos_ss$sample_size <- pcos_ss$sample_cases+pcos_ss$sample_controls
ckd_ss$sample_size <- ckd_ss$sample_cases+ckd_ss$sample_controls
hypertension_ss$sample_size <- hypertension_ss$sample_cases+hypertension_ss$sample_controls

#Let's get the percentage:

#t2d_suzuki$prevalence<- NA #we are not using it for anything else so we can actually do this!
t2d_ss$prevalence<- as.numeric(t2d_ss$prevalence)/100
chd_ss$prevalence<- as.numeric(chd_ss$prevalence)/100
nafld_ss$prevalence<- as.numeric(nafld_ss$prevalence)/100
pcos_ss$prevalence<- as.numeric(pcos_ss$prevalence)/100
ckd_ss$prevalence<- as.numeric(ckd_ss$prevalence)/100
hypertension_ss$prevalence<- as.numeric(hypertension_ss$prevalence)/100

##############
#More depots!#
##############

asat <- fread("output/1_curated_gwas/asat_curated_saaket.txt")
asatadjbmi <- fread("output/1_curated_gwas/asatadj_curated_saaket.txt")

asat$sample_size <- 38965
asatadjbmi$sample_size <- 37641 #data from supplementary info which is more accurate than overall in paper

gsat <- fread("output/1_curated_gwas/gfat_curated_saaket.txt")
gsatadjbmi <- fread("output/1_curated_gwas/gfatadj_curated_saaket.txt")

asat$sample_size <- 38965
asatadjbmi$sample_size <- 37641 #data from supplementary info which is more accurate than overall in paper

vat <- fread("output/1_curated_gwas/vat_curated_saaket.txt")
vatadjbmi <- fread("output/1_curated_gwas/vatadj_curated_saaket.txt")

vat$sample_size <- 38965
vatadjbmi$sample_size <- 37641 #data from supplementary info which is more accurate than overall in paper

liver_volume <- fread("output/1_curated_gwas/liver_volume_curated.txt")
liver_fat <- fread("output/1_curated_gwas/liver_fat_curated.txt")
pancreas_fat <- fread("output/1_curated_gwas/pancreas_fat_curated.txt")
anterior_thigh_fat <- fread("output/1_curated_gwas/anterior_thight_fat_curated.txt")
posterior_thigh_fat <- fread("output/1_curated_gwas/posterior_thight_fat_curated.txt")

######################################
#STEP 2: Let's load the data properly#
######################################

list_of_phenos <- list(fiadjbmi,
                       hdl,
                       tg,
                       isiadjbmi_match_ordered,
                       ifcadjbmi_match_ordered,
                       fgadjbmi_match_ordered,
                       thgadjbmi_match_ordered,
                       tg_hdl_ratio,
                       bmi_ss,
                       hc_ss,
                       hcadjbmi_ss,
                       wc_ss,
                       wcadjbmi_ss,
                       whr_ss,
                       whradjbmi_ss,
                       t2d_ss,
                       #t2d_suzuki,
                       chd_ss,
                       nafld_ss,
                       pcos_ss,
                       ckd_ss,
                       hypertension_ss,
                       asat,
                       asatadjbmi,
                       gsat,
                       gsatadjbmi, 
                       vat,
                       vatadjbmi,
                       liver_volume,
                       liver_fat,
                       pancreas_fat,
                       anterior_thigh_fat,
                       posterior_thigh_fat)

trait_names <- c("FIadjBMI",
                 "HDL",
                 "TG",
                 "ISIadjBMI",
                 "IFCadjBMI",
                 "FGadjBMI", 
                 "2hGadjBMI",
                 "tg_hdl_ratio",
                 "BMI",
                 "HC",
                 "HCadjBMI",
                 "WC",
                 "WCadjBMI",
                 "WHR",
                 "WHRadjBMI",
                 "T2D",
                 #"T2D (Suzuki)",
                 "CHD",
                 "NAFLD", 
                 "PCOS",
                 "CKD",
                 "Hypertension",
                 "ASAT",
                 "ASATadjBMI",
                 "GSAT",
                 "GSATadjBMI", 
                 "VAT", 
                 "VATadjBMI", 
                 "Liver volume",
                 "Liver fat %",
                 "Pancreas fat %",
                 "Anterior thigh fat",
                 "Posterior thigh fat")

type_of_trait <- c(rep("Insulin Resistance", 8), rep("Anthropometric", 7), rep("Cardiometabolic complication", 6), rep("MRI", 11))

########################################################
#STEP 1: let's align the cluster to the positive allele#
########################################################

data_df <- iradjbmi[which(iradjbmi$chr_pos%in%all_variants$chr_pos),]

#And align them to the positive allele:
  
new_a1 <- ifelse(as.numeric(data_df$beta) < 0, data_df$other_allele, data_df$effect_allele)
new_a2 <- ifelse(as.numeric(data_df$beta) < 0, data_df$effect_allele, data_df$other_allele)
new_beta <- ifelse(as.numeric(data_df$beta) < 0, as.numeric(data_df$beta)*(-1), as.numeric(data_df$beta))
#new_eaf <- ifelse(as.numeric(data_df$beta) < 0, 1-as.numeric(data_df$effect_allele_frequency), as.numeric(data_df$effect_allele_frequency))

data_df$effect_allele <- new_a1
data_df$other_allele <- new_a2
data_df$beta <- new_beta
#data_df$effect_allele_frequency <- new_eaf

#######################################
#STEP 2: recursively match ir variants#
#######################################

ir_df <- data_df[which(data_df$chr_pos%in%ir_variants$chr_pos),]

recursively_matching_data(exp_df = ir_df, list_of_phenos_ =list_of_phenos, trait_names_ = trait_names, output_path <- output_path, exposure = "ir_variants_")

