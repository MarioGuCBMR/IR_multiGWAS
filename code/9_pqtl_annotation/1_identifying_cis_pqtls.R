##############
#INTRODUCTION#
##############

#This code reads the data from Open platform and saves it:

###################
#Loading libraries#
###################

library(arrow)
library(dplyr)
library(tidyverse)
library(data.table)

###################
#Loading functions#
###################

merge_parquet_files <- function(folder_path) {
  # List all .parquet files in the folder
  parquet_files <- list.files(path = folder_path, pattern = "\\.parquet$", full.names = TRUE)
  
  if (length(parquet_files) == 0) {
    stop("No parquet files found in the specified directory.")
  }
  
  # Read and bind them together, with a message for each file
  merged_df <- parquet_files %>%
    lapply(function(file) {
      message("Reading: ", file)
      read_parquet(file)
    }) %>%
    bind_rows()
  
  return(merged_df)
}

protein_parser <- function(protein_id){
  
  #STEP 1: figure out what kind is it:
  
  check <- str_detect(protein_id, "UKB_PPP")
  
  if(check){
    
    prot <- unlist(str_split(protein_id, "UKB_PPP_EUR_"))[2]
    prot <- unlist(str_split(prot, "_"))[1]
    
    return(prot)
    
  } else {
    
    prot <- unlist(str_split(protein_id, "sun_2018_aptamer_plasma_"))[2]
    prot <- unlist(str_split(prot, "_"))[1]
    
    return(toupper(prot))
    
  }
  
}

other_allele_parser <- function(id_){
  #ID format: CHR_BP_REF_ALT

  allele <- unlist(str_split(id_, "_"))[3]
  
  return(allele)
  
}

effect_allele_parser <- function(id_){
  #ID format: CHR_BP_REF_ALT
  
  allele <- unlist(str_split(id_, "_"))[4]
  
  return(allele)
  
}


##############
#Loading data#
##############

fine_mapped_path <- "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/open_platform_data/raw_data/open_platform_fine_mapped_data/"

credible_sets_df <- as.data.frame(merge_parquet_files(fine_mapped_path))

################################
#STEP 1: let's do a quick match#
################################

#Let's get the variants from PCOS project and see what we can find through them:

path_wd <- "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/"

setwd(path_wd)

ir_new <- fread("output/4_ir_loci_discovery/3_novel_variants/3_novel_variants/282_ir_variants_with_ir_and_t2d_comparisons_w_closest_genes.txt")
proxies <- data.table::fread("output/4_ir_loci_discovery/4_proxies/proxies_build_37_and_38_4_282_variants.txt")

#The match is required to be done through OTG ID, otargen package can help us with this

proxies$id_otargen_1 <- paste(proxies$chr, "_", proxies$pos_hg38, "_", proxies$ref, "_", proxies$alt, sep = "")
proxies$id_otargen_2 <- paste(proxies$chr, "_", proxies$pos_hg38, "_", proxies$alt, "_", proxies$ref, sep = "")

######################
#Let's match the data#
######################

fine_mapped_match <- credible_sets_df[which(credible_sets_df$variantId%in%proxies$id_otargen_1 | credible_sets_df$variantId%in%proxies$id_otargen_2),]

#For the time being I am only interested in molecular QTL so let's remove the QTLs:

fine_mapped_match <- fine_mapped_match[which(fine_mapped_match$studyType != "gwas"),] #especially relevant for UKB-PPP aproaches. Way better than our previous approach
pqtl_match <- fine_mapped_match[which(fine_mapped_match$studyType == "pqtl"),]

pqtl_match$cis_trans <- ifelse(pqtl_match$isTransQtl == FALSE, "cis", "trans")
pqtl_match$protein <- as.character(unlist(sapply(pqtl_match$studyId, protein_parser)))

cis_pqtl <-pqtl_match[which(pqtl_match$cis_trans == "cis"),]
cis_pqtl$effect_allele <- as.character(unlist(sapply(cis_pqtl$variantId, effect_allele_parser)))
cis_pqtl$other_allele <- as.character(unlist(sapply(cis_pqtl$variantId, other_allele_parser)))

################################################################
#Let's limit the proxy data because it would just be too much!!#
################################################################

proxies <- proxies[which(proxies$id_otargen_1%in%cis_pqtl$variantId | proxies$id_otargen_2%in%cis_pqtl$variantId),]

##########################
#Let's first add the data#
##########################

proxies$protein <- NA
proxies$beta.pqtl <- NA
proxies$standard_error.pqtl <- NA
proxies$pval.pqtl <- NA
proxies$credible_set.pqtl <- NA 
proxies$logBF.pqtl <- NA 

for(index in seq(1, length(proxies$rsID))){
  
  #STEP 1: get the info we want:
  
  variant= proxies$rsID[index]
  query = proxies$query_snp_rsid[index]
  
  chr_38 <- proxies$chr[index]
  bp_38 = proxies$pos_hg38[index]
  
  #Then we query the data in pqtl_df
  
  pqtl_tmp <- cis_pqtl[which(cis_pqtl$chr == chr_38 & cis_pqtl$pos == bp_38),]
  
  if(is_empty(pqtl_tmp$variantId)){
    
    next()
    
  } else {
    
    #######################################################
    #Let's add the information from OTG to align the betas#
    #######################################################
    
    #STEP 1: retrieve info to do alignment (frequencies for lead and proxy)
    
    print(index)
    
    variant_otg <- tryCatch(otargen::variantInfo(variant), error=function(e) NA)
    lead_otg <- tryCatch(otargen::variantInfo(query), error=function(e) NA)
    
    ############################################################################################
    #Careful, if we have a mismatch for some reason, we will add a comment and jump on the next#
    ############################################################################################
    
    if(length(variant_otg) == 1| length(lead_otg) == 1){ #one of them have missing info:
      
      proxies$protein[index] <- paste(pqtl_tmp$protein, collapse = ";")

      proxies$ir_allele.pqtl[index] <- pqtl_tmp$effect_allele[1]
      proxies$is_allele.pqtl[index] <- pqtl_tmp$other_allele[1]
      
      proxies$beta.pqtl[index] <- paste(pqtl_tmp$beta, collapse = ";")
      proxies$beta.pqtl[index] <- paste("non_aligned", proxies$beta.pqtl[index], sep= "_") #adding that the allele might be wrong! Needs to be checked
      
      proxies$standard_error.pqtl[index] <- paste(pqtl_tmp$standardError, collapse = ";")
      proxies$p_value.pqtl[index] <- paste(pqtl_tmp$pValueExponent, collapse = ";")
      proxies$credible_set.pqtl[index] <- paste(pqtl_tmp$credibleSetIndex, collapse = ";")
      proxies$logBF.pqtl[index] <- paste(pqtl_tmp$credibleSetlog10BF, collapse = ";")
      
      next()
      
    }
    
    ############################################
    #If all good, let's extract the AF for both#
    ############################################
    
    #STEP 1: first for the lead:
    
    ir_allele <- ir_new$effect_allele[which(ir_new$variant == query)]
    
    af_ref <- ifelse(ir_allele == lead_otg$altAllele, as.numeric(lead_otg$gnomadNFE), 1-as.numeric(lead_otg$gnomadNFE))
    
    af_ref_res <- ifelse(af_ref < 0.5, "min", "max")
    
    #STEP 2: let's check the frequency for the proxy:
    
    pqtl_allele <- pqtl_tmp$effect_allele[1] # we just need one, all associations are aligned to the same allele. Hence, having only one of the alleles should be more than enough!
    
    af_qtl <- ifelse(pqtl_allele == variant_otg$altAllele, as.numeric(variant_otg$gnomadNFE), 1-as.numeric(variant_otg$gnomadNFE))
    
    af_qtl_res <- ifelse(af_qtl < 0.5, "min", "max")
    
    #####################################################################################################
    #IF af_qtl_res and af_ref_res are the same, we have the beta aligned to the same correlated alleles!#
    #####################################################################################################
    
    if(af_qtl_res == af_ref_res){
      
      proxies$protein[index] <- paste(pqtl_tmp$protein, collapse = ";")
      
      proxies$ir_allele.pqtl[index] <- pqtl_tmp$effect_allele[1]
      proxies$is_allele.pqtl[index] <- pqtl_tmp$other_allele[1]
      proxies$beta.pqtl[index] <- paste(pqtl_tmp$beta, collapse = ";")
      
      proxies$standard_error.pqtl[index] <- paste(pqtl_tmp$standardError, collapse = ";")
      proxies$p_value.pqtl[index] <- paste(pqtl_tmp$pValueExponent, collapse = ";")
      proxies$credible_set.pqtl[index] <- paste(pqtl_tmp$credibleSetIndex, collapse = ";")
      proxies$logBF.pqtl[index] <- paste(pqtl_tmp$credibleSetlog10BF, collapse = ";")
      
      
    } else {
      
      ###########################################
      #But if we don't then we need to change!!!#
      ###########################################
      
      print("exception")
      proxies$protein[index] <- paste(pqtl_tmp$protein, collapse = ";")
      
      proxies$ir_allele.pqtl[index] <- pqtl_tmp$other_allele[1]
      proxies$is_allele.pqtl[index] <- pqtl_tmp$effect_allele[1]
      proxies$beta.pqtl[index] <- paste(pqtl_tmp$beta*(-1), collapse = ";")
      
      proxies$standard_error.pqtl[index] <- paste(pqtl_tmp$standardError, collapse = ";")
      proxies$p_value.pqtl[index] <- paste(pqtl_tmp$pValueExponent, collapse = ";")
      proxies$credible_set.pqtl[index] <- paste(pqtl_tmp$credibleSetIndex, collapse = ";")
      proxies$logBF.pqtl[index] <- paste(pqtl_tmp$credibleSetlog10BF, collapse = ";")
      
    }
    
  }
  
}

###########################################################
#Let's divide the data into clusters and into novel or not#
###########################################################

novel_variants <- ir_new[which(ir_new$reported_ir == "no"),]
novel_variants <- novel_variants[which(novel_variants$ir_source != "/Oliveri_gw"),] #70, just in case

#Now let's add the data for each subset of variants:

#Let's load only the fat distribution variants:

bmi_ns <- fread("output/5_enrichment_analyses/3_go_shifter/ld_bmi_ns.txt")
bmi_ns <- ir_new[which(ir_new$variant%in%bmi_ns$RsIdA),]

bmi_neg <- fread("output/5_enrichment_analyses/3_go_shifter/ld_bmi_neg.txt")
bmi_neg <- ir_new[which(ir_new$variant%in%bmi_neg$RsIdA),]

bmi_pos <- fread("output/5_enrichment_analyses/3_go_shifter/ld_bmi_pos.txt")
bmi_pos <- ir_new[which(ir_new$variant%in%bmi_pos$RsIdA),]

bmi_outlier <- fread("output/5_enrichment_analyses/3_go_shifter/ld_outliers.txt")
bmi_outlier <- ir_new[which(ir_new$variant%in%bmi_outlier$RsIdA),]

###########################
#Let's try adding the data#
###########################

proxies$reported <- NA
proxies$source <- NA
proxies$cluster <- NA

for(index_proxy in seq(1, length(proxies$query_snp_rsid))){
  
  #Variant:
  
  rsid <- proxies$query_snp_rsid[index_proxy]
  
  ir_tmp <- ir_new[which(ir_new$variant == rsid),]
  
  proxies$reported[index_proxy] <- ir_tmp$reported_ir
  proxies$source[index_proxy] <- ir_tmp$ir_source
  
  if(rsid%in%bmi_ns$variant){
    
    proxies$cluster[index_proxy] <- "BMI (P>0.05)"
    
  }
  
  if(rsid%in%bmi_neg$variant){
    
    proxies$cluster[index_proxy] <- "BMI- (P<0.05)"
    
  }
  
  if(rsid%in%bmi_pos$variant){
    
    proxies$cluster[index_proxy] <- "BMI+ (P<0.05)"
    
  }
  
  if(rsid%in%bmi_outlier$variant){
    
    proxies$cluster[index_proxy] <- "BMI+;GSATadjBMI-"
    
  }
  
  
}


##################################################################
#Let's save this dataframe and get some numbers for the main text#
##################################################################

dir.create("output/9_pqtl_results")
dir.create("output/9_pqtl_results/1_cis_pqtl_matches")

fwrite(proxies, "output/9_pqtl_results/1_cis_pqtl_matches/cis_pqtl_matches.txt")

proxies <- fread("output/9_pqtl_results/1_cis_pqtl_matches/cis_pqtl_matches.txt")

#Let's check the results:

length(unique(proxies$protein)) #23
length(unique(proxies$query_snp_rsid)) #20

#############################################
#Let's try getting the instruments for MR!!!#
#############################################

matched_credible_sets_df <- credible_sets_df[which(credible_sets_df$studyId%in%cis_pqtl$studyId),]
matched_credible_sets_df <- matched_credible_sets_df[which(matched_credible_sets_df$isTransQtl == FALSE),]
matched_credible_sets_df$effect_allele <- as.character(unlist(sapply(matched_credible_sets_df$variantId, effect_allele_parser)))
matched_credible_sets_df$other_allele <- as.character(unlist(sapply(matched_credible_sets_df$variantId, other_allele_parser)))

#Let's remove indels:

yes_vect <- c("A", "G", "C", "T")
matched_credible_sets_df <- matched_credible_sets_df[which(matched_credible_sets_df$effect_allele%in%yes_vect & matched_credible_sets_df$other_allele%in%yes_vect),]

#Let's be careful and remove those that belong to the same credible set:

matched_credible_sets_df$prot_cred_set <- paste(matched_credible_sets_df$studyId, "_", matched_credible_sets_df$credibleSetIndex, sep = "")

#Let's order by protein and by pvalue:

matched_credible_sets_df <- matched_credible_sets_df[order(matched_credible_sets_df$credibleSetlog10BF, decreasing = TRUE),]
#matched_credible_sets_df <- matched_credible_sets_df[which(duplicated(matched_credible_sets_df$prot_cred_set) == FALSE),] #all of them are already independent it seems! There is a mistake in the OTG annotation...

matched_4_mr <- as.data.frame(table(matched_credible_sets_df$studyId))

#Let's filter for those that have <3:

matched_credible_sets_df <- matched_credible_sets_df[which(matched_credible_sets_df$studyId%in%matched_4_mr$Var1),] #233 variants!!

#We need to get proxies for each of these fellas. First let's find an rsID for each of them.
#We can do this with LDLink!

matched_credible_sets_df$variant <- NA
matched_credible_sets_df$chr_hg19 <- NA
matched_credible_sets_df$pos_hg19 <- NA
matched_credible_sets_df$af <- NA

#for(index in seq(1, 4)){
#Here is an index of variants that could not be found:
#17_42550742_C_T - monoallelic in European
#11_46311035_A_T - variant not in build 37
#19_50827100_G_A - this one either...
#7_150753877_C_T monoallelic too
#Same with 4_86977143_A_G
#17_41966448_A_G this one too
#7_150606818_A_G
#another one: 12_56876579_C_G
#This one got erased too: rs534242037 
#Out with this one too: 11_47357599_C_T
#12_56126190_C_A
#8_20584119_T_G
#11_47840798_T_C
#10_117519058_C_T
#17_41399455_C_T
#17_41368976_C_T

for(index in seq(233, length(matched_credible_sets_df$variantId))){
  
  print(index)
  print(matched_credible_sets_df$variantId[index])
  
  #STEP 1: get variant:
  
  variant= matched_credible_sets_df$variantId[index]
  variant=as.character(unlist(str_split(variant, "_")))
  variant=paste("chr", variant[1], ":", variant[2], sep = "")

  #STEP 2: get data from OTG
  
  skip_to_next <- FALSE
  # Let's make sure we retrieve the build37 data...
  variant_info <- tryCatch(LDlinkR::LDproxy(variant, pop = "EUR", genome_build = "grch38_high_coverage", token="04cad4ca4374"), error = function(e) { skip_to_next <<- TRUE})
  if(skip_to_next){
    skip_to_next <- FALSE
    next()
  }
  
  info=variant_info[which(variant_info$Coord==variant),]
  
  #STEP 3: retrive the info
  
  rsid <- info$RS_Number
  
  #We get the rsid, let's use the same approach to get the build 37:
  
  skip_to_next_2 <- FALSE
  rsid_info <- tryCatch(LDlinkR::LDproxy(rsid, pop = "EUR", token="04cad4ca4374"), error = function(e) { skip_to_next_2 <<- TRUE})
  if(skip_to_next_2){
    skip_to_next_2 <- FALSE
    next()
  }
  
  final_info = rsid_info[which(rsid_info$Distance==0 & rsid_info$R2==1),]
  
  chr_37 <- as.character(unlist(str_split(final_info$Coord, ":"))[1])
  chr_37 <- as.numeric(as.character(unlist(str_split(chr_37, "chr"))[2]))
  
  bp_37 <- as.numeric(as.character(unlist(str_split(final_info$Coord, ":"))[2]))
  af_ <- final_info$MAF
   
  #And add the data
  
  matched_credible_sets_df$variant[index] <- final_info$RS_Number
  matched_credible_sets_df$chr_hg19[index] <- chr_37
  matched_credible_sets_df$pos_hg19[index] <- bp_37
  matched_credible_sets_df$af[index] <- af_
  
  rm(variant)
  rm(variant_info)
  rm(info)
  rm(final_info)
  
  
}

#This took a long while, let's just save it and update it later if necessary

#matched_credible_sets_df_ <- matched_credible_sets_df %>%
#  select(-c("locus", "qualityControls", "ldSet"))

#fwrite(matched_credible_sets_df_, "output/9_pqtl_results/2_mr/0_data_4_input/all_associations_per_protein_18032026.txt")

######################################################################
#There are 10 that we need to retrieve manually it seems, let's do it#
######################################################################

matched_credible_sets_df <- matched_credible_sets_df %>%
  select(-c("locus", "qualityControls", "ldSet"))

completed_matched_credible_set_df <- matched_credible_sets_df[which(is.na(matched_credible_sets_df$variant) == FALSE),]

#NOTE: I HAVE DOUBLE CHECKED ALL CONVERSION IN DBSNP - ALL IS GOOD

#Let's save this data then, though we are going to run this just in case again

dir.create("output/9_pqtl_results/2_mr/")
dir.create("output/9_pqtl_results/2_mr/0_data_4_input")

fwrite(completed_matched_credible_set_df, "output/9_pqtl_results/2_mr/0_data_4_input/all_associations_per_protein_18032026.txt")

#check <- fread("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/output/9_pqtl_results/2_mr/0_data_4_input/all_associations_per_protein_18032026.txt")

###############################################################
#Now this is the turn to get those that have at least 3 SNPs!!#
###############################################################

#After filtering out rare alleles, that is.

matched_4_mr_clean <- completed_matched_credible_set_df[which(as.numeric(completed_matched_credible_set_df$af) > 0.01 &
                                                                as.numeric(completed_matched_credible_set_df$af < 0.99)),]


matched_4_mr_clean$protein <- as.character(unlist(sapply(matched_4_mr_clean$studyId, protein_parser)))

#Let's get only the UKBB to be able to get proxies:

matched_4_mr_clean_ukbb <- matched_4_mr_clean[which(str_detect(matched_4_mr_clean$studyId, "UKB")),]

#Now check those that have more than 3 signals:

freq_check <- as.data.frame(table(matched_4_mr_clean_ukbb$protein))
freq_check <- freq_check[which(freq_check$Freq > 2),] #18

View(freq_check)

matched_4_mr_clean_ukbb <- matched_4_mr_clean_ukbb[which(matched_4_mr_clean_ukbb$protein%in%freq_check$Var1),]

fwrite(matched_4_mr_clean_ukbb, "output/9_pqtl_results/2_mr/0_data_4_input/all_associations_per_protein_w_more_than_3_credible_sets_18032026.txt")

######################################################
#Let's get the proxies behind each independent signal#
######################################################

fwrite(as.data.frame(matched_4_mr_clean_ukbb$variant), "output/9_pqtl_results/2_mr/0_data_4_input/proxies_4_haploreg_18032026.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = " ")

check <- fread("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/output/9_pqtl_results/2_mr/0_data_4_input/all_associations_per_protein_w_more_than_3_credible_sets_18032026.txt.txt")
