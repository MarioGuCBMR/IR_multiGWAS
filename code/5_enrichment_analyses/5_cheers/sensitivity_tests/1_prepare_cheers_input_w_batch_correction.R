##############
#INTRODUCTION#
##############

#This code will prepare the input data for CHEERS. We are doing this for the other metabolic clusters from Lotta and Suzuki et al.
#We do not need to normalize the data, since that has already been done before!
#We just need to create new inputs!

###################
#Load the packages#
###################

library(tidyverse)
library(data.table)
library(DESeq2)
library(limma)
library(data.table)

###################
#Loading functions#
###################

adding_b37_otg <- function(proxies_){
  #This code takes the haploreg output, searches for proxies in OTG and retrieves their chromosome and base pair location in build 37
  
  #######################################################
  #STEP 1: create columns where to store the information#
  #######################################################
  
  proxies_$chr_37 <- NA
  proxies_$bp_37 <- NA
  
  #####################################
  #STEP 2: let's loop over the proxies#
  #####################################
  
  for(index in seq(1, length(proxies_$rsID))){
    
    print(proxies_$rsID[index])
    
    variant= proxies_$rsID[index]
    query = proxies_$query_snp_rsid[index]
    
    skip_to_next <- FALSE
    
    # Let's make sure we retrieve the build37 data...
    
    variant_info <- tryCatch(otargen::variantInfo(variant), error = function(e) { skip_to_next <<- TRUE})
    
    #If we skip, the chromosome and base_pair_location will remain = NA
    
    if(skip_to_next){
      
      skip_to_next <- FALSE
      next()
      
    } else {
      
      proxies_$chr_37[index] <- variant_info$chromosomeB37
      proxies_$bp_37[index] <- variant_info$positionB37
      
    }
    
    #And just in case let's remove this:
    
    rm(variant_info)
    
  }
  
  #Finally, let's get back the proxies dataframe:
  
  return(proxies_)
  
}

formatting_data_4_cheers <- function(data_df){
  
  #STEP 0: be careful we might have weird rsID;
  
  data_df <- data_df[which(data_df$rsID != ""),]
  
  #STEP 1: let's add the chromosome and base_pair_location to all the proxies:
  
  data_df_b37 <- adding_b37_otg(data_df)

  #STEP 2: let's clean this data from NAs
  
  data_df_b37_clean <- data_df_b37[which(is.na(data_df_b37$bp_37) == FALSE),]
  
  #sTEP 3: reformat the data:
  
  data_df_b37_clean <- data_df_b37_clean %>% dplyr::select(rsID, chr_37, bp_37)
  data_df_b37_clean$chr_37=paste0("chr",data_df_b37_clean$chr_37)
  
  #STEP 4: return the cleaned dataframe:
  
  return(data_df_b37_clean)
  
}

##########################################################################
#STEP 1: Let's prepare the data for the other lipodystrophy-like clusters#
##########################################################################

path_2_input <- "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/"

#Let's load the data:

lotta <- fread("output/5_enrichment_analyses/5_cheers/3_sensitivity_tests/0_proxies/lotta_haploreg_output.txt")
lipodystrophy <- fread("output/5_enrichment_analyses/5_cheers/3_sensitivity_tests/0_proxies/lipodystrophy_haploreg_output.txt")

beta_cell_neg <- fread("output/5_enrichment_analyses/5_cheers/3_sensitivity_tests/0_proxies/beta_cell_neg_haploreg_output.txt", fill=TRUE)
beta_cell_pos <- fread("output/5_enrichment_analyses/5_cheers/3_sensitivity_tests/0_proxies/beta_cell_pos_haploreg_output.txt", fill=TRUE)
residual <- fread("output/5_enrichment_analyses/5_cheers/3_sensitivity_tests/0_proxies/residual_haploreg_output.txt", fill=TRUE)

metabolic_syndrome <- fread("output/5_enrichment_analyses/5_cheers/3_sensitivity_tests/0_proxies/metabolic_syndrome_haploreg_output.txt", fill=TRUE)
body_fat <- fread("output/5_enrichment_analyses/5_cheers/3_sensitivity_tests/0_proxies/body_fat_haploreg_output.txt", fill=TRUE)

obesity <- fread("output/5_enrichment_analyses/5_cheers/3_sensitivity_tests/0_proxies/obesity_haploreg_output.txt", fill=TRUE)
liver <- fread("output/5_enrichment_analyses/5_cheers/3_sensitivity_tests/0_proxies/liver_haploreg_output.txt", fill=TRUE)

#Let's get the clean data and save it:

dir.create("output/5_enrichment_analyses/5_cheers/3_sensitivity_tests/1_cheers_input")
dir.create("output/5_enrichment_analyses/5_cheers/3_sensitivity_tests/1_cheers_input/variant_data")

setwd("output/5_enrichment_analyses/5_cheers/3_sensitivity_tests/1_cheers_input/variant_data/")

lotta_b37_clean <- formatting_data_4_cheers(lotta)
lipodystrophy_b37_clean <- formatting_data_4_cheers(lipodystrophy)

beta_cell_neg_b37_clean <- formatting_data_4_cheers(beta_cell_neg)
beta_cell_pos_b37_clean <- formatting_data_4_cheers(beta_cell_pos)
residual_b37_clean <- formatting_data_4_cheers(residual)

metabolic_syndrome_b37_clean <- formatting_data_4_cheers(metabolic_syndrome)
body_fat_b37_clean <- formatting_data_4_cheers(body_fat)

obesity_b37_clean <- formatting_data_4_cheers(obesity)
liver_b37_clean <- formatting_data_4_cheers(liver)

#######################################
#Let's save all of this fantastic data#
#######################################

write.table(lotta_b37_clean, file = "proxies_lotta.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

write.table(lipodystrophy_b37_clean, file = "proxies_lipodystrophy.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(metabolic_syndrome_b37_clean, file = "proxies_metabolic_syndrome.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(body_fat_b37_clean, file = "proxies_body_fat.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

write.table(beta_cell_neg_b37_clean, file = "proxies_beta_cell_neg.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(beta_cell_pos_b37_clean, file = "proxies_beta_cell_pos.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(residual_b37_clean, file = "proxies_residual.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")

write.table(obesity_b37_clean, file = "proxies_obesity.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
write.table(liver_b37_clean, file = "proxies_liver.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")


