##############
#INTRODUCTION#
##############

#This code annotates the fine-mapped variants that pass cS2G filtering.

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)
library(arrow)

###################
#Loading libraries#
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

###############
#Load the data#
###############

path_2_input <- "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/"

setwd(path_2_input)

#Let's load the IR variants:

ir_variants <- fread("output/8_fine_mapping/all_loci/seed_4_annotation.txt")

#Get the fine-mapped data

fine_mapped <- fread("output/8_fine_mapping/all_loci/gene_prioritization/input/variants_4_cS2G_dictionary.txt")

############################################################
#Let's load the prioritized variants and see what we can do#
############################################################

fine_mapped <- fine_mapped[which(as.numeric(fine_mapped$pip) > 0.8),] #43
fine_mapped$annotation <- NA

for(i in seq(1, length(fine_mapped$variant))){
  
  #STEP 1: get the rsID
  
  rsid <- fine_mapped$variant[i]
  
  #STEP 2:
  
  tmp_df  <- tryCatch({
    otargen::variantInfo(rsid)
  }, error = function(e) {
    return(NULL)  # Or return a tibble with NA fields for robustness
  })
  
  #STEP 3: get the most severe consequence:
  
  fine_mapped$annotation[i] <- ifelse(is.null(tmp_df), NA, tmp_df$mostSevereConsequence)
  
}

#Let's check what we have:

table(fine_mapped$annotation)

#3_prime_UTR_variant       5_prime_UTR_variant   downstream_gene_variant        intergenic_variant 
#3                         2                         1                         8 
#intron_variant          missense_variant regulatory_region_variant     upstream_gene_variant 
#18                         2                         3                         6 

#Let's focus on the missense variant, nothing more we can do:

yes_vect<-  c("missense_variant")
fine_mapped <- fine_mapped[which(fine_mapped$annotation%in%yes_vect),]

############################
#NEWER VERSION for OPT DATA# just checking annotations
############################

mapped_path <- "raw_data/opt_variant_data/"

variant_df <- as.data.frame(merge_parquet_files(mapped_path))
