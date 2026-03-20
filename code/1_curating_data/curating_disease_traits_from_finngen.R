##############
#INTRODUCTION#
##############

#This is a code to curate FinnGen data!
#We are going to add the cnotrols, cases and prevalences from Ristey R9!

###########
#libraries#
###########

library(data.table)
library(tidyverse)

###################
#Loading functions#
###################

add_build_37 <- function(og_df, dict_){
  
  #STEP 0: make dummy variables so that we can build the function:
  
  #og_df <- tmp_df
  #dict_ <- dict
  
  #STEP 1: match data with dictionary for non-duplicates
  
  og_df_non_dupl <- og_df[which(duplicated(og_df$variant) == FALSE),]

  og_df_non_dupl_match <- og_df_non_dupl[which(og_df_non_dupl$variant%in%dict_$variant),]
  dict_match <- dict_[which(dict_$variant%in%og_df_non_dupl$variant),]
  dict_ordered <- dict_match[order(match(dict_match$variant, og_df_non_dupl_match$variant)),]
  
  #Let's check the data:
  
  length(which(dict_ordered$variant == og_df_non_dupl_match$variant)) #all
  
  #Let's add the chr_pos:
  
  og_df_non_dupl_match$chromosome <- dict_ordered$chromosome
  og_df_non_dupl_match$base_pair_location <-  dict_ordered$base_pair_location_37
  og_df_non_dupl_match$chr_pos <- paste("chr", dict_ordered$chromosome, ":", dict_ordered$base_pair_location_37, sep = "")
  
  #STEP 2: let's work with the duplicates now:
  
  og_df_dupl <- og_df[which(duplicated(og_df$variant) == TRUE),]
  
  #Let's get the data from the dictionary
  
  dict_non_dupl <- dict_[which(duplicated(dict_$variant) == FALSE),]
  dict_match_dupl <- dict_non_dupl[which(dict_non_dupl$variant%in%og_df_dupl$variant),] 
  
  #And now match it again:
  
  og_df_dupl <- og_df_dupl[which(og_df_dupl$variant%in%dict_match_dupl$variant),]
  
  #We can have several duplicates due to multialleic, so let's make a quick loop:
  
  og_df_dupl$chr_pos <- NA
  og_df_dupl$chromosome <- NA
  og_df_dupl$base_pair_location <- NA
  
  for(index_variant in seq(1, length(og_df_dupl$variant))){
    
    #First let's find the data from the dictionary:
    
    dict_tmp <- dict_match_dupl[which(dict_match_dupl$variant == og_df_dupl$variant[index_variant]),]
    
    #the dictionary should not have any duplicates so this works out:
    
    og_df_dupl$chromosome[index_variant] <- dict_tmp$chromosome
    og_df_dupl$base_pair_location[index_variant] <-  dict_tmp$base_pair_location_37
    og_df_dupl$chr_pos[index_variant] <-  paste("chr", dict_tmp$chromosome, ":", dict_tmp$base_pair_location_37, sep = "")
    
  }
  
  #STEP 3: let's bind everything:
  
  final_df <- rbind(og_df_non_dupl_match, og_df_dupl)
  
  return(final_df)

}

recursively_curate <- function(list_of_files, final_names, dict, sample_size, sample_cases, sample_controls, prevalence){
  
  #STEP 0: start loop to recursively save data:
  
  for(files_index in seq(1, length(list_of_files))){
    
    tmp_df <- fread(paste("raw_data/", list_of_files[files_index], sep =""))
    
    colnames(tmp_df) <- c("chromosome", "base_pair_location", "other_allele", "effect_allele", "variant", "gene", "p_value", "logpp", "beta", "standard_error", "effect_allele_frequency", "alt_cases", "alt_controls")
    
    tmp_df <- tmp_df %>%
      select(variant, effect_allele, other_allele, effect_allele_frequency, beta, standard_error, p_value)
    
    tmp_df_w_chr_pos <- add_build_37(tmp_df, dict)
    
    tmp_df_w_chr_pos$sample_cases <- sample_cases[files_index]
    tmp_df_w_chr_pos$sample_controls <- sample_controls[files_index]
    tmp_df_w_chr_pos$prevalence <- prevalence[files_index]
    
    fwrite(tmp_df_w_chr_pos, paste("output/1_curated_data/", final_names[files_index], sep =""))
    
  }
  
}

#############################
#Preparing data for curation#
#############################

#We are gonna load the tg_hdl_ratio from 2024. 

project_path <- "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2023/" #change it with your own path.

setwd(project_path)

files_2_curate <- list.files("raw_data")
files_2_curate <- files_2_curate[which(str_detect(files_2_curate, "finngen"))]

#################################################################
#Let's add the dictionary so that we can add as much information#
#################################################################

dict <- fread("output/1_curated_data/build_38_2_37_for_finngen_only_ss.txt")

#######################################################
#Let's split the data between bmi and fat distribution#
#######################################################

#We are splitting the data like this because they come from different studies and, thus, have different columns.

files_ <- c("pcos_curated.txt", "chd_curated.txt", 
            "hypertension_curated.txt", "ckd_curated.txt",
            "nafld_curated.txt", "t2d_curated.txt") #Let's put them in the same order as they appear so that we can loop through indexes without an issue. 

cases_samples <- c(34388,
                   46959,
                   122996,
                   10039,
                   2568,
                   65085)

control_samples <- c(195922,
                     365222,
                     289117,
                     396706,
                     409613,
                     335112)

total_sample <- cases_samples+control_samples

prevalences <- c(15.05,
                 14.52,
                 31.11,
                 2.98,
                 0.66,
                 16.65)

recursively_curate(list_of_files=files_2_curate, final_names=files_, dict=dict, sample_size = total_sample, sample_cases = cases_samples, sample_controls = control_samples, prevalence = prevalences)

########################
#Let's do a quick check#
########################

ckd <- fread("output/1_curated_data/ckd_curated.txt")

head(ckd) #checked with original dataframe! Perfect

