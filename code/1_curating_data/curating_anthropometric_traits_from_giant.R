##############
#INTRODUCTION#
##############

#This is a code to curate GIANT data without UKBB for the following traits:

#BMI,
#WHR; WHRadjBMI
#WC, WCadjBMI,
#HC, HCadjBMI

#All of them in their sex-combined, male and female versions. 

#Importantly, we won't be removing either MAF nor MHC regions. Those were done in the exposure IR data already.
#MAF is not excluded because we have greater sample size in the exposure, making it a better approach for frequencies.
#MHC is already removed in exposures, so it won't affect us either way. 

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
  
  if(length(og_df_dupl$variant) == 0){
    
      return(og_df_non_dupl_match)    
    
  }
  
  #Let's get the data from the dictionary
  
  dict_non_dupl <- dict_[which(duplicated(dict_$variant) == FALSE),]
  dict_match_dupl <- dict_non_dupl[which(dict_non_dupl$variant%in%og_df_dupl$variant),] 
  
  #And now match it again:
  
  og_df_dupl <- og_df_dupl[which(og_df_dupl$variant%in%dict_match_dupl$variant),]
  
  #We can have several duplicates due to multialleic, so let's make a quick loop:
  
  og_df_dupl$chromosome <- NA
  og_df_dupl$base_pair_location <- NA
  og_df_dupl$chr_pos <- NA
  
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

recursively_curate_bmi <- function(list_of_files, final_names, dict){
  
  #STEP 0: start loop to recursively save data:
  
  for(files_index in seq(1, length(list_of_files))){
    
    tmp_df <- fread(paste("raw_data/anthropometric_traits_giant/", list_of_files[files_index], sep =""))
    
    #Format correctly...
    
    colnames(tmp_df) <- c("variant", "effect_allele", "other_allele", "effect_allele_frequency", "beta", "standard_error", "p_value", "sample_size")
    
    #Now let's add the chr_pos
    
    tmp_df_w_chr_pos <- add_build_37(tmp_df, dict)
    
    fwrite(tmp_df_w_chr_pos, paste("output/1_curated_data/", final_names[files_index], sep =""))
    
  }
  
}

recursively_curate_fat<- function(list_of_files, final_names, dict){
  
  #STEP 0: start loop to recursively save data:
  
  for(files_index in seq(1, length(list_of_files))){
    
    tmp_df <- fread(paste("raw_data/anthropometric_traits_giant/", list_of_files[files_index], sep =""))
    
    #Format correctly...
    
    if(length(colnames(tmp_df)) > 8){
      
      tmp_df <- tmp_df %>%
        select(-c("Chr", "Pos"))
      
    }
    
    colnames(tmp_df) <- c("variant", "effect_allele", "other_allele", "effect_allele_frequency", "beta", "standard_error", "p_value", "sample_size")
    
    tmp_df_w_chr_pos <- add_build_37(tmp_df, dict)
    
    fwrite(tmp_df_w_chr_pos, paste("output/1_curated_data/", final_names[files_index], sep =""))
    
  }
  
}

#############################
#Preparing data for curation#
#############################

#We are gonna load the tg_hdl_ratio from 2024. 

project_path <- "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2023/" #change it with your own path.

setwd(project_path)

files_2_curate <- list.files("raw_data/anthropometric_traits_giant/")
files_2_curate <- files_2_curate[which(files_2_curate != "raw")]

#################################################################
#Let's add the dictionary so that we can add as much information#
#################################################################

dict <- fread("output/1_curated_data/build_36_2_37_for_giant_only_ss.txt")

#######################################################
#Let's split the data between bmi and fat distribution#
#######################################################

#We are splitting the data like this because they come from different studies and, thus, have different columns.

bmi_2_curate <- files_2_curate[which(str_detect(files_2_curate, "mc_merge_nogc"))]

bmi_files <- c("bmi_giant_male_curated.txt", "bmi_giant_curated.txt", "bmi_giant_female_curated.txt") #Let's put them in the same order as they appear so that we can loop through indexes without an issue. 

recursively_curate_bmi(list_of_files=bmi_2_curate, final_names=bmi_files, dict=dict)

########################
#Let's do a quick check#
########################

bmi_combined <- fread("output/1_curated_data/bmi_giant_curated.txt")
bmi_female <- fread("output/1_curated_data/bmi_giant_female_curated.txt")
bmi_male <- fread("output/1_curated_data/bmi_giant_male_curated.txt")

head(bmi_combined) #checked with original dataframe! Perfect
head(bmi_male) #checked with original dataframe! Perfect
head(bmi_female) #checked with original dataframe! Perfect

################################################
#Let's do the same with fat distribution traits#
################################################

fat_2_curate <- files_2_curate[which(str_detect(files_2_curate, "mc_merge_nogc") == FALSE)]

fat_files <- c("hc_giant_curated.txt",       "hc_giant_female_curated.txt",        "hc_giant_male_curated.txt",  "hcadjbmi_giant_curated.txt",  "hcadjbmi_giant_female_curated.txt",  "hcadjbmi_giant_male_curated.txt",    "wc_giant_curated.txt",       
               "wc_giant_female_curated.txt",         "wc_giant_male_curated.txt",           "wcadjbmi_giant_curated.txt",  "wcadjbmi_giant_female_curated.txt",   "wcadjbmi_giant_male_curated.txt",     "whr_giant_curated.txt",      
               "whr_giant_female_curated.txt",        "whr_giant_male_curated.txt",          "whradjbmi_giant_curated.txt", "whradjbmi_giant_female_curated.txt",  "whradjbmi_giant_male_curated.txt") #Let's put them in the same order as they appear so that we can loop through indexes without an issue. 

recursively_curate_fat(list_of_files=fat_2_curate, final_names=fat_files, dict=dict)

########################
#Let's do a quick check#
########################

whr_combined <- fread("output/1_curated_data/whr_giant_curated.txt")
whr_female <- fread("output/1_curated_data/whr_giant_female_curated.txt")
whr_male <- fread("output/1_curated_data/whr_giant_male_curated.txt")

head(whr_combined) #checked with original dataframe! Perfect
head(whr_male) #checked with original dataframe! Perfect
head(whr_female) #checked with original dataframe! Perfect

#Funny p-vals, but they are like that in the original! Nothing wrong here!
