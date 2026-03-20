##############
#INTRODUCTION#
##############

#This is a code makes a dataframe with chr and positions in build 36 and 37 so that we can properly clean all the curated data. 

#HC sex-combined GWAS has all the postiions for the summary statistics, making it easier for us to retrieve the rest.

###########
#libraries#
###########

library(data.table)
library(tidyverse)

###################
#Loading functions#
###################

rsid_parser <- function(id){
  
  rsid <- str_split(id, ":")[[1]][1]
  
  return(rsid)
  
}

#################################################################################################
#STRATEGY: WE ARE GOING TO USE THE BUILD 37 DATA FROM THE META-ANALYSIS WITH UKBB AS FIRST CHECK#
#################################################################################################

project_path <- "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2023/" #change it with your own path.

setwd(project_path)

ref <- fread("../../mho_variants_common_info/raw_data/whradjbmi.giant-ukbb.meta-analysis.combined.23May2018.txt.gz")
hcadjbmi <- fread("raw_data/anthropometric_traits_giant/GIANT_2015_HIPadjBMI_COMBINED_EUR.txt")

######################################
#Let's clean the reference one second#
######################################

ref <- ref[which(is.na(ref$INFO) == FALSE)] #these have the chromosome and base_pair_location wrong
ref$variant <- as.character(unlist(sapply(ref$SNP, rsid_parser)))

###################################
#Let's prepare the dictionary here#
###################################

dict <- hcadjbmi %>%
  select(MarkerName, Chr, Pos)

colnames(dict) <- c("variant", "chromosome", "base_pair_location_18")

################################
#Let's perform a first matching#
################################

ref_2_dict <- ref[which(ref$variant%in%dict$variant),]
dict_2_ref <- dict[which(dict$variant%in%ref$variant),]

#Reference has multi-allelic SNPs. 
#We can remove them, we just want the positions.

ref_2_dict <- ref_2_dict[which(duplicated(ref_2_dict$variant) == FALSE),]

dict_2_ref <- dict_2_ref[order(match(dict_2_ref$variant, ref_2_dict$variant)),]

#Let's check that they are correctly aligned:

length(which(dict_2_ref$variant == ref_2_dict$variant)) #perfect match

#Finally add:

dict_2_ref$base_pair_location_37 <- ref_2_dict$POS

#########################################################################
#Let's work with the rest of the variants - we will put them in liftover#
#########################################################################

dict_miss <- dict[which(!(dict$variant%in%ref$variant)),] #55K are missing. That is not so much!

data_4_liftover <- dict_miss

data_4_liftover$start <-   as.numeric(data_4_liftover$base_pair_location_18)-1
data_4_liftover$end <-   as.numeric(data_4_liftover$base_pair_location_18)+1

#Final dataframe:

data_4_liftover <- data_4_liftover %>%
  select(chromosome, start, end)

colnames(data_4_liftover) <- c("chrom", "start", "end")

#Let's save the data:

fwrite(as.data.frame(data_4_liftover), "output/1_curated_data/input_4_liftover.txt", col.names = FALSE, row.names = FALSE, quote = FALSE, sep = " ")

#It seems that  54,676/54,685 do not even exist in build 37. They were deleted afterwards.

#With these numbers I think it is save to say that we can save the dict_2_ref dataframe and work with that without an issue.

###################
#Saving dictionary#
###################

fwrite(dict_2_ref, "output/1_curated_data/build_36_2_37_for_giant_only_ss.txt")
