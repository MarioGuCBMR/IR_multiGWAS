##############
#INTRODUCTION#
##############

#This code prepares the data for cS2G.

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)

##############
#Loading data#
##############

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/output")

list_of_files <- list.files("8_fine_mapping/all_loci/carma_res_clean/")

for(file_ in list_of_files){
  
  if(!(exists("final_df"))){
    
    final_df <- fread(paste("8_fine_mapping/all_loci/carma_res_clean/", file_, sep=""))
    
  } else {
    
    tmp_df <- fread(paste("8_fine_mapping/all_loci/carma_res_clean/", file_, sep=""))
    
    final_df <- rbind(final_df, tmp_df)
    
    
  }
  
}


####################################################
#Let's remove the duplicates and get the data ready#
####################################################

final_df <- final_df[order(final_df$pip, decreasing = TRUE),]

final_non_dupl_df <- final_df[which(duplicated(final_df$chr_pos) == FALSE),]

#Let's remove those that are not equal to 0:

final_non_dupl_df <- final_non_dupl_df[which(as.numeric(final_non_dupl_df$pip) != 0),]

##########################
#Selecting chr_pos column#
##########################

cs2g_df <- final_non_dupl_df %>%
  select(chr_pos, pip)

colnames(cs2g_df) <- c("Coord", "PIP")

#Let's add the folder:

dir.create("8_fine_mapping/all_loci/gene_prioritization")
dir.create("8_fine_mapping/all_loci/gene_prioritization/input")
dir.create("8_fine_mapping/all_loci/gene_prioritization/output")

#Let's save this data:

fwrite(cs2g_df, "8_fine_mapping/all_loci/gene_prioritization/input/variants_4_cS2G.txt", quote=FALSE, col.names = TRUE, row.names = FALSE, sep = " ")

########################################################
#Let''s mae a dictionary so that I don't go crazy later#
########################################################

fwrite(final_non_dupl_df, "8_fine_mapping/all_loci/gene_prioritization/input/variants_4_cS2G_dictionary.txt")

