##############
#INTRODUCTION#
##############

#This code reads 15-state chromatin data from all cell-types and tissues in ROADMAP and prepares so that it can be run in go shifter:

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)

###################
#Loading functions#
###################

parse_and_save <- function(file_){
  
  print(file_)
  
  #STEP 1: read the file:
  
  e_tmp <- data.table::fread(paste("raw_data/atac_seq_data/", file_, sep = ""))

  #STEP 2: parse for the annotation we want:
  
  promoter_and_enhancer_15 <- c("1_TssA", "2_TssAFlnk", "6_EnhG", "7_Enh", "10_TssBiv", "11_BivFlnk", "12_EnhBiv") #Currin et al (2019)
  
  e_tmp <- e_tmp[which(e_tmp$V4%in%promoter_and_enhancer_15),]
  
  colnames(e_tmp) <- c("#chr", "start", "end", "state", "0", ".","start_2", "end_2", "other")
  
  #STEP 3: save the data: 
  
  e_clean <- e_tmp %>%
    dplyr::select("#chr","start","end")
  
  fwrite(unique(e_clean), paste("output/5_enrichment_analyses/3_go_shifter/curated_atac_seq/", file_, sep = ""), sep = "\t", col.names = F) #fwrite is able to compress without a problem
  
}

#####################
#Let's read the data#
#####################

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/")

list_of_files <- list.files("raw_data/atac_seq_data/")
list_of_files <- list_of_files[which(str_detect(list_of_files, "15_coreMarks"))]

##########################
#Let's loop over the data#
##########################

for(i in seq(1, length(list_of_files))){
  
  parse_and_save(list_of_files[i])
  
}
