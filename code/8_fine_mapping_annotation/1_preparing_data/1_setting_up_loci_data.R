##############
#INTRODUCTION#
##############

#This code sets up the input for fine-mapping of several loci of interest.

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)

##############
#Loading data#
##############

path_2_files <-  "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/" #change it with your own path.

setwd(path_2_files)

ir_new <- fread("output/4_ir_loci_discovery/3_novel_variants/4_colocalization/352_variants_with_all_info_plus_fiadjbmi_and_tg_hdl_coloc.txt")
ir_new <- ir_new[which(ir_new$reported_ir != ""),]

#########################################
#Let's get the loci and see what happens#
#########################################

ir_new$start_ <- as.numeric(ir_new$base_pair_location)-500000
ir_new$end_ <- as.numeric(ir_new$base_pair_location)+500000

ir_new$start_ <- ifelse(as.numeric(ir_new$start_) < 0, 0, ir_new$start_)

################################################################
#Alright, let's go and set up the data needed for running CARMA#
################################################################

ir_new <- ir_new %>%
  select(chromosome, start_, end_, chr_pos) #names got changed. For this version we are gonna run it this way!

dir.create("output/8_fine_mapping")
dir.create("output/8_fine_mapping/all_loci")

fwrite(ir_new, "output/8_fine_mapping/all_loci/all_ir_loci.txt", quote=FALSE, sep = " ", col.names = FALSE, row.names = FALSE)
