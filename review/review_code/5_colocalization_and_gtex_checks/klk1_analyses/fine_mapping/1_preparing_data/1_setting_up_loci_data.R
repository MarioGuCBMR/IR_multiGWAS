##############
#INTRODUCTION#
##############

#This code sets up the input for fine-mapping of the KLK1 locus.

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)

##############
#Loading data#
##############

all_variants=readxl::read_xlsx("/projects/kilpelainen-AUDIT/people/zlc436/IR_GSEM_2025/review/manuscript/03062026/Supplementary Tables_20March2026.xlsx", sheet = 3)
all_variants = as.data.frame(all_variants)
all_variants = all_variants[2:nrow(all_variants),] #remove extra first row that explains stuff
colnames(all_variants) = all_variants[1,]  #change column names
all_variants = all_variants[2:nrow(all_variants),] #remove columns names in the first line now that we have used them

#We do not need to many columns AND we have functions ready for a certain subset of variables.

all_variants = all_variants[,1:12]
colnames(all_variants) = c("variant", "chromosome", "base_pair_location", "minimum_allele_frequency", "effect_allele", "other_allele", "novel_or_reported", "bmi_subgroup", "VEP", "beta", "standard_error", "pvalue")
all_variants$sample_size = 2e06 #an approximation, not needed for the actual analyses, just to make the functions work

#############################
#Let's get the KLK1 variant!#
#############################

klk1_lead = all_variants[which(all_variants$variant=="rs485186"),]

#########################################
#Let's get the loci and see what happens#
#########################################

klk1_lead$start_ <- as.numeric(klk1_lead$base_pair_location)-500000
klk1_lead$end_ <- as.numeric(klk1_lead$base_pair_location)+500000

klk1_lead$start_ <- ifelse(as.numeric(klk1_lead$start_) < 0, 0, klk1_lead$start_)

################################################################
#Alright, let's go and set up the data needed for running CARMA#
################################################################

klk1_lead$chr_pos = paste("chr", klk1_lead$chromosome, ":", klk1_lead$base_pair_location, sep = "")

klk1_lead <- klk1_lead %>%
  select(chromosome, start_, end_, chr_pos) #names got changed. For this version we are gonna run it this way!

dir.create("/maps/projects/kilpelainen-AUDIT/people/zlc436/IR_GSEM_2025/review/review_output/4_colocalization_and_gtex_checks/klk1_analyses")
dir.create("/maps/projects/kilpelainen-AUDIT/people/zlc436/IR_GSEM_2025/review/review_output/4_colocalization_and_gtex_checks/klk1_analyses/all_loci")

fwrite(klk1_lead, "/maps/projects/kilpelainen-AUDIT/people/zlc436/IR_GSEM_2025/review/review_output/4_colocalization_and_gtex_checks/klk1_analyses/all_loci/all_ir_loci.txt", quote=FALSE, sep = " ", col.names = FALSE, row.names = FALSE)
