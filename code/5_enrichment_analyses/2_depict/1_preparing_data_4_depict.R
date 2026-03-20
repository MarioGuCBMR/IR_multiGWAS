##############
#INTRODUCTION#
##############

#This code prepares the data to run DEPICT.

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)

###################
#Loading functions#
###################

depict_parser <- function(data_df, snps){
  
  #STEP 1: filter the SNPs:
  
  data_df <- data_df[which(data_df$variant%in%snps),]
  
  #STEP 2 filter for the columns:
  
  depict_df <- data_df %>%
    select(variant, chromosome, base_pair_location, p_value)
  
  colnames(depict_df) <- c("SNP", "Chr", "Pos", "P")
  
  return(depict_df)
  
   
}

###################
#STEP 1: LOAD DATA#
###################

#Get original data for exposure to obtain the samplesizes and good allele frequencies:

path_2_input <- "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/"

setwd(path_2_input)

iradjbmi <- fread("output/2_mv_gwas/fiadjbmi_hdl_tg_common_dwls_curated.txt")

#We are going to use the merged data with BMI since we are working with those.
#The first analyses are done: cluster_all takes all 282 IR variants. 
#Next we are going to run DEPICT with the subset of clusters according to BMI.

all_ir <- fread("output/5_enrichment_analyses/1_prs/1_expo_outcome_df/ir_variants_BMI.txt") #282
gsat_match <- fread("output/5_enrichment_analyses/1_prs/1_expo_outcome_df/ir_variants_GSATadjBMI.txt") #some weird outlier...

#Let's divide them accordingly:

bmi_neg <- all_ir[which(as.numeric(all_ir$beta.outcome) < 0 & as.numeric(all_ir$pval.outcome) < 0.05),]
bmi_pos <- all_ir[which(as.numeric(all_ir$beta.outcome) > 0 & as.numeric(all_ir$pval.outcome) < 0.05),]
bmi_no_sign <- all_ir[which(as.numeric(all_ir$pval.outcome) > 0.05),]

#Be careful with the pleiotropic effects:

gsat_match <- gsat_match[which(gsat_match$SNP%in%bmi_pos$SNP),]
gsat_match <- gsat_match[which(gsat_match$beta.outcome < 0 & gsat_match$pval.outcome < 0.05),] #15

bmi_pos_new <- bmi_pos[which(!(bmi_pos$SNP%in%gsat_match$SNP)),]

################
#Filtering data#
################

cluster_1 <- depict_parser(iradjbmi, bmi_neg$SNP)
cluster_2 <- depict_parser(iradjbmi, bmi_no_sign$SNP)
cluster_3 <- depict_parser(iradjbmi, bmi_pos$SNP)
pleiotropy <- depict_parser(iradjbmi, gsat_match$SNP)

all_combined <- depict_parser(iradjbmi, all_ir$SNP)

#####################
#Let's save the data#
#####################

#Let's create the file:

dir.create("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/output/5_enrichment_analyses/2_depict/")
dir.create("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/output/5_enrichment_analyses/2_depict/input")
dir.create("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/output/5_enrichment_analyses/2_depict/output")

fwrite(cluster_1, "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/output/5_enrichment_analyses/2_depict/input/cluster_bmi_neg.txt", col.names = TRUE, row.names = FALSE, quote=FALSE, sep = "\t")
fwrite(cluster_2, "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/output/5_enrichment_analyses/2_depict/input/cluster_bmi_ns.txt", col.names = TRUE, row.names = FALSE, quote=FALSE, sep = "\t")
fwrite(cluster_3, "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/output/5_enrichment_analyses/2_depict/input/cluster_bmi_pos.txt", col.names = TRUE, row.names = FALSE, quote=FALSE, sep = "\t")
fwrite(pleiotropy, "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/output/5_enrichment_analyses/2_depict/input/pleiotropy.txt", col.names = TRUE, row.names = FALSE, quote=FALSE, sep = "\t")
fwrite(all_combined, "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_functional_annotation/output/3_depict/input/cluster_all.txt", col.names = TRUE, row.names = FALSE, quote=FALSE, sep = "\t")

