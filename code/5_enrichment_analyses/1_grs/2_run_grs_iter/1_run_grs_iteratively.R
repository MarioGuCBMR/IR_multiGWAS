##############
#INTRODUCTION#
##############

#This code runs MR iteratively for all outcomes!!

###################
#LOADING LIBRARIES#
###################

library(tidyverse)
library(data.table)
library(TwoSampleMR)
library(gtx)
library(meta)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(rmarkdown)
library(data.table)
library(LDlinkR)

###################
#STEP 1: LOAD DATA#
###################

#Get original data for exposure to obtain the samplesizes and good allele frequencies:

project_path <- "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/" #change it with your own path.

setwd(project_path)

list_of_files <- list.files("output/5_enrichment_analyses/1_prs/1_expo_outcome_df/")

dir.create("output/5_enrichment_analyses/1_prs/2_prs")
output_path <- "output/5_enrichment_analyses/1_prs/2_prs/"

##############################
#STEP 1: LET's START THE LOOP#
##############################

for(index_file in seq(1, length(list_of_files))){
  
  #STEP 1.1: load the data:
  
  grs_df <- fread(paste("output/5_enrichment_analyses/1_prs/1_expo_outcome_df//", list_of_files[index_file], sep=""))
  
  #STEP 1.2: compute PRS. Weighted and unweighted:
  
  weighted_grs_res <- gtx::grs.summary(w=grs_df$beta.exposure, b=grs_df$beta.outcome, s = grs_df$se.outcome, grs_df$samplesize.outcome)
  unweighted_grs_res <- gtx::grs.summary(w=rep(1, length(grs_df$beta.exposure)), b=grs_df$beta.outcome, s = grs_df$se.outcome, grs_df$samplesize.outcome)
  
  saveRDS(weighted_grs_res, paste(output_path, "weighted_", gsub(".txt", ".RDS", list_of_files[index_file]), sep=""))
  saveRDS(unweighted_grs_res, paste(output_path, "unweighted_", gsub(".txt", ".RDS", list_of_files[index_file]), sep=""))

}

############################################################################################################
#STEP 2: let's run the same analyses but selecting the associations for BMI according to their effect sizes#
############################################################################################################

bmi_match <- fread("output/5_enrichment_analyses/1_prs/1_expo_outcome_df/ir_variants_BMI.txt")
gsat_match <- fread("output/5_enrichment_analyses/1_prs/1_expo_outcome_df/ir_variants_GSATadjBMI.txt") #some weird outlier...

#Let's divide them accordingly:

bmi_neg <- bmi_match[which(as.numeric(bmi_match$beta.outcome) < 0 & as.numeric(bmi_match$pval.outcome) < 0.05),]
bmi_pos <- bmi_match[which(as.numeric(bmi_match$beta.outcome) > 0 & as.numeric(bmi_match$pval.outcome) < 0.05),]
bmi_no_sign <- bmi_match[which(as.numeric(bmi_match$pval.outcome) > 0.05),]

##############################
#Let's loop the negative data#
##############################

dir.create("output/5_enrichment_analyses/1_prs/2_prs_bmi_neg")
output_path <- "output/5_enrichment_analyses/1_prs/2_prs_bmi_neg/"

for(index_file in seq(1, length(list_of_files))){
  
  #STEP 1.1: load the data:
  
  grs_df <- fread(paste("output/5_enrichment_analyses/1_prs/1_expo_outcome_df//", list_of_files[index_file], sep=""))
  grs_df <- grs_df[which(grs_df$SNP%in%bmi_neg$SNP),]
  
  #STEP 1.2: compute PRS. Weighted and unweighted:
  
  weighted_grs_res <- gtx::grs.summary(w=grs_df$beta.exposure, b=grs_df$beta.outcome, s = grs_df$se.outcome, grs_df$samplesize.outcome)
  unweighted_grs_res <- gtx::grs.summary(w=rep(1, length(grs_df$beta.exposure)), b=grs_df$beta.outcome, s = grs_df$se.outcome, grs_df$samplesize.outcome)
  
  saveRDS(weighted_grs_res, paste(output_path, "weighted_", gsub(".txt", ".RDS", list_of_files[index_file]), sep=""))
  saveRDS(unweighted_grs_res, paste(output_path, "unweighted_", gsub(".txt", ".RDS", list_of_files[index_file]), sep=""))
  
}

###################
#Let's do positive#
###################

dir.create("output/5_enrichment_analyses/1_prs/2_prs_bmi_pos")
output_path <- "output/5_enrichment_analyses/1_prs/2_prs_bmi_pos/"

for(index_file in seq(1, length(list_of_files))){
  
  #STEP 1.1: load the data:
  
  grs_df <- fread(paste("output/5_enrichment_analyses/1_prs/1_expo_outcome_df//", list_of_files[index_file], sep=""))
  grs_df <- grs_df[which(grs_df$SNP%in%bmi_pos$SNP),]
  
  #STEP 1.2: compute PRS. Weighted and unweighted:
  
  weighted_grs_res <- gtx::grs.summary(w=grs_df$beta.exposure, b=grs_df$beta.outcome, s = grs_df$se.outcome, grs_df$samplesize.outcome)
  unweighted_grs_res <- gtx::grs.summary(w=rep(1, length(grs_df$beta.exposure)), b=grs_df$beta.outcome, s = grs_df$se.outcome, grs_df$samplesize.outcome)
  
  saveRDS(weighted_grs_res, paste(output_path, "weighted_", gsub(".txt", ".RDS", list_of_files[index_file]), sep=""))
  saveRDS(unweighted_grs_res, paste(output_path, "unweighted_", gsub(".txt", ".RDS", list_of_files[index_file]), sep=""))
  
}

##############################
#Let's do positive wo outlier#
##############################

dir.create("output/5_enrichment_analyses/1_prs/2_prs_bmi_pos_wo_outlier")
output_path <- "output/5_enrichment_analyses/1_prs/2_prs_bmi_pos_wo_outlier/"

gsat_match <- gsat_match[which(gsat_match$SNP%in%bmi_pos$SNP),]
gsat_match <- gsat_match[which(gsat_match$beta.outcome < 0 & gsat_match$pval.outcome < 0.05),] #15

bmi_pos_new <- bmi_pos[which(!(bmi_pos$SNP%in%gsat_match$SNP)),]

for(index_file in seq(1, length(list_of_files))){
  
  #STEP 1.1: load the data:
  
  grs_df <- fread(paste("output/5_enrichment_analyses/1_prs/1_expo_outcome_df//", list_of_files[index_file], sep=""))
  grs_df <- grs_df[which(grs_df$SNP%in%bmi_pos_new$SNP),]
  
  #STEP 1.2: compute PRS. Weighted and unweighted:
  
  weighted_grs_res <- gtx::grs.summary(w=grs_df$beta.exposure, b=grs_df$beta.outcome, s = grs_df$se.outcome, grs_df$samplesize.outcome)
  unweighted_grs_res <- gtx::grs.summary(w=rep(1, length(grs_df$beta.exposure)), b=grs_df$beta.outcome, s = grs_df$se.outcome, grs_df$samplesize.outcome)
  
  saveRDS(weighted_grs_res, paste(output_path, "weighted_", gsub(".txt", ".RDS", list_of_files[index_file]), sep=""))
  saveRDS(unweighted_grs_res, paste(output_path, "unweighted_", gsub(".txt", ".RDS", list_of_files[index_file]), sep=""))
  
}

######################################
#Let's do positive with only outliers#
######################################

dir.create("output/5_enrichment_analyses/1_prs/2_prs_bmi_pos_outlier")
output_path <- "output/5_enrichment_analyses/1_prs/2_prs_bmi_pos_outlier/"

bmi_pos_outlier <- bmi_pos[which(bmi_pos$SNP%in%gsat_match$SNP),]

for(index_file in seq(1, length(list_of_files))){
  
  #STEP 1.1: load the data:
  
  grs_df <- fread(paste("output/5_enrichment_analyses/1_prs/1_expo_outcome_df//", list_of_files[index_file], sep=""))
  grs_df <- grs_df[which(grs_df$SNP%in%bmi_pos_outlier$SNP),]
  
  #STEP 1.2: compute PRS. Weighted and unweighted:
  
  weighted_grs_res <- gtx::grs.summary(w=grs_df$beta.exposure, b=grs_df$beta.outcome, s = grs_df$se.outcome, grs_df$samplesize.outcome)
  unweighted_grs_res <- gtx::grs.summary(w=rep(1, length(grs_df$beta.exposure)), b=grs_df$beta.outcome, s = grs_df$se.outcome, grs_df$samplesize.outcome)
  
  saveRDS(weighted_grs_res, paste(output_path, "weighted_", gsub(".txt", ".RDS", list_of_files[index_file]), sep=""))
  saveRDS(unweighted_grs_res, paste(output_path, "unweighted_", gsub(".txt", ".RDS", list_of_files[index_file]), sep=""))
  
}

##################################################
#Finally, let's do those that are not significant#
##################################################

dir.create("output/5_enrichment_analyses/1_prs/2_prs_bmi_ns")
output_path <- "output/5_enrichment_analyses/1_prs/2_prs_bmi_ns/"

for(index_file in seq(1, length(list_of_files))){
  
  #STEP 1.1: load the data:
  
  grs_df <- fread(paste("output/5_enrichment_analyses/1_prs/1_expo_outcome_df//", list_of_files[index_file], sep=""))
  grs_df <- grs_df[which(grs_df$SNP%in%bmi_no_sign$SNP),]
  
  #STEP 1.2: compute PRS. Weighted and unweighted:
  
  weighted_grs_res <- gtx::grs.summary(w=grs_df$beta.exposure, b=grs_df$beta.outcome, s = grs_df$se.outcome, grs_df$samplesize.outcome)
  unweighted_grs_res <- gtx::grs.summary(w=rep(1, length(grs_df$beta.exposure)), b=grs_df$beta.outcome, s = grs_df$se.outcome, grs_df$samplesize.outcome)
  
  saveRDS(weighted_grs_res, paste(output_path, "weighted_", gsub(".txt", ".RDS", list_of_files[index_file]), sep=""))
  saveRDS(unweighted_grs_res, paste(output_path, "unweighted_", gsub(".txt", ".RDS", list_of_files[index_file]), sep=""))
  
}
