##############
#INTRODUCTION#
##############

#This code matches our SNPs data with the traits of interest
#harmonizes them so that they are properly aligned.
#And finally performs GRS, which are saved ina dataframe.

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)

###################
#Loading functions#
###################

source("N:/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/Mediation_analysis/code/0_functions/functions_4_mediation.R")

load_grs_data <- function(list_of_files){
  
  for(index_file in seq(1, length(list_of_files))){
    
    tmp_df <- readRDS(list_of_files[index_file])
    
    if(!(exists("grs_df"))){
      
      #Let's make columns...
      
      exposure <- "IR"
      outcome <- str_split(list_of_files[index_file], ".RDS")[[1]][1]
      outcome <- str_split(outcome, "_")[[1]]
      outcome <- outcome[length(outcome)]
      
      betas <- tmp_df$ahat
      ses <- tmp_df$aSE
      pvals <- tmp_df$pval
      lowerci <- lower_ci(betas, ses)
      upperci <- upper_ci(betas, ses)
      nsnps <- tmp_df$m
      
      final_grs <- c(exposure, outcome, betas, ses, pvals, lowerci, upperci, nsnps)
      final_grs <- as.data.frame(t(final_grs))
      colnames(final_grs) <- c("exposure", "outcome", "betas", "ses", "pvals", "lower_ci", "upper_ci", "nsnps")
      
      grs_df <- final_grs
      
    } else {
      
      exposure <- "IR"
      outcome <- str_split(list_of_files[index_file], ".RDS")[[1]][1]
      outcome <- str_split(outcome, "_")[[1]]
      outcome <- outcome[length(outcome)]
      
      betas <- tmp_df$ahat
      ses <- tmp_df$aSE
      pvals <- tmp_df$pval
      lowerci <- lower_ci(betas, ses)
      upperci <- upper_ci(betas, ses)
      nsnps <- tmp_df$m
      
      final_grs <- c(exposure, outcome, betas, ses, pvals, lowerci, upperci, nsnps)
      final_grs <- as.data.frame(t(final_grs))
      colnames(final_grs) <- c("exposure", "outcome", "betas", "ses", "pvals", "lower_ci", "upper_ci", "nsnps")
      
      grs_df <- rbind(grs_df, final_grs)
      
    }
    
  }
  
  return(grs_df)
  
}

grs_plotter_first <- function(grs_df, title, weighted){
  
  # Format betas and CIs
  betas_rounded     <- format(as.numeric(grs_df$betas), scientific = TRUE, digits = 4)
  lower_ci_rounded  <- format(as.numeric(grs_df$lower_ci), scientific = TRUE, digits = 4)
  upper_ci_rounded  <- format(as.numeric(grs_df$upper_ci), scientific = TRUE, digits = 4)
  
  grs_df$full_effect <- paste(betas_rounded, " (", lower_ci_rounded, ",", upper_ci_rounded, ")", sep = "")
  
  grs_df$Trait <- factor(grs_df$outcome, levels = c("FIadjBMI", "HDL", "TG",  "TG/HDL", "ISIadjBMI", "IFCadjBMI", "FGadjBMI", "2hGadjBMI", "BMI", "HC", "HCadjBMI", "WC", "WCadjBMI", "WHR", "WHRadjBMI", "CHD", "CKD", "Hypertension", "NAFLD", "PCOS", "T2D", "ASAT", "ASATadjBMI", "GSAT", "GSATadjBMI", "VAT", "VATadjBMI", "Anterior thigh fat", "Posterior thigh fat", "Liver volume", "Liver fat %", "Pancreas fat %"))
  grs_df$type  <- factor(grs_df$type, levels = c("Hallmark IR", "Glycemic", "Anthropometric", "IR-derived disease", "Fat deposition (MRI)"))
  
  # Ensure pvals numeric and format
  grs_df$pvals_num <- as.numeric(grs_df$pvals)
  grs_df$pval_formatted <- paste0("p=", format(grs_df$pvals_num, scientific = TRUE, digits = 2))
  
  # Create significance stars column based on p-value thresholds
  grs_df$signif_stars <- ""
  grs_df$signif_stars[grs_df$pvals_num < 0.05] <- "*"
  grs_df$signif_stars[grs_df$pvals_num < 5e-8] <- "**"
  
  # Position to place stars slightly above upper_ci
  star_y_pos <- as.numeric(grs_df$upper_ci) + ifelse(weighted == TRUE, 0.1, 0.001)  # tweak this offset if needed
  
  #Now we can make our plot!
  
  plotio <- ggplot(grs_df, aes(x = fct_rev(Trait), y = as.numeric(betas))) +
    geom_point(aes(color=type), size = 5, position = position_dodge(width = 0.75)) +
    geom_errorbar(aes(ymin = as.numeric(lower_ci), ymax = as.numeric(upper_ci)),
                  width = 0.15, position = position_dodge(width = 0.75),
                  color = "black", linewidth = 1) +
    geom_text(aes(y = star_y_pos, label = signif_stars, group = type),
              position = position_dodge(width = 0.75),
              size = 6, color = "black", fontface = "bold") +
    facet_wrap(~type, strip.position = "left", nrow = 7, scales = "free_y") +
    geom_hline(yintercept = 0, color = "red", size = 0.5, linetype = 2) +
    
    coord_flip() +
    xlab("Trait") +
    ylab("effect sizes") +
    ggtitle(title) +
    theme_bw() +
    theme(
      legend.position = "none",
      legend.text = element_text(size = 11),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black"),
      axis.text = element_text(size = 11),
      axis.title = element_text(size = 10, face = "bold"),
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
    )
  
  #Let's do one last switch to make all PRS comparable:
  
  if (weighted == TRUE) {
    plotio <- plotio + scale_y_continuous(limits = c(-6, 6))
  } else {
    
    plotio <- plotio + scale_y_continuous(limits = c(-0.1, 0.1))
    
  }
  
  return(plotio)
}

grs_plotter_mid <- function(grs_df, title, weighted){
  
  # Format betas and CIs
  betas_rounded     <- format(as.numeric(grs_df$betas), scientific = TRUE, digits = 4)
  lower_ci_rounded  <- format(as.numeric(grs_df$lower_ci), scientific = TRUE, digits = 4)
  upper_ci_rounded  <- format(as.numeric(grs_df$upper_ci), scientific = TRUE, digits = 4)
  
  grs_df$full_effect <- paste(betas_rounded, " (", lower_ci_rounded, ",", upper_ci_rounded, ")", sep = "")
  
  grs_df$Trait <- factor(grs_df$outcome, levels = c("FIadjBMI", "HDL", "TG",  "TG/HDL", "ISIadjBMI", "IFCadjBMI", "FGadjBMI", "2hGadjBMI", "BMI", "HC", "HCadjBMI", "WC", "WCadjBMI", "WHR", "WHRadjBMI", "CHD", "CKD", "Hypertension", "NAFLD", "PCOS", "T2D", "ASAT", "ASATadjBMI", "GSAT", "GSATadjBMI", "VAT", "VATadjBMI", "Anterior thigh fat", "Posterior thigh fat", "Liver volume", "Liver fat %", "Pancreas fat %"))
  grs_df$type  <- factor(grs_df$type, levels = c("Hallmark IR", "Glycemic", "Anthropometric", "IR-derived disease", "Fat deposition (MRI)"))
  
  # Ensure pvals numeric and format
  grs_df$pvals_num <- as.numeric(grs_df$pvals)
  grs_df$pval_formatted <- paste0("p=", format(grs_df$pvals_num, scientific = TRUE, digits = 2))
  
  # Create significance stars column based on p-value thresholds
  grs_df$signif_stars <- ""
  grs_df$signif_stars[grs_df$pvals_num < 0.05] <- "*"
  grs_df$signif_stars[grs_df$pvals_num < 5e-8] <- "**"
  
  # Position to place stars slightly above upper_ci
  star_y_pos <- as.numeric(grs_df$upper_ci) + ifelse(weighted == TRUE, 0.1, 0.001)  # tweak this offset if needed
  
  plotio <- ggplot(grs_df, aes(x = fct_rev(Trait), y = as.numeric(betas))) +
    geom_point(aes(color=type), size = 5, position = position_dodge(width = 0.75)) +
    geom_errorbar(aes(ymin = as.numeric(lower_ci), ymax = as.numeric(upper_ci)),
                  width = 0.15, position = position_dodge(width = 0.75),
                  color = "black", linewidth = 1) +
    geom_text(aes(y = star_y_pos, label = signif_stars, group = type),
              position = position_dodge(width = 0.75),
              size = 6, color = "black", fontface = "bold") +
    facet_wrap(~type, strip.position = "left", nrow = 7, scales = "free_y") +
    geom_hline(yintercept = 0, color = "red", size = 0.5, linetype = 2) +
    
    coord_flip() +
    xlab("Trait") +
    ylab("effect sizes") +
    ggtitle(title) +
    theme_bw() +
    theme(
      legend.position = "none",
      legend.text = element_text(size = 11),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_blank(),
      
      # Keep axis text for "x" (i.e., traits), but strip it for "y"
      axis.text.x = element_text(size = 11),       # Trait labels (flipped to y)
      axis.title.x = element_text(size = 10),      # If you want a label for Trait
      
      axis.text.y = element_blank(),               # No effect size ticks
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
    )
  
  #Let's do one last switch to make all PRS comparable:
  
  if (weighted == TRUE) {
    plotio <- plotio + scale_y_continuous(limits = c(-6, 6))
  } else {
    
    plotio <- plotio + scale_y_continuous(limits = c(-0.1, 0.1))
    
  }
  
  return(plotio)
}

grs_plotter_last <- function(grs_df, title, weighted){
  
  # Format betas and CIs
  betas_rounded     <- format(as.numeric(grs_df$betas), scientific = TRUE, digits = 4)
  lower_ci_rounded  <- format(as.numeric(grs_df$lower_ci), scientific = TRUE, digits = 4)
  upper_ci_rounded  <- format(as.numeric(grs_df$upper_ci), scientific = TRUE, digits = 4)
  
  grs_df$full_effect <- paste(betas_rounded, " (", lower_ci_rounded, ",", upper_ci_rounded, ")", sep = "")
  
  grs_df$Trait <- factor(grs_df$outcome, levels = c("FIadjBMI", "HDL", "TG",  "TG/HDL", "ISIadjBMI", "IFCadjBMI", "FGadjBMI", "2hGadjBMI", "BMI", "HC", "HCadjBMI", "WC", "WCadjBMI", "WHR", "WHRadjBMI", "CHD", "CKD", "Hypertension", "NAFLD", "PCOS", "T2D", "ASAT", "ASATadjBMI", "GSAT", "GSATadjBMI", "VAT", "VATadjBMI", "Anterior thigh fat", "Posterior thigh fat", "Liver volume", "Liver fat %", "Pancreas fat %"))
  grs_df$type  <- factor(grs_df$type, levels = c("Hallmark IR", "Glycemic", "Anthropometric", "IR-derived disease", "Fat deposition (MRI)"))
  
  # Ensure pvals numeric and format
  grs_df$pvals_num <- as.numeric(grs_df$pvals)
  grs_df$pval_formatted <- paste0("p=", format(grs_df$pvals_num, scientific = TRUE, digits = 2))
  
  # Create significance stars column based on p-value thresholds
  grs_df$signif_stars <- ""
  grs_df$signif_stars[grs_df$pvals_num < 0.05] <- "*"
  grs_df$signif_stars[grs_df$pvals_num < 5e-8] <- "**"
  
  # Position to place stars slightly above upper_ci
  star_y_pos <- as.numeric(grs_df$upper_ci) + ifelse(weighted == TRUE, 0.1, 0.001)  # tweak this offset if needed
  
  plotio <- ggplot(grs_df, aes(x = fct_rev(Trait), y = as.numeric(betas))) +
    geom_point(aes(color=type), size = 5, position = position_dodge(width = 0.75), show.legend = TRUE) +
    geom_errorbar(aes(ymin = as.numeric(lower_ci), ymax = as.numeric(upper_ci)),
                  width = 0.15, position = position_dodge(width = 0.75),
                  color = "black", linewidth = 1) +
    geom_text(aes(y = star_y_pos, label = signif_stars, group = type),
              position = position_dodge(width = 0.75),
              size = 6, color = "black", fontface = "bold") +
    facet_wrap(~type, strip.position = "left", nrow = 7, scales = "free_y") +
    geom_hline(yintercept = 0, color = "red", size = 0.5, linetype = 2) +
    
    coord_flip() +
    xlab("Trait") +
    ylab("effect sizes") +
    ggtitle(title) +
    theme_bw() +
    theme(
      legend.position = "right",
      legend.text = element_text(size = 11),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_blank(),
      
      # Keep axis text for "x" (i.e., traits), but strip it for "y"
      axis.text.x = element_text(size = 11),       # Trait labels (flipped to y)
      axis.title.x = element_text(size = 10),      # If you want a label for Trait
      
      axis.text.y = element_blank(),               # No effect size ticks
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
    )
  
  #Let's make sure we add the same axis:
  
  if (weighted == TRUE) {
    plotio <- plotio + scale_y_continuous(limits = c(-6, 6))
  } else {
    
    plotio <- plotio + scale_y_continuous(limits = c(-0.1, 0.1))
    
  }
  
  return(plotio)
}

compiling_and_cleaning_scores<- function(path_, get_unweighted=TRUE, analysis_){
  #An easier compiler that uses load_grs to get the data we want
  
  #STEP 0: setting up the working directory:
  
  #setwd(path_)
  
  #STEP 1: get the files
  
  list_of_files <- list.files(path_, full.names = TRUE)
  
  #STEP 2: get either weighted or unweighted and obtain only data for the ir_variants
  
  list_of_files <- list_of_files[which(str_detect(list_of_files, "unweighted") == get_unweighted)]
  list_of_files <- list_of_files[which(str_detect(list_of_files, "ir_variants"))]
  list_of_files <- list_of_files[which(str_detect(list_of_files, "non") == FALSE)]
  
  #STEP 3: load analyses for each trait:
  
  ir_282_variants <- load_grs_data(list_of_files)
  
  #STEP 4: let's clean the weird results:
  
  ir_282_variants$outcome <- ifelse(ir_282_variants$outcome == "ratio", "TG/HDL", ir_282_variants$outcome)
  ir_282_variants$outcome <- ifelse(ir_282_variants$outcome == "Anterior thight fat", "Anterior thigh fat", ir_282_variants$outcome)
  ir_282_variants$outcome <- ifelse(ir_282_variants$outcome == "Poster thigh fat", "Posterior thigh fat", ir_282_variants$outcome)
  
  #STEP 5: remove duplicates:Careful with duplicates due to different runs:
  
  ir_282_variants <- ir_282_variants[which(duplicated(ir_282_variants) == FALSE),]
  
  #STEP 6: return data with analyses data
  
  ir_282_variants$analysis <- analysis_
  
  return(ir_282_variants)
  
}

#######################
#Loading weighted data#
#######################

path_all <- "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/output/5_enrichment_analyses/1_prs/2_prs/"
path_bmi_ns <- "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/output/5_enrichment_analyses/1_prs/2_prs_bmi_ns/"
path_bmi_neg <- "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/output/5_enrichment_analyses/1_prs/2_prs_bmi_neg/"
path_bmi_pos <- "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/output/5_enrichment_analyses/1_prs/2_prs_bmi_pos/"
#path_bmi_pos_out <- "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/output/5_enrichment_analyses/1_prs/2_prs_bmi_pos_outlier/"
#path_bmi_pos_wo_out <- "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/output/5_enrichment_analyses/1_prs/2_prs_bmi_pos_wo_outlier/"

#Now let' load the data for the weighted ones:

all_prs <- compiling_and_cleaning_scores(path_all, get_unweighted=FALSE, analysis_ = "282 IR SNPs")
bmi_ns_prs <- compiling_and_cleaning_scores(path_bmi_ns, get_unweighted=FALSE, analysis_ = "BMI-neutral IR SNPs")
bmi_neg_prs <- compiling_and_cleaning_scores(path_bmi_neg, get_unweighted=FALSE, analysis_ = "BMI-decreasing IR SNPs")
bmi_pos_prs <- compiling_and_cleaning_scores(path_bmi_pos, get_unweighted=FALSE, analysis_ = "BMI-increasing IR SNPs")
#bmi_pos_w_outlier_prs <- compiling_and_cleaning_scores(path_bmi_pos_out, get_unweighted=FALSE, analysis_ = "BMI-increasing (only outliers)")
#bmi_pos_wo_outlier_prs <- compiling_and_cleaning_scores(path_bmi_pos_wo_out, get_unweighted=FALSE, analysis_ = "BMI-increasing (without outliers)")

weighted_prs <- rbind(all_prs, bmi_ns_prs, bmi_neg_prs, bmi_pos_prs)

weighted_prs$beta <- round(as.numeric(weighted_prs$betas), 2)
weighted_prs$confidence_intervals <- paste("[", as.character(round(as.numeric(weighted_prs$lower_ci), 2)), ",", as.character(round(as.numeric(weighted_prs$upper_ci), 2)), "]", sep = "")
weighted_prs$ses <- round(as.numeric(weighted_prs$ses), 4)
#weighted_prs$pvals <- format(as.numeric(weighted_prs$pvals), scientific = TRUE, digits = 3)

#Let's add a type so that we can order the traits:

weighted_prs$type = NA

IR <- c("FIadjBMI", "HDL", "TG")
glycemic <- c("TG/HDL", "ISIadjBMI", "IFCadjBMI", "FGadjBMI", "2hGadjBMI")
anthropometric <- c("BMI", "WHR", "WHRadjBMI", "WC", "WCadjBMI", "HC", "HCadjBMI")
ir_disease <- c("T2D", "CHD", "CKD", "Hypertension", "NAFLD", "PCOS")
mri_fat <- c("ASAT", "ASATadjBMI", "GSAT", "GSATadjBMI", "VAT", "VATadjBMI")
ectopic_fat <- c("Anterior thigh fat", "Liver fat %", "Liver volume", "Pancreas fat %", "Posterior thigh fat")

weighted_prs <- weighted_prs[which(weighted_prs$outcome%in%c(IR, glycemic, anthropometric, ir_disease, mri_fat, ectopic_fat)),]

weighted_prs$type <- ifelse(weighted_prs$outcome%in%IR, "Hallmark IR", weighted_prs$type)
weighted_prs$type <- ifelse(weighted_prs$outcome%in%glycemic, "Glycemic", weighted_prs$type)
weighted_prs$type <- ifelse(weighted_prs$outcome%in%anthropometric, "Anthropometric", weighted_prs$type)
weighted_prs$type <- ifelse(weighted_prs$outcome%in%ir_disease, "IR-derived disease", weighted_prs$type)
weighted_prs$type <- ifelse(weighted_prs$outcome%in%mri_fat, "Fat depot volumes (MRI)", weighted_prs$type)
weighted_prs$type <- ifelse(weighted_prs$outcome%in%ectopic_fat, "Ectopic fat deposition (MRI)", weighted_prs$type)

#Let's filter the data so that we can merge it easier in the supplementary tables:

weighted_prs_filter = weighted_prs %>%
  dplyr::select(type, analysis, outcome, nsnps, beta, confidence_intervals, ses, pvals)

#Let's try to order the data:

library(dplyr)
library(forcats)

# Analysis order (based on nsnps or name)
analysis_order <- c(
  "282 IR SNPs",
  "BMI-neutral IR SNPs",
  "BMI-decreasing IR SNPs",
  "BMI-increasing IR SNPs"
)

outcome_order <- c(
  "Hallmark IR",
  "Glycemic",
  "IR-derived disease",
  "Anthropometric",
  "Fat depot volumes (MRI)",
  "Ectopic fat deposition (MRI)"
)

weighted_prs_filter <- weighted_prs_filter %>%
  mutate(
    # Apply custom factor levels
    analysis = factor(analysis, levels = analysis_order),
    type  = factor(type,  levels = outcome_order)
  ) %>%
  arrange(analysis, type, outcome)


fwrite(weighted_prs_filter, "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/output/5_enrichment_analyses/1_prs/3_clean_results/all_clusters_weighted_prs.csv")

#########################
#Loading unweighted data#
#########################

#Now let' load the data for the weighted ones:

all_prs <- compiling_and_cleaning_scores(path_all, get_unweighted=TRUE, analysis_ = "282 IR SNPs")
bmi_ns_prs <- compiling_and_cleaning_scores(path_bmi_ns, get_unweighted=TRUE, analysis_ = "BMI-neutral IR SNPs")
bmi_neg_prs <- compiling_and_cleaning_scores(path_bmi_neg, get_unweighted=TRUE, analysis_ = "BMI-decreasing IR SNPs")
bmi_pos_prs <- compiling_and_cleaning_scores(path_bmi_pos, get_unweighted=TRUE, analysis_ = "BMI-increasing IR SNPs")
#bmi_pos_w_outlier_prs <- compiling_and_cleaning_scores(path_bmi_pos_out, get_unweighted=FALSE, analysis_ = "BMI-increasing (only outliers)")
#bmi_pos_wo_outlier_prs <- compiling_and_cleaning_scores(path_bmi_pos_wo_out, get_unweighted=FALSE, analysis_ = "BMI-increasing (without outliers)")

unweighted_prs <- rbind(all_prs, bmi_ns_prs, bmi_neg_prs, bmi_pos_prs)

#unweighted_prs$beta <- round(as.numeric(unweighted_prs$betas), 4)
unweighted_prs$confidence_intervals <- paste("[", format(as.numeric(unweighted_prs$lower_ci), scientific=TRUE, digits=3),
                                             ",", 
                                             format(as.numeric(unweighted_prs$upper_ci), scientific=TRUE, digits=3), "]", sep = "")
#unweighted_prs$ses <- round(as.numeric(unweighted_prs$ses)
#weighted_prs$pvals <- format(as.numeric(weighted_prs$pvals), scientific = TRUE, digits = 3)

#Let's add a type so that we can order the traits:

unweighted_prs$type = NA

IR <- c("FIadjBMI", "HDL", "TG")
glycemic <- c("TG/HDL", "ISIadjBMI", "IFCadjBMI", "FGadjBMI", "2hGadjBMI")
anthropometric <- c("BMI", "WHR", "WHRadjBMI", "WC", "WCadjBMI", "HC", "HCadjBMI")
ir_disease <- c("T2D", "CHD", "CKD", "Hypertension", "NAFLD", "PCOS")
mri_fat <- c("ASAT", "ASATadjBMI", "GSAT", "GSATadjBMI", "VAT", "VATadjBMI")
ectopic_fat <- c("Anterior thigh fat", "Liver fat %", "Liver volume", "Pancreas fat %", "Posterior thigh fat")

unweighted_prs <- unweighted_prs[which(unweighted_prs$outcome%in%c(IR, glycemic, anthropometric, ir_disease, mri_fat, ectopic_fat)),]

unweighted_prs$type <- ifelse(unweighted_prs$outcome%in%IR, "Hallmark IR", unweighted_prs$type)
unweighted_prs$type <- ifelse(unweighted_prs$outcome%in%glycemic, "Glycemic", unweighted_prs$type)
unweighted_prs$type <- ifelse(unweighted_prs$outcome%in%anthropometric, "Anthropometric", unweighted_prs$type)
unweighted_prs$type <- ifelse(unweighted_prs$outcome%in%ir_disease, "IR-derived disease", unweighted_prs$type)
unweighted_prs$type <- ifelse(unweighted_prs$outcome%in%mri_fat, "Fat depot volumes (MRI)", unweighted_prs$type)
unweighted_prs$type <- ifelse(unweighted_prs$outcome%in%ectopic_fat, "Ectopic fat deposition (MRI)", unweighted_prs$type)

#Let's filter the data so that we can merge it easier in the supplementary tables:

unweighted_prs_filter = unweighted_prs %>%
  dplyr::select(type, analysis, outcome, nsnps, betas, confidence_intervals, ses, pvals)

#Let's try to order the data:

library(dplyr)
library(forcats)

# Analysis order (based on nsnps or name)
analysis_order <- c(
  "282 IR SNPs",
  "BMI-neutral IR SNPs",
  "BMI-decreasing IR SNPs",
  "BMI-increasing IR SNPs"
)

outcome_order <- c(
  "Hallmark IR",
  "Glycemic",
  "IR-derived disease",
  "Anthropometric",
  "Fat depot volumes (MRI)",
  "Ectopic fat deposition (MRI)"
)

unweighted_prs_filter <- unweighted_prs_filter %>%
  mutate(
    # Apply custom factor levels
    analysis = factor(analysis, levels = analysis_order),
    type  = factor(type,  levels = outcome_order)
  ) %>%
  arrange(analysis, type, outcome)

fwrite(unweighted_prs_filter, "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/output/5_enrichment_analyses/1_prs/3_clean_results/all_clusters_unweighted_prs.csv")

