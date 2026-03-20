##############
#INTRODUCTION#
##############

#This code plots the weighted PRS for 282 IR variants for Figure 2A of the manuscript.
#NOTE: I updated the code to also plots CHEERS results down below for Figure 3.

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
  
  grs_df$Trait <- factor(grs_df$outcome, levels = c("FIadjBMI", "HDL", "TG",  "ISIadjBMI", "IFCadjBMI", "TG/HDL", "2hGadjBMI", "FGadjBMI", "BMI", "HC", "HCadjBMI", "WC", "WCadjBMI", "WHR", "WHRadjBMI", "CHD", "CKD", "Hypertension", "NAFLD", "PCOS", "T2D", "ASAT", "ASATadjBMI", "GSAT", "GSATadjBMI", "VAT", "VATadjBMI", "Anterior thigh fat", "Posterior thigh fat", "Liver volume", "Liver fat %", "Pancreas fat %"))
  grs_df$type  <- factor(grs_df$type, levels = c("Hallmark IR", "Glycemic", "IR-derived disease", "Anthropometric", "Fat depot volumes (MRI)", "Ectopic fat deposition (MRI)"))
  
  # Ensure pvals numeric and format
  grs_df$pvals_num <- as.numeric(grs_df$pvals)
  grs_df$pval_formatted <- paste0("p=", format(grs_df$pvals_num, scientific = TRUE, digits = 2))
  
  # Create significance stars column based on p-value thresholds
  grs_df$signif_stars <- ""
  grs_df$signif_stars[grs_df$pvals_num < 0.05] <- "*"
  grs_df$signif_stars[grs_df$pvals_num < 5e-8] <- "**"
  
  # Position to place stars slightly above upper_ci
  star_y_pos <- as.numeric(grs_df$upper_ci) + ifelse(weighted == TRUE, 0.3, 0.01)  # tweak this offset if needed
  
  #Now we can make our plot!
  
  plotio <- ggplot(grs_df, aes(x = fct_rev(Trait), y = as.numeric(betas))) +
    geom_point(aes(color=type), size = 6, position = position_dodge(width = 0.75)) +
    geom_errorbar(aes(ymin = as.numeric(lower_ci), ymax = as.numeric(upper_ci)),
                  width = 0.15, position = position_dodge(width = 0.75),
                  color = "black", linewidth = 1) +
    geom_text(aes(y = star_y_pos, label = signif_stars, group = type),
              position = position_dodge(width = 0.75),
              size = 6, color = "black", fontface = "bold") +
    facet_wrap(~type, strip.position = "left", 
               nrow = 7, scales = "free_y") +
    geom_hline(yintercept = 0, color = "red", size = 0.5, linetype = 2) +
    
    coord_flip() +
    xlab("Trait") +
    ylab("effect sizes") +
    ggtitle(title) +
    theme_bw() +
    theme(
      legend.position = "none",
      legend.text = element_text(size = 11),
      panel.grid.major.x = element_line(color = "grey80", size = 0.1, linetype = 3),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_blank(),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 14),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
      
    )
  
  #Let's do one last switch to make all PRS comparable:
  
  if (weighted == TRUE) {
    plotio <- plotio + scale_y_continuous(limits = c(-4, 4),  breaks = seq(-4, 4, by = 1))
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
  
  grs_df$Trait <- factor(grs_df$outcome, levels = c("FIadjBMI", "HDL", "TG",  "ISIadjBMI", "IFCadjBMI", "TG/HDL", "2hGadjBMI", "FGadjBMI", "BMI", "HC", "HCadjBMI", "WC", "WCadjBMI", "WHR", "WHRadjBMI", "CHD", "CKD", "Hypertension", "NAFLD", "PCOS", "T2D", "ASAT", "ASATadjBMI", "GSAT", "GSATadjBMI", "VAT", "VATadjBMI", "Anterior thigh fat", "Posterior thigh fat", "Liver volume", "Liver fat %", "Pancreas fat %"))
  grs_df$type  <- factor(grs_df$type, levels = c("Hallmark IR", "Glycemic", "IR-derived disease", "Anthropometric", "Fat depot volumes (MRI)", "Ectopic fat deposition (MRI)"))
  
  # Ensure pvals numeric and format
  grs_df$pvals_num <- as.numeric(grs_df$pvals)
  grs_df$pval_formatted <- paste0("p=", format(grs_df$pvals_num, scientific = TRUE, digits = 2))
  
  # Create significance stars column based on p-value thresholds
  grs_df$signif_stars <- ""
  grs_df$signif_stars[grs_df$pvals_num < 0.05] <- "*"
  grs_df$signif_stars[grs_df$pvals_num < 5e-8] <- "**"
  
  # Position to place stars slightly above upper_ci
  star_y_pos <- as.numeric(grs_df$upper_ci) + ifelse(weighted == TRUE, 0.3, 0.01)  # tweak this offset if needed
  
  plotio <- ggplot(grs_df, aes(x = fct_rev(Trait), y = as.numeric(betas))) +
    geom_point(aes(color=type), size = 6, position = position_dodge(width = 0.75), show.legend = TRUE) +
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
      panel.grid.major.x = element_line(color = "grey80", size = 0.1, linetype = 3),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_blank(),
      
      # Keep axis text for "x" (i.e., traits), but strip it for "y"
      axis.text.x = element_text(size = 12),       # Trait labels (flipped to y)
      axis.title.x = element_text(size = 14),      # If you want a label for Trait
      
      axis.text.y = element_blank(),               # No effect size ticks
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      strip.text = element_blank(), #removes facet in this plot to avoid redundancy
      strip.background = element_rect(color = "black", fill = NA, linewidth = 0.5),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    )
  
  #Let's do one last switch to make all PRS comparable:
  
  if (weighted == TRUE) {
    plotio <- plotio + scale_y_continuous(limits = c(-4, 4),  breaks = seq(-4, 4, by = 1))
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
  
  grs_df$Trait <- factor(grs_df$outcome, levels = c("FIadjBMI", "HDL", "TG",  "ISIadjBMI", "IFCadjBMI", "TG/HDL", "2hGadjBMI", "FGadjBMI", "BMI", "HC", "HCadjBMI", "WC", "WCadjBMI", "WHR", "WHRadjBMI", "CHD", "CKD", "Hypertension", "NAFLD", "PCOS", "T2D", "ASAT", "ASATadjBMI", "GSAT", "GSATadjBMI", "VAT", "VATadjBMI", "Anterior thigh fat", "Posterior thigh fat", "Liver volume", "Liver fat %", "Pancreas fat %"))
  grs_df$type  <- factor(grs_df$type, levels = c("Hallmark IR", "Glycemic", "IR-derived disease", "Anthropometric", "Fat depot volumes (MRI)", "Ectopic fat deposition (MRI)"))
  
  # Ensure pvals numeric and format
  grs_df$pvals_num <- as.numeric(grs_df$pvals)
  grs_df$pval_formatted <- paste0("p=", format(grs_df$pvals_num, scientific = TRUE, digits = 2))
  
  # Create significance stars column based on p-value thresholds
  grs_df$signif_stars <- ""
  grs_df$signif_stars[grs_df$pvals_num < 0.05] <- "*"
  grs_df$signif_stars[grs_df$pvals_num < 5e-8] <- "**"
  
  # Position to place stars slightly above upper_ci
  star_y_pos <- as.numeric(grs_df$upper_ci) + ifelse(weighted == TRUE, 0.3, 0.01)  # tweak this offset if needed
  
  plotio <- ggplot(grs_df, aes(x = fct_rev(Trait), y = as.numeric(betas))) +
    geom_point(aes(color=type), size = 6, position = position_dodge(width = 0.75)) +
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
      legend.text = element_blank(),
      panel.grid.major.x = element_line(color = "grey80", size = 0.1, linetype = 3),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      axis.line = element_blank(),
      
      # Keep axis text for "x" (i.e., traits), but strip it for "y"
      axis.text.x = element_text(size = 12),       # Trait labels (flipped to y)
      axis.title.x = element_text(size = 14),      # If you want a label for Trait
      
      axis.text.y = element_blank(),               # No effect size ticks
      axis.title.y = element_blank(),
      axis.ticks.y = element_blank(),
      strip.text = element_blank(), #removes facet in this plot to avoid redundancy
      strip.background = element_rect(color = "black", fill = NA, linewidth = 0.5),
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5)
    )
  
  #Let's make sure we add the same axis:
  
  if (weighted == TRUE) {
    plotio <- plotio + scale_y_continuous(limits = c(-4, 4),  breaks = seq(-4, 4, by = 1))
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
bmi_ns_prs <- compiling_and_cleaning_scores(path_bmi_ns, get_unweighted=FALSE, analysis_ = "141 BMI-neutral")
bmi_neg_prs <- compiling_and_cleaning_scores(path_bmi_neg, get_unweighted=FALSE, analysis_ = "63 BMI-decreasing")
bmi_pos_prs <- compiling_and_cleaning_scores(path_bmi_pos, get_unweighted=FALSE, analysis_ = "78 BMI-increasing")
#bmi_pos_w_outlier_prs <- compiling_and_cleaning_scores(path_bmi_pos_out, get_unweighted=FALSE, analysis_ = "BMI-increasing (only outliers)")
#bmi_pos_wo_outlier_prs <- compiling_and_cleaning_scores(path_bmi_pos_wo_out, get_unweighted=FALSE, analysis_ = "63 BMI-increasing")

weighted_prs <- rbind(all_prs, bmi_ns_prs, bmi_neg_prs, bmi_pos_prs)

##########################################################################################
#Let's now plot the data for the weighted variants, we are including all traits this time#
##########################################################################################

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

#Get the data for the first panel

weighted_prs_panel_a <- weighted_prs[which(weighted_prs$type%in%c("Hallmark IR", "Glycemic", "IR-derived disease")),]
weighted_prs_panel_b <- weighted_prs[which(weighted_prs$type%in%c("Anthropometric", "Fat depot volumes (MRI)", "Ectopic fat deposition (MRI)")),]

#Let's save this so that we can fit it in the final table:

weigthed_all_a <- grs_plotter_first(weighted_prs_panel_a[which(weighted_prs_panel_a$analysis == "282 IR SNPs"),], "282 IR SNPs", weighted= TRUE)
weigthed_bmi_ns_a <- grs_plotter_mid(weighted_prs_panel_a[which(weighted_prs_panel_a$analysis == "141 BMI-neutral"),], "141 BMI-neutral IR SNPs", weighted= TRUE)
weigthed_bmi_neg_a <- grs_plotter_mid(weighted_prs_panel_a[which(weighted_prs_panel_a$analysis == "63 BMI-decreasing"),], "63 BMI-decreasing IR SNPs", weighted= TRUE)
weigthed_bmi_pos_a<- grs_plotter_last(weighted_prs_panel_a[which(weighted_prs_panel_a$analysis == "78 BMI-increasing"),], "78 BMI-increasing IR SNPs", weighted= TRUE)

weigthed_all_b <- grs_plotter_first(weighted_prs_panel_b[which(weighted_prs_panel_b$analysis == "282 IR SNPs"),], "282 IR SNPs", weighted= TRUE)
weigthed_bmi_ns_b <- grs_plotter_mid(weighted_prs_panel_b[which(weighted_prs_panel_b$analysis == "141 BMI-neutral"),], "141 BMI-neutral IR SNPs", weighted= TRUE)
weigthed_bmi_neg_b <- grs_plotter_mid(weighted_prs_panel_b[which(weighted_prs_panel_b$analysis == "63 BMI-decreasing"),], "63 BMI-decreasing IR SNPs", weighted= TRUE)
weigthed_bmi_pos_b <- grs_plotter_last(weighted_prs_panel_b[which(weighted_prs_panel_b$analysis == "78 BMI-increasing"),], "78 BMI-increasing IR SNPs", weighted= TRUE)

#####################
#Let's mix the data!#
#####################

# #Let's finally make the figure!
# 
# library(patchwork)
# 
# p_all_a <- weigthed_all_a + weigthed_bmi_ns_a + weigthed_bmi_neg_a + weigthed_bmi_pos_a + plot_layout(ncol = 4)
# p_all_b <- weigthed_all_b + weigthed_bmi_ns_b + weigthed_bmi_neg_b + weigthed_bmi_pos_b + plot_layout(ncol = 4)
# 
# ggsave(
#   filename = "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/manuscript/figures/drafts/prs_clusters_weighted_main_figure_panel_a.svg",
#   plot = p_all_a,
#   width = 14,
#   height = 9,
#   units = "in",  # or "in" depending on your preference, but cm is common
#   device = "svg"
# )
# 
# ggsave(
#   filename = "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/manuscript/figures/drafts/prs_clusters_weighted_main_figure_panel_b.svg",
#   plot = p_all_b,
#   width = 14,
#   height = 9,
#   units = "in",  # or "in" depending on your preference, but cm is common
#   device = "svg"
# )

###########################################################################################################
#Let's make a version of panel A that has only the 282 IR variants! Updated figure 1 after Tuomas comments#
###########################################################################################################

weighted_prs <- all_prs

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

#Get the data for the first panel

weighted_prs_panel_a <- weighted_prs[which(weighted_prs$type%in%c("Hallmark IR", "Glycemic", "IR-derived disease", "Anthropometric")),]

#Let's save this so that we can fit it in the final table:

weigthed_all_a <- grs_plotter_first(weighted_prs_panel_a[which(weighted_prs_panel_a$analysis == "282 IR SNPs"),], "282 IR SNPs", weighted= TRUE)

ggsave(
  filename = "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/manuscript/figures/drafts/prs_282_panel_A_figure_1_31102025.svg",
  plot = weigthed_all_a,
  width = 6,
  height = 9,
  units = "in",  # or "in" depending on your preference, but cm is common
  device = "svg"
)

#################################################################################################################################
#Let's do another special editing here! We are going to make a new figure out for MRI-derived PRS and CHEERS results per cluster#
#################################################################################################################################

weighted_prs <- rbind(bmi_ns_prs, bmi_neg_prs, bmi_pos_prs)

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

#Get the data for the first panel

weighted_prs_panel_a <- weighted_prs[which(weighted_prs$type%in%c("Fat depot volumes (MRI)", "Ectopic fat deposition (MRI)")),]

#Let's save this so that we can fit it in the final table:

weigthed_bmi_ns_a <- grs_plotter_first(weighted_prs_panel_a[which(weighted_prs_panel_a$analysis == "141 BMI-neutral"),], "141 BMI-neutral IR SNPs", weighted= TRUE)
weigthed_bmi_neg_a <- grs_plotter_mid(weighted_prs_panel_a[which(weighted_prs_panel_a$analysis == "63 BMI-decreasing"),], "63 BMI-decreasing IR SNPs", weighted= TRUE)
weigthed_bmi_pos_a<- grs_plotter_last(weighted_prs_panel_a[which(weighted_prs_panel_a$analysis == "78 BMI-increasing"),], "78 BMI-increasing IR SNPs", weighted= TRUE)

###########################################################
#Now let's add the CHEERS data so that we can combine it! #
###########################################################

Sys.setenv(LANG = "en")
library(tidyverse)
library(openxlsx)
library(data.table)
library(stringr)
library(GenomicRanges)

library(ggplot2)
library(dplyr)
library(forcats)
library(stringr)
library(scales)

#functions

cleaning_samples <- function(batch_){
  
  #STEP 1: remove combat_corrected_
  
  batch_clean <- unlist(str_split(batch_, "_combat_corrected"))[1]
  
  batch_data <- unlist(str_split(batch_clean, "_"))
  
  cell <- batch_data[1]
  day <- batch_data[2]
  batch <- batch_data[3]
  rep <- batch_data[4]
  
  #Let's clean the day:
  
  day <- unlist(str_split(day, "day"))[2]
  day <- paste("Day", day)
  
  #Same with batch:
  
  batch <- unlist(str_split(batch, "batch"))[2]
  batch <- paste("Batch", batch)
  
  #Now with replicates:
  
  rep <- unlist(str_split(rep, "rep"))[2]
  rep <- paste("Replicate", rep) 
  
  #Final string:
  
  batch_end <- paste(day, " - ", batch, " (", rep, ")", sep = "")
  #batch_end <- paste(batch, " (", rep, ")", sep = "")
  
  return(batch_end)
  
}

day_parser <- function(batch_){
  
  #STEP 1: remove combat_corrected_
  
  batch_clean <- unlist(str_split(batch_, "_combat_corrected"))[1]
  
  batch_data <- unlist(str_split(batch_clean, "_"))
  
  day <- batch_data[2]
  
  return(day)
  
}

# Define a consistent color palette for all plots
cell_stage_palette <- c(
  "Day 0"  = "#E69F00",  # orange
  "Day 4"  = "#56B4E9",  # sky blue
  "Day 14"  = "#009E73"  # bluish green
)

# Core plot function with customizable theme tweaks
cheers_plot_base <- function(cheers_enrichment, title, strip_y = FALSE, show_legend = FALSE) {
  
  colnames(cheers_enrichment) <- c("Samples", "Pval", "Day")
  
  cheers_enrichment <- cheers_enrichment %>%
    mutate(Day = case_when(
      str_detect(Day, "day0") ~ "Day 0",
      str_detect(Day, "day2") ~ "Day 2",
      str_detect(Day, "day4") ~ "Day 4",
      TRUE ~ "Day 14"
    )) %>%
    mutate(Day = factor(Day, levels = c("Day 0", "Day 4", "Day 14")),
           Samples = factor(Samples, levels = unique(Samples[order(Day)])))
  
  p <- ggplot(cheers_enrichment, aes(x = fct_rev(Samples), y = -log10(Pval), fill = Day)) +
    geom_col(width = 0.8) +
    coord_flip() +
    geom_hline(yintercept = -log10(0.05 / 25), linetype = "dashed", color = "firebrick", size = 1) +
    scale_fill_manual(values = cell_stage_palette) +
    labs(
      title = title,
      y = expression(-log[10](italic("P"))),
      x = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      axis.text.y = if (strip_y) element_blank() else element_text(size = 12),
      axis.text.x = element_text(size = 12),
      axis.ticks.y = if (strip_y) element_blank() else element_line(),
      axis.title.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = if (show_legend) "right" else "none",
      legend.title = element_text(size = 11),
      legend.text = element_text(size = 10),
      axis.line = element_line(colour = "black"),
      plot.background = element_rect(fill = "white", color = NA)
    )
  
  return(p)
}

# Specific versions using the base with tweaks

cheers_plotting_unique <- function(cheers_enrichment, title) {
  cheers_plot_base(cheers_enrichment, title, strip_y = FALSE, show_legend = TRUE)
}


cheers_plotting_first <- function(cheers_enrichment, title) {
  cheers_plot_base(cheers_enrichment, title, strip_y = FALSE, show_legend = FALSE)
}

cheers_plotting_mid <- function(cheers_enrichment, title) {
  cheers_plot_base(cheers_enrichment, title, strip_y = TRUE, show_legend = FALSE)
}

cheers_plotting_last <- function(cheers_enrichment, title) {
  cheers_plot_base(cheers_enrichment, title, strip_y = TRUE, show_legend = TRUE)
}


#####################
#Let's load the data#
#####################

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/")

cheers_enrichment_bmi_ns = read.table("output/5_enrichment_analyses/5_cheers/2_cheers_output/proxies_bmi_ns_disease_enrichment_pValues.txt")
cheers_enrichment_bmi_neg = read.table("output/5_enrichment_analyses/5_cheers/2_cheers_output/proxies_bmi_neg_disease_enrichment_pValues.txt")
cheers_enrichment_bmi_pos = read.table("output/5_enrichment_analyses/5_cheers/2_cheers_output/proxies_bmi_pos_disease_enrichment_pValues.txt")

#Let's clean a bit the data:

cheers_enrichment_bmi_ns$day <- sapply(cheers_enrichment_bmi_ns$V1, day_parser)
cheers_enrichment_bmi_ns$V1 <- sapply(cheers_enrichment_bmi_ns$V1, cleaning_samples) #the function will take care of the naming conventions

cheers_enrichment_bmi_neg$day <- sapply(cheers_enrichment_bmi_neg$V1, day_parser)
cheers_enrichment_bmi_neg$V1 <- sapply(cheers_enrichment_bmi_neg$V1, cleaning_samples) #the function will take care of the naming conventions

cheers_enrichment_bmi_pos$day <- sapply(cheers_enrichment_bmi_pos$V1, day_parser)
cheers_enrichment_bmi_pos$V1 <- sapply(cheers_enrichment_bmi_pos$V1, cleaning_samples) #the function will take care of the naming conventions

bmi_ns_plot <- cheers_plotting_first(cheers_enrichment_bmi_ns, title="141 BMI-neutral IR loci")
bmi_neg_plot <- cheers_plotting_mid(cheers_enrichment_bmi_neg, title="63 BMI-decreasing IR loci")
bmi_pos_plot <- cheers_plotting_last(cheers_enrichment_bmi_pos, title="78 BMI+ IR loci")

#Let's combine the plots:

library(patchwork)

p_combined=weigthed_bmi_ns_a + weigthed_bmi_neg_a + weigthed_bmi_pos_a + bmi_ns_plot + bmi_neg_plot + bmi_pos_plot + plot_layout(ncol=3)

ggsave(
  filename = "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/manuscript/figures/drafts/prs_cheers_for_clusters_figure_2_31102025.svg",
  plot = p_combined,
  width = 12,
  height = 12,
  units = "in",  # or "in" depending on your preference, but cm is common
  device = "svg"
)
