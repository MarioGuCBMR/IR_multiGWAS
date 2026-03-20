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
  
  grs_df$Trait <- factor(grs_df$outcome, levels = c("FIadjBMI", "HDL", "TG",  "TG/HDL", "ISIadjBMI", "IFCadjBMI", "FGadjBMI", "2hGadjBMI", "BMI", "HC", "HCadjBMI", "WC", "WCadjBMI", "WHR", "WHRadjBMI", "CHD", "CKD", "Hypertension", "NAFLD", "PCOS", "T2D", "ASAT", "ASATadjBMI", "GSAT", "GSATadjBMI", "VAT", "VATadjBMI", "Anterior thigh fat", "Posterior thigh fat", "Liver fat %", "Pancreas fat %"))
  grs_df$type  <- factor(grs_df$type, levels = c("Hallmark IR", "Glycemic", "Anthropometric", "IR-derived disease"))
  
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
      axis.line = element_line(colour = "black"),
      axis.text = element_text(size = 11),
      axis.title = element_text(size = 10, face = "bold"),
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
    )
  
  return(plotio)
}

grs_plotter_mid <- function(grs_df, title, weighted){
  
  # Format betas and CIs
  betas_rounded     <- format(as.numeric(grs_df$betas), scientific = TRUE, digits = 4)
  lower_ci_rounded  <- format(as.numeric(grs_df$lower_ci), scientific = TRUE, digits = 4)
  upper_ci_rounded  <- format(as.numeric(grs_df$upper_ci), scientific = TRUE, digits = 4)
  
  grs_df$full_effect <- paste(betas_rounded, " (", lower_ci_rounded, ",", upper_ci_rounded, ")", sep = "")
  
  grs_df$Trait <- factor(grs_df$outcome, levels = c("FIadjBMI", "HDL", "TG",  "TG/HDL", "ISIadjBMI", "IFCadjBMI", "FGadjBMI", "2hGadjBMI", "BMI", "HC", "HCadjBMI", "WC", "WCadjBMI", "WHR", "WHRadjBMI", "CHD", "CKD", "Hypertension", "NAFLD", "PCOS", "T2D", "ASAT", "ASATadjBMI", "GSAT", "GSATadjBMI", "VAT", "VATadjBMI", "Anterior thigh fat", "Posterior thigh fat", "Liver fat %", "Pancreas fat %"))
  grs_df$type  <- factor(grs_df$type, levels = c("Hallmark IR", "Glycemic", "Anthropometric", "IR-derived disease"))
  
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
  
  return(plotio)
}

grs_plotter_last <- function(grs_df, title, weighted){
  
  # Format betas and CIs
  betas_rounded     <- format(as.numeric(grs_df$betas), scientific = TRUE, digits = 4)
  lower_ci_rounded  <- format(as.numeric(grs_df$lower_ci), scientific = TRUE, digits = 4)
  upper_ci_rounded  <- format(as.numeric(grs_df$upper_ci), scientific = TRUE, digits = 4)
  
  grs_df$full_effect <- paste(betas_rounded, " (", lower_ci_rounded, ",", upper_ci_rounded, ")", sep = "")
  
  grs_df$Trait <- factor(grs_df$outcome, levels = c("FIadjBMI", "HDL", "TG",  "TG/HDL", "ISIadjBMI", "IFCadjBMI", "FGadjBMI", "2hGadjBMI", "BMI", "HC", "HCadjBMI", "WC", "WCadjBMI", "WHR", "WHRadjBMI", "CHD", "CKD", "Hypertension", "NAFLD", "PCOS", "T2D", "ASAT", "ASATadjBMI", "GSAT", "GSATadjBMI", "VAT", "VATadjBMI", "Anterior thigh fat", "Posterior thigh fat", "Liver fat %", "Pancreas fat %"))
  grs_df$type  <- factor(grs_df$type, levels = c("Hallmark IR", "Glycemic", "Anthropometric", "IR-derived disease"))
  
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
  
  return(plotio)
}

##############
#Loading data#
##############

path_4_grs <- "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/output/5_enrichment_analyses/1_prs/2_prs/"

setwd(path_4_grs)

#######################################################################
#STEP 1: Let's generate the plot for all the variants for the Weighted#
#######################################################################

list_of_files <- list.files(path_4_grs)

#Let's do the unweighted first:

list_of_files <- list_of_files[which(str_detect(list_of_files, "unweighted") == FALSE)]

ir_282_variants <- list_of_files[which(str_detect(list_of_files, "ir_variants"))]
ir_282_variants <- ir_282_variants[which(str_detect(ir_282_variants, "non") == FALSE)]

#Let's get the files for each trait:

ir_282_variants <- load_grs_data(ir_282_variants)

#Let's clean the weird one:

ir_282_variants$outcome <- ifelse(ir_282_variants$outcome == "ratio", "TG/HDL", ir_282_variants$outcome)
ir_282_variants$outcome <- ifelse(ir_282_variants$outcome == "ratio", "TG/HDL", ir_282_variants$outcome)
ir_282_variants$outcome <- ifelse(ir_282_variants$outcome == "Anterior thight fat", "Anterior thigh fat", ir_282_variants$outcome)
ir_282_variants$outcome <- ifelse(ir_282_variants$outcome == "Poster thigh fat", "Posterior thigh fat", ir_282_variants$outcome)

#Careful with duplicates due to different runs:

ir_282_variants <- ir_282_variants[which(duplicated(ir_282_variants) == FALSE),]

#Great, let's make confidence intervals here and save the whole dataframe. 
#Later we can modify it to make the plots.

dir.create("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/output/5_enrichment_analyses/1_prs/3_clean_results/")

fwrite(ir_282_variants, "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/output/5_enrichment_analyses/1_prs/3_clean_results/282_weighted_prs.csv")

#Let's plot the data for the genetic correlations

ir_282_variants$type = NA

IR <- c("FIadjBMI", "HDL", "TG")
glycemic <- c("TG/HDL", "ISIadjBMI", "IFCadjBMI", "FGadjBMI", "2hGadjBMI")
anthropometric <- c("BMI", "WHR", "WHRadjBMI", "WC", "WCadjBMI", "HC", "HCadjBMI")
ir_disease <- c("T2D", "CHD", "CKD", "Hypertension", "NAFLD", "PCOS")

ir_282_variants <- ir_282_variants[which(ir_282_variants$outcome%in%c(IR, glycemic, anthropometric, ir_disease)),]

ir_282_variants$type <- ifelse(ir_282_variants$outcome%in%IR, "Hallmark IR", ir_282_variants$type)
ir_282_variants$type <- ifelse(ir_282_variants$outcome%in%glycemic, "Glycemic", ir_282_variants$type)
ir_282_variants$type <- ifelse(ir_282_variants$outcome%in%anthropometric, "Anthropometric", ir_282_variants$type)
ir_282_variants$type <- ifelse(ir_282_variants$outcome%in%ir_disease, "IR-derived disease", ir_282_variants$type)

#Let's save this so that we can fit it in the final table:

p_1 <- grs_plotter_first(ir_282_variants, "Weighted PRS with 282 IR variants", weighted= TRUE)

######################################
#Let's do now those that decrease BMI#
######################################

list_of_files <- list.files(path_4_grs)

list_of_files <- list_of_files[which(str_detect(list_of_files, "unweighted") == TRUE)]

ir_282_variants <- list_of_files[which(str_detect(list_of_files, "ir_variants"))]
ir_282_variants <- ir_282_variants[which(str_detect(ir_282_variants, "non") == FALSE)]

#Let's get the files for each trait:

ir_282_variants <- load_grs_data(ir_282_variants)

#Let's clean the weird one:

ir_282_variants$outcome <- ifelse(ir_282_variants$outcome == "ratio", "TG/HDL", ir_282_variants$outcome)
ir_282_variants$outcome <- ifelse(ir_282_variants$outcome == "ratio", "TG/HDL", ir_282_variants$outcome)
ir_282_variants$outcome <- ifelse(ir_282_variants$outcome == "Anterior thight fat", "Anterior thigh fat", ir_282_variants$outcome)
ir_282_variants$outcome <- ifelse(ir_282_variants$outcome == "Poster thigh fat", "Posterior thigh fat", ir_282_variants$outcome)

#Careful with duplicates due to different runs:

ir_282_variants <- ir_282_variants[which(duplicated(ir_282_variants) == FALSE),]

#Great, let's make confidence intervals here and save the whole dataframe. 
#Later we can modify it to make the plots.

dir.create("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/output/5_enrichment_analyses/1_prs/3_clean_results/")

fwrite(ir_282_variants, "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/output/5_enrichment_analyses/1_prs/3_clean_results/282_unweighted_prs.csv")

#Let's plot the data for the genetic correlations

ir_282_variants$type = NA

IR <- c("FIadjBMI", "HDL", "TG")
glycemic <- c("TG/HDL", "ISIadjBMI", "IFCadjBMI", "FGadjBMI", "2hGadjBMI")
anthropometric <- c("BMI", "WHR", "WHRadjBMI", "WC", "WCadjBMI", "HC", "HCadjBMI")
ir_disease <- c("T2D", "CHD", "CKD", "Hypertension", "NAFLD", "PCOS")

ir_282_variants <- ir_282_variants[which(ir_282_variants$outcome%in%c(IR, glycemic, anthropometric, ir_disease)),]

ir_282_variants$type <- ifelse(ir_282_variants$outcome%in%IR, "Hallmark IR", ir_282_variants$type)
ir_282_variants$type <- ifelse(ir_282_variants$outcome%in%glycemic, "Glycemic", ir_282_variants$type)
ir_282_variants$type <- ifelse(ir_282_variants$outcome%in%anthropometric, "Anthropometric", ir_282_variants$type)
ir_282_variants$type <- ifelse(ir_282_variants$outcome%in%ir_disease, "IR-derived disease", ir_282_variants$type)

#Let's save this so that we can fit it in the final table:

p_2 <- grs_plotter_last(ir_282_variants, "Unweighted PRS with 282 IR variants", weighted= FALSE)

######################
#Wrapping up the plot#
######################

library(patchwork)

p_all <- p_1 + p_2 + plot_layout(ncol = 2)

ggsave(
  filename = "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/manuscript/figures/prs_282_weighted_vs_unweigthed_supplementary_figure.svg",
  plot = p_all,
  width = 14,
  height = 10,
  units = "in",  # or "in" depending on your preference, but cm is common
  device = "svg"
)
