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



grs_plotter_first <- function(grs_df, title){
  
  # Format betas and CIs
  betas_rounded     <- format(as.numeric(grs_df$betas), scientific = TRUE, digits = 4)
  lower_ci_rounded  <- format(as.numeric(grs_df$lower_ci), scientific = TRUE, digits = 4)
  upper_ci_rounded  <- format(as.numeric(grs_df$upper_ci), scientific = TRUE, digits = 4)
  
  grs_df$full_effect <- paste(betas_rounded, " (", lower_ci_rounded, ",", upper_ci_rounded, ")", sep = "")
  
  grs_df$Trait <- factor(grs_df$outcome, levels = c("FIadjBMI", "HDL", "ISIadjBMI", "TG",  "TG/HDL", "BMI", "HC", "HCadjBMI", "WC", "WCadjBMI", "WHR", "WHRadjBMI", "CHD", "CKD", "Hypertension", "NAFLD", "PCOS", "T2D", "ASAT", "ASATadjBMI", "GSAT", "GSATadjBMI", "VAT", "VATadjBMI", "Anterior thigh fat", "Posterior thigh fat", "Liver fat %", "Pancreas fat %"))
  grs_df$type  <- factor(grs_df$type, levels = c("Anthropometric", "Fat accumulation (MRI)", "IR-derived disease"))
  
  # Ensure that pvals is numeric and format it
  grs_df$pval_formatted <- paste0("p=", format(as.numeric(grs_df$pvals), scientific = TRUE, digits = 2))
  
  # Get global rightmost point for consistent label placement
  label_x_pos <- max(as.numeric(grs_df$upper_ci), na.rm = TRUE) + 0.02  # you can tweak the 0.02 spacing
  
  plotio <- ggplot(grs_df, aes(x = fct_rev(Trait), y = as.numeric(betas))) +
    geom_point(aes(color=type), size = 5, position = position_dodge(width = 0.75)) +
    geom_errorbar(aes(ymin = as.numeric(lower_ci), ymax = as.numeric(upper_ci)),
                  width = 0.15, position = position_dodge(width = 0.75),
                  color = "black", linewidth = 1) +
    facet_wrap(~type, strip.position = "left", nrow = 7, scales = "free_y") +
    geom_hline(yintercept = 0, color = "red", size = 0.5, linetype = 2) +
    
    # Add aligned p-values
    #geom_text(aes(y = label_x_pos, label = pval_formatted),
    #          hjust = 0, size = 3.5, color = "black") +
    
    # Ensure enough space to show p-value labels
    #expand_limits(y = label_x_pos + 0.02) +
    
    coord_flip() +
    xlab("Trait") +
    ylab("effect sizes") +
    ggtitle(paste0(title)) +
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

grs_plotter_mid <- function(grs_df, title){
  
  # Format betas and CIs
  betas_rounded     <- format(as.numeric(grs_df$betas), scientific = TRUE, digits = 4)
  lower_ci_rounded  <- format(as.numeric(grs_df$lower_ci), scientific = TRUE, digits = 4)
  upper_ci_rounded  <- format(as.numeric(grs_df$upper_ci), scientific = TRUE, digits = 4)
  
  grs_df$full_effect <- paste(betas_rounded, " (", lower_ci_rounded, ",", upper_ci_rounded, ")", sep = "")
  
  grs_df$Trait <- factor(grs_df$outcome, levels = c("FIadjBMI", "HDL", "ISIadjBMI", "TG",  "TG/HDL", "BMI", "HC", "HCadjBMI", "WC", "WCadjBMI", "WHR", "WHRadjBMI", "CHD", "CKD", "Hypertension", "NAFLD", "PCOS", "T2D", "ASAT", "ASATadjBMI", "GSAT", "GSATadjBMI", "VAT", "VATadjBMI", "Anterior thigh fat", "Posterior thigh fat", "Liver fat %", "Pancreas fat %"))
  grs_df$type  <- factor(grs_df$type, levels = c("Anthropometric", "Fat accumulation (MRI)", "IR-derived disease"))
  
  # Ensure that pvals is numeric and format it
  grs_df$pval_formatted <- paste0("p=", format(as.numeric(grs_df$pvals), scientific = TRUE, digits = 2))
  
  # Get global rightmost point for consistent label placement
  label_x_pos <- max(as.numeric(grs_df$upper_ci), na.rm = TRUE) + 0.02  # you can tweak the 0.02 spacing
  
  plotio <- ggplot(grs_df, aes(x = fct_rev(Trait), y = as.numeric(betas))) +
    geom_point(aes(color=type), size = 5, position = position_dodge(width = 0.75)) +
    geom_errorbar(aes(ymin = as.numeric(lower_ci), ymax = as.numeric(upper_ci)),
                  width = 0.15, position = position_dodge(width = 0.75),
                  color = "black", linewidth = 1) +
    facet_wrap(~type, strip.position = "left", nrow = 7, scales = "free_y") +
    geom_hline(yintercept = 0, color = "red", size = 0.5, linetype = 2) +
    
    # Add aligned p-values
    #geom_text(aes(y = label_x_pos, label = pval_formatted),
    #          hjust = 0, size = 3.5, color = "black") +
    
    # Ensure enough space to show p-value labels
    #expand_limits(y = label_x_pos + 0.02) +
    
    coord_flip() +
    xlab("Trait") +
    ylab("effect sizes") +
    ggtitle(paste0(title)) +
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

######################################
#Let's do now those that decrease BMI#
######################################

path_4_grs <- "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/output/5_enrichment_analyses/1_prs/2_prs_bmi_neg/"

setwd(path_4_grs)

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

#Let's try making some types so that we can make it easier:

ir_282_variants$type = NA

fat_depot <- c("ASAT", "ASATadjBMI", "GSAT", "GSATadjBMI", "VAT", "VATadjBMI", "Liver fat %", "Pancreas fat %", "Anterior thigh fat", "Posterior thigh fat")
anthropometric <- c("BMI", "WHR", "WHRadjBMI", "WC", "WCadjBMI", "HC", "HCadjBMI")
ir_disease <- c("T2D", "CHD", "CKD", "Hypertension", "NAFLD", "PCOS")

ir_282_variants <- ir_282_variants[which(ir_282_variants$outcome%in%c(fat_depot, anthropometric, ir_disease)),]

ir_282_variants$type <- ifelse(ir_282_variants$outcome%in%fat_depot, "Fat accumulation (MRI)", ir_282_variants$type)
ir_282_variants$type <- ifelse(ir_282_variants$outcome%in%anthropometric, "Anthropometric", ir_282_variants$type)
ir_282_variants$type <- ifelse(ir_282_variants$outcome%in%ir_disease, "IR-derived disease", ir_282_variants$type)

#Let's save this so that we can fit it in the final table:

p_2 <- grs_plotter_mid(ir_282_variants, "63 IR variants with BMI-decreasing effects")

#####################################################
#Let's do now those that are not associated with BMI#
#####################################################

path_4_grs <- "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/output/5_enrichment_analyses/1_prs/2_prs_bmi_ns/"

setwd(path_4_grs)

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

#Let's try making some types so that we can make it easier:

ir_282_variants$type = NA

fat_depot <- c("ASAT", "ASATadjBMI", "GSAT", "GSATadjBMI", "VAT", "VATadjBMI", "Liver fat %", "Pancreas fat %", "Anterior thigh fat", "Posterior thigh fat")
anthropometric <- c("BMI", "WHR", "WHRadjBMI", "WC", "WCadjBMI", "HC", "HCadjBMI")
ir_disease <- c("T2D", "CHD", "CKD", "Hypertension", "NAFLD", "PCOS")

ir_282_variants <- ir_282_variants[which(ir_282_variants$outcome%in%c(fat_depot, anthropometric, ir_disease)),]

ir_282_variants$type <- ifelse(ir_282_variants$outcome%in%fat_depot, "Fat accumulation (MRI)", ir_282_variants$type)
ir_282_variants$type <- ifelse(ir_282_variants$outcome%in%anthropometric, "Anthropometric", ir_282_variants$type)
ir_282_variants$type <- ifelse(ir_282_variants$outcome%in%ir_disease, "IR-derived disease", ir_282_variants$type)

#Let's save this so that we can fit it in the final table:

p_1 <- grs_plotter_first(ir_282_variants, "141 IR variants with BMI non-significant effects")

################################################################
#Let's do now those that are associated with an increase of BMI#
################################################################

path_4_grs <- "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/output/5_enrichment_analyses/1_prs/2_prs_bmi_pos_wo_outlier/"

setwd(path_4_grs)

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

#Let's try making some types so that we can make it easier:

ir_282_variants$type = NA

fat_depot <- c("ASAT", "ASATadjBMI", "GSAT", "GSATadjBMI", "VAT", "VATadjBMI", "Liver fat %", "Pancreas fat %", "Anterior thigh fat", "Posterior thigh fat")
anthropometric <- c("BMI", "WHR", "WHRadjBMI", "WC", "WCadjBMI", "HC", "HCadjBMI")
ir_disease <- c("T2D", "CHD", "CKD", "Hypertension", "NAFLD", "PCOS")

ir_282_variants <- ir_282_variants[which(ir_282_variants$outcome%in%c(fat_depot, anthropometric, ir_disease)),]

ir_282_variants$type <- ifelse(ir_282_variants$outcome%in%fat_depot, "Fat accumulation (MRI)", ir_282_variants$type)
ir_282_variants$type <- ifelse(ir_282_variants$outcome%in%anthropometric, "Anthropometric", ir_282_variants$type)
ir_282_variants$type <- ifelse(ir_282_variants$outcome%in%ir_disease, "IR-derived disease", ir_282_variants$type)

#Let's save this so that we can fit it in the final table:

p_3 <- grs_plotter_mid(ir_282_variants, "63 IR variants with BMI-increasing effects")

################################################################
#Let's do now those that are associated with an increase of BMI#
################################################################

path_4_grs <- "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/output/5_enrichment_analyses/1_prs/2_prs_bmi_pos_outlier/"

setwd(path_4_grs)

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

#Let's try making some types so that we can make it easier:

ir_282_variants$type = NA

fat_depot <- c("ASAT", "ASATadjBMI", "GSAT", "GSATadjBMI", "VAT", "VATadjBMI", "Liver fat %", "Pancreas fat %", "Anterior thigh fat", "Posterior thigh fat")
anthropometric <- c("BMI", "WHR", "WHRadjBMI", "WC", "WCadjBMI", "HC", "HCadjBMI")
ir_disease <- c("T2D", "CHD", "CKD", "Hypertension", "NAFLD", "PCOS")

ir_282_variants <- ir_282_variants[which(ir_282_variants$outcome%in%c(fat_depot, anthropometric, ir_disease)),]

ir_282_variants$type <- ifelse(ir_282_variants$outcome%in%fat_depot, "Fat accumulation (MRI)", ir_282_variants$type)
ir_282_variants$type <- ifelse(ir_282_variants$outcome%in%anthropometric, "Anthropometric", ir_282_variants$type)
ir_282_variants$type <- ifelse(ir_282_variants$outcome%in%ir_disease, "IR-derived disease", ir_282_variants$type)

#Let's save this so that we can fit it in the final table:

p_4 <- grs_plotter_mid(ir_282_variants, "15 IR variants with BMI-increasing, pleiotropic effects")

######################
#Wrapping up the plot#
######################

library(patchwork)

p_all <- p_1 + p_2 + p_3 + plot_layout(ncol = 3)

#Let's save the plot and see...

tiff("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/output/5_enrichment_analyses/1_prs/3_plots/prs_weigthed_figure_3.tiff", width = 5000, height = 3000, res=300)
p_all
dev.off()
