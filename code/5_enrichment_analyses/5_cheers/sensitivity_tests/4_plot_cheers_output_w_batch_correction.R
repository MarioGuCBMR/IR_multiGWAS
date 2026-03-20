##############
#INTRODUCTION#
##############

#This code plots the results for CHEERS enrichment!

###################
#Loading libraries#
###################

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

###################
#Loading functions#
###################

# cleaning_samples <- function(batch_){
#   
#   #STEP 1: remove combat_corrected_
#   
#   batch_clean <- unlist(str_split(batch_, "_combat_corrected"))[1]
#   
#   batch_data <- unlist(str_split(batch_clean, "_"))
#   
#   cell <- batch_data[1]
#   day <- batch_data[2]
#   batch <- batch_data[3]
#   rep <- batch_data[4]
#   
#   #Let's clean the day:
#   
#   day <- unlist(str_split(day, "day"))[2]
#   day <- paste("Day", day)
#   
#   #Same with batch:
#   
#   batch <- unlist(str_split(batch, "batch"))[2]
#   batch <- paste("Batch", batch)
#   
#   #Now with replicates:
#   
#   rep <- unlist(str_split(rep, "rep"))[2]
#   rep <- paste("Replicate", rep) 
#   
#   #Final string:
#   
#   batch_end <- paste(day, " (", batch, ";", rep, ")", sep = "")
#   
#   return(batch_end)
#   
# }
# 
# 
# # Define a consistent color palette for all plots
# cell_stage_palette <- c(
#   "Day 0"  = "#E69F00",  # orange
#   "Day 4"  = "#56B4E9",  # sky blue
#   "Day 14"  = "#009E73"  # bluish green
# )
# 
# # Core plot function with customizable theme tweaks
# cheers_plot_base <- function(cheers_enrichment, title, strip_y = FALSE, show_legend = FALSE) {
#   
#   colnames(cheers_enrichment) <- c("Samples", "Pval")
#   
#   cheers_enrichment <- cheers_enrichment %>%
#     mutate(Cell_stage = case_when(
#       str_detect(Samples, "Day 0") ~ "Day 0",
#       str_detect(Samples, "Day 2") ~ "Day 2",
#       str_detect(Samples, "Day 4") ~ "Day 4",
#       TRUE ~ "Day 14"
#     )) %>%
#     mutate(Cell_stage = factor(Cell_stage, levels = c("Day 0", "Day 4", "Day 14")),
#            Samples = factor(Samples, levels = unique(Samples[order(Cell_stage)])))
#   
#   p <- ggplot(cheers_enrichment, aes(x = fct_rev(Samples), y = -log10(Pval), fill = Cell_stage)) +
#     geom_col(width = 0.8) +
#     coord_flip() +
#     geom_hline(yintercept = -log10(0.05 / 25), linetype = "dashed", color = "firebrick", size = 1) +
#     scale_fill_manual(values = cell_stage_palette) +
#     labs(
#       title = title,
#       y = expression(-log[10](italic("P"))),
#       x = NULL
#     ) +
#     theme_minimal(base_size = 12) +
#     theme(
#       plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
#       axis.text.y = if (strip_y) element_blank() else element_text(size = 10),
#       axis.ticks.y = if (strip_y) element_blank() else element_line(),
#       axis.title.y = element_blank(),
#       panel.grid.major.y = element_blank(),
#       panel.grid.minor = element_blank(),
#       legend.position = if (show_legend) "right" else "none",
#       legend.title = element_blank(),
#       legend.text = element_text(size = 10),
#       axis.line = element_line(colour = "black"),
#       plot.background = element_rect(fill = "white", color = NA)
#     )
#   
#   return(p)
# }
# 
# # Specific versions using the base with tweaks
# cheers_plotting_first <- function(cheers_enrichment, title) {
#   cheers_plot_base(cheers_enrichment, title, strip_y = FALSE, show_legend = FALSE)
# }
# 
# cheers_plotting_mid <- function(cheers_enrichment, title) {
#   cheers_plot_base(cheers_enrichment, title, strip_y = TRUE, show_legend = FALSE)
# }
# 
# cheers_plotting_last <- function(cheers_enrichment, title) {
#   cheers_plot_base(cheers_enrichment, title, strip_y = TRUE, show_legend = TRUE)
# }

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

prepping_data <- function(cheers_enrichment){
  #Uses the parser and cleaning functions to clean the data:
  
  cheers_enrichment$Day <- sapply(cheers_enrichment$V1, day_parser)
  cheers_enrichment$V1 <- sapply(cheers_enrichment$V1, cleaning_samples)
  
  return(cheers_enrichment)
  
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
  
  p <- p + scale_y_continuous(limits = c(0, 5),  breaks = seq(0, 5, by = 1))
  
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
cheers_enrichment_metabolic_syndrome = read.table("output/5_enrichment_analyses/5_cheers/3_sensitivity_tests/2_cheers_output//IR_proxies_cluster_metabolic_syndrome_disease_enrichment_pValues.txt")
cheers_enrichment_lipodystrophy = read.table("output/5_enrichment_analyses/5_cheers/3_sensitivity_tests/2_cheers_output/IR_proxies_cluster_lipodystrophy_disease_enrichment_pValues.txt")
cheers_enrichment_obesity = read.table("output/5_enrichment_analyses/5_cheers/3_sensitivity_tests/2_cheers_output/IR_proxies_cluster_obesity_disease_enrichment_pValues.txt")

cheers_enrichment_beta_cell_pi_pos = read.table("output/5_enrichment_analyses/5_cheers/3_sensitivity_tests/2_cheers_output//IR_proxies_cluster_beta_cell_pos_disease_enrichment_pValues.txt")
cheers_enrichment_beta_cell_pi_neg = read.table("output/5_enrichment_analyses/5_cheers/3_sensitivity_tests/2_cheers_output/IR_proxies_cluster_beta_cell_neg_disease_enrichment_pValues.txt")
cheers_enrichment_body_fat = read.table("output/5_enrichment_analyses/5_cheers/3_sensitivity_tests/2_cheers_output/IR_proxies_cluster_body_fat_disease_enrichment_pValues.txt")
cheers_enrichment_residual = read.table("output/5_enrichment_analyses/5_cheers/3_sensitivity_tests/2_cheers_output/IR_proxies_cluster_residual_disease_enrichment_pValues.txt")
cheers_enrichment_lipid = read.table("output/5_enrichment_analyses/5_cheers/3_sensitivity_tests/2_cheers_output/IR_proxies_cluster_liver_disease_enrichment_pValues.txt")

#Let's clean a bit the data:

cheers_enrichment_bmi_ns = prepping_data(cheers_enrichment_bmi_ns)
cheers_enrichment_bmi_neg =  prepping_data(cheers_enrichment_bmi_neg)
cheers_enrichment_bmi_pos =  prepping_data(cheers_enrichment_bmi_pos)
cheers_enrichment_metabolic_syndrome =  prepping_data(cheers_enrichment_metabolic_syndrome)
cheers_enrichment_lipodystrophy = prepping_data(cheers_enrichment_lipodystrophy)
cheers_enrichment_obesity =  prepping_data(cheers_enrichment_obesity)

cheers_enrichment_beta_cell_pi_pos =  prepping_data(cheers_enrichment_beta_cell_pi_pos)
cheers_enrichment_beta_cell_pi_neg =  prepping_data(cheers_enrichment_beta_cell_pi_neg)
cheers_enrichment_body_fat =  prepping_data(cheers_enrichment_body_fat)
cheers_enrichment_residual =  prepping_data(cheers_enrichment_residual)
cheers_enrichment_lipid =  prepping_data(cheers_enrichment_lipid)

#Let's plot the data:

#all_proxies_plot <- cheers_plotting_first(cheers_enrichment_full, title="")
bmi_ns_plot <- cheers_plotting_first(cheers_enrichment_bmi_ns, title="141 BMI-neutral IR loci")
bmi_neg_plot <- cheers_plotting_mid(cheers_enrichment_bmi_neg, title="63 BMI-decreasing  IR loci")
bmi_pos_plot <- cheers_plotting_last(cheers_enrichment_bmi_pos, title="63 BMI-increasing IR loci")

metabolic_syndrome_plot <- cheers_plotting_first(cheers_enrichment_metabolic_syndrome, title="166 metabolic syndrome T2D loci")
lipodystrophy_plot <- cheers_plotting_mid(cheers_enrichment_lipodystrophy, title="45 lipodystrophy T2D loci")
body_fat_plot <- cheers_plotting_last(cheers_enrichment_body_fat, title="273 body fat T2D loci")

# beta_cell_pos_plot <- cheers_plotting_first(cheers_enrichment_beta_cell_pi_pos, title="166 metabolic syndrome T2D SNPs")
# beta_cell_neg_plot <- cheers_plotting_mid(cheers_enrichment_beta_cell_pi_neg, title="45 lipodystrophy T2D SNPs")
# residual_plot <- cheers_plotting_mid(cheers_enrichment_residual, title="45 lipodystrophy T2D SNPs")
# 
# obesity_plot <- cheers_plotting_last(cheers_enrichment_obesity, title="53 obesity IR SNPs")
# lipid_plot <- cheers_plotting_last(cheers_enrichment_lipid, title="53 obesity IR SNPs")

#Let's do the patchwork to make all plots together:

library(patchwork)

ir_decreasing_clusters <- bmi_ns_plot + bmi_neg_plot + bmi_pos_plot + metabolic_syndrome_plot + lipodystrophy_plot + body_fat_plot  +
  #beta_cell_pos_plot + beta_cell_neg_plot + residual_plot + obesity_plot + lipid_plot 
   plot_layout(ncol = 3)

ggsave(
  filename = "manuscript/figures/cheers_clusters_and_lotta_suzuki_ir_enrichment.svg",
  plot = ir_decreasing_clusters,
  width = 13,
  height = 11,
  units = "in",  # or "in" depending on your preference, but cm is common
  device = "svg"
)

