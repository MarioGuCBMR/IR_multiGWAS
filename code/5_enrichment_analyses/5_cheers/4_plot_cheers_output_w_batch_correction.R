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
    scale_y_continuous(limits = c(0, 6)) +  # <--- sets range of x-axis (flipped)
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

cheers_enrichment_full = read.table("output/5_enrichment_analyses/5_cheers/2_cheers_output/proxies_all_disease_enrichment_pValues.txt")
cheers_enrichment_bmi_ns = read.table("output/5_enrichment_analyses/5_cheers/2_cheers_output/proxies_bmi_ns_disease_enrichment_pValues.txt")
cheers_enrichment_bmi_neg = read.table("output/5_enrichment_analyses/5_cheers/2_cheers_output/proxies_bmi_neg_disease_enrichment_pValues.txt")
cheers_enrichment_bmi_pos = read.table("output/5_enrichment_analyses/5_cheers/2_cheers_output/proxies_bmi_pos_disease_enrichment_pValues.txt")

#Let's clean a bit the data:

cheers_enrichment_full$day <- sapply(cheers_enrichment_full$V1, day_parser)
cheers_enrichment_full$V1 <- sapply(cheers_enrichment_full$V1, cleaning_samples) #the function will take care of the naming conventions

cheers_enrichment_bmi_ns$day <- sapply(cheers_enrichment_bmi_ns$V1, day_parser)
cheers_enrichment_bmi_ns$V1 <- sapply(cheers_enrichment_bmi_ns$V1, cleaning_samples) #the function will take care of the naming conventions

cheers_enrichment_bmi_neg$day <- sapply(cheers_enrichment_bmi_neg$V1, day_parser)
cheers_enrichment_bmi_neg$V1 <- sapply(cheers_enrichment_bmi_neg$V1, cleaning_samples) #the function will take care of the naming conventions

cheers_enrichment_bmi_pos$day <- sapply(cheers_enrichment_bmi_pos$V1, day_parser)
cheers_enrichment_bmi_pos$V1 <- sapply(cheers_enrichment_bmi_pos$V1, cleaning_samples) #the function will take care of the naming conventions

#Let's plot the data:

all_proxies_plot <- cheers_plotting_unique(cheers_enrichment_full, title="")
bmi_ns_plot <- cheers_plotting_first(cheers_enrichment_bmi_ns, title="141 BMI (P>0.05) IR loci")
bmi_neg_plot <- cheers_plotting_mid(cheers_enrichment_bmi_neg, title="63 BMI- (P<0.05) IR loci")
bmi_pos_plot <- cheers_plotting_last(cheers_enrichment_bmi_pos, title="78 BMI+ (P<0.05) IR loci")

#Let's do the patchwork to make all plots together:

library(patchwork)

ir_clusters <- bmi_ns_plot + bmi_neg_plot + bmi_pos_plot + plot_layout(ncol = 3)

#####################
#Let's save the data#
#####################

ggsave(
  filename = "manuscript/figures/cheers_all_ir_enrichment.svg",
  plot = all_proxies_plot,
  width = 6,
  height = 8,
  units = "in",  # or "in" depending on your preference, but cm is common
  device = "svg"
)

ggsave(
 filename = "manuscript/figures/cheers_clusters_ir_enrichment.svg",
 plot = ir_clusters,
 width = 14,
 height = 8,
 units = "in",  # or "in" depending on your preference, but cm is common
 device = "svg"
)
