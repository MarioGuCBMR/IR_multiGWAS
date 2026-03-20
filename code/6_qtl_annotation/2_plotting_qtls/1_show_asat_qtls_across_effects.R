##############
#INTRODUCTION#
##############

#This code takes the summary statistics of several traits and plots them so that we can see the pattern of associations 
#of the IR variants with certain fat depot traits

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)
library(TwoSampleMR)
library(jsonlite)
library(dplyr)
library(tidyr)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggrepel)

###################
#Loading functions#
###################

retrieving_qtls <- function(ir_set, qtl_set){
  
  #Let's update the variants:
  
  ir_set$asat_qtl <- NA
  
  for(index in seq(1, length(ir_set$variant))){
    
    #STEP 1: get row with genes:
    
    snp_row <- ir_set[index,]
    
    qtl_info_match <- qtl_set[which(qtl_set$query_snp_rsid%in%snp_row$variant),]
    
    #Let's split the data:
    
    name_asat <- qtl_info_match$gene_name.asat[which(str_detect(qtl_info_match$gene_name.asat, "ENSG") == FALSE)]
    name_asat_gtex <- qtl_info_match$gene_name.asat_gtex[which(str_detect(qtl_info_match$gene_name.asat_gtex, "ENSG") == FALSE)]
    name_asat_sqtl <- qtl_info_match$gene_name.asat_sqtl[which(str_detect(qtl_info_match$gene_name.asat_sqtl, "ENSG") == FALSE)]
    
    gene_info <- c(name_asat, name_asat_gtex, name_asat_sqtl)
    gene_info <- gene_info[which(is.na(gene_info) == FALSE)]
    gene_info <- as.character(unlist(str_split(gene_info, ";")))
    gene_info <- unique(gene_info)
    gene_info <- gene_info[order(gene_info)]
    
    gene_info <- paste(gene_info, collapse = "/")
    
    #Let's add it to the variant:
    
    ir_set$asat_qtl[index] <- gene_info
    
    
  }
  
  return(ir_set)
  
}

preparing_data <- function(ir_set, pheno_df){
  
  ir_classes <- ir_set[, c("variant", "reported_ir", "asat_qtl")]
  data_classes <- merge(pheno_df, ir_classes, by.x = "SNP", by.y = "variant", all.x = T) ; nrow(data_classes) ; head(data_classes, 3)
  
  #Now let's add a phenotype_group variable:
  
  data_all <- data_classes
  
  #Let's do some changes so that this is visually appealing:
  
  data_all$phenotype_name <- ifelse(data_all$phenotype_name == "tg", "TG/HDL", data_all$phenotype_name)
  data_all$phenotype_name <- ifelse(data_all$phenotype_name == "Anterior thight fat", "Anterior thigh fat", data_all$phenotype_name)
  data_all$phenotype_name <- ifelse(data_all$phenotype_name == "Poster thigh fat", "Posterior thigh fat", data_all$phenotype_name)
  
  #Let's only take the following traits:
  
  data_all <- data_all[which(data_all$phenotype_name%in%c("BMI", "WHRadjBMI", "HC", "GSAT", "GSATadjBMI", "VAT", "VATadjBMI", "Liver fat %", "Liver volume", "Pancreas fat %", "Posterior thigh fat", "Anterior thigh fat")),]
  
  anthropometric <- c("BMI", "WHR", "WHRadjBMI", "HC", "WC", "GSAT", "GSATadjBMI")
  ectopic_fat <- c("VAT", "VATadjBMI", "Liver fat %", "Liver volume", "Pancreas fat %", "Posterior thigh fat", "Anterior thigh fat")
  
  data_all$phenotype_group <- NA
  data_all$phenotype_group <- ifelse(data_all$phenotype_name%in%anthropometric, "Fat distribution", data_all$phenotype_group)
  data_all$phenotype_group <- ifelse(data_all$phenotype_name%in%ectopic_fat, "Ectopic fat", data_all$phenotype_group)
  
  #We need a couple of additional columns:
  
  data_all$pValue <- data_all$pval.outcome 
  data_all$beta <- data_all$beta.outcome 
  data_all$cluster <- data_all$reported_ir 
  
  return(data_all)
  
}

plotting_info_first <- function(combined_df){
  #This functions performs a set of transformations and makes the first plot to be generated in a row of a subset of variants.
  
  # Define p-value bins
  combined_df <- combined_df %>%
    mutate(p_bin = factor(case_when(
      pValue < 5e-8  ~ "<5E-08",
      pValue < 0.05 & pValue > 5e-08 ~ "<0.05",
      TRUE ~ ">=0.05"
    ), levels = c("<5E-08", "<0.05", ">=0.05"))) # Ensure correct order
  
  # Assign triangle shape based on beta sign
  combined_df <- combined_df %>%
    mutate(shape = ifelse(beta > 0, "▲", "▼"))
  
  combined_df <- combined_df %>%
    mutate(cluster = recode(cluster,
                            "no" = "Novel",
                            "yes" = "Previously known"))
  combined_df <- combined_df %>%
    mutate(cluster = factor(cluster, levels = c("Novel", "Previously known")))
  
  # combined_df_sel <- rbind(head(combined_df, 50), tail(combined_df, 50))
  
  # Make sure phenotype_name is a factor to control order
  combined_df$phenotype_name <- factor(combined_df$phenotype_name, levels = unique(combined_df$phenotype_name))
  
  #Let's remove the asat_qtl with either NAs or with ""
  
  combined_df <- combined_df[which(combined_df$asat_qtl != "" & is.na(combined_df$asat_qtl) == FALSE),]
  
  # Make sure asat_qtl (or gene/SNP label column) is a factor with appropriate order
  combined_df$asat_qtl <- factor(combined_df$asat_qtl, levels = sort(unique(combined_df$asat_qtl)))
  
  # Plot
  plot_transposed <- ggplot(combined_df, aes(x = phenotype_name, y = fct_rev(asat_qtl))) +
    geom_tile(aes(fill = shape, alpha = p_bin), color = "white", size = 0.5) +
    scale_alpha_manual(values = c("<5E-08" = 1, "<0.05" = 0.6, ">=0.05" = 0.1)) +
    geom_text(data = subset(combined_df, p_bin == "<5E-08"), aes(label = "*"), size = 5) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),  # make phenotype labels readable
      axis.text.y = element_text(size = 8),
      legend.position = "top"
    ) +
    facet_grid(cluster ~ phenotype_group, scales = "free", space = "free", switch = "y") +
    labs(title = "",
         x = "",
         y = "",
         alpha = "P-value",
         fill = "Beta")
  
}

plotting_info_first_no_x <- function(combined_df){
  
  combined_df <- combined_df %>%
    mutate(p_bin = factor(case_when(
      pValue < 5e-8 ~ "<5E-08",
      pValue < 0.05 & pValue > 5e-08 ~ "<0.05",
      TRUE ~ ">=0.05"
    ), levels = c("<5E-08", "<0.05", ">=0.05")),
    shape = ifelse(beta > 0, "▲", "▼"),
    cluster = recode(cluster, "no" = "Novel", "yes" = "Previously known"),
    cluster = factor(cluster, levels = c("Novel", "Previously known"))
    )
  
  combined_df$phenotype_name <- factor(combined_df$phenotype_name, levels = unique(combined_df$phenotype_name))
  combined_df <- combined_df[which(combined_df$asat_qtl != "" & !is.na(combined_df$asat_qtl)), ]
  combined_df$asat_qtl <- factor(combined_df$asat_qtl, levels = sort(unique(combined_df$asat_qtl)))
  
  ggplot(combined_df, aes(x = phenotype_name, y = fct_rev(asat_qtl))) +
    geom_tile(aes(fill = shape, alpha = p_bin), color = "white", size = 0.5) +
    scale_alpha_manual(values = c("<5E-08" = 1, "<0.05" = 0.6, ">=0.05" = 0.1)) +
    geom_text(data = subset(combined_df, p_bin == "<5E-08"), aes(label = "*"), size = 5) +
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 8),
      legend.position = "top"
    ) +
    facet_grid(cluster ~ phenotype_group, scales = "free", space = "free", switch = "y") +
    labs(title = "", x = "", y = "", alpha = "P-value", fill = "Beta")
}

plotting_info_first_no_x_no_legend <- function(combined_df){
  
  combined_df <- combined_df %>%
    mutate(p_bin = factor(case_when(
      pValue < 5e-8 ~ "<5E-08",
      pValue < 0.05 & pValue > 5e-08 ~ "<0.05",
      TRUE ~ ">=0.05"
    ), levels = c("<5E-08", "<0.05", ">=0.05")),
    shape = ifelse(beta > 0, "▲", "▼"),
    cluster = recode(cluster, "no" = "Novel", "yes" = "Previously known"),
    cluster = factor(cluster, levels = c("Novel", "Previously known"))
    )
  
  combined_df$phenotype_name <- factor(combined_df$phenotype_name, levels = unique(combined_df$phenotype_name))
  combined_df <- combined_df[which(combined_df$asat_qtl != "" & !is.na(combined_df$asat_qtl)), ]
  combined_df$asat_qtl <- factor(combined_df$asat_qtl, levels = sort(unique(combined_df$asat_qtl)))
  
  ggplot(combined_df, aes(x = phenotype_name, y = fct_rev(asat_qtl))) +
    geom_tile(aes(fill = shape, alpha = p_bin), color = "white", size = 0.5) +
    scale_alpha_manual(values = c("<5E-08" = 1, "<0.05" = 0.6, ">=0.05" = 0.1)) +
    geom_text(data = subset(combined_df, p_bin == "<5E-08"), aes(label = "*"), size = 5) +
    theme_classic() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 8),
      legend.position = "none"
    ) +
    facet_grid(cluster ~ phenotype_group, scales = "free", space = "free", switch = "y") +
    labs(title = "", x = "", y = "")
}

plotting_info_first_no_legend <- function(combined_df){
  
  combined_df <- combined_df %>%
    mutate(p_bin = factor(case_when(
      pValue < 5e-8 ~ "<5E-08",
      pValue < 0.05 & pValue > 5e-08 ~ "<0.05",
      TRUE ~ ">=0.05"
    ), levels = c("<5E-08", "<0.05", ">=0.05")),
    shape = ifelse(beta > 0, "▲", "▼"),
    cluster = recode(cluster, "no" = "Novel", "yes" = "Previously known"),
    cluster = factor(cluster, levels = c("Novel", "Previously known"))
    )
  
  combined_df$phenotype_name <- factor(combined_df$phenotype_name, levels = unique(combined_df$phenotype_name))
  combined_df <- combined_df[which(combined_df$asat_qtl != "" & !is.na(combined_df$asat_qtl)), ]
  combined_df$asat_qtl <- factor(combined_df$asat_qtl, levels = sort(unique(combined_df$asat_qtl)))
  
  ggplot(combined_df, aes(x = phenotype_name, y = fct_rev(asat_qtl))) +
    geom_tile(aes(fill = shape, alpha = p_bin), color = "white", size = 0.5) +
    scale_alpha_manual(values = c("<5E-08" = 1, "<0.05" = 0.6, ">=0.05" = 0.1)) +
    geom_text(data = subset(combined_df, p_bin == "<5E-08"), aes(label = "*"), size = 5) +
    theme_classic() +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
      axis.text.y = element_text(size = 8),
      legend.position = "none"
    ) +
    facet_grid(cluster ~ phenotype_group, scales = "free", space = "free", switch = "y") +
    labs(title = "", x = "", y = "")
}

##############
#Loading data#
##############

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/")

ir_variants <- fread("output/4_ir_loci_discovery/3_novel_variants/3_novel_variants/282_ir_variants_with_ir_and_t2d_comparisons_w_closest_genes.txt")
qtl_df <- readRDS("output/6_qtl_annotation/conditional_and_fine_mapped_asat_and_vat_eqtls_sqtls.txt")

#Let's retrieve the ASAT QTLs and add them as a column:

ir_variants <- retrieving_qtls(ir_variants, qtl_df)

#Now let's add the data for each subset of variants:

#Let's load only the fat distribution variants:

bmi_ns <- fread("output/5_enrichment_analyses/3_go_shifter/ld_bmi_ns.txt")
bmi_ns <- ir_variants[which(ir_variants$variant%in%bmi_ns$RsIdA),]

bmi_neg <- fread("output/5_enrichment_analyses/3_go_shifter/ld_bmi_neg.txt")
bmi_neg <- ir_variants[which(ir_variants$variant%in%bmi_neg$RsIdA),]

bmi_pos <- fread("output/5_enrichment_analyses/3_go_shifter/ld_bmi_pos.txt")
bmi_pos <- ir_variants[which(ir_variants$variant%in%bmi_pos$RsIdA),]

bmi_outlier <- fread("output/5_enrichment_analyses/3_go_shifter/ld_outliers.txt")
bmi_outlier <- ir_variants[which(ir_variants$variant%in%bmi_outlier$RsIdA),]

###################################
#Let's load now the phenotype data#
###################################

#Let's set up the data we want:

list_of_files <- list.files("output/5_enrichment_analyses/1_prs/1_expo_outcome_df/")
list_of_files <- list_of_files[which(str_detect(list_of_files, "all") == FALSE)]
list_of_files <- list_of_files[which(str_detect(list_of_files, "no") == FALSE)]

#Let's filter for the phenotypes from UKBB that we know we are going to find data:

list_of_files <- list_of_files[which(str_detect(list_of_files, "thigh") == TRUE |
                                     str_detect(list_of_files, "Liver") == TRUE | 
                                     str_detect(list_of_files, "Pancreas") == TRUE | 
                                     str_detect(list_of_files, "tg_hdl") == TRUE |
                                       str_detect(list_of_files, "WHR") == TRUE |
                                       str_detect(list_of_files, "WC") == TRUE |
                                       str_detect(list_of_files, "HC") == TRUE |
                                       str_detect(list_of_files, "BMI") == TRUE |
                                       str_detect(list_of_files, "BMI") == TRUE |
                                       str_detect(list_of_files, "AT") == TRUE)]

no_vect <- c("ir_variants_FIadjBMI.txt", "ir_variants_ISIadjBMI.txt", "ir_variants_IFCadjBMI.txt", "ir_variants_HCadjBMI.txt", "ir_variants_WCadjBMI.txt", "ir_variants_T2D (Suzuki).txt")

list_of_files <- list_of_files[which(!(list_of_files%in%no_vect))]
                                       
for(i in seq(length(list_of_files))){
  
  #STEP 1: load the file
  
  tmp_df <- fread(paste("output/5_enrichment_analyses/1_prs/1_expo_outcome_df/", list_of_files[i], sep = ""))
  
  print(length(tmp_df$SNP))
  
  #STEP 2: add phenotype name
  
  pheno <- unlist(str_split(list_of_files[i], ".txt"))[1]
  pheno <- unlist(str_split(pheno, "_"))[3]
  
  tmp_df$phenotype_name <- pheno
  
  #STEP 3: make input dataframe:
  
  if(!(exists("final_df"))){
    
    final_df <- tmp_df
    
  } else {
    
    final_df <- rbind(final_df, tmp_df)
  }
  
}

####################################################################
#Let's format the data so that we can actually see the associations#
####################################################################

bmi_ns_4_plotting <- preparing_data(bmi_ns, final_df)
bmi_neg_4_plotting <- preparing_data(bmi_neg, final_df)
bmi_pos_4_plotting <- preparing_data(bmi_pos, final_df)
bmi_out_4_plotting <- preparing_data(bmi_outlier, final_df)

#And now let's plot:

bmi_ns_plot <- plotting_info_first_no_x(bmi_ns_4_plotting)
bmi_neg_plot <- plotting_info_first_no_legend(bmi_neg_4_plotting)
bmi_pos_plot <- plotting_info_first_no_x(bmi_pos_4_plotting)
bmi_outlier_plot <- plotting_info_first_no_legend(bmi_out_4_plotting)

###############################################
#Let's add info for beautiful plotting options#
###############################################

library(patchwork)

plot_all <- bmi_ns_plot + bmi_pos_plot +  bmi_neg_plot + bmi_outlier_plot + patchwork::plot_layout(ncol = 2)

tiff("output/6_qtl_annotation/asat_qtls_with_effects_on_fat.tiff", height = 9000, width = 9000, res=300)
plot_all
dev.off()
