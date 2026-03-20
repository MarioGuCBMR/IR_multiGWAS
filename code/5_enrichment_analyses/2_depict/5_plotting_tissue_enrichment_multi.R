##############
#INTRODUCTION#
##############

#This code matches our SNPs data with the traits of interest
#harmonizes them so that they are properly aligned.
#And finally performs GRS, which are saved ina dataframe.

###################
#Loading libraries#
###################

library(ggplot2)
library(dplyr)
library(scales)
library(data.table)
library(tidyverse)

###################
#Loading functions#
###################

source("N:/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/Mediation_analysis/code/0_functions/functions_4_mediation.R")

# create_enrichment_plot_grouped <- function(enrichment_data, tissue_col = "Name", logP_col = "logP", FDR_col = "FDR", 
#                                            mesh_col = "MESH", title = "Enrichment Analysis by Tissue", 
#                                            color_significant = "red", color_less_significant="blue", color_nonsignificant = "gray") {
#   
#   # Ensure that FDR and MESH are factors for proper coloring and grouping
#   enrichment_data[[FDR_col]] <- factor(enrichment_data[[FDR_col]], levels = c("<0.01", "<0.05", ">0.20"))
#   enrichment_data[[mesh_col]] <- factor(enrichment_data[[mesh_col]])  # MESH term as a factor
#   
#   # Create the plot
#   ggplot(enrichment_data, aes_string(x = tissue_col, y = logP_col, fill = FDR_col)) +
#     geom_bar(stat = "identity", show.legend = FALSE, width = 0.7) +  # Fixed bar width
#     scale_fill_manual(values = c("<0.01" = color_significant, "<0.05"=color_less_significant, ">0.05" = color_nonsignificant)) +
#     theme_minimal(base_size = 14) +  # Clean, minimal theme with larger font
#     theme(
#       axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12, face = "bold"),  # Angled labels for tissue names
#       axis.title = element_text(face = "bold", size = 14),
#       axis.text.y = element_text(size = 12),
#       panel.grid = element_blank(),  # Clean grid lines
#       plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
#       strip.text = element_text(size = 14, face = "bold")  # Bold MESH labels
#     ) +
#     labs(
#       x = "Tissue",  # Label for x-axis
#       y = expression(-log[10](p)),  # Y-axis label with math expression
#       title = title  # Customizable plot title
#     ) +
#     scale_y_continuous(limits = c(0, 12)) +  # Fix Y-axis range across all facets
#     facet_wrap(~ get(mesh_col), scales = "free_x", ncol = length(unique(enrichment_data$MESH)))  # Facet by MESH term with free y-scales
# }

create_enrichment_plot_grouped <- function(enrichment_data, tissue_col = "Name", logP_col = "logP", FDR_col = "FDR", 
                                           mesh_col = "MESH", title = "Enrichment Analysis by Tissue", 
                                           color_significant = "red", color_less_significant="blue", color_nonsignificant = "gray") {
  
  # Ensure that FDR and MESH are factors for proper coloring and grouping
  enrichment_data[[FDR_col]] <- factor(enrichment_data[[FDR_col]], levels = c("<0.01", "<0.05", ">0.20"))
  enrichment_data[[mesh_col]] <- factor(enrichment_data[[mesh_col]])  # MESH term as a factor
  
  # Get the maximum value of logP for the y-axis limit
  max_logP <- max(enrichment_data[[logP_col]], na.rm = TRUE)
  
  # Get the top 5 hits based on logP per MESH term
  top_hits <- enrichment_data %>%
    group_by(get(mesh_col)) %>%
    top_n(5, !!sym(logP_col)) %>%  # Get top 5 based on logP
    ungroup()
  
  # Create a new data frame to store the final enriched data with padding for each MESH term
  padded_data <- data.frame()
  
  # Loop over each unique MESH term
  for (mesh_term in unique(enrichment_data[[mesh_col]])) {
    
    # Subset the top hits for this MESH term
    mesh_data <- top_hits %>% filter(get(mesh_col) == mesh_term)
    
    # Calculate how many rows to add (if there are fewer than 5, add dummy rows)
    rows_to_add <- 5 - nrow(mesh_data)
    
    # Create dummy rows if needed
    if (rows_to_add > 0) {
      # Create dummy rows with unique identifiers (z_dummy_1, z_dummy_2, etc.)
      dummy_rows <- mesh_data[1, ]  # Copy the first row to create dummy rows
      dummy_rows[[logP_col]] <- 0  # Set logP to 0 for dummy rows (invisible)
      dummy_rows[[FDR_col]] <- ">0.20"  # Set FDR to ">0.20" for dummy rows (invisible)
      
      # Create dummy tissue names (z_dummy_1, z_dummy_2, etc.)
      for (i in 1:rows_to_add) {
        dummy_rows[[tissue_col]] <- paste0("z_dummy_", i)
        
        # Repeat the dummy row for the required number of rows
        mesh_data <- rbind(mesh_data, dummy_rows)
      }
    }
    
    # Add this data to the padded_data dataframe
    padded_data <- rbind(padded_data, mesh_data)
  }
  
  # Create the plot using the padded data (which includes dummy rows)
  ggplot(padded_data, aes_string(x = tissue_col, y = logP_col, fill = FDR_col)) +
    geom_bar(stat = "identity", show.legend = FALSE, width = 0.7) +  # Fixed bar width
    scale_fill_manual(values = c("<0.01" = color_significant, "<0.05" = color_less_significant, ">0.20" = color_nonsignificant)) +
    theme_minimal(base_size = 14) +  # Clean, minimal theme with larger font
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12, face = "bold"),  # Angled labels for tissue names
      axis.title = element_text(face = "bold", size = 14),
      axis.text.y = element_text(size = 12),
      axis.line = element_line(colour = "black", size = 1),  # Add axis lines
      axis.ticks = element_line(colour = "black"),  # Add axis ticks
      panel.grid = element_blank(),  # Clean grid lines
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      strip.text = element_text(size = 14, face = "bold"),  # Bold MESH labels
      strip.background = element_rect(fill = "white", color = "white"),  # Remove background of facet labels
      panel.spacing = unit(1, "lines"),  # Control spacing between facets
      axis.ticks.x = element_blank()  # Remove axis ticks for a cleaner look
    ) +
    labs(
      x = "Tissue",  # Label for x-axis
      y = expression(-log[10](p)),  # Y-axis label with math expression
      #title = title  # Customizable plot title
    ) +
    scale_y_continuous(limits = c(0, max_logP)) +  # Dynamic Y-axis range
    facet_wrap(~ get(mesh_col), scales = "free_x", ncol = length(unique(enrichment_data[[mesh_col]]))) +  # Facet by MESH term with free y-scales
    scale_x_discrete(labels = function(x) ifelse(grepl("^z_dummy_", x), "", x))  # Hide dummy labels
}

create_enrichment_plot_grouped_no_x_labels <- function(enrichment_data, tissue_col = "Name", logP_col = "logP", FDR_col = "FDR", 
                                                       mesh_col = "MESH", title = "Enrichment Analysis by Tissue", 
                                                       color_significant = "red", color_less_significant="blue", color_nonsignificant = "gray") {
  
  # Ensure that FDR and MESH are factors for proper coloring and grouping
  enrichment_data[[FDR_col]] <- factor(enrichment_data[[FDR_col]], levels = c("<0.01", "<0.05", ">0.20"))
  enrichment_data[[mesh_col]] <- factor(enrichment_data[[mesh_col]])  # MESH term as a factor
  
  # Get the maximum value of logP for the y-axis limit
  max_logP <- max(enrichment_data[[logP_col]], na.rm = TRUE)
  
  # Get the top 5 hits based on logP per MESH term
  top_hits <- enrichment_data %>%
    group_by(get(mesh_col)) %>%
    top_n(5, !!sym(logP_col)) %>%  # Get top 5 based on logP
    ungroup()
  
  # Create a new data frame to store the final enriched data with padding for each MESH term
  padded_data <- data.frame()
  
  # Loop over each unique MESH term
  for (mesh_term in unique(enrichment_data[[mesh_col]])) {
    
    # Subset the top hits for this MESH term
    mesh_data <- top_hits %>% filter(get(mesh_col) == mesh_term)
    
    # Calculate how many rows to add (if there are fewer than 5, add dummy rows)
    rows_to_add <- 5 - nrow(mesh_data)
    
    # Create dummy rows if needed
    if (rows_to_add > 0) {
      # Create dummy rows with unique identifiers (z_dummy_1, z_dummy_2, etc.)
      dummy_rows <- mesh_data[1, ]  # Copy the first row to create dummy rows
      dummy_rows[[logP_col]] <- 0  # Set logP to 0 for dummy rows (invisible)
      dummy_rows[[FDR_col]] <- ">0.20"  # Set FDR to ">0.20" for dummy rows (invisible)
      
      # Create dummy tissue names (z_dummy_1, z_dummy_2, etc.)
      for (i in 1:rows_to_add) {
        dummy_rows[[tissue_col]] <- paste0("z_dummy_", i)
        
        # Repeat the dummy row for the required number of rows
        mesh_data <- rbind(mesh_data, dummy_rows)
      }
    }
    
    # Add this data to the padded_data dataframe
    padded_data <- rbind(padded_data, mesh_data)
  }
  
  # Create the plot using the padded data (which includes dummy rows)
  ggplot(padded_data, aes_string(x = tissue_col, y = logP_col, fill = FDR_col)) +
    geom_bar(stat = "identity", show.legend = FALSE, width = 0.7) +  # Fixed bar width
    scale_fill_manual(values = c("<0.01" = color_significant, "<0.05" = color_less_significant, ">0.20" = color_nonsignificant)) +
    theme_minimal(base_size = 14) +  # Clean, minimal theme with larger font
    theme(
      axis.text.x = element_blank(),  # Remove x-axis labels
      axis.title = element_blank(),   # Remove x-axis title
      axis.text.y = element_text(size = 12),
      axis.line = element_line(colour = "black", size = 1),  # Add axis lines
      axis.ticks = element_line(colour = "black"),  # Add axis ticks
      panel.grid = element_blank(),  # Clean grid lines
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      strip.text = element_text(size = 14, face = "bold"),  # Bold MESH labels
      strip.background = element_rect(fill = "white", color = "white"),  # Remove background of facet labels
      panel.spacing = unit(1, "lines"),  # Control spacing between facets
      axis.ticks.x = element_blank()  # Remove axis ticks for a cleaner look
    ) +
    labs(
      y = expression(-log[10](p)),  # Y-axis label with math expression
      #title = title  # Customizable plot title
    ) +
    scale_y_continuous(limits = c(0, max_logP)) +  # Dynamic Y-axis range
    facet_wrap(~ get(mesh_col), scales = "free_x", ncol = length(unique(enrichment_data[[mesh_col]]))) +  # Facet by MESH term with free y-scales
    scale_x_discrete(labels = function(x) ifelse(grepl("^z_dummy_", x), "", x))  # Hide dummy labels
}

create_enrichment_plot_grouped_no_x_labels_no_facet_titles <- function(enrichment_data, tissue_col = "Name", logP_col = "logP", FDR_col = "FDR", 
                                                                       mesh_col = "MESH", title = "Enrichment Analysis by Tissue", 
                                                                       color_significant = "red", color_less_significant="blue", color_nonsignificant = "gray") {
  
  # Ensure that FDR and MESH are factors for proper coloring and grouping
  enrichment_data[[FDR_col]] <- factor(enrichment_data[[FDR_col]], levels = c("<0.01", "<0.05", ">0.20"))
  enrichment_data[[mesh_col]] <- factor(enrichment_data[[mesh_col]])  # MESH term as a factor
  
  # Get the maximum value of logP for the y-axis limit
  max_logP <- max(enrichment_data[[logP_col]], na.rm = TRUE)
  
  # Get the top 5 hits based on logP per MESH term
  top_hits <- enrichment_data %>%
    group_by(get(mesh_col)) %>%
    top_n(5, !!sym(logP_col)) %>%  # Get top 5 based on logP
    ungroup()
  
  # Create a new data frame to store the final enriched data with padding for each MESH term
  padded_data <- data.frame()
  
  # Loop over each unique MESH term
  for (mesh_term in unique(enrichment_data[[mesh_col]])) {
    
    # Subset the top hits for this MESH term
    mesh_data <- top_hits %>% filter(get(mesh_col) == mesh_term)
    
    # Calculate how many rows to add (if there are fewer than 5, add dummy rows)
    rows_to_add <- 5 - nrow(mesh_data)
    
    # Create dummy rows if needed
    if (rows_to_add > 0) {
      # Create dummy rows with unique identifiers (z_dummy_1, z_dummy_2, etc.)
      dummy_rows <- mesh_data[1, ]  # Copy the first row to create dummy rows
      dummy_rows[[logP_col]] <- 0  # Set logP to 0 for dummy rows (invisible)
      dummy_rows[[FDR_col]] <- ">0.20"  # Set FDR to ">0.20" for dummy rows (invisible)
      
      # Create dummy tissue names (z_dummy_1, z_dummy_2, etc.)
      for (i in 1:rows_to_add) {
        dummy_rows[[tissue_col]] <- paste0("z_dummy_", i)
        
        # Repeat the dummy row for the required number of rows
        mesh_data <- rbind(mesh_data, dummy_rows)
      }
    }
    
    # Add this data to the padded_data dataframe
    padded_data <- rbind(padded_data, mesh_data)
  }
  
  # Create the plot using the padded data (which includes dummy rows)
  ggplot(padded_data, aes_string(x = tissue_col, y = logP_col, fill = FDR_col)) +
    geom_bar(stat = "identity", show.legend = FALSE, width = 0.7) +  # Fixed bar width
    scale_fill_manual(values = c("<0.01" = color_significant, "<0.05" = color_less_significant, ">0.20" = color_nonsignificant)) +
    theme_minimal(base_size = 14) +  # Clean, minimal theme with larger font
    theme(
      axis.text.x = element_blank(),  # Remove x-axis labels
      axis.title.x = element_blank(), # Remove x-axis title
      axis.title.y = element_text(face = "bold", size = 14),  # Keep y-axis title
      axis.text.y = element_text(size = 12),
      axis.line = element_line(colour = "black", size = 1),  # Add axis lines
      axis.ticks = element_line(colour = "black"),  # Add axis ticks
      panel.grid = element_blank(),  # Clean grid lines
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      strip.text = element_blank(),  # Remove facet titles (MESH categories)
      strip.background = element_rect(fill = "white", color = "white"),  # Remove background of facet labels
      panel.spacing = unit(1, "lines"),  # Control spacing between facets
      axis.ticks.x = element_blank()  # Remove axis ticks for a cleaner look
    ) +
    labs(
      y = expression(-log[10](p)),  # Y-axis label with math expression
      #title = title  # Customizable plot title
    ) +
    scale_y_continuous(limits = c(0, max_logP)) +  # Dynamic Y-axis range
    facet_wrap(~ get(mesh_col), scales = "free_x", ncol = length(unique(enrichment_data[[mesh_col]]))) +  # Facet by MESH term with free y-scales
    scale_x_discrete(labels = function(x) ifelse(grepl("^z_dummy_", x), "", x))  # Hide dummy labels
}

create_enrichment_plot_grouped_no_facet_no_y_axis_title <- function(enrichment_data, tissue_col = "Name", logP_col = "logP", FDR_col = "FDR", 
                                                                    mesh_col = "MESH", title = "Enrichment Analysis by Tissue", 
                                                                    color_significant = "red", color_less_significant="blue", color_nonsignificant = "gray") {
  
  # Ensure that FDR and MESH are factors for proper coloring and grouping
  enrichment_data[[FDR_col]] <- factor(enrichment_data[[FDR_col]], levels = c("<0.01", "<0.05", ">0.20"))
  enrichment_data[[mesh_col]] <- factor(enrichment_data[[mesh_col]])  # MESH term as a factor
  
  # Get the maximum value of logP for the y-axis limit
  max_logP <- max(enrichment_data[[logP_col]], na.rm = TRUE)
  
  # Get the top 5 hits based on logP per MESH term
  top_hits <- enrichment_data %>%
    group_by(get(mesh_col)) %>%
    top_n(5, !!sym(logP_col)) %>%  # Get top 5 based on logP
    ungroup()
  
  # Create a new data frame to store the final enriched data with padding for each MESH term
  padded_data <- data.frame()
  
  # Loop over each unique MESH term
  for (mesh_term in unique(enrichment_data[[mesh_col]])) {
    
    # Subset the top hits for this MESH term
    mesh_data <- top_hits %>% filter(get(mesh_col) == mesh_term)
    
    # Calculate how many rows to add (if there are fewer than 5, add dummy rows)
    rows_to_add <- 5 - nrow(mesh_data)
    
    # Create dummy rows if needed
    if (rows_to_add > 0) {
      # Create dummy rows with unique identifiers (z_dummy_1, z_dummy_2, etc.)
      dummy_rows <- mesh_data[1, ]  # Copy the first row to create dummy rows
      dummy_rows[[logP_col]] <- 0  # Set logP to 0 for dummy rows (invisible)
      dummy_rows[[FDR_col]] <- ">0.20"  # Set FDR to ">0.20" for dummy rows (invisible)
      
      # Create dummy tissue names (z_dummy_1, z_dummy_2, etc.)
      for (i in 1:rows_to_add) {
        dummy_rows[[tissue_col]] <- paste0("z_dummy_", i)
        
        # Repeat the dummy row for the required number of rows
        mesh_data <- rbind(mesh_data, dummy_rows)
      }
    }
    
    # Add this data to the padded_data dataframe
    padded_data <- rbind(padded_data, mesh_data)
  }
  
  # Create the plot using the padded data (which includes dummy rows)
  ggplot(padded_data, aes_string(x = tissue_col, y = logP_col, fill = FDR_col)) +
    geom_bar(stat = "identity", show.legend = FALSE, width = 0.7) +  # Fixed bar width
    scale_fill_manual(values = c("<0.01" = color_significant, "<0.05" = color_less_significant, ">0.20" = color_nonsignificant)) +
    theme_minimal(base_size = 14) +  # Clean, minimal theme with larger font
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12, face = "bold"),  # Angled labels for tissue names
      axis.title.x = element_text(face = "bold", size = 14),  # Keep x-axis title
      axis.text.y = element_text(size = 12),
      axis.line = element_line(colour = "black", size = 1),  # Add axis lines
      axis.ticks = element_line(colour = "black"),  # Add axis ticks
      panel.grid = element_blank(),  # Clean grid lines
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      strip.text = element_blank(),  # Remove facet titles (MESH categories)
      strip.background = element_rect(fill = "white", color = "white"),  # Remove background of facet labels
      panel.spacing = unit(1, "lines"),  # Control spacing between facets
      axis.ticks.x = element_blank()  # Remove axis ticks for a cleaner look
    ) +
    labs(
      x = "Tissue",  # Label for x-axis
      y = NULL,  # Remove y-axis title
      #title = title  # Customizable plot title
    ) +
    scale_y_continuous(limits = c(0, max_logP)) +  # Dynamic Y-axis range
    facet_wrap(~ get(mesh_col), scales = "free_x", ncol = length(unique(enrichment_data[[mesh_col]]))) +  # Facet by MESH term with free y-scales
    scale_x_discrete(labels = function(x) ifelse(grepl("^z_dummy_", x), "", x))  # Hide dummy labels
}


create_enrichment_plot_grouped_no_x_labels_no_facet_titles_no_y_axis_title <- function(enrichment_data, tissue_col = "Name", logP_col = "logP", FDR_col = "FDR", 
                                                                                       mesh_col = "MESH", title = "Enrichment Analysis by Tissue", 
                                                                                       color_significant = "red", color_less_significant="blue", color_nonsignificant = "gray") {
  
  # Ensure that FDR and MESH are factors for proper coloring and grouping
  enrichment_data[[FDR_col]] <- factor(enrichment_data[[FDR_col]], levels = c("<0.01", "<0.05", ">0.20"))
  enrichment_data[[mesh_col]] <- factor(enrichment_data[[mesh_col]])  # MESH term as a factor
  
  # Get the maximum value of logP for the y-axis limit
  max_logP <- max(enrichment_data[[logP_col]], na.rm = TRUE)
  
  # Get the top 5 hits based on logP per MESH term
  top_hits <- enrichment_data %>%
    group_by(get(mesh_col)) %>%
    top_n(5, !!sym(logP_col)) %>%  # Get top 5 based on logP
    ungroup()
  
  # Create a new data frame to store the final enriched data with padding for each MESH term
  padded_data <- data.frame()
  
  # Loop over each unique MESH term
  for (mesh_term in unique(enrichment_data[[mesh_col]])) {
    
    # Subset the top hits for this MESH term
    mesh_data <- top_hits %>% filter(get(mesh_col) == mesh_term)
    
    # Calculate how many rows to add (if there are fewer than 5, add dummy rows)
    rows_to_add <- 5 - nrow(mesh_data)
    
    # Create dummy rows if needed
    if (rows_to_add > 0) {
      # Create dummy rows with unique identifiers (z_dummy_1, z_dummy_2, etc.)
      dummy_rows <- mesh_data[1, ]  # Copy the first row to create dummy rows
      dummy_rows[[logP_col]] <- 0  # Set logP to 0 for dummy rows (invisible)
      dummy_rows[[FDR_col]] <- ">0.20"  # Set FDR to ">0.20" for dummy rows (invisible)
      
      # Create dummy tissue names (z_dummy_1, z_dummy_2, etc.)
      for (i in 1:rows_to_add) {
        dummy_rows[[tissue_col]] <- paste0("z_dummy_", i)
        
        # Repeat the dummy row for the required number of rows
        mesh_data <- rbind(mesh_data, dummy_rows)
      }
    }
    
    # Add this data to the padded_data dataframe
    padded_data <- rbind(padded_data, mesh_data)
  }
  
  # Create the plot using the padded data (which includes dummy rows)
  ggplot(padded_data, aes_string(x = tissue_col, y = logP_col, fill = FDR_col)) +
    geom_bar(stat = "identity", show.legend = FALSE, width = 0.7) +  # Fixed bar width
    scale_fill_manual(values = c("<0.01" = color_significant, "<0.05" = color_less_significant, ">0.20" = color_nonsignificant)) +
    theme_minimal(base_size = 14) +  # Clean, minimal theme with larger font
    theme(
      axis.text.x = element_blank(),  # Remove x-axis labels
      axis.title.x = element_blank(),  # Remove x-axis title
      axis.text.y = element_text(size = 12),
      axis.line = element_line(colour = "black", size = 1),  # Add axis lines
      axis.ticks = element_line(colour = "black"),  # Add axis ticks
      panel.grid = element_blank(),  # Clean grid lines
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      strip.text = element_blank(),  # Remove facet titles (MESH categories)
      strip.background = element_rect(fill = "white", color = "white"),  # Remove background of facet labels
      panel.spacing = unit(1, "lines"),  # Control spacing between facets
      axis.ticks.x = element_blank()  # Remove axis ticks for a cleaner look
    ) +
    labs(
      y = NULL,  # Remove y-axis title
      #title = title  # Customizable plot title
    ) +
    scale_y_continuous(limits = c(0, max_logP)) +  # Dynamic Y-axis range
    facet_wrap(~ get(mesh_col), scales = "free_x", ncol = length(unique(enrichment_data[[mesh_col]]))) +  # Facet by MESH term with free y-scales
    scale_x_discrete(labels = function(x) ifelse(grepl("^z_dummy_", x), "", x))  # Hide dummy labels
}

##############
#Loading data#
##############

path_4_grs <- "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025//"

setwd(path_4_grs)

###################################################
#STEP 1: Let's load all the tissue enrichment data#
###################################################

bmi_ns <- fread("output/5_enrichment_analyses/2_depict/output/cluster_bmi_ns/cluster_bmi_ns_tissueenrichment.txt")
bmi_neg <- fread("output/5_enrichment_analyses/2_depict/output/cluster_bmi_neg/cluster_bmi_neg_tissueenrichment.txt")
bmi_pos <- fread("output/5_enrichment_analyses/2_depict/output/cluster_bmi_pos/cluster_bmi_pos_tissueenrichment.txt")

#Let's change the name of the dataframe:

bmi_ns <- bmi_ns %>%
  select(Name, `MeSH first level term`, `Nominal P value`, `False discovery rate`)

bmi_neg <- bmi_neg %>%
  select(Name, `MeSH first level term`, `Nominal P value`, `False discovery rate`)

bmi_pos <- bmi_pos %>%
  select(Name, `MeSH first level term`, `Nominal P value`, `False discovery rate`)

#Let's switch the columns:

colnames(bmi_ns) <- c("Name", "MESH", "pvalue", "FDR")
colnames(bmi_neg) <- c("Name", "MESH", "pvalue", "FDR")
colnames(bmi_pos) <- c("Name", "MESH", "pvalue", "FDR")

#Let's do some small changes:

bmi_ns$logP <- -log10(as.numeric(bmi_ns$pvalue))
bmi_neg$logP <- -log10(as.numeric(bmi_neg$pvalue))
bmi_pos$logP <- -log10(as.numeric(bmi_pos$pvalue))

#Let's make a list of MESH terms that are significant for at least one!

yes_vect <- c("<0.01", "<0.05")

terms <- c(bmi_ns$MESH[which(bmi_ns$FDR%in%yes_vect)], bmi_neg$MESH[which(bmi_neg$FDR%in%yes_vect)], 
           bmi_pos$MESH[which(bmi_pos$FDR%in%yes_vect)])

tissues <- c(bmi_ns$Name[which(bmi_ns$FDR%in%yes_vect)], bmi_neg$Name[which(bmi_neg$FDR%in%yes_vect)], 
             bmi_pos$Name[which(bmi_pos$FDR%in%yes_vect)])

terms <- unique(terms)
tissues <- unique(tissues)

#Let's change those that have a different color, don't care about those:

bmi_ns$FDR <- ifelse(bmi_ns$FDR%in%yes_vect, bmi_ns$FDR, ">0.05")
bmi_neg$FDR <- ifelse(bmi_neg$FDR%in%yes_vect, bmi_neg$FDR, ">0.05")
bmi_pos$FDR <- ifelse(bmi_pos$FDR%in%yes_vect, bmi_pos$FDR, ">0.05")

#Let's get those that are in terms:

bmi_ns <- bmi_ns[which(bmi_ns$MESH%in%terms),]
bmi_neg <- bmi_neg[which(bmi_neg$MESH%in%terms),]
bmi_pos <- bmi_pos[which(bmi_pos$MESH%in%terms),]

bmi_ns <- bmi_ns[which(bmi_ns$Name%in%tissues),]
bmi_neg <- bmi_neg[which(bmi_neg$Name%in%tissues),]
bmi_pos <- bmi_pos[which(bmi_pos$Name%in%tissues),]

#Let's limit the data to best hits across all datasets because it is too difficult to compare:

bmi_ns$data <- "BMI (P>0.05)"
bmi_neg$data <- "BMI- (P<0.05)"
bmi_pos$data <- "BMI+ (P<0.05)"

all_data_sets <- rbind(bmi_ns, bmi_neg, bmi_pos)

all_data_sets <- all_data_sets[order(all_data_sets$pvalue),]

for(term in terms){
  
  best_terms <- all_data_sets[which(all_data_sets$MESH == term),]
  
  #let's get the best one
  
  best_terms <- best_terms[which(duplicated(best_terms$Name)==FALSE),]
  
  best_terms <- best_terms[1:5,] #top 5
  
  best_terms <- best_terms[which(is.na(best_terms$Name) == FALSE),]
  
  if(!(exists("final_terms"))){
    
    final_terms <- best_terms
    
  } else {
    
    final_terms <- rbind(final_terms, best_terms)
    
  }
  
}

#Finally let's change the data:

bmi_ns <- bmi_ns[which(bmi_ns$Name%in%final_terms$Name),]
bmi_neg <- bmi_neg[which(bmi_neg$Name%in%final_terms$Name),]
bmi_pos <- bmi_pos[which(bmi_pos$Name%in%final_terms$Name),]

#Let's plot the expected data:

p_1 <- create_enrichment_plot_grouped_no_x_labels(bmi_ns)
p_2 <- create_enrichment_plot_grouped_no_x_labels_no_facet_titles_no_y_axis_title(bmi_neg)
p_3 <- create_enrichment_plot_grouped_no_facet_no_y_axis_title(bmi_pos)

######################################
#Let's combine all the plots together#
######################################

library(patchwork)

p_all <- p_1 + p_2 + p_3 + plot_layout(ncol = 1)

ggsave(
  filename = "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/manuscript/figures/tissue_enrichment_depict_clusters.svg",
  plot = p_all,
  width = 15,
  height = 10,
  units = "in",  # or "in" depending on your preference, but cm is common
  device = "svg"
)
