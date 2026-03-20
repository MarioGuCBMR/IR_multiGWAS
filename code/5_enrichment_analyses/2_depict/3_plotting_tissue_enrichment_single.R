##############
#INTRODUCTION#
##############

#This makes a figure for DEPICT results.

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
     # theme(
     #   axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12, face = "bold"),  # Angled labels for tissue names
     #   axis.title = element_text(face = "bold", size = 14),
     #   axis.text.y = element_text(size = 12),
     #   axis.line = element_line(colour = "black", size = 1),  # Add axis lines
     #   axis.ticks = element_line(colour = "black"),  # Add axis ticks
     #   panel.grid = element_blank(),  # Clean grid lines
     #   plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
     #   strip.text = element_text(size = 14, face = "bold"),  # Bold MESH labels
     #   strip.background = element_rect(fill = "white", color = "white"),  # Remove background of facet labels
     #   panel.spacing = unit(1, "lines"),  # Control spacing between facets
     #   axis.ticks.x = element_blank()  # Remove axis ticks for a cleaner look
     # ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 14, face = "bold"),
      axis.text.y = element_text(size = 14),
      axis.title = element_text(face = "bold", size = 16),
      axis.line = element_line(colour = "black", size = 1),
      axis.ticks = element_line(colour = "black"),
      panel.grid = element_blank(),
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      strip.text = element_text(size = 16, face = "bold"),
      strip.background = element_rect(fill = "white", color = "white"),
      panel.spacing = unit(1, "lines"),
      axis.ticks.x = element_blank()
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



##############
#Loading data#
##############

path_4_grs <- "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/"

setwd(path_4_grs)

###################################################
#STEP 1: Let's load all the tissue enrichment data#
###################################################

ir_res <- fread("output/5_enrichment_analyses/2_depict/output/cluster_all/cluster_all_tissueenrichment.txt")

#Let's change the name of the dataframe:

ir_res <- ir_res %>%
  select(Name, `MeSH first level term`, `Nominal P value`, `False discovery rate`)

#Let's switch the columns:

colnames(ir_res) <- c("Name", "MESH", "pvalue", "FDR")

#Let's do some small changes:

ir_res$logP <- -log10(as.numeric(ir_res$pvalue))

#Let's make a list of MESH terms that are significant for at least one! 
#The code was done to compare several tissues, hence why it looks a bit weird, but does the trick

yes_vect <- c("<0.01", "<0.05")

terms <- ir_res$MESH[which(ir_res$FDR%in%yes_vect)]

tissues <- ir_res$Name[which(ir_res$FDR%in%yes_vect)]

terms <- unique(terms)
tissues <- unique(tissues)

#Let's change those that have a different color, don't care about those:

ir_res$FDR <- ifelse(ir_res$FDR%in%yes_vect, ir_res$FDR, ">0.05")

#Let's get those that are in terms:

ir_res <- ir_res[which(ir_res$MESH%in%terms),]

#Let's limit the data to best hits across all datasets because it is too difficult to compare:

ir_res$data <- "282 IR variants"

all_data_sets <- ir_res

all_data_sets <- all_data_sets[order(all_data_sets$pvalue),]

#Let's plot the ir_res data:

p_1 <- create_enrichment_plot_grouped(all_data_sets)

#Let's save the data:

ggsave(
  filename = "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/manuscript/figures/tissue_enrichment_depict_all_variants.svg",
  plot = p_1,
  width = 20,
  height = 8,
  units = "in",  # or "in" depending on your preference, but cm is common
  device = "svg"
)
