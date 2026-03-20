##############
#INTRODUCTION#
##############

#This code reads 15-state chromatin data from all cell-types and tissues in ROADMAP and prepares so that it can be run in go shifter.
#The code makes functions to load data from different clusters more easily.

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)
library(ggrepel)

###################
#Loading functions#
###################

parse_and_update <- function(list_file_, tissue_df, path_){
  
  for(file_ in list_file_){
  
    print(file_)
    
    #STEP 1: read the file:
    
    e_tmp <- data.table::fread(paste(path_, file_, sep = ""))
    
    #STEP 2: calculate the enrichment:
    
    original_overlap <- e_tmp$nSnpOverlap[1]
    original_enrich = e_tmp$enrichment[1]
    
    e_test <- e_tmp[2:length(e_tmp$nperm),]
    
    pvalue= as.numeric(length(which(e_test$enrichment >= original_enrich)))/10000
    
    #STEP 2: update the tissue_df:
    
    id_ <- unlist(str_split(file_, "_"))[1]
    
    tissue_df$enrichment_pval[which(tissue_df$ID == id_)] <- pvalue
    tissue_df$overlap[which(tissue_df$ID == id_)] <- original_overlap
    
    
  }
  
  #STEP 3: return the data:
  
  return(tissue_df)

}

updating_df_w_adipose <- function(roadmap_update, path_){
  #A file that load the results for METSIM data
  #Adds it to the roadmap_update results
  #and returns it!
  
  
  adipose <- fread(paste(path_, "/adipose.nperm10000.enrich", sep = ""))
  sbgs_day0 <-fread(paste(path_, "/sgbs_day0.nperm10000.enrich", sep= ""))
  sbgs_day4 <- fread(paste(path_, "/sgbs_day4.nperm10000.enrich", sep = ""))
  sbgs_day14 <- fread(paste(path_, "/sgbs_day14.nperm10000.enrich", sep=""))
  
  #Let's compute the data for adipose and add it:
  
  adipose_overlap <- adipose$nSnpOverlap[1]
  adipose_enrich = adipose$enrichment[1]
  
  adipose_test <- adipose[2:length(adipose$nperm),]
  
  adipose_pvalue= as.numeric(length(which(adipose_test$enrichment >= adipose_enrich)))/10000
  
  adipose_df <- as.data.frame(t(c("METSIM_001", "Adipose tissue consensus", adipose_pvalue, adipose_overlap)))
  colnames(adipose_df) <- colnames(roadmap_update)
  
  roadmap_update <- rbind(roadmap_update, adipose_df)
  
  #Let's compute the data for day 0 and add it:
  
  sbgs_day0_overlap <- sbgs_day0$nSnpOverlap[1]
  sbgs_day0_enrich = sbgs_day0$enrichment[1]
  
  sbgs_day0_test <- sbgs_day0[2:length(sbgs_day0$nperm),]
  
  sbgs_day0_pvalue= as.numeric(length(which(sbgs_day0_test$enrichment >= sbgs_day0_enrich)))/10000
  
  sbgs_day0_df <- as.data.frame(t(c("METSIM_002", "SGBS day 0", sbgs_day0_pvalue, sbgs_day0_overlap)))
  colnames(sbgs_day0_df) <- colnames(roadmap_update)
  
  roadmap_update <- rbind(roadmap_update, sbgs_day0_df)
  
  #Let's compute the data for adipose and add it:
  
  sbgs_day4_overlap <- sbgs_day4$nSnpOverlap[1]
  sbgs_day4_enrich = sbgs_day4$enrichment[1]
  
  sbgs_day4_test <- sbgs_day4[2:length(sbgs_day4$nperm),]
  
  sbgs_day4_pvalue= as.numeric(length(which(sbgs_day4_test$enrichment >= sbgs_day4_enrich)))/10000
  
  sbgs_day4_df <- as.data.frame(t(c("METSIM_003", "SGBS day 4", sbgs_day4_pvalue, sbgs_day4_overlap)))
  colnames(sbgs_day4_df) <- colnames(roadmap_update)
  
  roadmap_update <- rbind(roadmap_update, sbgs_day4_df)
  
  #Let's compute the data for adipose and add it:
  
  sbgs_day14_overlap <- sbgs_day14$nSnpOverlap[1]
  sbgs_day14_enrich = sbgs_day14$enrichment[1]
  
  sbgs_day14_test <- sbgs_day14[2:length(sbgs_day14$nperm),]
  
  sbgs_day14_pvalue= as.numeric(length(which(sbgs_day14_test$enrichment >= sbgs_day14_enrich)))/10000
  
  sbgs_day14_df <- as.data.frame(t(c("METSIM_004", "SGBS day 14", sbgs_day14_pvalue, sbgs_day14_overlap)))
  colnames(sbgs_day14_df) <- colnames(roadmap_update)
  
  roadmap_update <- rbind(roadmap_update, sbgs_day14_df)
  
  #Finally, let's return the data:
  
  return(roadmap_update)
  
}

formatting_4_graph <- function(roadmap_update){
  #A function that formats the data so that it is ready to be printed
  
  #First let's add the logP
  
  roadmap_update$logP <- -log10(as.numeric(roadmap_update$enrichment_pval))
  
  roadmap_update$logP[is.infinite(roadmap_update$logP)] <- 7
  roadmap_update$overlap <- as.numeric(roadmap_update$overlap)
  
  #Add threshold settings for the graph
  
  threshold_0.05 <- -log10(0.05)
  threshold_0.01 <- -log10(1e-04)
  
  #We identify points to highlight -significant or adipose-related-
  
  roadmap_update$highlight <- ifelse(
    roadmap_update$logP == 7 | 
      str_detect(roadmap_update$sample, "(?i)adipo") |
      str_detect(roadmap_update$sample, "(?i)SGBS") |
      as.numeric(roadmap_update$enrichment_pval) <= 1e-04,
    #roadmap_update$overlap >= 170, 
    TRUE, 
    FALSE
  )
  
  #And finally return:
  
  return(roadmap_update)
  
  
}

plotting_data <- function(roadmap_update){
  #This function makes a plot in poppish styles.
  
  # Set significance threshold lines
  
  threshold_0.05 <- -log10(0.05) #nominal
  threshold_0.01 <- -log10(1e-04) #permutation FDR
  
  #Let's make the plot!!
  
  plotio <- ggplot(roadmap_update, aes(x = overlap, y = logP)) +

    # Set significance threshold lines
    
    geom_hline(yintercept = threshold_0.05, linetype = "dashed", color = "red", size = 0.8) +
    geom_hline(yintercept = threshold_0.01, linetype = "dashed", color = "darkred", size = 0.8) +
    
    # Points with different colors for highlighted ones
    geom_point(aes(color = highlight), size = 1.75, alpha = 0.55) +
    scale_color_manual(values = c("FALSE" = "#505050", "TRUE" = "#E69F00")) +  # Grey & Orange
    
    # Labels for highlighted points
    geom_text_repel(
      aes(label = ifelse(highlight, sample, '')),
      size = 9, 
      box.padding = 0.7, 
      point.padding = 0.5, 
      max.overlaps = Inf,
      family = "Times New Roman"
    ) +
    
    # Labels and theme
    labs(
      x = "Overlap", 
      y = expression(-log[10](P)),  # Proper math notation for logP
      title = ""
    ) +
    
    # Nature-style theme
    theme_classic(base_size = 18, base_family = "Times New Roman") +  
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
      axis.title = element_text(face = "bold", size = 18),
      axis.text = element_text(size = 16),
      axis.line = element_line(size = 1.2),  # Thicker axis lines
      axis.ticks = element_line(size = 1.2),  # Thicker ticks
      axis.ticks.length = unit(0.3, "cm"),  # Inward ticks
      panel.grid.major = element_line(color = "grey80", size = 0.4),  # Subtle major grid
      panel.grid.minor = element_blank(),  # No minor grid
      legend.position = "none"  # Remove legend for clean look
    )
  
  #And return it:
  
  return(plotio)
  
}

#####################
#Let's read the data#
#####################

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/")

bmi_ns <- list.files("output/5_enrichment_analyses/3_go_shifter/output/bmi_ns/roadmap_systematic/")
bmi_ns <- bmi_ns[which(str_detect(bmi_ns, ".enrich"))]

bmi_neg <- list.files("output/5_enrichment_analyses/3_go_shifter/output/bmi_neg/roadmap_systematic/")
bmi_neg <- bmi_neg[which(str_detect(bmi_neg, ".enrich"))]

bmi_pos <- list.files("output/5_enrichment_analyses/3_go_shifter/output/bmi_pos/roadmap_systematic/")
bmi_pos <- bmi_pos[which(str_detect(bmi_pos, ".enrich"))]

roadmap_df <- fread("raw_data/atac_seq_data/nar-00860-z-2017-File021.txt")

###############################################################################
#Let's update the tissue data with the results!! Let's clean a bit and do this#
###############################################################################

roadmap_df <- roadmap_df %>%
  dplyr::select(EpigenomeID, StdEpigenomeName)

colnames(roadmap_df) <- c("ID", "sample")

roadmap_df$enrichment_pval <- NA
roadmap_df$overlap <- NA

#And update:

bmi_ns_update <- parse_and_update(bmi_ns, roadmap_df, path_ = "output/5_enrichment_analyses/3_go_shifter/output/bmi_ns/roadmap_systematic/")
bmi_neg_update <- parse_and_update(bmi_neg, roadmap_df, path_ = "output/5_enrichment_analyses/3_go_shifter/output/bmi_neg/roadmap_systematic/")
bmi_pos_update <- parse_and_update(bmi_pos, roadmap_df, path_ = "output/5_enrichment_analyses/3_go_shifter/output/bmi_pos/roadmap_systematic/")

#######################################################
#Let's add the SGBS cells too, to see what is going on#
#######################################################

bmi_ns_updated <- updating_df_w_adipose(bmi_ns_update, "output/5_enrichment_analyses/3_go_shifter/output/bmi_ns/metsim/")
bmi_neg_updated <- updating_df_w_adipose(bmi_neg_update, "output/5_enrichment_analyses/3_go_shifter/output/bmi_neg/metsim/")
bmi_pos_updated <- updating_df_w_adipose(bmi_pos_update, "output/5_enrichment_analyses/3_go_shifter/output/bmi_pos/metsim/")

#################################################################
#Let's add -logP and add the necessary data to compute the graph#
#################################################################

bmi_ns_updated <- formatting_4_graph(bmi_ns_updated)
bmi_neg_updated <- formatting_4_graph(bmi_neg_updated)
bmi_pos_updated <- formatting_4_graph(bmi_pos_updated)

#Let's save the data first..., before running stuff...

fwrite(bmi_ns_updated, "output/5_enrichment_analyses/3_go_shifter/output/enrichment_bmi_ns_df.txt", sep = "\t")
fwrite(bmi_neg_updated, "output/5_enrichment_analyses/3_go_shifter/output/enrichment_bmi_neg_df.txt", sep = "\t")
fwrite(bmi_pos_updated, "output/5_enrichment_analyses/3_go_shifter/output/enrichment_bmi_pos_df.txt", sep = "\t")

#########################
#Let's plot the results!#
#########################

bmi_ns_plot <- plotting_data(bmi_ns_updated)
bmi_neg_plot <- plotting_data(bmi_neg_updated)
bmi_pos_plot <- plotting_data(bmi_pos_updated)

#Let's combine the data:

library(patchwork)

p_all <- bmi_ns_plot + bmi_neg_plot + bmi_pos_plot + plot_layout(ncol = 3)

ggsave(
  filename = "manuscript/figures/go_shifter_enrichment_4_clusters_plot.svg",
  plot = p_all,
  width = 35,
  height = 12,
  units = "in",  # or "in" depending on your preference, but cm is common
  device = "svg"
)
