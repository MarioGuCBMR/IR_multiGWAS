##############
#INTRODUCTION#
##############

#This code reads 15-state chromatin data from all cell-types and tissues in ROADMAP and prepares so that it can be run in go shifter:

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)
library(ggrepel)

###################
#Loading functions#
###################

parse_and_update <- function(list_file_, tissue_df){
  
  for(file_ in list_file_){
  
    print(file_)
    
    #STEP 1: read the file:
    
    e_tmp <- data.table::fread(paste("output/5_enrichment_analyses/3_go_shifter/output/bmi_ns/roadmap_systematic/", file_, sep = ""))
    
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

#####################
#Let's read the data#
#####################

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/")

list_of_files <- list.files("output/5_enrichment_analyses/3_go_shifter/output/bmi_ns/roadmap_systematic/")
list_of_files <- list_of_files[which(str_detect(list_of_files, ".enrich"))]

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

roadmap_update <- parse_and_update(list_of_files, roadmap_df)

#######################################################
#Let's add the SGBS cells too, to see what is going on#
#######################################################

adipose <- fread("output/5_enrichment_analyses/3_go_shifter/output/bmi_ns/metsim/adipose.nperm10000.enrich")
sbgs_day0 <- fread("output/5_enrichment_analyses/3_go_shifter/output/bmi_ns/metsim/sgbs_day0.nperm10000.enrich")
sbgs_day4 <- fread("output/5_enrichment_analyses/3_go_shifter/output/bmi_ns/metsim/sgbs_day4.nperm10000.enrich")
sbgs_day14 <- fread("output/5_enrichment_analyses/3_go_shifter/output/bmi_ns/metsim/sgbs_day14.nperm10000.enrich")

#Let's compute the data for adipose and add it:

adipose_overlap <- adipose$nSnpOverlap[1]
adipose_enrich = adipose$enrichment[1]

adipose_test <- adipose[2:length(adipose$nperm),]

adipose_pvalue= as.numeric(length(which(adipose_test$enrichment >= adipose_enrich)))/10000

adipose_df <- as.data.frame(t(c("METSIM_001", "adipose_tissue", adipose_pvalue, adipose_overlap)))
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

##################################
#Let's add -logP and make a graph#
##################################

roadmap_update$logP <- -log10(as.numeric(roadmap_update$enrichment_pval))

roadmap_update$logP[is.infinite(roadmap_update$logP)] <- 7
roadmap_update$overlap <- as.numeric(roadmap_update$overlap)

threshold_0.05 <- -log10(0.05)
threshold_0.01 <- -log10(1e-04)

# Identify points to highlight
roadmap_update$highlight <- ifelse(
  roadmap_update$logP == 7 | 
    str_detect(roadmap_update$sample, "(?i)adipo") |
    str_detect(roadmap_update$sample, "(?i)SGBS") |
    as.numeric(roadmap_update$enrichment_pval) <= 1e-04 |
    roadmap_update$overlap >= 170, 
  TRUE, 
  FALSE
)

#Let's save the data first...#

fwrite(roadmap_update, "output/5_enrichment_analyses/3_go_shifter/output/enrichment_bmi_ns_df.txt")

#########################
#Let's plot the results!#
#########################

plot <- ggplot(roadmap_update, aes(x = overlap, y = logP)) +
  # Significance threshold lines
  geom_hline(yintercept = threshold_0.05, linetype = "dashed", color = "red", size = 0.8) +
  geom_hline(yintercept = threshold_0.01, linetype = "dashed", color = "darkred", size = 0.8) +
  
  # Points with different colors for highlighted ones
  geom_point(aes(color = highlight), size = 3.5, alpha = 0.85) +
  scale_color_manual(values = c("FALSE" = "#505050", "TRUE" = "#E69F00")) +  # Grey & Orange
  
  # Labels for highlighted points
  geom_text_repel(
    aes(label = ifelse(highlight, sample, '')),
    size = 4, 
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

tiff("output/5_enrichment_analyses/3_go_shifter/output/enrichment_plot_bmi_ns.tiff", width = 4000, height = 3500, res=300)
plot
dev.off()
