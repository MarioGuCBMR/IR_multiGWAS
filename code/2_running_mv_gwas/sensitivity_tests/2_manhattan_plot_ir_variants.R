##############
#INTRODUCTION#
##############

#For simplicity in the Manhattan plot we are going to highlight the 282 first and then the novel ones

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)
library(GenomicSEM)
library(vautils)
library(vroom)
library(tidyverse)
library(ggplot2)
library(ggrepel)

##############################
#Let's obtain the LDSC output#
##############################

path_2_input <- "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/"

setwd(path_2_input)

all_sumstats <- fread("output/2_mv_gwas/fiadjbmi_hdl_tg_common_dwls_curated.txt")

#Let's prepare several datasets with specific information so that we can generate an interesting Manhattan plot:

ir_dir <- fread("output/4_ir_loci_discovery/3_novel_variants/3_novel_variants/282_ir_variants_with_ir_and_t2d_comparisons.txt")

ir_new <- ir_dir[which(ir_dir$reported_ir == "no" & str_detect(ir_dir$ir_source, "gw") == FALSE),]

ir_lotta <- ir_dir[which(str_detect(ir_dir$ir_source, "Lotta")),]

ir_oliveri <- ir_dir[which(str_detect(ir_dir$ir_source, "Oliveri")),]

##################################################################
#Let's make the regions here and then we have the closest gene...#
##################################################################

start_ <- ifelse(as.numeric(ir_dir$base_pair_location)-500000<0, 0, as.numeric(ir_dir$base_pair_location)-500000)
end_ <- as.numeric(ir_dir$base_pair_location)+500000

ir_dir$loci <- paste(ir_dir$chromosome, ":", start_, "_", end_, sep ="")

close_gene <- find_nearest_gene(ir_dir, flanking = 100, build = "hg19",
                                   collapse = FALSE, snp = "variant", chr = "chromosome",
                                   bp = "base_pair_location")

for (rsid in ir_dir$variant) {
  # Subset the genes_nearest for the current rsid
  
  rsid_genes_nearest <- close_gene[close_gene$rsid == rsid, ]
  
  if (!is.na(rsid_genes_nearest$GENE[1])) {
  
  rsid_genes_nearest$distance_start <- abs(rsid_genes_nearest$position - rsid_genes_nearest$geneSTART)
  rsid_genes_nearest$distance_stop <- abs(rsid_genes_nearest$position - rsid_genes_nearest$geneSTOP)

  # Find the row with the minimum distance_start or distance_stop
  min_start_row <- rsid_genes_nearest[which.min(rsid_genes_nearest$distance_start), ]
  min_stop_row <- rsid_genes_nearest[which.min(rsid_genes_nearest$distance_stop), ]

  # Determine the closest gene based on the minimum distance
  closest_gene <- ifelse(min_start_row$distance_start <= min_stop_row$distance_stop, min_start_row$GENE, min_stop_row$GENE)
  
  close_gene[close_gene$rsid == rsid, "closest_gene"] <- closest_gene
  
  }

}

#Now let's loop and get it all...

ir_dir$nearest_gene <-ir_dir$variant

for(rsid in ir_dir$variant){
  
  gene_sel <- close_gene[which(close_gene$rsid == rsid),]$closest_gene[1]
  
  ir_dir$nearest_gene[which(ir_dir$variant == rsid)] <- gene_sel
  
}

##############################
#Let's do some quick cleaning#
##############################

#Let's see how many are genome-wide significant...

fiadjbmi_hdl_tg_gw <- all_sumstats[which(all_sumstats$p_value < 0.05),]

#And which of them are heterogeneic and which of them are truly shared:

significant_variants <- fiadjbmi_hdl_tg_gw[which(fiadjbmi_hdl_tg_gw$variant %in% ir_dir$variant),] #31

#And now preparing some dataframe here too:

variants_sel <- fiadjbmi_hdl_tg_gw[, c("variant", "chromosome", "base_pair_location", "p_value")]
colnames(variants_sel) <- c("SNP", "CHR", "BP", "P")
snpsOfInterest_list <- ir_dir$variant
snpsOfInterest_list <- snpsOfInterest_list[which(duplicated(snpsOfInterest_list) == FALSE)]

############################
#Prepare GW dataframe -don-#
############################

# Prepare the dataset
don <- variants_sel %>% 
  
  # Compute chromosome size
  group_by(CHR) %>% 
  summarise(chr_len=max(BP)) %>% 
  
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-as.numeric(chr_len)) %>%
  select(-chr_len) %>%
  
  # Add this info to the initial dataset
  left_join(variants_sel, ., by=c("CHR"="CHR")) %>%
  
  # Add a cumulative position of each SNP
  arrange(CHR, BP) %>%
  mutate(BPcum=BP+tot) %>%
  
  # Add highlight and annotation information
  mutate(is_highlight=ifelse(SNP %in% snpsOfInterest_list, "yes", "no")) %>%
  mutate(is_annotate=ifelse(P<1e-6, "yes", "no")) 

# Prepare X axis
axisdf <- don %>% group_by(CHR) %>% summarize(center=(max(BPcum) + min(BPcum)) / 2 )

####################################################
#Let's make manhattan plots for the common variants#
####################################################

# Subset don for ggrepel
genes_nearest_sel <- ir_dir[, c("variant", "nearest_gene")]
colnames(genes_nearest_sel) <- c("SNP", "closest_gene")
new_don <- merge(don, genes_nearest_sel, by = "SNP", all.x = TRUE)

# Subset don for ggrepel for Lotta!
new_don_subset <- new_don[which(new_don$SNP%in%ir_lotta$variant), ]
# Let's split new_don with Lotta signals

# new_don_lotta <- new_don[which(new_don$SNP%in%ir_lotta$variant),]
# new_don_oliveri <- new_don[which(new_don$SNP%in%ir_oliveri$variant),]
# new_don_novel <- new_don[which(new_don$SNP%in%ir_new$variant),]

# Define genome-wide significance threshold
gw_sig <- -log10(0.005)

# Optional: CUD color scheme for chromosomes
cud_colors_chr <- rep(c("#999999", "#0072B2"), 22)  # grey and blue alternating

# Add a logical column indicating whether the variant is in ir_lotta
new_don$in_lotta <- new_don$SNP %in% ir_lotta$variant

# Build the plot
p_1 <- ggplot(new_don, aes(x = BPcum, y = -log10(P))) +
  
  # Main point layer with conditional alpha for in_lotta
  geom_point(
    aes(color = as.factor(CHR), alpha = in_lotta), 
    size = 1.2
  ) +
  scale_color_manual(values = cud_colors_chr) +
  #scale_alpha_manual(values = c("TRUE" = 0.9, "FALSE" = 0.01)) +
  
  # Genome-wide significant hits, colored by whether gene is known
  geom_point(
    data = subset(new_don, is_highlight == "yes" & isTRUE(in_lotta) & !is.na(closest_gene)),
    color = "#D55E00", size = 2.5
  ) +
  geom_point(
    data = subset(new_don, is_highlight == "yes" & isTRUE(in_lotta) & is.na(closest_gene)),
    color = "black", size = 2.5
  ) +
  
  # Add labels (gene or SNP name)
  geom_label_repel(
    data = new_don_subset,
    aes(label = ifelse(is.na(closest_gene), SNP, closest_gene)),
    size = 5,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.2,
    segment.size = 0.3,
    min.segment.length = 0,
    segment.color = "gray40",
    nudge_y = 1.5
  ) +
  
  # Customize axes
  scale_x_continuous(
    label = axisdf$CHR,
    breaks = axisdf$center,
    expand = c(0.01, 0)
  ) +
  scale_y_continuous(
    limits = c(gw_sig, 35),
    breaks = seq(5, 35, 5),
    expand = c(0, 0)
  ) +
  
  # Add horizontal line for genome-wide significance
  geom_hline(yintercept = 6, linetype = "dashed", color = "gray40", linewidth = 0.4) +
  
  # Labels and theme
  labs(
    x = "Chromosome",
    y = expression(-log[10](italic(p))),
    title = "53 IR variants reported by Lotta et al"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
    axis.text.y = element_text(size = 13),
    axis.title = element_text(size = 15),
    plot.margin = margin(10, 10, 10, 10),
    panel.grid.major.y = element_line(color = "gray90", size = 0.2),
    panel.grid.minor = element_blank()
  )


#############################################################
#Let's do the same for novel, so that we can highlight those#
#############################################################

#Highlight the vairants:

# Subset don for ggrepel for Lotta!
new_don_subset <- new_don[which(new_don$SNP%in%ir_new$variant), ]

#And add in the original dataframe a binary variable saying wihch ones you want
new_don$in_novel <- new_don$SNP %in% ir_new$variant

# Build the plot
p_2 <- ggplot(new_don, aes(x = BPcum, y = -log10(P))) +
  
  # Main point layer with conditional alpha for in_lotta
  geom_point(
    aes(color = as.factor(CHR), alpha = in_novel), 
    size = 1.2
  ) +
  scale_color_manual(values = cud_colors_chr) +
  #scale_alpha_manual(values = c("TRUE" = 0.9, "FALSE" = 0.01)) +
  
  # Genome-wide significant hits, colored by whether gene is known
  geom_point(
    data = subset(new_don, is_highlight == "yes" & isTRUE(in_novel) & !is.na(closest_gene)),
    color = "#D55E00", size = 2.5
  ) +
  geom_point(
    data = subset(new_don, is_highlight == "yes" & isTRUE(in_novel) & is.na(closest_gene)),
    color = "black", size = 2.5
  ) +
  
  # Add labels (gene or SNP name)
  geom_label_repel(
    data = new_don_subset,
    aes(label = ifelse(is.na(closest_gene), SNP, closest_gene)),
    size = 5,
    max.overlaps = Inf,
    box.padding = 0.4,
    point.padding = 0.2,
    segment.size = 0.3,
    min.segment.length = 0,
    segment.color = "gray40",
    nudge_y = 1.5
  ) +
  
  # Customize axes
  scale_x_continuous(
    label = axisdf$CHR,
    breaks = axisdf$center,
    expand = c(0.01, 0)
  ) +
  scale_y_continuous(
    limits = c(gw_sig, 15),
    breaks = seq(5, 15, 5),
    expand = c(0, 0)
  ) +
  
  # Add horizontal line for genome-wide significance
  geom_hline(yintercept = 6, linetype = "dashed", color = "gray40", linewidth = 0.4) +
  
  # Labels and theme
  labs(
    x = "Chromosome",
    y = expression(-log[10](italic(p))),
    title = "70 novel IR variants"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 13),
    axis.text.y = element_text(size = 13),
    axis.title = element_text(size = 15),
    plot.margin = margin(10, 10, 10, 10),
    panel.grid.major.y = element_line(color = "gray90", size = 0.2),
    panel.grid.minor = element_blank()
  )


###################################################
#Let's write the ir_dir data with the closest gene#
###################################################

library(patchwork)

p_all <- p_1 + p_2 + plot_layout(ncol = 1)

ggsave(
  filename = "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/manuscript/figures/manhattan_plots_lotta_vs_novel.svg",
  plot = p_all,
  width = 25,
  height = 20,
  units = "in",  # or "in" depending on your preference, but cm is common
  device = "svg"
)

ggsave(
  filename = "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/manuscript/figures/manhattan_plots_lotta_vs_novel.tif",
  plot = p_all,
  width = 25,
  height = 20,
  units = "in",  # or "in" depending on your preference, but cm is common
  device = "tif"
)
