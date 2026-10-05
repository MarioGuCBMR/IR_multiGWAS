##############
#INTRODUCTION#
##############

#This code runs an R script that sets up the data for to replicate the  
#Lotta et al. 2017 Fig 3d Replication: Matched-Locus Binomial Enrichment
#Most 


# Method (Lotta et al. Online Methods):
#   - 1M permutations to estimate expected overlap proportion
#   - Binomial test per epigenome
#   - Annotations: ROADMAP 15-state model, 98 H3K27ac epigenomes
#   - Active enhancers: 6_EnhG + 7_Enh

###################
#Loading libraries#
###################

library(data.table)
library(GenomicRanges)
library(rtracklayer)
library(ggplot2)
library(svglite)
library(tidyverse)
library(ggrepel)

###################
#Loading functions#
###################

plot_unique=function(cluster_df, title_added){
  
  #STEP 0: make dummy df:
  
  #cluster_df=all_variants
  
  #STEP 1: get the -log10P to make lines
  
  cluster_df$nominal_lp = -log10(as.numeric(cluster_df$binomial_p))
  cluster_df$fdr_lp = -log10(as.numeric(cluster_df$bh_fdr))
  
  nominal_thr = -log10(0.05)
  bonferroni_thr = -log10(0.05/104)
  
  #STEP 2: add the most relevant tissues here: top 5 Bonferroni-passing + always 3 reference tissues ---
  #bonf_labels = c("Adipose Nuclei", "Liver", "Skeletal Muscle Female")  # Adipose Nuclei, Liver, Skeletal Muscle
  
  #And get the top 5:
  
  cluster_df = cluster_df[order(as.numeric(cluster_df$bh_fdr)),]
  check_df = cluster_df[1:10,]
  check_df = check_df[which(as.numeric(check_df$bh_fdr) < 0.05),]
  
  ref_labels = check_df$tissue
  
  #label_data = unique(c(bonf_labels, ref_labels)) #merging and removing duplicates
  #ref_labels = unique(c(bonf_labels, ref_labels)) #merging and removing duplicates
  label_data = ref_labels
  
  cluster_df$label_data = ifelse(cluster_df$tissue%in%label_data, TRUE, FALSE)
  cluster_df$tissue_clean = ifelse(cluster_df$label_data==TRUE, cluster_df$tissue, "")
  
  
  # --- Build plot ---
  p = ggplot(cluster_df, aes(x = observed_prop, y = nominal_lp)) +
    
    # Base points
    geom_point(aes(color = "BH > 0.05"), alpha = 0.4, size = 1.4) +
    
    # Significant tissues
    geom_point(data = subset(cluster_df, bh_fdr < 0.05), aes(color = "BH < 0.05"), alpha = 0.7, size = 2.4) +
    
    scale_color_manual(
      name = "BH-adjusted p-value",
      values = c(
        "BH > 0.05" = "grey60",
        "BH < 0.05" = "#0279EE"
      )
    ) +
    
    # Labels
    geom_text_repel(
      data = cluster_df,
      aes(label = tissue_clean,
          fontface = ifelse(label_data, "bold", "plain")),
      size = 4,
      max.overlaps = 30,
      force = 5,
      segment.size = 0.25,
      segment.alpha = 0.4,
      box.padding = 0.4,
      min.segment.length = 0.05,
      show.legend = FALSE
    ) +
    
    # Significance thresholds
    geom_hline(yintercept = nominal_thr, linetype = "dashed", color = "orange", linewidth = 0.5) +
    geom_hline(yintercept = bonferroni_thr, linetype = "dotted", color = "red", linewidth = 1.0) +
    
    annotate(
      "text",
      x = 0.01,
      y = nominal_thr + 0.15,
      label = "Nominal (p=0.05)",
      size = 4.5,
      hjust = 0,
      color = "black"
    ) +
    annotate(
      "text",
      x = 0.01,
      y = bonferroni_thr + 0.15,
      label = "Bonferroni (0.05/104)",
      size = 4.5,
      hjust = 0,
      color = "black"
    ) +
    
    # Axes
    labs(
      title = title_added,
      x = "Observed proportion (overlap / nSNPs)",
      y = "-log10P"
    ) +
    
    theme_minimal(base_family = "Liberation Sans") +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust=0.5),
      strip.text = element_text(size = 10, face = "bold"),
      axis.text = element_text(size = 12),
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      panel.grid.minor = element_blank(),
      legend.title = element_text(size = 12, face = "bold"),
      legend.position = "right",
      legend.text = element_text(size = 11)
    ) +
    
    scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2)
    ) +
    scale_y_continuous(
      limits = c(0, 12),
      breaks = seq(0, 12, by = 2)
    )
  
}

plot_subclusters=function(cluster_df, title_){
  
  threshold_0.05 = -log10(0.05)
  threshold_0.01 = -log10(0.05/102)
  cluster_df = cluster_df[order(as.numeric(cluster_df$bonferroni), decreasing = FALSE),]
  cluster_df$highlight = FALSE
  cluster_df$highlight[1:5] = TRUE
  #adipose_vect=c("Adipose Nuclei", "Adipose (Perrin et al)", "SGBS day 0 (Perrin et al)",  "SGBS day 4 (Perrin et al)", "SGBS day 14 (Perrin et al)")
  #cluster_df$highlight = ifelse(cluster_df$tissue%in%adipose_vect, TRUE, cluster_df$highlight)
  
  cluster_df$tissue[which(cluster_df$tissue == "Adipose (Perrin et al)")] = "Adipose Tissue"
  cluster_df$tissue[which(cluster_df$tissue == "SGBS day 0 (Perrin et al)")] = "SGBS day 0"
  cluster_df$tissue[which(cluster_df$tissue == "SGBS day 4 (Perrin et al)")] = "SGBS day 4"
  cluster_df$tissue[which(cluster_df$tissue == "SGBS day 14 (Perrin et al)")] = "SGBS day 14"
  cluster_df$tissue[which(cluster_df$tissue == "Bone Marrow Derived Cultured Mesenchymal Stem Cells")] = "BM Mesenchymal Stem Cells"
  cluster_df$tissue[which(cluster_df$tissue == "Mesenchymal Stem Cell Derived Chondrocyte Cultured Cells")] = "MSC-derived Chondrocytes"
  cluster_df$tissue[which(cluster_df$tissue == "Foreskin Fibroblast Primary Cells skin02")] = "Foreskin Fibroblasts"
  cluster_df$tissue[which(cluster_df$tissue == "Monocytes-CD14+ RO01746 Primary Cells")] = "Monocytes-CD14+"
  
  
  plot_all <- ggplot(cluster_df, aes(x = observed_prop*100, y = -log10(binomial_p))) +
    # Significance threshold lines
    geom_hline(yintercept = threshold_0.05, linetype = "dashed", color = "red", size = 0.8) +
    geom_hline(yintercept = threshold_0.01, linetype = "dashed", color = "darkred", size = 0.8) +
    
    # Points with different colors for highlighted ones
    geom_point(aes(color = highlight), size = 1.75, alpha = 0.55) +
    scale_color_manual(values = c("FALSE" = "#505050", "TRUE" = "#E69F00")) +  # Grey & Orange
    
    # Labels for highlighted points
    geom_text_repel(
      aes(label = ifelse(highlight, tissue, '')),
      size = 6, 
      box.padding = 0.7, 
      point.padding = 0.5, 
      max.overlaps = Inf,
      family = "Arial", 
      force=3
    ) +
    
    # Labels and theme
    labs(
      x = "Overlap %", 
      y = expression(-log[10](P)),  # Proper math notation for logP
      title = title_
    ) +
    
    #Scales 
    
    scale_x_continuous(limits = c(0, 100),
                       breaks = seq(0, 100, by = 10)) +
    
    scale_y_continuous( limits = c(0, 10),
                        breaks = seq(0, 10, by = 2)) +
    
    # Nature-style theme
    theme_classic(base_size = 18, base_family = "Arial") +  
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
      axis.title = element_text(face = "bold", size = 18),
      axis.text = element_text(size = 16),
      axis.line = element_line(size = 1.2),  # Thicker axis lines
      axis.ticks = element_line(size = 1.2),  # Thicker ticks
      axis.ticks.length = unit(0.3, "cm"),  # Inward ticks
      panel.grid.major = element_blank(),  # No major grid
      panel.grid.minor = element_blank(),  # No minor grid
      legend.position = "none"  # Remove legend for clean look
    )
  
  
  
}

#######################
#Setting up parameters#
#######################

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025//")

#First let's run the SNPs:

all_variants = fread("review/review_output/5_enrichment_analyses/2_binomial_lotta/output/all_282_ir/results_binomial.tsv")
all_variants[, bh_fdr := p.adjust(binomial_p, method = "BH")]
all_variants[, bonferroni := p.adjust(binomial_p, method = "bonferroni")]

#Now the lipodystrophy ones:

lotta = fread("review/review_output/5_enrichment_analyses/2_binomial_lotta/output/lotta_test/results_binomial.tsv")
lotta[, bh_fdr := p.adjust(binomial_p, method = "BH")]
lotta[, bonferroni := p.adjust(binomial_p, method = "bonferroni")]

#Now the lipodystrophy ones:

bmi_neutral = fread("review/review_output/5_enrichment_analyses/2_binomial_lotta/output/ns_141_ir/results_binomial.tsv")
bmi_neutral[, bh_fdr := p.adjust(binomial_p, method = "BH")]
bmi_neutral[, bonferroni := p.adjust(binomial_p, method = "bonferroni")]

#Now the protective ones:

bmi_decreasing = fread("review/review_output/5_enrichment_analyses/2_binomial_lotta/output/neg_63_ir/results_binomial.tsv")
bmi_decreasing[, bh_fdr := p.adjust(binomial_p, method = "BH")]
bmi_decreasing[, bonferroni := p.adjust(binomial_p, method = "bonferroni")]

#Now the thigh_neg ones:

bmi_increasing = fread("review/review_output/5_enrichment_analyses/2_binomial_lotta/output/pos_78_ir/results_binomial.tsv")
bmi_increasing[, bh_fdr := p.adjust(binomial_p, method = "BH")]
bmi_increasing[, bonferroni := p.adjust(binomial_p, method = "bonferroni")]

dir.create("review/review_output/5_enrichment_analyses/2_binomial_lotta/output/summarized_results")

fwrite(all_variants, "review/review_output/5_enrichment_analyses/2_binomial_lotta/output/summarized_results/all_282_ir_loci.txt")
fwrite(lotta, "review/review_output/5_enrichment_analyses/2_binomial_lotta/output/summarized_results/lotta_test.txt")
fwrite(bmi_neutral, "review/review_output/5_enrichment_analyses/2_binomial_lotta/output/summarized_results/ns_141.txt")
fwrite(bmi_decreasing, "review/review_output/5_enrichment_analyses/2_binomial_lotta/output/summarized_results/neg_63.txt")
fwrite(bmi_increasing, "review/review_output/5_enrichment_analyses/2_binomial_lotta/output/summarized_results/pos_78.txt")

#############################################
#Let's make a table out of this shit already#
#############################################

final_df = all_variants %>%
  dplyr::select("eid", "tissue")

final_df$'n_observed (282 IR loci)' = all_variants$observed
final_df$'n_expected (282 IR loci)' = all_variants$expected
final_df$'prop_observed (282 IR loci)' = all_variants$observed_prop
final_df$'prop_expected (282 IR loci)' = all_variants$expected_prop
final_df$'p_value_binomial_test (282 IR loci)' = all_variants$binomial_p
final_df$'bonferroni_fdr_binomial_test (282 IR loci)' = all_variants$bonferroni

#Now let's order the rest:

lotta = lotta[order(match(lotta$eid, final_df$eid)),]
print(length(which(lotta$eid == final_df$eid)))

final_df$'n_observed (53 Lotta IR loci)' = lotta$observed
final_df$'n_expected (53 Lotta IR loci)' = lotta$expected
final_df$'prop_observed (53 Lotta IR loci)' = lotta$observed_prop
final_df$'prop_expected (53 Lotta IR loci)' = lotta$expected_prop
final_df$'p_value_binomial_test (53 Lotta IR loci)' = lotta$binomial_p
final_df$'bonferroni_fdr_binomial_test (53 Lotta IR loci)' = lotta$bonferroni

#BMI neutral:

bmi_neutral = bmi_neutral[order(match(bmi_neutral$eid, final_df$eid)),]
print(length(which(bmi_neutral$eid == final_df$eid)))

final_df$'n_observed (141 BMI-neutral IR loci)' = bmi_neutral$observed
final_df$'n_expected (141 BMI-neutral IR loci)' = bmi_neutral$expected
final_df$'prop_observed (141 BMI-neutral IR loci)' = bmi_neutral$observed_prop
final_df$'prop_expected (141 BMI-neutral IR loci)' = bmi_neutral$expected_prop
final_df$'p_value_binomial_test (141 BMI-neutral IR loci)' = bmi_neutral$binomial_p
final_df$'bonferroni_fdr_binomial_test (141 BMI-neutral IR loci)' = bmi_neutral$bonferroni

#And now we do the rest of the clusters:

bmi_decreasing = bmi_decreasing[order(match(bmi_decreasing$eid, final_df$eid)),]
print(length(which(bmi_decreasing$eid == final_df$eid)))

final_df$'n_observed (63 BMI-decreasing IR)' = bmi_decreasing$observed
final_df$'n_expected (63 BMI-decreasing IR)' = bmi_decreasing$expected
final_df$'prop_observed (63 BMI-decreasing IR loci)' = bmi_decreasing$observed_prop
final_df$'prop_expected (63 BMI-decreasing IR loci)' = bmi_decreasing$expected_prop
final_df$'p_value_binomial_test (63 BMI-decreasing IR loci)' = bmi_decreasing$binomial_p
final_df$'bonferroni_binomial_test (63 BMI-decreasing IR loci)' = bmi_decreasing$bonferroni

#And now we do the rest of the clusters:

bmi_increasing = bmi_increasing[order(match(bmi_increasing$eid, final_df$eid)),]
print(length(which(bmi_increasing$eid == final_df$eid)))

final_df$'n_observed (78 BMI-increasing IR)' = bmi_increasing$observed
final_df$'n_expected (78 BMI-increasing IR)' = bmi_increasing$expected
final_df$'prop_observed (78 BMI-increasing IR loci)' = bmi_increasing$observed_prop
final_df$'prop_expected (78 BMI-increasing IR loci)' = bmi_increasing$expected_prop
final_df$'p_value_binomial_test (78 BMI-increasing IR loci)' = bmi_increasing$binomial_p
final_df$'bonferroni_binomial_test (78 BMI-increasing IR loci)' = bmi_increasing$bonferroni

#And we can save this too:

fwrite(final_df, "review/review_output/5_enrichment_analyses/2_binomial_lotta/output/summarized_results/all_clusters_combined.txt")

##################
#Let's make plots#
##################

threshold_0.05 = -log10(0.05)
threshold_0.01 = -log10(0.05/102)
all_variants = all_variants[order(as.numeric(all_variants$bonferroni), decreasing = FALSE),]
all_variants$highlight = FALSE
all_variants$highlight[1:5] = TRUE
adipose_vect=c("Adipose Nuclei", "Adipose (Perrin et al)", "SGBS day 0 (Perrin et al)",  "SGBS day 4 (Perrin et al)", "SGBS day 14 (Perrin et al)")
all_variants$highlight = ifelse(all_variants$tissue%in%adipose_vect, TRUE, all_variants$highlight)

all_variants$tissue[which(all_variants$tissue == "Adipose (Perrin et al)")] = "Adipose Tissue"
all_variants$tissue[which(all_variants$tissue == "SGBS day 0 (Perrin et al)")] = "SGBS day 0"
all_variants$tissue[which(all_variants$tissue == "SGBS day 4 (Perrin et al)")] = "SGBS day 4"
all_variants$tissue[which(all_variants$tissue == "SGBS day 14 (Perrin et al)")] = "SGBS day 14"
all_variants$tissue[which(all_variants$tissue == "Bone Marrow Derived Cultured Mesenchymal Stem Cells")] = "BM Mesenchymal Stem Cells"

plot_all <- ggplot(all_variants, aes(x = observed_prop*100, y = -log10(binomial_p))) +
  # Significance threshold lines
  geom_hline(yintercept = threshold_0.05, linetype = "dashed", color = "red", size = 0.8) +
  geom_hline(yintercept = threshold_0.01, linetype = "dashed", color = "darkred", size = 0.8) +
  
  # Points with different colors for highlighted ones
  geom_point(aes(color = highlight), size = 1.75, alpha = 0.55) +
  scale_color_manual(values = c("FALSE" = "#505050", "TRUE" = "#E69F00")) +  # Grey & Orange
  
  # Labels for highlighted points
  geom_text_repel(
    aes(label = ifelse(highlight, tissue, '')),
    size = 6, 
    box.padding = 0.7, 
    point.padding = 0.5, 
    max.overlaps = Inf,
    family = "Arial", 
    force=3
  ) +
  
  # Labels and theme
  labs(
    x = "Overlap %", 
    y = expression(-log[10](P)),  # Proper math notation for logP
    title = ""
  ) +
  
  #Scales 
  
  scale_x_continuous(limits = c(0, 100),
                     breaks = seq(0, 100, by = 10)) +
  
  scale_y_continuous( limits = c(0, 14),
                      breaks = seq(0, 14, by = 2)) +
  
  # Nature-style theme
  theme_classic(base_size = 18, base_family = "Arial") +  
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
    axis.title = element_text(face = "bold", size = 18),
    axis.text = element_text(size = 16),
    axis.line = element_line(size = 1.2),  # Thicker axis lines
    axis.ticks = element_line(size = 1.2),  # Thicker ticks
    axis.ticks.length = unit(0.3, "cm"),  # Inward ticks
    panel.grid.major = element_blank(),  # No major grid
    panel.grid.minor = element_blank(),  # No minor grid
    legend.position = "none"  # Remove legend for clean look
  )

dir.create("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/review/manuscript/figures")
dir.create("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/review/manuscript/figures/raw")
dir.create("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/review/manuscript/figures/clean")

ggsave(
  filename = "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/review/manuscript/figures/raw/282_ir_loci_binomial_enrichment_plot.svg",
  plot = plot_all,
  width = 12,
  height = 10,
  units = "in",  # or "in" depending on your preference, but cm is common
  device = "svg"
)

#######################################################################################################################
#Let's add the other plots together with the same format utilizing a function, but adding a title to make thing easier#
#######################################################################################################################

lotta_plot = plot_subclusters(lotta, "53 IR loci (Lotta et al)")
decreasing_plot = plot_subclusters(bmi_decreasing, "63 BMI-decreasing IR loci")
neutral_plot = plot_subclusters(bmi_neutral, "141 BMI-neutral IR loci")
increasing_plot = plot_subclusters(bmi_increasing, "78 BMI-increasing IR loci")

library(patchwork)

combined_plot=lotta_plot + neutral_plot + decreasing_plot + increasing_plot  + plot_layout(ncol=2)

ggsave(
  filename = "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/review/manuscript/figures/raw/cluster_binomial_enrichment_plot.svg",
  plot = combined_plot,
  width = 16,
  height = 16,
  units = "in",  # or "in" depending on your preference, but cm is common
  device = "svg"
)
