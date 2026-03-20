##############
#INTRODUCTION#
##############

#This code generates the Manhattan and QQ-plot for the IRadjBMI GWAS (IR mvGWAS)

#GWASTools requires of Matrix1.6.0 which had to be downloaded from CRAN!
#https://cran.r-project.org/src/contrib/Archive/Matrix/Matrix_1.6-0.tar.gz

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)
library(GWASTools)
library(Cairo)
library(ggplot2)
library(ggsci)

######################################
#Let's get the lambda for all of them#
######################################


lambda_compute= function(pvals){
  
  chisq_vals <- qchisq(1 - pvals, df = 1)
  lambda_gc <- median(chisq_vals) / qchisq(0.5, df = 1)
  return(lambda_gc)
  
}


##############
#Loading data#
##############

path_2_input <- "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/"

setwd(path_2_input)

#Loading IR GWAS.

fiadjbmi <- fread("output/1_curated_gwas/fiadjbmi_curated.txt")
hdl <- fread("output/1_curated_gwas/hdl_curated.txt")
tg <- fread("output/1_curated_gwas/tg_curated.txt")
tg_hdl_ratio <- fread("output/1_curated_gwas/tg_hdl_ratio_curated.txt") #just a test
iradjbmi <- fread("output/2_mv_gwas/fiadjbmi_hdl_tg_common_dwls_curated.txt")

tg_hdl_ratio$variant <- tg_hdl_ratio$rs_id

#Let's compute the lambda:

lambda_fiadjbmi = lambda_compute(as.numeric(fiadjbmi$p_value)) #1.06
lambda_hdl = lambda_compute(as.numeric(hdl$p_value)) #1.06
lambda_tg = lambda_compute(as.numeric(tg$p_value)) #1.09
lambda_tg_hdl = lambda_compute(as.numeric(tg_hdl_ratio$p_value)) #1.32

####################
#Generating QQ-plot#
####################

sample_size <- 10000
fiadjbmi_pvals <- sample(fiadjbmi$p_value, size = sample_size, replace = FALSE)

CairoSVG("manuscript/figures/qq_plot_fiadjbmi.svg", width = 6, height = 6)
qqPlot(fiadjbmi_pvals)
dev.off()

tg_pvals <- sample(tg$p_value, size = sample_size, replace = FALSE)

CairoSVG("manuscript/figures/qq_plot_tg.svg", width = 6, height = 6)
qqPlot(as.numeric(tg_pvals))
dev.off()

hdl_pvals <- sample(hdl$p_value, size = sample_size, replace = FALSE)

CairoSVG("manuscript/figures/qq_plot_hdl.svg", width = 6, height = 6)
qqPlot(as.numeric(hdl_pvals))
dev.off()

common_pvals <- sample(iradjbmi$p_value, size = sample_size, replace = FALSE)

CairoSVG("manuscript/figures/qq_plot_common_gwas.svg", width = 6, height = 6)
qqPlot(common_pvals)
dev.off()

########################################
#Let's try to see what we can work with#
########################################

# Step 1: Filter each dataset first (keep only -log10(p) < 20)
filter_low_significance <- function(df) {
  df %>% filter(-log10(as.numeric(p_value)) < 20)
}

fiadjbmi_filt <- filter_low_significance(fiadjbmi)
hdl_filt <- filter_low_significance(hdl)
tg_filt <- filter_low_significance(tg)
tg_hdl_ratio_filt <- filter_low_significance(tg_hdl_ratio)
iradjbmi_filt <- filter_low_significance(iradjbmi)

# Step 2: Get common variants after filtering
common_variants <- Reduce(intersect, list(
  fiadjbmi_filt$variant,
  hdl_filt$variant,
  tg_filt$variant,
  tg_hdl_ratio_filt$variant,
  iradjbmi_filt$variant
))

# Step 3: Sample 10,000 variants
set.seed(123)
subset_variants <- sample(common_variants, 10000)

# Step 4: Extract p-values for these variants
get_pvals <- function(df, name) {
  df %>%
    filter(variant %in% subset_variants) %>%
    select(variant, p_value) %>%
    mutate(trait = name)
}

fiadjbmi_sub <- get_pvals(fiadjbmi, "FIadjBMI")
hdl_sub <- get_pvals(hdl, "HDL")
tg_sub <- get_pvals(tg, "TG")
tg_hdl_ratio_sub <- get_pvals(tg_hdl_ratio, "TG/HDL ratio")
iradjbmi_sub <- get_pvals(iradjbmi, "Common factor GWAS")

# Step 5: Combine all into long format
all_pvals <- rbind(
  fiadjbmi_sub,
  hdl_sub,
  tg_sub,
  tg_hdl_ratio_sub,
  iradjbmi_sub
)

# Step 6: Compute expected vs observed
all_pvals <- all_pvals %>%
  group_by(trait) %>%
  arrange(as.numeric(p_value)) %>%
  mutate(
    observed = -log10(as.numeric(p_value)),
    expected = -log10(ppoints(n()))
  )

# Manually defined high-contrast, colorblind-safe palette
plot <- ggplot(all_pvals, aes(x = expected, y = observed, color = trait)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.7) +
  geom_point(alpha = 0.75, size = 1.8) +
  scale_color_manual(
    values = cud_colors,
    name = "Trait",
    guide = guide_legend(override.aes = list(size = 5))  # BIGGER legend points
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 15, face = "bold"),
    legend.text = element_text(size = 14),
    legend.key.size = unit(1.8, "lines"),  # Bigger legend keys
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 13, hjust = 0.5),
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 12)
  ) +
  labs(
    #title = "Overlapping QQ Plot Across Traits",
    #subtitle = "10,000 variants with -log10(p) < 20 in all traits",
    x = expression(Expected~~-log[10](italic(p))),
    y = expression(Observed~~-log[10](italic(p)))
  )

##############################
#Finally, let's save the plot#
##############################

ggsave(
  filename = "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/manuscript/figures/qqplots_combined_10000_shared_LOGP_max_20.svg",
  plot = plot,
  width = 12,
  height = 9,
  units = "in",  # or "in" depending on your preference, but cm is common
  device = "svg"
)
