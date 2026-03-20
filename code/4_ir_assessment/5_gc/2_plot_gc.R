##############
#INTRODUCTION#
##############

#This code performs genetic correlations

###################
#Loading libraries#
###################

library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(forcats)
library(reshape2)

###################
#Loading functions#
###################

curated_2_munging_mv <- function(curated_df){
  
  #STEP 0: let's remove those that have not rsID, cuz LDSC gets confused:
  
  curated_df <- curated_df[which(str_detect(as.character(curated_df$variant), "rs") == TRUE)]
  
  #STEP 1: get hte MAF:
  
  curated_df$MAF <- curated_df$minimum_allele_frequency
  
  #STEP 2: get the right columns:
  
  curated_df_4_ldsc <- curated_df %>%
    select(variant, effect_allele, other_allele, beta, standard_error, p_value, MAF, sample_size)
  
  colnames(curated_df_4_ldsc) <- c("SNP", "A1", "A2", "BETA", "SE", "P", "MAF", "N")
  
  #STEP 3: Save the data:
  
  return(curated_df_4_ldsc)
  
}

curated_2_munging <- function(curated_df){
  
  #STEP 0: let's remove those that have not rsID, cuz LDSC gets confused:
  
  curated_df <- curated_df[which(str_detect(as.character(curated_df$variant), "rs") == TRUE)]
  
  #STEP 1: get hte MAF:
  
  curated_df$MAF <- ifelse(as.numeric(curated_df$effect_allele_frequency) > 0.50, 1-as.numeric(curated_df$effect_allele_frequency), as.numeric(curated_df$effect_allele_frequency))
  
  #STEP 2: get the right columns:
  
  curated_df_4_ldsc <- curated_df %>%
    select(variant, effect_allele, other_allele, beta, standard_error, p_value, MAF, sample_size)
  
  colnames(curated_df_4_ldsc) <- c("SNP", "A1", "A2", "BETA", "SE", "P", "MAF", "N")
  
  #STEP 3: Save the data:
  
  return(curated_df_4_ldsc)
  
}

curated_2_binary <- function(curated_df){
  
  #STEP 0: let's remove those that have not rsID, cuz LDSC gets confused:
  
  curated_df <- curated_df[which(str_detect(as.character(curated_df$variant), "rs") == TRUE)]
  
  #STEP 1: get hte MAF:
  
  curated_df$MAF <- ifelse(as.numeric(curated_df$effect_allele_frequency) > 0.50, 1-as.numeric(curated_df$effect_allele_frequency), as.numeric(curated_df$effect_allele_frequency))
  
  #STEP 2: get the effect sample size:
  
  curated_df$Neff<-4/((2*curated_df$MAF*(1-curated_df$MAF))*curated_df$standard_error^2)
  
  #STEP 2: get the right columns:
  
  curated_df_4_ldsc <- curated_df %>%
    select(variant, effect_allele, other_allele, beta, standard_error, p_value, MAF, Neff)
  
  colnames(curated_df_4_ldsc) <- c("SNP", "A1", "A2", "BETA", "SE", "P", "MAF", "N")
  
  #STEP 3: Save the data:
  
  return(curated_df_4_ldsc)
  
}

get_lower_tri_index <- function(i, j) {
  if (i < j) {
    tmp <- i
    i <- j
    j <- tmp
  }
  return(i * (i - 1) / 2 + j)
}

##############
#Loading data#
##############

ir_fiadjbmi_fiadjbmi=c(0.4679,0.0406)

ir_fiadjbmi_hdl=c(-0.8467,0.0388)

ir_fiadjbmi_tg=c(0.844,0.0431)

ir_fiadjbmi_tg_hdl_ratio=c(0.9322,0.0427)

ir_fiadjbmi_isiadjbmi=c(-0.3526,0.0403)

ir_fiadjbmi_ifcadjbmi=c(0.1629,0.0483)

ir_fiadjbmi_fgadjbmi=c(0.1593,0.0279)

ir_fiadjbmi_thgadjbmi=c(0.2243,0.0403)

ir_fiadjbmi_bmi=c(0.467,0.0198)

ir_fiadjbmi_hc=c(0.334,0.0234)

ir_fiadjbmi_hcadjbmi=c(-0.1276,0.0259)

ir_fiadjbmi_wc=c(0.5444,0.0203)

ir_fiadjbmi_wcadjbmi=c(0.2384,0.0284)

ir_fiadjbmi_whr=c(0.6098,0.0258)

ir_fiadjbmi_whradjbmi=c(0.396,0.031)

ir_fiadjbmi_t2d=c(0.5617,0.028)

ir_fiadjbmi_nafld=c(0.657,0.0738)

ir_fiadjbmi_chd=c(0.3062,0.0301)

ir_fiadjbmi_ckd=c(0.3665,0.0491)

ir_fiadjbmi_pcos=c(0.2353,0.0341)

ir_fiadjbmi_hypertension=c(0.3378,0.0231)

# List of all variable names (can automate this if all start with ir_fiadjbmi_)
var_names <- ls(pattern = "^ir_fiadjbmi_")

# Extract correlation and SE, build data frame
cor_list <- lapply(var_names, function(vname) {
  vals <- get(vname)
  trait2 <- sub("^ir_fiadjbmi_", "", vname)
  data.frame(
    Trait1 = "ir_fiadjbmi",
    Trait2 = trait2,
    Correlation = vals[1],
    SE = vals[2]
  )
})

cor_df <- do.call(rbind, cor_list)

# Calculate additional stats
cor_df <- cor_df %>%
  mutate(
    Z = Correlation / SE,
    P = 2 * pnorm(-abs(Z)),
    CI_lower = Correlation - 1.96 * SE,
    CI_upper = Correlation + 1.96 * SE,
    CI = sprintf("[%.3f, %.3f]", CI_lower, CI_upper),
    P_formatted = signif(P, 2),
    Partner = Trait2
  )

# Save to CSV
library(data.table)
fwrite(
  cor_df %>% select(Partner, Correlation, SE, CI, P_formatted),
  "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/output/4_ir_loci_discovery/5_gc/3_gc/genetic_correlations_res_clean_manual.csv"
)

##########################
#Now let's plot the data!#
##########################

# Your data prep (mapping trait names, filtering, stars, CIs) stays the same
cor_df <- cor_df %>%
  mutate(
    Trait1 = "Common factor GWAS",
    Pvalue = P
  )

pretty_names <- c(
  "fiadjbmi" = "FIadjBMI",
  "hdl" = "HDL",
  "tg" = "TG",
  "isiadjbmi" = "ISIadjBMI",
  "ifcadjbmi" = "IFCadjBMI",
  "tg_hdl_ratio" = "TG/HDL",
  "thgadjbmi" = "2hGadjBMI",
  "fgadjbmi" = "FGadjBMI",
  "chd" = "CHD",
  "ckd" = "CKD",
  "hypertension" = "Hypertension",
  "nafld"="NAFLD",
  "pcos" = "PCOS",
  "t2d" = "T2D",
  "bmi" = "BMI", 
  "hc" = "HC",
  "hcadjbmi" = "HCadjBMI",
  "wc" = "WC",
  "wcadjbmi" = "WCadjBMI",
  "whr" = "WHR",
  "whradjbmi" = "WHRadjBMI"
  
)

trait_order <- c("Common factor GWAS",    
                 "FIadjBMI",
                  "HDL",
                 "TG",
                  "ISIadjBMI",
                 "IFCadjBMI",
                 "TG/HDL",
                 "2hGadjBMI",
                 "FGadjBMI",
                 "CHD",
                 "CKD",
                 "Hypertension",
                 "NAFLD",
                 "PCOS",
                 "T2D",
                 "BMI", 
                 "HC",
                 "HCadjBMI",
                 "WC",
                 "WCadjBMI",
                 "WHR",
                 "WHRadjBMI")

cor_df <- cor_df %>%
  mutate(
    Trait2 = recode(Trait2, !!!pretty_names)
  ) %>%
  filter(Trait2 %in% trait_order) %>%
  filter(Trait1 %in% trait_order) %>%
  mutate(
    stars = case_when(
      Pvalue < 5e-8 ~ "**",
      Pvalue < 0.05 ~ "*",
      TRUE ~ ""
    ),
    label = paste0(sprintf("%.2f", Correlation), stars),
    CI_lower = Correlation - 1.96 * SE,
    CI_upper = Correlation + 1.96 * SE
  )

# Order Trait2 by increasing Correlation
cor_df$Trait2 <- factor(cor_df$Trait2, levels = cor_df$Trait2[order(cor_df$Correlation)])

#Let's add a column for type:

cor_df$TraitType <- NA
cor_df$TraitType <- ifelse(cor_df$Trait2%in%c("FIadjBMI", "HDL", "TG"), "Hallmark IR", cor_df$TraitType)
cor_df$TraitType <- ifelse(cor_df$Trait2%in%c("ISIadjBMI",
                                              "IFCadjBMI",
                                              "TG/HDL",
                                              "2hGadjBMI",
                                              "FGadjBMI"), "Glycemic", cor_df$TraitType)

cor_df$TraitType <- ifelse(cor_df$Trait2%in%c("CHD",
                                              "CKD",
                                              "Hypertension",
                                              "NAFLD",
                                              "PCOS",
                                              "T2D"), "IR-derived disease", cor_df$TraitType)

cor_df$TraitType <- ifelse(cor_df$Trait2%in%c("BMI", 
                                              "HC",
                                               "HCadjBMI",
                                               "WC",
                                               "WCadjBMI",
                                               "WHR",
                                               "WHRadjBMI"), "Anthropometric", cor_df$TraitType)

cor_df$TraitType <- factor(cor_df$TraitType, levels = c("Hallmark IR",
                                                        "Glycemic",
                                                        "IR-derived disease",
                                                        "Anthropometric"))

#Let's change TG/HDL so that it makes sense:

cor_df$CI_upper[which(cor_df$Trait2 == "TG/HDL")] <- 1.00

plotio <- ggplot(cor_df, aes(x = Correlation, y = Trait2)) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), 
                 height = 0.3, color = "black", size = 1) +  # black error bars
  geom_point(aes(color = TraitType, fill=TraitType), size = 5, stroke = 1.5, shape = 21) +  # earthy points
  geom_text(aes(label = stars), nudge_x = 0.12, size = 7, color = "black", fontface = "bold") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black", size = 0.8) +
  scale_x_continuous(
    limits = c(-1.1, 1.1),
    breaks = seq(-1, 1, 0.25),
    expand = c(0, 0)
  ) +
  facet_grid(rows = vars(TraitType), scales = "free_y", space = "free_y") +   
  theme_minimal(base_size = 15, base_family = "Times") +
  labs(
    x = "Genetic Correlation",
    y = NULL
  ) +
  theme(
    legend.position = "none",
    axis.text.y = element_text(face = "bold", size = 14, color = "black"),
    axis.text.x = element_text(size = 13, color = "black", face = "bold"),
    axis.title.x = element_text(face = "bold", size = 15, color = "black"),
    strip.text.y = element_text(face = "bold", size = 14, color = "black"), 
    panel.grid.major.y = element_blank(),
    panel.grid.major.x = element_line(color = "grey70", size = 0.5, linetype = "dashed"),  # 🔥 thinner & dashed
    panel.grid.minor = element_blank(),
    plot.margin = margin(15, 15, 15, 15),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)  # 🔥 adds boxes
  )



ggsave(
  filename = "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/manuscript/figures/gc_glycemic_traits_forest_plot.svg",
  plot = plotio,
  width = 10,
  height = 10,
  units = "in",  # or "in" depending on your preference, but cm is common
  device = "svg"
)



