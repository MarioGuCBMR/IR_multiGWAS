##############
#INTRODUCTION#
##############

#This code generates the LDSC forest plots for the ir_mvgwas<-

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)

###################
#Loading functions#
###################

source("N:/SUN-CBMR-Kilpelainen-Group/Team projects/Mario&German/Smoking_WHR/CODE_4_Replication/Mediation_analysis/code/0_functions/functions_4_mediation.R")

ldsc_plot_producer <- function(grs_df, title){
  
  #STEP 0: let's transform in factors the traits that we columns that we need:
  
  grs_df$Trait <- factor(grs_df$trait, levels = c("FIadjBMI", "HDL", "ISIadjBMI", "TG",  "TG/HDL", "BMI", "HC", "HCadjBMI", "WC", "WCadjBMI", "WHR", "WHRadjBMI", "CHD", "CKD", "Hypertension", "NAFLD", "PCOS", "T2D"))
  
  grs_df$type <- c("IR", "IR", "IR", "IR", "IR", "Anthropometric", "Anthropometric", "Anthropometric", "Anthropometric", "Anthropometric", "Anthropometric", "Anthropometric", "IR-derived disease", "IR-derived disease", "IR-derived disease", "IR-derived disease", "IR-derived disease", "IR-derived disease")
  
  grs_df$type = factor(grs_df$type, levels=c("IR",
                                           "Anthropometric",
                                           "IR-derived disease"))
  
  grs_df$LDSC = factor(grs_df$LDSC, levels=c("IRadjBMI"))
  
  plotio <- 
    
    #SETTING the info where we are all gonna work
    
    ggplot(grs_df, aes(x= fct_rev(Trait),y = as.numeric(gc), ymin=-1, ymax=1, shape=LDSC)) +
    
    #Generating geom_point:
    
    geom_point(aes(shape=LDSC, color=Trait), size = 5,  position = position_dodge(width = 0.75)) +
    scale_shape_manual(name= "LDSC", values=c("IR" = 17,
                                             "IRadjBMI" = 16)) +
    
    scale_y_continuous(breaks = scales::pretty_breaks(n = 5))+
    
    #Adding the error bars:
    
    geom_errorbar(aes(ymin=as.numeric(lower_ci), ymax=as.numeric(upper_ci)), width = 0.15, cex= 1, position = position_dodge(width = 0.75), color="black", linewidth=1) +
    
    #Generating a facet to distribute the data beautifully:
    
    facet_wrap(~type, strip.position="left", nrow=7, scales = "free_y") +
    geom_hline(yintercept = 0, color = "red", size = 0.5, linetype=2) +
    
    #First geom_text:
    
    #geom_text(aes(label=grs_df$nsnps),
    #          color = "black", position = position_dodge(width = 0.75), vjust = 1.6, size=3.5) +
    
    #Second geom_text:
    
    #geom_text(aes(label=as.character(grs_df$full_effect)),
    #          color = "black", position = position_dodge(width = 0.75), vjust = 1.9, size=3.5) +
    
  
  coord_flip() + 
    
    
    #SETTING AXIS OPTIONS:
    
    xlab("Trait") +
    ylab("Genetic correlation") +
    
    #SETTING TITLE OPTIONS:
    
    ggtitle(paste0(title)) +
    theme_bw() +
    theme(legend.position="none") +
    #theme(legend.text=element_text(size=11)) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),panel.background = element_blank(), axis.line = element_line(colour = "black"))+
    theme(axis.text=element_text(size=11), axis.title=element_text(size=10,face="bold"), plot.title = element_text(size = 12, face = "bold", hjust=0.5)) 
  
  return(plotio)
  
}


###################################
#Loading data from the log of LDSC# J:\CBMR\SUN-CBMR-Kilpelainen-Group\Mario_Tools\IR_GSEM_2023\output\5_comparing_fiadjbmi_fi\3_data_4_munging\4_munged_data\ir_fiadjbmi.sumstats.gz_ir_fi.sumstats.gz_fiadjbmi.sumstats.gz_fi.sumstats.gz_hdl.sumstats.gz_tg.sum_ldsc
###################################

ir_fiadjbmi_fiadjbmi<- c(0.4679 ,0.0429)
ir_fiadjbmi_hdl<- c(-0.8467 ,0.0395)
ir_fiadjbmi_tg<- c(0.844 ,0.044)
ir_fiadjbmi_isiadjbmi<- c(-0.3527, 0.0405)
ir_fiadjbmi_tg_hdl_ratio<- c(0.9322, 0.0445)

ir_fiadjbmi_bmi<- c(0.467 ,0.0198)
ir_fiadjbmi_hc<- c(0.334, 0.0229)
ir_fiadjbmi_hcadjbmi <- c(-0.1276 ,0.0259)
ir_fiadjbmi_wc<- c(0.5438, 0.0203)
ir_fiadjbmi_wcadjbmi<- c(0.2384 ,0.0286)
ir_fiadjbmi_whr<- c(0.6074, 0.0266)
ir_fiadjbmi_whradjbmi<- c(0.3949 ,0.0314)

ir_fiadjbmi_t2d<- c(0.559 ,0.0286)
ir_fiadjbmi_nafld<- c(0.6551 ,0.0699)
ir_fiadjbmi_chd<- c(0.3009 ,0.0316)
ir_fiadjbmi_ckd<- c(0.3603 ,0.0482)
ir_fiadjbmi_pcos<- c(0.2341 ,0.0333)
ir_fiadjbmi_hypertension<- c(0.3382 ,0.025)

#############################
#Merge the data for IRadjBMI#
#############################

iradjbmi <- rbind(as.data.frame(t(ir_fiadjbmi_fiadjbmi)), as.data.frame(t(ir_fiadjbmi_hdl)), as.data.frame(t(ir_fiadjbmi_isiadjbmi)), as.data.frame(t(ir_fiadjbmi_tg)), as.data.frame(t(ir_fiadjbmi_tg_hdl_ratio)),
                  as.data.frame(t(ir_fiadjbmi_bmi)), as.data.frame(t(ir_fiadjbmi_hc)),  as.data.frame(t(ir_fiadjbmi_hcadjbmi)), as.data.frame(t(ir_fiadjbmi_wc)), as.data.frame(t(ir_fiadjbmi_wcadjbmi)), as.data.frame(t(ir_fiadjbmi_whr)), as.data.frame(t(ir_fiadjbmi_whradjbmi)),
                  as.data.frame(t(ir_fiadjbmi_chd)), as.data.frame(t(ir_fiadjbmi_ckd)), as.data.frame(t(ir_fiadjbmi_hypertension)), as.data.frame(t(ir_fiadjbmi_nafld)), as.data.frame(t(ir_fiadjbmi_pcos)), as.data.frame(t(ir_fiadjbmi_t2d))
                  )

colnames(iradjbmi) <- c("gc", "se")

iradjbmi$trait <- c("FIadjBMI", "HDL", "ISIadjBMI", "TG",  "TG/HDL", "BMI", "HC", "HCadjBMI", "WC", "WCadjBMI", "WHR", "WHRadjBMI", "CHD", "CKD", "Hypertension", "NAFLD", "PCOS", "T2D")
iradjbmi$type <- c("IR", "IR", "IR", "IR", "IR", "Anthropometric", "Anthropometric", "Anthropometric", "Anthropometric", "Anthropometric", "Anthropometric", "Anthropometric", "IR-derived disease", "IR-derived disease", "IR-derived disease", "IR-derived disease", "IR-derived disease", "IR-derived disease")

########################################################################
#Now let's get the LDSC column and set up the LDSC confidence intervals#
########################################################################

iradjbmi$LDSC <- "IRadjBMI"

ldsc_df <- iradjbmi

ldsc_df$lower_ci <- lower_ci(ldsc_df$gc, ldsc_df$se)
ldsc_df$upper_ci <- upper_ci(ldsc_df$gc, ldsc_df$se)

betas_rounded <- format(ldsc_df$gc, scientific = FALSE, digits = 4)
lower_ci_rounded <- format(ldsc_df$lower_ci,scientific = FALSE, digits = 4)
upper_ci_rounded <- format(ldsc_df$upper_ci, scientific = FALSE, digits = 4)

ldsc_df$full_effect <- paste(betas_rounded,
                            " (",
                            lower_ci_rounded,
                            ",",
                            upper_ci_rounded,
                            ")",
                            sep=""
)

#Let's save this table for the supplementary information:

data.table::fwrite(ldsc_df, "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2023/report/tables/separated_tables/genetic_correlations.csv")

#Let's make the plot:

p <- ldsc_plot_producer(ldsc_df, "Genetic Correlations for IR mvGWAS")

#Save it...

saveRDS(p, "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2023/output/5_comparing_fiadjbmi_fi/3_gc/3_gc/plot_genetic_correlation_clean.RDS")
