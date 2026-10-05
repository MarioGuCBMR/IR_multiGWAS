##############
#INTRODUCTION#
##############

#Let's try providing a better false discovery strategy for our genetic correlations and PRS.

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)

#######################################
#Loading data for genetic correlations#
#######################################

suppl3=readxl::read_xlsx("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/review/manuscript/10082026/Supplementary Tables_11August2026.xlsx", sheet=4)
suppl3= as.data.frame(suppl3)

#Genetic correlations is 21 rows:

gc=suppl3[3:23,] #6 columns are p-values:
gc=gc[,1:6]
colnames(gc) = c("Type", "Trait", "Genetic Correlation", "95% Confidence intervals", "Standard error", "P-value")
raw_pvals=as.numeric(gc[,6]) #correct
gc$bh_pvals=p.adjust(raw_pvals,method = "BH")
gc$bonferroni_pvals = p.adjust(raw_pvals,method = "bonferroni")

dir.create("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/review/manuscript/tables")

fwrite(gc, "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/review/manuscript/tables/gc_with_new_fdr.csv")

####################################################
#Let's do the same for weighted PRS for 282 IR loci#
####################################################

weigthed_282 = suppl3[27:58,]
weigthed_282=weigthed_282[,1:7]
colnames(weigthed_282) = c("Type", "Trait", "Number of SNPs", "PRS", "95% Confidence intervals", "Standard error", "P-value")
raw_pvals=as.numeric(weigthed_282[,7]) #correct
weigthed_282$bh_pvals=p.adjust(raw_pvals,method = "BH")
weigthed_282$bonferroni_pvals = p.adjust(raw_pvals,method = "bonferroni")

fwrite(weigthed_282, "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/review/manuscript/tables/weighted_282_prs_with_new_fdr.csv")

######################################################
#Let's do the same for unweighted PRS for 282 IR loci#
######################################################

unweigthed_282 = suppl3[62:93,] #start in 59, but we are skipping column names
unweigthed_282=unweigthed_282[,1:7]
colnames(unweigthed_282) = c("Type", "Trait", "Number of SNPs", "PRS", "95% Confidence intervals", "Standard error", "P-value")
raw_pvals=as.numeric(unweigthed_282[,7]) #correct
unweigthed_282$bh_pvals=p.adjust(raw_pvals,method = "BH")
unweigthed_282$bonferroni_pvals = p.adjust(raw_pvals,method = "bonferroni")

fwrite(unweigthed_282, "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/review/manuscript/tables/unweighted_282_prs_with_new_fdr.csv")

####################################################
#Let's do the same for weighted PRS for 141 IR loci#
####################################################

weigthed_141 = suppl3[27:58,]
weigthed_141 = weigthed_141[,c(1:2,8:12)]
colnames(weigthed_141) = c("Type", "Trait", "Number of SNPs", "PRS", "95% Confidence intervals", "Standard error", "P-value")
raw_pvals=as.numeric(weigthed_141[,7]) #correct
weigthed_141$bh_pvals=p.adjust(raw_pvals,method = "BH")
weigthed_141$bonferroni_pvals = p.adjust(raw_pvals,method = "bonferroni")

fwrite(weigthed_141, "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/review/manuscript/tables/weighted_141_prs_with_new_fdr.csv")

####################################################
#Let's do the same for weighted PRS for 141 IR loci#
####################################################

unweigthed_141 = suppl3[62:93,] #start in 59, but we are skipping column names
unweigthed_141 = unweigthed_141[,c(1:2,8:12)]
colnames(unweigthed_141) = c("Type", "Trait", "Number of SNPs", "PRS", "95% Confidence intervals", "Standard error", "P-value")
raw_pvals=as.numeric(unweigthed_141[,7]) #correct
unweigthed_141$bh_pvals=p.adjust(raw_pvals,method = "BH")
unweigthed_141$bonferroni_pvals = p.adjust(raw_pvals,method = "bonferroni")

fwrite(unweigthed_141, "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/review/manuscript/tables/unweighted_141_prs_with_new_fdr.csv")

####################################################
#Let's do the same for weighted PRS for 63 IR loci#
####################################################

weigthed_63 = suppl3[27:58,]
weigthed_63 = weigthed_63[,c(1:2,13:17)]
colnames(weigthed_63) = c("Type", "Trait", "Number of SNPs", "PRS", "95% Confidence intervals", "Standard error", "P-value")
raw_pvals=as.numeric(weigthed_63[,7]) #correct
weigthed_63$bh_pvals=p.adjust(raw_pvals,method = "BH")
weigthed_63$bonferroni_pvals = p.adjust(raw_pvals,method = "bonferroni")

fwrite(weigthed_63, "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/review/manuscript/tables/weighted_63_prs_with_new_fdr.csv")

####################################################
#Let's do the same for weighted PRS for 141 IR loci#
####################################################

unweigthed_63 = suppl3[62:93,] #start in 59, but we are skipping column names
unweigthed_63 = unweigthed_63[,c(1:2,13:17)]
colnames(unweigthed_63) = c("Type", "Trait", "Number of SNPs", "PRS", "95% Confidence intervals", "Standard error", "P-value")
raw_pvals=as.numeric(unweigthed_63[,7]) #correct
unweigthed_63$bh_pvals=p.adjust(raw_pvals,method = "BH")
unweigthed_63$bonferroni_pvals = p.adjust(raw_pvals,method = "bonferroni")

fwrite(unweigthed_63, "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/review/manuscript/tables/unweighted_63_prs_with_new_fdr.csv")

####################################################
#Let's do the same for weighted PRS for 78 IR loci#
####################################################

weigthed_78 = suppl3[27:58,]
weigthed_78 = weigthed_78[,c(1:2,18:22)]
colnames(weigthed_78) = c("Type", "Trait", "Number of SNPs", "PRS", "95% Confidence intervals", "Standard error", "P-value")
weigthed_78$`P-value`[18] = 6.50e-321 #a comma that could not be properly converted fuck this up. Solved!
raw_pvals=as.numeric(weigthed_78[,7]) #correct
weigthed_78$bh_pvals=p.adjust(raw_pvals,method = "BH")
weigthed_78$bonferroni_pvals = p.adjust(raw_pvals,method = "bonferroni")

fwrite(weigthed_78, "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/review/manuscript/tables/weighted_78_prs_with_new_fdr.csv")

####################################################
#Let's do the same for weighted PRS for 141 IR loci#
####################################################

unweigthed_78 = suppl3[62:93,] #start in 59, but we are skipping column names
unweigthed_78 = unweigthed_78[,c(1:2,18:22)]
colnames(unweigthed_78) = c("Type", "Trait", "Number of SNPs", "PRS", "95% Confidence intervals", "Standard error", "P-value")
raw_pvals=as.numeric(unweigthed_78[,7]) #correct
unweigthed_78$bh_pvals=p.adjust(raw_pvals,method = "BH")
unweigthed_78$bonferroni_pvals = p.adjust(raw_pvals,method = "bonferroni")

fwrite(unweigthed_78, "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/review/manuscript/tables/unweighted_78_prs_with_new_fdr.csv")
