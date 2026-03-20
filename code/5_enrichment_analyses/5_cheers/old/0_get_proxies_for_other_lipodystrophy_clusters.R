##############
#INTRODUCTION#
##############

#This code performs CPASSOC with our three traits!!

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)

##############
#Loading data#
##############

#Let's change the working directory:

path_2_files <-  "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/" #change it with your own path.

#########################################################
#Let's load all the previous lipodystrophy-like clusters#
#########################################################

lotta <- as.data.frame(readxl::read_excel("raw_data/previous_loci/lotta_et_al.xlsx"))

suzuki <- as.data.frame(readxl::read_excel("raw_data/previous_loci/suzuki_et_al.xlsx"))

#Let's take Suzuki's hits that are lipodystrophy-like:

beta_cell_pi_neg <- suzuki[which(suzuki$`Cluster assignment` == "Beta cell -PI"),]
beta_cell_pi_pos <- suzuki[which(suzuki$`Cluster assignment` == "Beta cell +PI"),]
body_fat <- suzuki[which(suzuki$`Cluster assignment` == "Body fat"),]
lipodystrophy <- suzuki[which(suzuki$`Cluster assignment` == "Lipodystrophy"),]
liver <- suzuki[which(suzuki$`Cluster assignment` == "Liver/lipid metabolism"),]
metabolic_syndrome <- suzuki[which(suzuki$`Cluster assignment` == "Metabolic syndrome"),]
obesity <- suzuki[which(suzuki$`Cluster assignment` == "Obesity"),]
residual <- suzuki[which(suzuki$`Cluster assignment` == "Residual glycaemic"),]

######################################################
#Let's format each of the data so that we can combine#
######################################################

dir.create("output/3_cheers")
dir.create("output/3_cheers/0_proxies")

#Let's save thje data

fwrite(as.data.frame(lotta$SNP), "output/3_cheers/0_proxies/lotta_4_haploreg.txt", col.names = TRUE, row.names = FALSE, quote=FALSE, sep = "\t")
fwrite(as.data.frame(beta_cell_pi_neg$`Index SNV`), "output/3_cheers/0_proxies/suzuki_beta_cell_neg_4_haploreg.txt", col.names = TRUE, row.names = FALSE, quote=FALSE, sep = "\t")
fwrite(as.data.frame(beta_cell_pi_pos$`Index SNV`), "output/3_cheers/0_proxies/suzuki_beta_cell_pos_4_haploreg.txt", col.names = TRUE, row.names = FALSE, quote=FALSE, sep = "\t")
fwrite(as.data.frame(body_fat$`Index SNV`), "output/3_cheers/0_proxies/suzuki_body_fat_4_haploreg.txt", col.names = TRUE, row.names = FALSE, quote=FALSE, sep = "\t")
fwrite(as.data.frame(lipodystrophy$`Index SNV`), "output/3_cheers/0_proxies/suzuki_lipodystrophy_4_haploreg.txt", col.names = TRUE, row.names = FALSE, quote=FALSE, sep = "\t")
fwrite(as.data.frame(liver$`Index SNV`), "output/3_cheers/0_proxies/suzuki_liver_4_haploreg.txt", col.names = TRUE, row.names = FALSE, quote=FALSE, sep = "\t")
fwrite(as.data.frame(metabolic_syndrome$`Index SNV`), "output/3_cheers/0_proxies/suzuki_metabolic_syndrome_4_haploreg.txt", col.names = TRUE, row.names = FALSE, quote=FALSE, sep = "\t")
fwrite(as.data.frame(obesity$`Index SNV`), "output/3_cheers/0_proxies/suzuki_obesity_4_haploreg.txt", col.names = TRUE, row.names = FALSE, quote=FALSE, sep = "\t")
fwrite(as.data.frame(residual$`Index SNV`), "output/3_cheers/0_proxies/suzuki_residual_glycameic_4_haploreg.txt", col.names = TRUE, row.names = FALSE, quote=FALSE, sep = "\t")
