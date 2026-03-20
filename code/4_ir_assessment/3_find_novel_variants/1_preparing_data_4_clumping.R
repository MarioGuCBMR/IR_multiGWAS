##############
#INTRODUCTION#
##############

#This code performs CPASSOC with our three traits!!

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)
library(LDlinkR)

###################
#Loading functions#
###################

formatting_data_4_clumping <- function(ss_gw){
  #This function formats the data for the clumping.
  
  #STEP 2: now make the SNP column which is SNP:A1:A2
  
  rsid <- ifelse(str_detect(ss_gw$variant, "rs") == TRUE, ss_gw$variant, ss_gw$chr_pos)
  SNP <- paste(rsid, ":", ss_gw$other_allele, ":", ss_gw$effect_allele, sep = "")
  
  #Now we select the columns that we want:
  
  ss_gw$SNP <- SNP
  ss_gw$rsid <- rsid
  
  final_df <- ss_gw %>%
    select(SNP, p_value, chr_pos, effect_allele, other_allele, rsid)
  
  colnames(final_df) <- c("SNP", "pval", "chr_pos", "effect_allele", "other_allele", "rsid")
  
  return(final_df)
  
}

##############
#Loading data#
##############

#Let's change the working directory:

path_2_files <-  "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_functional_annotation//" #change it with your own path.

setwd(path_2_files)

ir_variants <- fread("../IR_GSEM_2023/output/5_comparing_fiadjbmi_fi/1_data_following_ir_dir/282_ir_fiadjbmi_hdl_tg_dwls_snps.txt")
ir_gwas <- fread("../IR_GSEM_2023/output/3_mv_gwas/fiadjbmi_hdl_tg/")

##################################
#Let's load all the previous loci#
##################################

lotta <- as.data.frame(readxl::read_excel("raw_data/previous_loci/lotta_et_al.xlsx"))

magic <- as.data.frame(readxl::read_excel("raw_data/previous_loci/magic.xlsx"))
magic <- magic[which(duplicated(magic$rsID) == FALSE),]

oliveri <- as.data.frame(readxl::read_excel("raw_data/previous_loci/oliveri_et_al.xlsx"))
oliveri$SNP <- paste(oliveri$rsID, ":", oliveri$OA, ":", oliveri$EA, sep ="")
oliveri$chr_pos <- paste("chr", oliveri$CHR, ":", oliveri$`POS (hg19/b37)`, sep ="")

deforest <- as.data.frame(readxl::read_excel("raw_data/previous_loci/deforest_et_al.xlsx"))
deforest$SNP <- paste(deforest$`Lead SNP`, ":", deforest$`Other Allele`, ":", deforest$`Effect Allele`, sep ="")
deforest$`Chr:Pos`[which(is.na(deforest$`Chr:Pos`))] <- "1:27021913"
deforest$chr_pos <- paste("chr", deforest$`Chr:Pos`, sep ="")

suzuki <- as.data.frame(readxl::read_excel("raw_data/previous_loci/suzuki_et_al.xlsx"))
suzuki$chr_pos <- paste("chr", suzuki$Chromosome, ":", suzuki$`Position     (bp, b37)`, sep ="")

######################################################
#Let's format each of the data so that we can combine#
######################################################

ir_4_clumping <- formatting_data_4_clumping(ir_variants)

#We will need the alleles, which might be missing from some datasets.
#To properly perform the matching we will:

#1) load the data from FIadjBMI for Lotta
#2) Oliveri has all data, that is no problem
#3) Deforest has all data - no problem.
#4) Load T2D from Suzuki:

fiadjbmi <- fread("../IR_GSEM_2023/output/1_curated_data/fiadjbmi_curated.txt")
t2d <- fread("../../Team projects/MariaJose&Mario&Raquel/retrieving_m6A_T2D_variants/output/1_curated_data/t2d_curated.txt")

#Let's do the MAGIC ones first

fiadjbmi_lotta <- fiadjbmi[which(fiadjbmi$variant%in%lotta$SNP),]
fiadjbmi_magic <- fiadjbmi[which(fiadjbmi$variant%in%magic$rsID),] #127, two in HLA another one is rare. SKIP.

lotta_4_clumping <- formatting_data_4_clumping(fiadjbmi_lotta) #this is good
magic_4_clumping <- formatting_data_4_clumping(fiadjbmi_magic) #this is good

#Let's only format Oliveri and DeForest

oliveri_4_clumping <- oliveri %>%
  select(SNP, P, chr_pos, EA, OA, rsID)

colnames(oliveri_4_clumping) <- colnames(lotta_4_clumping)

deforest_4_clumping <- deforest %>%
  select(SNP, `P-value (corrected)`, chr_pos, `Effect Allele`, `Other Allele`, `Lead SNP`)

colnames(deforest_4_clumping) <- colnames(lotta_4_clumping)

#Finally, let's get Suzuki...

t2d_match <- t2d[which(t2d$chr_pos%in%suzuki$chr_pos),] #1197
t2d_mismatch <- suzuki[which(!(suzuki$chr_pos%in%t2d$chr_pos)),] #removed those rare in EUR - amazing

suzuki_match <- suzuki[which(suzuki$chr_pos%in%t2d$chr_pos),]
suzuki_ordered <- suzuki_match[order(match(suzuki_match$chr_pos, t2d_match$chr_pos)),]

length(which(suzuki_ordered$chr_pos == t2d_match$chr_pos)) #perfect match

suzuki_ordered$effect_allele <- t2d_match$effect_allele
suzuki_ordered$other_allele <- t2d_match$other_allele
suzuki_ordered$p_value <- t2d_match$other_allele
suzuki_ordered$SNP <- paste(suzuki_ordered$`Index SNV`, ":", suzuki_ordered$other_allele, ":", suzuki_ordered$effect_allele, sep = "")

t2d_4_clumping <- suzuki_ordered %>%
  select(SNP, p_value, chr_pos, effect_allele, other_allele, `Index SNV`)

colnames(t2d_4_clumping) <- colnames(lotta_4_clumping)

#######################################
#Let's format the data for all of them#
#######################################

ir_4_clumping$source <- "mvGWAS"
lotta_4_clumping$source <- "Lotta"
magic_4_clumping$source <- "MAGIC"
oliveri_4_clumping$source <- "Oliveri"
deforest_4_clumping$source <- "DeForest"
t2d_4_clumping$source <- "Suzuki_t2d"

#We are gonna set the p-values of all the rest of the variants to a very low threshold so that they are always chosen. 
#What we care if our SNPs are not correlated with our new variants at all

lotta_4_clumping$pval <- 5e-300
magic_4_clumping$pval <- 5e-300
oliveri_4_clumping$pval <- 5e-300
deforest_4_clumping$pval <- 5e-300
t2d_4_clumping$pval <- 5e-300

##########################################################
#Let's combine the results in pairs and assess duplicates#
##########################################################

ir_lotta <- rbind(ir_4_clumping, lotta_4_clumping)
ir_magic <- rbind(ir_4_clumping, magic_4_clumping)
ir_oliveri <- rbind(ir_4_clumping, oliveri_4_clumping)
ir_deforest <- rbind(ir_4_clumping, deforest_4_clumping)
ir_t2d <- rbind(ir_4_clumping, t2d_4_clumping)

#Let's order by p-value and remove duplicates:

ir_lotta_ordered <- ir_lotta[order(as.numeric(ir_lotta$pval)),]
ir_magic_ordered <- ir_magic[order(as.numeric(ir_magic$pval)),]
ir_oliveri_ordered <- ir_oliveri[order(as.numeric(ir_oliveri$pval)),]
ir_deforest_ordered <- ir_deforest[order(as.numeric(ir_deforest$pval)),]
ir_t2d_ordered <- ir_t2d[order(as.numeric(ir_t2d$pval)),]

ir_lotta_non_dupl <- ir_lotta_ordered[which(duplicated(ir_lotta_ordered$rsid) == FALSE)] #done with RSID to remove tri-allelics. We choose the best one.
ir_magic_non_dupl <- ir_magic_ordered[which(duplicated(ir_magic_ordered$rsid) == FALSE)] #done with RSID to remove tri-allelics. We choose the best one.

ir_oliveri_non_dupl <- ir_oliveri_ordered[which(duplicated(ir_oliveri_ordered$rsid) == FALSE)] #done with RSID to remove tri-allelics. We choose the best one.
ir_deforest_non_dupl <- ir_deforest_ordered[which(duplicated(ir_deforest_ordered$rsid) == FALSE)] #done with RSID to remove tri-allelics. We choose the best one.

ir_t2d_non_dupl <- ir_t2d_ordered[which(duplicated(ir_t2d_ordered$rsid) == FALSE)] #done with RSID to remove tri-allelics. We choose the best one.

################################
#Now let's load the data for IR#
################################

dir.create("output/1_novel_ir")
dir.create("output/1_novel_ir/1_data_4_clumping")
dir.create("output/1_novel_ir/2_clumped_data")

fwrite(ir_lotta_non_dupl, "output/1_novel_ir/1_data_4_clumping/ir_vs_lotta.txt", col.names = TRUE, row.names = FALSE, quote=FALSE, sep = "\t")
fwrite(ir_magic_non_dupl, "output/1_novel_ir/1_data_4_clumping/ir_vs_magic.txt", col.names = TRUE, row.names = FALSE, quote=FALSE, sep = "\t")
fwrite(ir_oliveri_non_dupl, "output/1_novel_ir/1_data_4_clumping/ir_vs_oliveri.txt", col.names = TRUE, row.names = FALSE, quote=FALSE, sep = "\t")
fwrite(ir_deforest_non_dupl, "output/1_novel_ir/1_data_4_clumping/ir_vs_deforest.txt", col.names = TRUE, row.names = FALSE, quote=FALSE, sep = "\t")
fwrite(ir_t2d_non_dupl, "output/1_novel_ir/1_data_4_clumping/ir_vs_t2d.txt", col.names = TRUE, row.names = FALSE, quote=FALSE, sep = "\t")

