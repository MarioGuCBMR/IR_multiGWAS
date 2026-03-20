##############
#INTRODUCTION#
##############

#This code will prepare FIadjBMI, HDL and TG data necessary to run CPASSOC.

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)
library(bigsnpr)

#memory.limit(size=80000000)

##############
#Loading data#
##############

path_2_files <- "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/"

setwd(path_2_files)

fi_curated <- fread("output/1_curated_gwas/fiadjbmi_curated_his.txt")
hdl_curated <- fread("output/1_curated_gwas/hdl_curated_his.txt")
tg_curated <- fread("output/1_curated_gwas/tg_curated_his.txt")

###############################################################
#Let's merge the data aligning all to FIadjBMI-increasing SNPs#
###############################################################

#Now let's get the common SNPs:

common_snps <- Reduce(intersect, list(fi_curated$chr_pos, hdl_curated$chr_pos, tg_curated$chr_pos)) 

fi_common <- fi_curated[which(fi_curated$chr_pos%in%common_snps),]

#Let's align these fellas:

new_a1 <- ifelse(as.numeric(fi_common$beta) < 0, fi_common$other_allele, fi_common$effect_allele)
new_a2 <- ifelse(as.numeric(fi_common$beta) < 0, fi_common$effect_allele, fi_common$other_allele)
new_beta <- ifelse(as.numeric(fi_common$beta) < 0, as.numeric(fi_common$beta)*(-1), as.numeric(fi_common$beta))
new_eaf <- ifelse(as.numeric(fi_common$beta) < 0, 1-as.numeric(fi_common$effect_allele_frequency), as.numeric(fi_common$effect_allele_frequency))

fi_common$final_effect_allele <- new_a1
fi_common$final_other_allele <- new_a2
fi_common$final_beta <- new_beta
fi_common$final_effect_allele_frequency <- new_eaf

head(fi_common)

summary(fi_common$final_beta) #all positive! Great.

fi_common_pos <- fi_common %>%
  select("chromosome","base_pair_location", "final_effect_allele", "final_other_allele", "final_effect_allele_frequency", "final_beta", "standard_error", "p_value", "sample_size", "chr_pos")

colnames(fi_common_pos) <- c("chr","pos", "a1", "a0", "effect_allele_frequency", "beta", "standard_error", "p_value", "sample_size", "ID")

############################
#MATCHING DATA with bigsnpr#
############################

#we are going to use bigsnpr for this:

colnames(hdl_curated) <- c("variant", "chr","pos", "a1", "a0", "effect_allele_frequency", "beta", "standard_error", "p_value", "sample_size", "ID")
colnames(tg_curated) <- c("variant", "chr","pos", "a1", "a0", "effect_allele_frequency", "beta", "standard_error", "p_value", "sample_size", "ID")

#Careful chr and pos should be as.character()

fi_common_pos$chr <- as.character(fi_common_pos$chr)
fi_common_pos$pos <- as.character(fi_common_pos$pos)

hdl_curated$chr <- as.character(hdl_curated$chr)
hdl_curated$pos <- as.character(hdl_curated$pos)

tg_curated$chr <- as.character(tg_curated$chr)
tg_curated$pos <- as.character(tg_curated$pos)

#Let's get to do some matching ;)

fi_hdl_1 <- bigsnpr::snp_match(hdl_curated, fi_common_pos, strand_flip = FALSE, join_by_pos = TRUE, remove_dups = FALSE, match.min.prop=0.00)
fi_hdl_2 <- bigsnpr::snp_match(hdl_curated, fi_common_pos, strand_flip = TRUE, join_by_pos = TRUE, remove_dups = FALSE, match.min.prop=0.00)
fi_hdl_2 <- fi_hdl_2[which(!(fi_hdl_2$ID%in%fi_hdl_1$ID)),]

fi_hdl_df <- rbind(fi_hdl_1, fi_hdl_2)

#Let's now change the name of the columns and apply the same, but for TG:

fi_hdl_df_clean <- fi_hdl_df %>%
  select("variant", "ID.ss",  "chr", "pos", "a0", "a1", "effect_allele_frequency.ss",
         "beta.ss", "standard_error.ss", "p_value.ss",  "sample_size.ss",
         "beta", "standard_error", "p_value",  "sample_size")


colnames(fi_hdl_df_clean) <- c("variant", "ID",  "chr", "pos", "a0", "a1", "effect_allele_frequency.hdl",
                                     "beta.hdl", "standard_error.hdl", "p_value.hdl",  "sample_size.hdl",
                                     "beta", "standard_error", "p_value",  "sample_size")

#I think now we can do the same, but with TG:

fi_hdl_tg_1 <- bigsnpr::snp_match(tg_curated, fi_hdl_df_clean, strand_flip = FALSE, join_by_pos = TRUE, remove_dups = FALSE, match.min.prop=0.00)
fi_hdl_tg_2 <- bigsnpr::snp_match(tg_curated, fi_hdl_df_clean, strand_flip = TRUE, join_by_pos = TRUE, remove_dups = FALSE, match.min.prop=0.00)
fi_hdl_tg_2 <- fi_hdl_tg_2[which(!(fi_hdl_tg_2$ID%in%fi_hdl_tg_1$ID)),]

fi_hdl_tg_df <- rbind(fi_hdl_tg_1, fi_hdl_tg_2)

####################################################################################
#Let's clean the data just a little bit, so that we can have the data to work later#
####################################################################################

final_sumstats_df <- fi_hdl_tg_df %>%
  select("variant.ss", "ID.ss",  "chr", "pos", "a1", "a0", "effect_allele_frequency",
         "beta", "standard_error", "p_value",  "sample_size",
         "beta.hdl", "standard_error.hdl", "p_value.hdl",  "sample_size.hdl",
         "beta.ss", "standard_error.ss", "p_value.ss",  "sample_size.ss")

colnames(final_sumstats_df) <- c("variant", "chr_pos",  "chr", "pos", "effect_allele", "other_allele", "effect_allele_frequency",
                                 "beta.fiadjbmi", "standard_error.fiadjbmi", "p_value.fiadjbmi",  "sample_size.fiadjbmi",
                                 "beta.hdl", "standard_error.hdl", "p_value.hdl",  "sample_size.hdl",
                                 "beta.tg", "standard_error.tg", "p_value.tg",  "sample_size.tg")

#I tripled checked directions - they all make sense
#Take best SNP from this df to see:

#check <- final_sumstats_df[which(final_sumstats_df$p_value.hdl < 0.005 & final_sumstats_df$p_value.fiadjbmi < 0.005 & final_sumstats_df$p_value.tg < 0.005),]

dir.create("output/4_ir_loci_discovery/1_cpassoc/multi_ancestry")
dir.create("output/4_ir_loci_discovery/1_cpassoc/multi_ancestry/HIS")

fwrite(final_sumstats_df, "output/4_ir_loci_discovery/1_cpassoc/multi_ancestry/HIS/merged_fiadjbmi_hdl_tg_sumstats_4_cpassoc.txt")

#########################################################
#Now let's perform some pruning to get this data rolling#
#########################################################

#CPASSOC cannot have causal SNPs so let's calculate the z-scores for each trait and remove those that have >1.96 z-score in absolute values:

final_sumstats_df$z_score.fiadjbmi <- as.numeric(final_sumstats_df$beta.fiadjbmi)/as.numeric(final_sumstats_df$standard_error.fiadjbmi)
final_sumstats_df$z_score.hdl <- as.numeric(final_sumstats_df$beta.hdl)/as.numeric(final_sumstats_df$standard_error.hdl)
final_sumstats_df$z_score.tg <- as.numeric(final_sumstats_df$beta.tg)/as.numeric(final_sumstats_df$standard_error.tg)

sumstats_4_pruning <- final_sumstats_df[which(abs(final_sumstats_df$z_score.fiadjbmi) < 1.96 & abs(final_sumstats_df$z_score.hdl) < 1.96 & abs(final_sumstats_df$z_score.tg) < 1.96),]

#Let's see:

summary(abs(sumstats_4_pruning$z_score.fiadjbmi))
summary(abs(sumstats_4_pruning$z_score.hdl))
summary(abs(sumstats_4_pruning$z_score.tg))

fwrite(sumstats_4_pruning, "output/4_ir_loci_discovery/1_cpassoc/multi_ancestry/HIS/merged_fiadjbmi_hdl_tg_sumstats_w_maxZ_1.96.txt")

#sumstats_4_pruning <- fread("output/2_variant_discovery/1_triangulation/4_data_4_cpassoc_and_curated_results/merged_fiadjbmi_hdl_tg_sumstats_w_maxZ_1.96.txt")

#And now let's finally prune the shit out of them.
#If we do ld-clumping with these set of varaiants it is going to be OK.

sumstats_4_pruning$rsid <- sumstats_4_pruning$variant
sumstats_4_pruning$pval <- 1 #this allows us to prune the data!!

#And now we can run the ld-clumping with 0.2 and now p-value threshold:

#First version for supercomputer at CBMR
sumstats_pruned <- ieugwasr::ld_clump_local(sumstats_4_pruning, bfile = "C:/Users/zlc436/Desktop/1kg.v3/AMR", plink_bin ="C:/Users/zlc436/Documents/R/R-4.1.3/library/plinkbinr/bin/plink_Windows.exe", clump_kb = 500, clump_r2 = 0.2, clump_p = 1) 

#Second version for laptop 
#sumstats_pruned <- ieugwasr::ld_clump_local(sumstats_4_pruning, bfile = "C:/Users/zlc436/Desktop/1000G_ieugwasr/1kg.v3/1kg.v3/EUR", plink_bin ="C:/Users/zlc436/Documents/R/R-4.1.2/library/plinkbinr/bin/plink_Windows.exe", clump_kb = 500, clump_r2 = 0.2, clump_p = 1) 

fwrite(sumstats_pruned, "output/4_ir_loci_discovery/1_cpassoc/multi_ancestry/HIS/pruned_fiadjbmi_hdl_tg_sumstats_w_maxZ_1.96.txt") #800k

######################################################
#I can calculate the correlation without an issue now#
######################################################

Dtcorr <- sumstats_pruned %>%
  select(z_score.fiadjbmi, z_score.hdl, z_score.tg)

fiadjbmi_hdl <- cor(Dtcorr[,1], Dtcorr[,2])
fiadjbmi_tg <- cor(Dtcorr[,1], Dtcorr[,3])
hdl_tg <- cor(Dtcorr[,2], Dtcorr[,3])

corr_vect <- 
  c(1,fiadjbmi_hdl, fiadjbmi_tg,
  fiadjbmi_hdl, 1, hdl_tg,
  fiadjbmi_tg, hdl_tg, 1)

corr_matrix <- matrix(corr_vect,nrow =  3, ncol =3)

saveRDS(corr_matrix, "output/4_ir_loci_discovery/1_cpassoc/multi_ancestry/HIS/correlation_matrix_between_traits.RDS")

