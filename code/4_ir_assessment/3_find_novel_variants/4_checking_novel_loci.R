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

path_2_files <-  "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025//" #change it with your own path.

setwd(path_2_files)

ir_variants <- fread("output/4_ir_loci_discovery/2_ir_variants/282_ir_fiadjbmi_hdl_tg_dwls_snps_w_cpassoc.txt")

non_ir_variants <- ir_variants[which(ir_variants$ir_direction == "no"),]
ir_variants <- ir_variants[which(ir_variants$ir_direction == "yes"),]

##################################
#Let's load all the previous loci#
##################################

lotta <- fread("output/4_ir_loci_discovery/3_novel_variants/2_clumped_data/clumped_data/ir_vs_lotta_clumped.txt")
magic <- fread("output/4_ir_loci_discovery/3_novel_variants/2_clumped_data/clumped_data/ir_vs_magic_clumped.txt")
oliveri <- fread("output/4_ir_loci_discovery/3_novel_variants/2_clumped_data/clumped_data/ir_vs_oliveri_clumped.txt")
deforest <- fread("output/4_ir_loci_discovery/3_novel_variants/2_clumped_data/clumped_data/ir_vs_deforest_clumped.txt")
suzuki <- fread("output/4_ir_loci_discovery/3_novel_variants/2_clumped_data/clumped_data/ir_vs_t2d_clumped.txt")

lotta_mvgwas <- lotta[which(lotta$source == "mvGWAS"),]
magic_mvgwas <- magic[which(magic$source == "mvGWAS"),]
oliveri_mvgwas <- oliveri[which(oliveri$source == "mvGWAS"),]
deforest_mvgwas <- deforest[which(deforest$source == "mvGWAS"),]
t2d_mvgwas <- suzuki[which(suzuki$source == "mvGWAS"),]

ir_variants$reported_ir <- "no"

ir_variants$reported_ir <- ifelse(!(ir_variants$variant%in%lotta_mvgwas$rsid), "yes", ir_variants$reported_ir)
ir_variants$reported_ir <- ifelse(!(ir_variants$variant%in%magic_mvgwas$rsid), "yes", ir_variants$reported_ir)
ir_variants$reported_ir <- ifelse(!(ir_variants$variant%in%oliveri_mvgwas$rsid), "yes", ir_variants$reported_ir)
ir_variants$reported_ir <- ifelse(!(ir_variants$variant%in%deforest_mvgwas$rsid), "yes", ir_variants$reported_ir)

ir_variants$ir_source <- ""

ir_variants$ir_source <- ifelse(!(ir_variants$variant%in%lotta_mvgwas$rsid), paste(ir_variants$ir_source, "/", "Lotta", sep = ""), ir_variants$ir_source)
ir_variants$ir_source <- ifelse(!(ir_variants$variant%in%magic_mvgwas$rsid), paste(ir_variants$ir_source, "/", "MAGIC", sep = ""), ir_variants$ir_source)
ir_variants$ir_source <- ifelse(!(ir_variants$variant%in%oliveri_mvgwas$rsid), paste(ir_variants$ir_source, "/", "Oliveri", sep = ""), ir_variants$ir_source)
ir_variants$ir_source <- ifelse(ir_variants$p_value.tg_hdl < 5e-08, paste(ir_variants$ir_source, "/", "Oliveri_gw", sep = ""), ir_variants$ir_source)
ir_variants$ir_source <- ifelse(!(ir_variants$variant%in%deforest_mvgwas$rsid), paste(ir_variants$ir_source, "/", "DeForest", sep = ""), ir_variants$ir_source)

table(ir_variants$ir_source) #73

ir_variants$t2d <- ifelse(!(ir_variants$variant%in%t2d_mvgwas$rsid), "yes", "no")

table(ir_variants$t2d) #

#Let's add one more 

######################################################
#Let's check whether we can retrieve the cluster data# 
######################################################

for(chr_ in seq(1,22)){
  
  tmp_df <- fread(paste("output/4_ir_loci_discovery/3_novel_variants/2_clumped_data/clumps_of_t2d/independent_chr", chr_, ".clumped", sep = ""))
  
  if(!(exists("t2d_clusters_df"))){
    
    t2d_clusters_df <- tmp_df
    
  } else {
    
    t2d_clusters_df <- rbind(t2d_clusters_df, tmp_df)
    
  }
  
} #1076

#This data's format is a bit weird, but we can retrieve the T2D loci through the p-value:

t2d_hits <- t2d_clusters_df[which(t2d_clusters_df$P == 5e-300),] #961 - perfect

t2d_hits$chr_pos <- paste("chr", t2d_hits$CHR, ":", t2d_hits$BP, sep = "")

###################################################
#Let's get the original data to obtain the proxies#
###################################################

#STEP 1: retrieve the data

suzuki <- as.data.frame(readxl::read_excel("raw_data/previous_loci/suzuki_et_al.xlsx"))
suzuki$chr_pos <- paste("chr", suzuki$Chromosome, ":", suzuki$`Position     (bp, b37)`, sep ="")

#STEP 2: obtain all summary statistics:

t2d <- fread("../../Team projects/MariaJose&Mario&Raquel/retrieving_m6A_T2D_variants/output/1_curated_data/t2d_curated.txt")

t2d_match <- t2d[which(t2d$chr_pos%in%suzuki$chr_pos),] #1197
t2d_mismatch <- suzuki[which(!(suzuki$chr_pos%in%t2d$chr_pos)),] #removed those rare in EUR - amazing

#STEP 3: add the info:

suzuki_match <- suzuki[which(suzuki$chr_pos%in%t2d$chr_pos),]
suzuki_ordered <- suzuki_match[order(match(suzuki_match$chr_pos, t2d_match$chr_pos)),]

length(which(suzuki_ordered$chr_pos == t2d_match$chr_pos)) #perfect match

suzuki_ordered$effect_allele <- t2d_match$effect_allele
suzuki_ordered$other_allele <- t2d_match$other_allele
suzuki_ordered$p_value <- t2d_match$other_allele

suzuki_clumped <- suzuki_ordered[which(suzuki_ordered$chr_pos%in%t2d_hits$chr_pos),] #perfect match

suzuki_clumped_ordered <- suzuki_clumped[order(match(suzuki_clumped$chr_pos, t2d_hits$chr_pos)),]

length(which(suzuki_clumped_ordered$chr_pos == t2d_hits$chr_pos))

t2d_hits$cluster <- suzuki_clumped_ordered$`Cluster assignment`

##########################################################
#Finally, let's loop our variants according to each hit!!#
##########################################################

ir_variants$t2d_cluster <- NA
ir_variants$t2d_snp <- NA
ir_variants$t2d_snp_r2 <- NA

for(index in seq(1, length(ir_variants$variant))){
  
  #STEP 1: retreive rsID
  
  rsid <- ir_variants$variant[index]
  
  #STEP 2: search for data:
  
  index_match <- which(str_detect(t2d_hits$SP2, rsid) | str_detect(t2d_hits$SNP, rsid))
  
  if(is_empty(index_match)){
    
    next()
    
  } else {
    
    print(index)
    
    if(t2d_hits$SP2[index_match] == "NONE"){ #the variant is a lead T2D SNP
      
      print(t2d_hits[index_match,])
      
      ir_variants$t2d_snp[index] = rsid
      ir_variants$t2d_snp_r2[index] = 1
      ir_variants$t2d_cluster[index] <- t2d_hits$cluster[index_match]
      
    } else {
      
      proxy <- unlist(str_split(t2d_hits$SNP[index_match],":"))[1]
      
      r2 <- LDlinkR::LDpair(rsid, proxy, pop = "EUR", token="04cad4ca4374")
      
      ir_variants$t2d_snp[index] = proxy
      ir_variants$t2d_snp_r2[index] = r2$r2
      ir_variants$t2d_cluster[index] <- t2d_hits$cluster[index_match]
      
    }
    
  }
  
} #we recovered them all!!!

######################################################
#Let's save the data, we are gonna use it right now!!#
######################################################

fwrite(ir_variants, "output/4_ir_loci_discovery/3_novel_variants/3_novel_variants/282_ir_variants_with_ir_and_t2d_comparisons.txt")

