##############
#INTRODUCTION#
##############

#This is from consensus adipose tissue ATAC-seq data from Perrin et al.

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)

#####################
#Let's read the data#
#####################

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/")

ir_variants <- fread("output/4_ir_loci_discovery/3_novel_variants/3_novel_variants/282_ir_variants_with_ir_and_t2d_comparisons_w_closest_genes.txt")
proxies <- fread("output/4_ir_loci_discovery/4_proxies/proxies_build_37_and_38_4_282_variants.txt")

#########################################
#Let's format and save the lead variants#
#########################################

leads <- ir_variants %>%
  select(variant, chromosome, base_pair_location)

colnames(leads) <- c("SNP","Chrom","BP")
leads$Chrom <- paste0("chr",leads$Chrom)
fwrite(leads, "output/5_enrichment_analyses/3_go_shifter/input.txt", sep = "\t")

#############################################
#Let's filter the data with different values#
#############################################

proxies_lead_info_37 <- left_join(proxies, leads, by=c("query_snp_rsid"="SNP"))
proxies_lead_info_37$chr_37 <- paste0("Chr",proxies_lead_info_37$chr)
proxies_lead_info_37 <- proxies_lead_info_37 %>%
  mutate(Chrom = str_replace(Chrom, "^chr", "Chr"))
proxies_lead_info_37 <- proxies_lead_info_37 %>%
  mutate(distance = abs(BP - pos_hg19))
proxies_lead_info_37$chr <- paste("Chr", proxies_lead_info_37$chr, sep = "")

###########################################################################
#Get the D', which I forgot to add cuz we do not need it anywhere but here#
###########################################################################

proxies_og <- fread("output/4_ir_loci_discovery/4_proxies/282_proxies_core_15.txt", fill =TRUE)

proxies_og_match <- proxies_og[which(proxies_og$rsID%in%proxies_lead_info_37$rsID),]

proxies_og_match <- proxies_og_match[order(match(proxies_og_match$rsID, proxies_lead_info_37$rsID)),]

length(which(proxies_lead_info_37$rsID == proxies_og_match$rsID)) #perfect mtch

proxies_lead_info_37$"D'" <- proxies_og_match$`D'`

proxies_lead_info_37 <- proxies_lead_info_37 %>% dplyr::select("Chrom","BP","query_snp_rsid","chr","pos_hg19","rsID","distance","r2","D'")
colnames(proxies_lead_info_37) <- c("ChromA","PosA","RsIdA","ChromB","PosB","RsIdB","Distance","RSquared","Dprime") 
fwrite(proxies_lead_info_37, "output/5_enrichment_analyses/3_go_shifter/ld.txt", sep = "\t")

#####################################################################################
#Let's divide the data according to the clusters and see if we have any differences!#
#####################################################################################

bmi_match <- fread("output/5_enrichment_analyses/1_prs/1_expo_outcome_df/ir_variants_BMI.txt")
#gsat_match <- fread("output/5_enrichment_analyses/1_prs/1_expo_outcome_df/ir_variants_GSATadjBMI.txt") #some weird outlier...

#Let's divide them accordingly:

bmi_neg <- bmi_match[which(as.numeric(bmi_match$beta.outcome) < 0 & as.numeric(bmi_match$pval.outcome) < 0.05),]
bmi_pos <- bmi_match[which(as.numeric(bmi_match$beta.outcome) > 0 & as.numeric(bmi_match$pval.outcome) < 0.05),]
bmi_no_sign <- bmi_match[which(as.numeric(bmi_match$pval.outcome) > 0.05),]

###########################################
#Let's first divide the lead dataframes...#
###########################################

lead_bmi_ns <- leads[which(leads$SNP%in%bmi_no_sign$SNP),] #141
lead_bmi_neg <- leads[which(leads$SNP%in%bmi_neg$SNP),] #63
lead_bmi_pos <- leads[which(leads$SNP%in%bmi_pos$SNP),] #78

#Let's save them:

fwrite(lead_bmi_ns, "output/5_enrichment_analyses/3_go_shifter/input_bmi_ns.txt", sep = "\t")
fwrite(lead_bmi_neg, "output/5_enrichment_analyses/3_go_shifter/input_bmi_neg.txt", sep = "\t")
fwrite(lead_bmi_pos, "output/5_enrichment_analyses/3_go_shifter/input_bmi_pos.txt", sep = "\t")

###############################
#Let's next divide the proxies#
###############################

#RsIDA has the lead info

proxies_lead_bmi_ns <- proxies_lead_info_37[which(proxies_lead_info_37$RsIdA%in%lead_bmi_ns$SNP),] 
proxies_lead_bmi_neg <- proxies_lead_info_37[which(proxies_lead_info_37$RsIdA%in%lead_bmi_neg$SNP),] 
proxies_lead_bmi_pos <- proxies_lead_info_37[which(proxies_lead_info_37$RsIdA%in%lead_bmi_pos$SNP),] 

#Let's check

length(unique(proxies_lead_bmi_ns$RsIdA)) #141 leads
length(unique(proxies_lead_bmi_neg$RsIdA)) #63 leads
length(unique(proxies_lead_bmi_pos$RsIdA)) #63 leads

#Works out!

fwrite(proxies_lead_bmi_ns, "output/5_enrichment_analyses/3_go_shifter/ld_bmi_ns.txt", sep = "\t")
fwrite(proxies_lead_bmi_neg, "output/5_enrichment_analyses/3_go_shifter/ld_bmi_neg.txt", sep = "\t")
fwrite(proxies_lead_bmi_pos, "output/5_enrichment_analyses/3_go_shifter/ld_bmi_pos.txt", sep = "\t")

