##############
#INTRODUCTION#
##############

#This code take the variants in LAMB1, PLAUR, KLK1 and INPP5a loci and checks their associations in 3 different GWAS.

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)

###################
#Loading functions#
###################

##############
#Loading data#
##############

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/")

ir_variants = fread("manuscript/supplementary_material/supplementary_tables/drafts/supplementary_table_2.csv")

#Let's recover them with proxies:

proxies <- fread("output/4_ir_loci_discovery/4_proxies/proxies_build_37_and_38_4_282_variants.txt")

yes_vect <- c("A", "G", "C", "T")

proxies <- proxies[which(proxies$query_snp_rsid%in%ir_variants$variant),] #3188
proxies <- proxies[which(proxies$ref%in%yes_vect & proxies$alt%in%yes_vect),]

length(which(duplicated(proxies$query_snp_rsid) == FALSE)) #282

#Loading the GWAS:

clamp=fread("review/review_raw_data/meta_1000g_risc_ulsam_eugene2_stanford_invnorm_mac5_stderr_bmi_280213.1.txt")
isiadjbmi=fread("output/1_curated_gwas/isiadjbmi_curated.txt")
ifcadjbmi=fread("output/1_curated_gwas/ifcadjbmi_curated.txt")

##############################################
#Let's get the SNPs in the format of interest#
##############################################

test_variants = ir_variants[which(ir_variants$variant%in%c("rs4760", "rs1133400", "rs485186", "rs4727695")),]

clamp_match = clamp[which(clamp$MarkerName%in%test_variants$chr_pos),]
isi_match = isiadjbmi[which(isiadjbmi$chr_pos%in%test_variants$chr_pos),]
ifc_match = ifcadjbmi[which(ifcadjbmi$chr_pos%in%test_variants$chr_pos),]

#Let's align and add the data:

clamp_match = clamp_match[order(match(clamp_match$MarkerName, test_variants$chr_pos)),]
print(length(which(clamp_match$MarkerName == test_variants$chr_pos)))

clamp_match$Allele1 = toupper(clamp_match$Allele1)
clamp_match$Allele2 = toupper(clamp_match$Allele2)

clamp_match$effect_allele = ifelse(clamp_match$Allele1 != test_variants$effect_allele, clamp_match$Allele2, clamp_match$Allele1)
clamp_match$other_allele = ifelse(clamp_match$Allele1 != test_variants$effect_allele, clamp_match$Allele1, clamp_match$Allele2)
clamp_match$beta = ifelse(clamp_match$Allele1 != test_variants$effect_allele, -1*as.numeric(clamp_match$Effect), as.numeric(clamp_match$Effect))

test_variants$beta.m_value = clamp_match$beta
test_variants$standard_error.m_value = clamp_match$StdErr
test_variants$p_value.m_value = clamp_match$`P-value`
test_variants$sample_size.m_value = clamp_match$N

#Now for ISIadjBMI:

isi_match = isi_match[order(match(isi_match$chr_pos, test_variants$chr_pos)),]
print(length(which(isi_match$chr_pos == test_variants$chr_pos)))

isi_match$effect_allele_ = ifelse(isi_match$effect_allele != test_variants$effect_allele, isi_match$other_allele, isi_match$effect_allele)
isi_match$other_allele_ = ifelse(isi_match$effect_allele != test_variants$effect_allele, isi_match$effect_allele, isi_match$other_allele)
isi_match$beta_ = ifelse(isi_match$effect_allele != test_variants$effect_allele, -1*as.numeric(isi_match$beta), as.numeric(isi_match$beta))

test_variants$beta.isiadjbmi = isi_match$beta_
test_variants$standard_error.isiadjbmi = isi_match$standard_error
test_variants$p_value.isiadjbmi = isi_match$p_value
test_variants$sample_size.isiadjbmi = isi_match$sample_size

#Now for IFCadjBMI:

ifc_match = ifc_match[order(match(ifc_match$chr_pos, test_variants$chr_pos)),]
print(length(which(ifc_match$chr_pos == test_variants$chr_pos)))

ifc_match$effect_allele_ = ifelse(ifc_match$effect_allele != test_variants$effect_allele, ifc_match$other_allele, ifc_match$effect_allele)
ifc_match$other_allele_ = ifelse(ifc_match$effect_allele != test_variants$effect_allele, ifc_match$effect_allele, ifc_match$other_allele)
ifc_match$beta_ = ifelse(ifc_match$effect_allele != test_variants$effect_allele, -1*as.numeric(ifc_match$beta), as.numeric(ifc_match$beta))

test_variants$beta.ifcadjbmi = ifc_match$beta_
test_variants$standard_error.ifcadjbmi = ifc_match$standard_error
test_variants$p_value.ifcadjbmi = ifc_match$p_value
test_variants$sample_size.ifcadjbmi = ifc_match$sample_size

###########################################
#Let's filter fo the SNPs and for the data#
###########################################

final_df = test_variants %>%
  dplyr::select(variant, chromosome, base_pair_location, effect_allele, other_allele, minimum_allele_frequency, beta.fiadjbmi, standard_error.fiadjbmi, p_value.fiadjbmi, sample_size.fiadjbmi, 
                beta.m_value, standard_error.m_value, p_value.m_value, sample_size.m_value,
                beta.isiadjbmi, standard_error.isiadjbmi, p_value.isiadjbmi, sample_size.isiadjbmi,
                beta.ifcadjbmi, standard_error.ifcadjbmi, p_value.ifcadjbmi, sample_size.ifcadjbmi)

fwrite(final_df, "review/manuscript/tables/isi_clamp_based_associations_4_main_targets.csv")                
