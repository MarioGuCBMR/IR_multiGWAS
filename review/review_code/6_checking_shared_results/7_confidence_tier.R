##############
#INTRODUCTION#
##############

#This code computes the confidence tier list for the multi-trait IR associations taking into account several of the approaches
#suggest by the Reviewers.
#Additionally, we will wrap up a table that combines the latest T2D results!



###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)

##############
#Loading data#
##############

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/")

ir_variants = fread("manuscript/supplementary_material/supplementary_tables/drafts/supplementary_table_2.csv")
even_sample_data = fread("review/manuscript/tables/282_ir_with_gsem_cpassoc_even_sample_pvals.txt")
t2d_data = fread("review/manuscript/tables/282_ir_with_t2d_data.txt")

#######################################################################
#Let's merge the dataframes so that things are easier for use later on#
#######################################################################

even_sample_data = even_sample_data[order(match(even_sample_data$variant, ir_variants$variant)),]
t2d_data = t2d_data[order(match(t2d_data$variant, ir_variants$variant)),]

print(length(which(even_sample_data$variant == ir_variants$variant))) #282
print(length(which(t2d_data$variant == ir_variants$variant))) #282

#Amazing! Let's add the data:

ir_variants$pvalue.gsem_even_n = even_sample_data$pvalue.gsem_even_n
ir_variants$pvalue.cpassoc_even_n = even_sample_data$pvalue.cpassoc_even_n
ir_variants$t2d_signal = t2d_data$t2d_signal
ir_variants$r2_t2d_signal = t2d_data$r2_t2d_signal

############################################################
#We have to add the Suzuki et al sumstats, so let's do that#
############################################################

suzuki=fread("output/1_curated_gwas/t2d_suzuki_curated.txt")

suzuki_match = suzuki[which(suzuki$chr_pos%in%ir_variants$chr_pos),]

suzuki_match = suzuki_match[order(match(suzuki_match$chr_pos, ir_variants$chr_pos)),]

print(length(which(suzuki_match$chr_pos == ir_variants$chr_pos))) #perfect match

suzuki_match$beta_aligned = ifelse(suzuki_match$effect_allele != ir_variants$effect_allele, as.numeric(suzuki_match$beta)*(-1), as.numeric(suzuki_match$beta))

summary(suzuki_match$beta_aligned)
summary(as.numeric(suzuki_match$p_value))

#Let's take a look at those that are significant:

length(which(as.numeric(suzuki_match$p_value) < 0.05)) #212/282 - quite good

#Let's add the info:

ir_variants$beta.t2d = as.numeric(suzuki_match$beta_aligned)
ir_variants$standard_error.t2d = as.numeric(suzuki_match$standard_error)
ir_variants$p_value.t2d = as.numeric(suzuki_match$p_value)
ir_variants$sample_size.t2d = as.numeric(suzuki_match$sample_size)

fwrite(ir_variants, "review/manuscript/tables/ir_leads_with_even_samples_and_new_t2d_data.txt")

###########################################################################################
#Let's add the tier list, of course!! that is what we wanted to do from the very beginning#
###########################################################################################

#How many are nominal significant for the hallmark traits?

ir_variants$nomimal_hallmark = ifelse(as.numeric(ir_variants$p_value.fiadjbmi) < 0.05 &
                                        as.numeric(ir_variants$p_value.hdl) < 0.05 &
                                        as.numeric(ir_variants$p_value.tg) < 0.05,
                                        "yes", "no")

table(ir_variants$nomimal_hallmark) #143 yes, 139 no

#How many are genome-wide significant for both CPASSOC and GSEM?

ir_variants$genome_wide_discovery = ifelse(as.numeric(ir_variants$p_value) < 5e-08 &
                                        as.numeric(ir_variants$pval.shet) < 5e-08,
                                      "yes", "no")

table(ir_variants$genome_wide_discovery) #169 yes, 113 no

#How many nominally associated for both CPASSOC and GSEM in replication?

ir_variants$nominal_replication = ifelse(as.numeric(ir_variants$pvalue.gsem_even_n) < 5e-02 &
                                             as.numeric(ir_variants$pvalue.cpassoc_even_n) < 5e-02,
                                           "yes", "no")

table(ir_variants$nominal_replication) #217 yes, 47 no

#How many are associated with cross-ancestry association after FDR?

ir_variants$cross_ancestry = ifelse(as.numeric(ir_variants$pval.shet_afr) < 5e-02/282 |
                                           as.numeric(ir_variants$pval.shet_eas) < 5e-02/282 |
                                           as.numeric(ir_variants$pval.shet_his) < 5e-02/282 |
                                           as.numeric(ir_variants$pval.shet_sas) < 5e-02/282,
                                         "yes", "no")

table(ir_variants$cross_ancestry) #62 yes, 185 no

#Careful we have NA in some of them. Those were not found in some of the analyses. 
#Let's change NAs with NOs.

ir_variants$nomimal_hallmark = ifelse(is.na(ir_variants$nomimal_hallmark), "no", ir_variants$nomimal_hallmark)
ir_variants$genome_wide_discovery = ifelse(is.na(ir_variants$genome_wide_discovery), "no", ir_variants$genome_wide_discovery)
ir_variants$nominal_replication = ifelse(is.na(ir_variants$nominal_replication), "no", ir_variants$nominal_replication)
ir_variants$cross_ancestry = ifelse(is.na(ir_variants$cross_ancestry), "no", ir_variants$cross_ancestry)

####################################################
#All numbers makes sense with what we had before!!!#
####################################################

ir_variants$confidence_level = 1
ir_variants$confidence_level = ifelse(ir_variants$nomimal_hallmark == "yes", ir_variants$confidence_level+1, ir_variants$confidence_level)
ir_variants$confidence_level = ifelse(ir_variants$genome_wide_discovery == "yes", ir_variants$confidence_level+1, ir_variants$confidence_level)
ir_variants$confidence_level = ifelse(ir_variants$nominal_replication == "yes", ir_variants$confidence_level+1, ir_variants$confidence_level)
ir_variants$confidence_level = ifelse(ir_variants$cross_ancestry == "yes", ir_variants$confidence_level+1, ir_variants$confidence_level)

table(ir_variants$confidence_level)

fwrite(ir_variants, "review/manuscript/tables/ir_leads_with_even_samples_and_new_t2d_data_plus_confidence_tiers.txt")

################################################################
#Let's check associations with P<0.05 for T2D in novel variants#
################################################################

check = fread("review/manuscript/tables/ir_leads_with_even_samples_and_new_t2d_data_plus_confidence_tiers.txt")
check = check[which(check$ir_source == ""),]
length(which(check$p_value.t2d < 5e-08)) #21
length(which(check$p_value.t2d < 5e-02)) #56

