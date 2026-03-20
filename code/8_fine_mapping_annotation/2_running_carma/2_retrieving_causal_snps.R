##############
#INTRODUCTION#
##############

#This code retrieves which proxies from our lead SNPs present PIP>0.1.
#These proxies will be used to make gene prioritization of the novel signals.

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)

##############
#Loading data#
##############

path_2_input <- "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/"

setwd(path_2_input)

ir_variants <- fread("output/4_ir_loci_discovery/3_novel_variants/4_colocalization/352_variants_with_all_info_plus_fiadjbmi_and_tg_hdl_coloc.txt")
ir_new <- ir_variants[which(ir_variants$reported_ir != ""),] #282 IR variants, that is all we need

#########################################
#Let's get the loci and see what happens#
#########################################

ir_new$start_ <- as.numeric(ir_new$base_pair_location)-500000
ir_new$end_ <- as.numeric(ir_new$base_pair_location)+500000

ir_new$start_ <- ifelse(as.numeric(ir_new$start_) <0, 0, ir_new$start_)

#Let's just check...

ir_new$loci <- paste(ir_new$chromosome, "_", ir_new$start_, "_", ir_new$end_, sep="")

###################################################
#Let's get the SNPs with high PIP for each loci!!!#
###################################################

dir.create("output/8_fine_mapping/all_loci/carma_res_clean/")

for(loci_index in seq(1, length(ir_new$loci))){
  
  print(loci_index)
  
  #STEP 1: let's get the variant info:
  
  loci <- ir_new$loci[loci_index]
  lead <- ir_new$variant[loci_index]
  lead_chrpos <- ir_new$chr_pos[loci_index]
  
  #STEP 2: let's get data needed for fine-mapping loaded:
  
  path_2_res <- paste("output/8_fine_mapping/all_loci/carma_res/", loci, ".RDS", sep="")
  path_2_ss <- paste("output/8_fine_mapping/all_loci/loci_ss_aligned/", loci, ".txt", sep="")
  path_2_ld <- paste("output/8_fine_mapping/all_loci/ld_matrices/", loci, ".RDS", sep="")
  
  carma_res <- readRDS(path_2_res)
  ss_df <- fread(path_2_ss)
  ld_df <- as.data.frame(readRDS(path_2_ld))
  
  #STEP 3: get the variants with PIP > 0.1
  
  #index_pips <- which(as.numeric(unlist(carma_res[[1]][1])) > 0.1)
  
  pips <- as.numeric(unlist(carma_res[[1]][1]))
  causal_df <- ss_df
  
  causal_df$pip <- pips 
  
  #STEP 4: get the r2 between these variants and your SNP: 
  
  #-> the basics of this is that LD-df and ss-df have the same row order! Meaning that we can use the indexes of ss_df to extract the LD info from the LD matrix.
  
  index_variants_of_interest <- which(ss_df$variant%in%causal_df$variant)
  ld_row_lead <- ld_df[which(ss_df$variant == lead),]
  
  ld_cols <- ld_row_lead %>%
    dplyr::select(all_of(index_variants_of_interest)) #all of the variants are actually in order from smallest to largest base pair location, which makes this comparison valid. 
  
  r2_vect <- as.numeric(t(ld_cols))^2 #they were just correlations, we need r2 so we are gonna square this one up
  
  causal_df$r2_with_lead <- r2_vect #test in LDLink shows perfect replication of LD between variants. All steps above worked out
  
  #STEP 5: let's filter those with r2>0.8 and with PIP > 0.1:
  
  causal_clean_df <- causal_df[which(as.numeric(causal_df$r2_with_lead) > 0.8),]
  #causal_clean_df <- causal_clean_df[which(as.numeric(causal_clean_df$pip) > 0.1),]
  causal_clean_df$lead_variant <- lead
  
  #STEP 6 add this info the the data set:
  
  ir_new$cs[loci_index] <- paste(causal_clean_df$variant, collapse = ";")
  ir_new$cs_variant_pips[loci_index] <- paste(causal_clean_df$pip, collapse = ";")
  
  #STEP 7: let's save the raw data because this will help us running cS2G for the variants 

  output_path <- paste("output/8_fine_mapping/all_loci/carma_res_clean/", loci, ".txt", sep="")
  
  fwrite(causal_clean_df, output_path)
  
  #And now get the proxies:
  
  # causal_df_proxies_lead <- causal_df[which(causal_df$variant%in%proxies$rsID),]
  # 
  # print(causal_df_proxies_lead)
  # 
  # if(length(causal_df_proxies_lead$variant) == 0){
  #   
  #   next()
  #   
  # } else {
  #   
  #   output_path <- paste("output/5_comparing_fiadjbmi_fi/7_fine_mapping/carma_res/carma_res_clean/", loci, "_lead_proxies.txt", sep="")
  #   fwrite(causal_df_proxies_lead, output_path)
  #   
  #   
  # }
  
}

#####################################
#Let's save the annotated variants!!#
#####################################

dir.create("output/8_fine_mapping/all_loci/annotated_data")

fwrite(ir_new, "output/8_fine_mapping/all_loci/seed_4_annotation.txt")
