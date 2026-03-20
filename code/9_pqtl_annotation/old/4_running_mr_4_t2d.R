##############
#INTRODUCTION#
##############

#This code performs MR  for all cis-pQTL instruments - t2d systematically.

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)
library(TwoSampleMR)

###################
#Loading functions#
###################

protein_parser <- function(protein_id){
  
  #STEP 1: figure out what kind is it:
  
  check <- str_detect(protein_id, "UKB_PPP")
  
  if(check){
    
    prot <- unlist(str_split(protein_id, "UKB_PPP_EUR_"))[2]
    prot <- unlist(str_split(prot, "_"))[1]
    
    return(prot)
    
  } else {
    
    prot <- unlist(str_split(protein_id, "sun_2018_aptamer_plasma_"))[2]
    prot <- unlist(str_split(prot, "_"))[1]
    
    return(toupper(prot))
    
  }
  
}


chr_pos_parser <- function(id){
  
  #id = "12:104000470:T:C:imp:v1"
  
  tmp <- unlist(str_split(id, ":"))
  
  chr_pos <- paste("chr", tmp[1], ":", tmp[2], sep = "")
  
  return(chr_pos)
  
}


load_prot <- function(clumped_data_, protein_name){
  
  #STEP 1: load the protein data:
  
  prot_data <- clumped_data_[which(clumped_data_$protein == protein_name),]
  
  #Let's get the chromosomes:
  
  chr_vect <- unique(prot_data$chr_hg19)

  #STEP 3: get file of your proteins:
  
  files_ <- list.files("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/raw_data/pqtl_data/")
  files_ <- files_[which(str_detect(files_, ".tar") == FALSE)]
  file_ <- files_[str_detect(files_, paste0("^", protein_name, "_"))]
  
  for(chr_ in chr_vect){
    
    print(chr_)
    
    pepd_tmp <- fread(paste("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/raw_data/pqtl_data/", file_, "/discovery_", paste("chr", chr_, sep = ""), "_", file_, ".gz", sep = ""))
    
    if(!(exists("pepd_end"))){
      
      pepd_end <- pepd_tmp
      
      rm(pepd_tmp)
      
    } else {
      
      pepd_end <- rbind(pepd_end, pepd_tmp)  
      
      rm(pepd_tmp)
      
    }
    
  }
  
  #STEP 4: do some slight modificaitons:
  
  pepd_end$chr_pos <- as.character(unlist(sapply(pepd_end$ID, chr_pos_parser)))
  pepd_end$p_value <- 10^-as.numeric(pepd_end$LOG10P)
  
  #STEP 5: print the data
  
  return(pepd_end)
  
}

chr_cleaner <- function(chr_pos){
  
  chr_ <- str_split(chr_pos, ":")[[1]][1]
  chr_ <- str_split(chr_, "chr")[[1]][2]
  
  return(as.numeric(chr_))
}

pos_cleaner <- function(chr_pos){
  
  pos_ <- str_split(chr_pos, ":")[[1]][2]
  
  return(as.numeric(pos_))
  
}

data_aligner <- function(query_ss, other_ss){
  
  ################################################################################
  #This code uses: exposure and proxy dataframes that should be loaded beforehand#
  ################################################################################
  
  #########################################################################################
  #STEP 0: let's run the matching with TwoSampleMR so we need the data in a certain format#
  #########################################################################################
  
  #Let's first check if we have the effect_allele_frequency column:
  
  check <- which(colnames(query_ss) == "effect_allele_frequency")
  
  if(is_empty(check)){
    
    query_ss$effect_allele_frequency <- NA
    
  }
  
  #Now we can proceed
  
  exposure <- query_ss %>%
    select(variant, chr_19, pos_19, ALLELE1, ALLELE0, A1FREQ, BETA, SE, p_value, N, chr_pos)
  
  colnames(exposure) <- c("SNP", "chr.exposure", "pos.exposure", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure", "beta.exposure", "se.exposure", "pval.exposure", "samplesize.exposure", "chr_pos.exposure")
  
  exposure$exposure <- "exposure"
  exposure$id.exposure <- "exposure"
  
  #Now with the outcome:
  
  check <- which(colnames(other_ss) == "effect_allele_frequency")
  
  if(is_empty(check)){
    
    other_ss$effect_allele_frequency <- NA
    
  }
  
  #Now we can proceed
  
  if("sample_cases"%in%colnames(other_ss)){
    
    outcome <- other_ss %>%
      select(variant, chromosome, base_pair_location, effect_allele, other_allele, effect_allele_frequency, beta, standard_error, p_value, sample_size, sample_cases, sample_controls, prevalence, chr_pos)
    
    colnames(outcome) <- c("SNP", "chr.outcome", "pos.outcome", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome", "pval.outcome", "samplesize.outcome", "ncase.outcome", "ncontrol.outcome", "prevalence.outcome", "chr_pos.outcome")
    
  } else {
    
    outcome <- other_ss %>%
      select(variant, chromosome, base_pair_location, effect_allele, other_allele, effect_allele_frequency, beta, standard_error, p_value, sample_size, chr_pos)
    
    colnames(outcome) <- c("SNP", "chr.outcome", "pos.outcome", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome", "pval.outcome", "samplesize.outcome", "chr_pos.outcome")
    
  }
  
  outcome$outcome <- "outcome"
  outcome$id.outcome <- "outcome"
  
  ############################################################################################################
  #STEP 1: match the data. This will probably fail with tri-allelic SNPs. Here I think we are gonna be OK... #
  ############################################################################################################
  
  merged_df <- TwoSampleMR::harmonise_data(exposure, outcome, action=3)
  print(dim(merged_df))
  
  merged_df <- merged_df[which(merged_df$remove == FALSE),] #removing incompatible alleles
  
  if(nrow(merged_df) == 0){ #a weird behaviour... we need to take into account
    
    merged_df$remove <- NA
    merged_df$palindromic <- NA
    merged_df$ambiguous <- NA
    merged_df$action <- NA
    merged_df$SNP_index <- NA
    merged_df$mr_keep <- NA
    
  }
  
  ######################################
  #STEP 2: do we need to query proxies?#
  ######################################
  
  missing_snps <- query_ss[which(!(query_ss$chr_pos%in%merged_df$chr_pos.exposure)),] 
  
  if(length(missing_snps$chr_pos) == 0){
    
    return(merged_df)
    
  } else {
    
    #We have to deal with proxies!!
    
    #1: check proxies for missing variants:
    
    proxies_4_missing <- proxies[which(proxies$query_snp_rsid%in%missing_snps$variant),]
    
    #2: check which are available in outcome:
    
    proxies_in_outcome <- proxies_4_missing[which(proxies_4_missing$rsID%in%other_ss$variant),]
    
    #3. Let's order the proxies by exposure data:
    
    exposure_match <- exposure_df[which(exposure_df$chr_pos_38%in%proxies_in_outcome$chr_pos_38),]
    
    #4. Let's add lead SNP so that we can properly add data:
    
    proxies_match <- proxies_in_outcome[which(proxies_in_outcome$chr_pos_38%in%exposure_match$chr_pos_38),]
    proxies_ordered <- proxies_match[order(match(proxies_match$chr_pos_38, exposure_match$chr_pos_38)),]
    
    length(which(proxies_ordered$chr_pos_38 == exposure_match$chr_pos_38))
    
    exposure_match$lead_snp <- proxies_ordered$query_snp_rsid
    exposure_match$variant <- proxies_ordered$rsID #adding RSID to have ALL INFO. Maybe the summary statistics only has CHR:POS.
    exposure_match$rsq <- proxies_ordered$r2 #adding RSID to have ALL INFO. Maybe the summary statistics only has CHR:POS.
    
    #5. Let's get only one hit per signal, retaining only the best:
    
    exposure_best <- exposure_match[order(as.numeric(exposure_match$rsq), decreasing = TRUE),]
    exposure_best <- exposure_best[which(duplicated(exposure_best$lead_snp) == FALSE),] #7!
    
    #6. Let's get the data
    
    tmp_df <- exposure_best
    
    #7. Let's get the data ready for merging...
    
    exposure_proxies <- tmp_df %>%
      select(variant, chr_19, pos_19, ALLELE1, ALLELE0, A1FREQ, BETA, SE, p_value, N, chr_pos)
    
    colnames(exposure_proxies) <- c("SNP", "chr.exposure", "pos.exposure", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure", "beta.exposure", "se.exposure", "pval.exposure", "samplesize.exposure", "chr_pos.exposure")
    
    exposure_proxies$exposure <- "exposure"
    exposure_proxies$id.exposure <- "exposure"
    
    #And let's merge with the outcome, that should have ALL the data:
    
    merged_proxies_df <- TwoSampleMR::harmonise_data(exposure_proxies, outcome, action=2)
    
    print(dim(merged_proxies_df))
    
    merged_proxies_df <- merged_proxies_df[which(merged_proxies_df$remove == FALSE),] #removing incompatible alleles
    
    if(nrow(merged_proxies_df) == 0){ #a weird behaviour... we need to take into account
      
      merged_proxies_df$remove <- NA
      merged_proxies_df$palindromic <- NA
      merged_proxies_df$ambiguous <- NA
      merged_proxies_df$action <- NA
      merged_proxies_df$SNP_index <- NA
      merged_proxies_df$mr_keep <- NA
      
    }
    
    #8. Let's add this info to the previous version:
    
    merged_df <- rbind(merged_df, merged_proxies_df)
    
    #9. Let's double-check for independence:
    
    #merged_df <- filt_ld(merged_df)
    
    return(merged_df)
    
  }
  
}


mr_plots <- function(dat)
{
  require(TwoSampleMR)
  require(ggplot2)
  require(ggrepel)
  require(dplyr)
  
  temp <- subset(dat, outcome == outcome[1] & exposure == exposure[1])
  exposure_name <- temp$exposure[1]
  outcome_name <- temp$outcome[1]
  
  if(! "labels" %in% names(dat)) dat$labels <- NA
  
  exposure_units <- temp$units.exposure[1]
  outcome_units <- temp$units.outcome[1]
  
  mrs <- mr_singlesnp(temp, all_method=c("mr_egger_regression", "mr_ivw", "mr_weighted_median", "mr_weighted_mode"))
  mrl <- mr_leaveoneout(temp)
  
  mrs$index <- 1:nrow(mrs)
  mrl$index <- 1:nrow(mrl)
  
  mrs <- dplyr::arrange(merge(mrs, select(temp, SNP, labels), all.x=TRUE), index)
  mrl <- dplyr::arrange(merge(mrl, select(temp, SNP, labels), all.x=TRUE), index)
  
  mrres <- mr(temp, method_list=c("mr_egger_regression", "mr_ivw", "mr_weighted_median", "mr_weighted_mode"))
  
  p <- gridExtra::grid.arrange(
    mr_forest_plot(mrs)[[1]] +
      ggplot2::labs(
        title="a)",
        x=paste0("MR estimate (", outcome_units, " per ", exposure_units, ")")
      ),
    mr_scatter_plot(mrres, temp)[[1]] +
      ggplot2::labs(
        title="b)",
        x=paste0("SNP effect on ", exposure_name),
        y=paste0("SNP effect on ", outcome_name)
      ) +
      geom_label_repel(ggplot2::aes(label=labels), point.padding = unit(0.5, "lines")),
    mr_leaveoneout_plot(mrl)[[1]] +
      ggplot2::labs(
        title="c)",
        x=paste0("MR estimate (", outcome_units, " per ", exposure_units, ")"),
        y="Excluded variant"
      ),
    mr_funnel_plot(mrs)[[1]] +
      ggplot2::labs(title="d)") +
      ggplot2::theme(legend.position="none") +
      ggrepel::geom_label_repel(ggplot2::aes(label=labels), point.padding = unit(0.5, "lines")),
    ncol=2
  )
  
  return(p)
}


remove_outlier <- function(mr_tmp_df, rucker){
  #OUT.0: let's make vectors for the models:
  ivw_vect <- c("A", "B")
  egger_vect <- c("C", "D")
  #OUT.1: Let's check the model that fits the data better: IVW (A, B) or Egger (C, D).
  rucker_model <- rucker[[1]]$res
  #Let's format the data:
  df_4_radial <- RadialMR::format_radial(BXG = mr_tmp_df$beta.exposure, BYG = mr_tmp_df$beta.outcome,
                                         seBXG = mr_tmp_df$se.exposure, seBYG = mr_tmp_df$se.outcome,
                                         RSID = mr_tmp_df$SNP)
  if(rucker_model%in%ivw_vect){
    radial_output <- tryCatch(RadialMR::ivw_radial(r_input = df_4_radial, alpha = 0.05, weights = 3),  error = function(e){
      return(NA)
    })
  } else {
    radial_output <- tryCatch(RadialMR::egger_radial(r_input = df_4_radial, alpha = 0.05, weights = 3), error = function(e){
      return(NA)
    })
  }
  if(length(radial_output) == 1){
    return(NA)
  }
  #IF WE DON'T HAVE OUTLIERS RETURN ORIGINAL
  if(length(radial_output$outlier) ==  1){
    mr_wo_outliers <- mr_tmp_df
  } else {
    outliers <- radial_output$outliers$SNP
    mr_wo_outliers <- mr_tmp_df[which(!(mr_tmp_df$SNP%in%outliers)),]
  }
  return(mr_wo_outliers)
}


compute_strict_het <- function(dat_1_filt_1){
  
  dat_1_filt_1$mr <- dat_1_filt_1$beta.outcome/dat_1_filt_1$beta.exposure
  dat_1_filt_1$mr_se <- ((dat_1_filt_1$mr*((dat_1_filt_1$se.exposure/dat_1_filt_1$beta.exposure)^2+(dat_1_filt_1$se.outcome/dat_1_filt_1$beta.outcome)^2)^0.5)^2)^0.5
  het_Q_Isq <- meta::metagen(dat_1_filt_1$mr, dat_1_filt_1$mr_se)
  
  return(het_Q_Isq)
  
}

#####################################################
#STEP 1: let's obtain the data for each associations#
#####################################################

#Set pathway to find data:

path_2_input <- "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/"

setwd(path_2_input)

#Let's load the original protein data with all instruments:

prot_og_df <- fread("output/9_pqtl_results/2_mr/0_data_4_input/all_associations_per_protein_w_more_than_3_credible_sets.txt") #131
prot_og_df$protein <- as.character(unlist(sapply(prot_og_df$studyId, protein_parser)))
prot_og_df$chr_pos <- paste("chr", prot_og_df$chr_hg19, ":", prot_og_df$pos_hg19, sep = "") #in build 37 to match associations from original GWAS for each protein

pqtl_matched <- fread("output/9_pqtl_results/1_cis_pqtl_matches/cis_pqtl_matches.txt") #131

prot_list <- prot_og_df$prot_cred_set
prot_list <- prot_list[which(str_detect(prot_list, "UK"))]
prot_list <- as.character(unlist(sapply(prot_list, protein_parser)))
prot_list <- unique(prot_list) #20! Amazing
prot_list <- prot_list[order(prot_list)]

#Let's load the outcome too:

outcome <- fread("output/1_curated_gwas/t2d_curated.txt")
outcome$sample_size <- outcome$sample_cases+outcome$sample_controls

#Finally, let's load the proxies: 

proxies <- fread("output/9_pqtl_results/2_mr/0_data_4_input/haploreg_output.txt", fill = TRUE)

#We need to clean the chromosome cuz it might lead to issues...

proxies <- proxies[which(proxies$chr != ""),]
proxies <- proxies[which(proxies$chr != "Array"),]

chr_clean <- as.character(unlist(str_split(proxies$chr, "Array")))
chr_clean <- chr_clean[which(chr_clean != "")]
proxies$chr <- chr_clean

proxies$chr_pos_38 <- paste("chr", proxies$chr, ":", proxies$pos_hg38, sep = "")

length(unique(proxies$query_snp_rsid)) #128/141 - mostly rare variant end up being removed

##########################
#Let's loop over the data#
##########################

for(protein_ in prot_list){
  
  #STEP 0: print protein to check for dramas
  
  print(protein_)
  
  #STEP 1: load the data needed for that protein
  
  exposure_df <- load_prot(clumped_data_ = prot_og_df, protein_name = protein_)
  exposure_df$chr_19 <- as.character(unlist(sapply(exposure_df$chr_pos, chr_cleaner)))
  exposure_df$pos_19 <- as.character(unlist(sapply(exposure_df$chr_pos, pos_cleaner)))
  exposure_df$chr_pos_38 <- paste("chr", exposure_df$CHROM, ":", exposure_df$GENPOS, sep = "")
  
  #STEP 2: find the SNPs in the exposure data:
  
  #Find instruments for that protein
  
  prot_tmp <- prot_og_df[which(prot_og_df$protein == protein_),]
  exposure_match <- exposure_df[which(exposure_df$chr_pos%in%prot_tmp$chr_pos),]
  prot_tmp <- prot_tmp[which(prot_tmp$chr_pos%in%exposure_match$chr_pos),]
  prot_tmp <- prot_tmp[order(match(prot_tmp$chr_pos, exposure_match$chr_pos)),]
  
  print(length(which(prot_tmp$chr_pos == exposure_match$chr_pos)))
  exposure_match$variant <- prot_tmp$variant
  
  #Careful!! Let's remove the rare allele frequency variants:
  
  #exposure_match <- exposure_match[which(as.numeric(exposure_match$A1FREQ) > 0.01 & as.numeric(exposure_match$A1FREQ) < 0.99),]
  
  #############################
  #STEP 3: perform the match!!#
  #############################
  
  mr_df <- data_aligner(exposure_match, outcome)
  #mr_df <- mr_df[which(as.numeric(mr_df$eaf.outcome) > 0.01 & as.numeric(mr_df$eaf.outcome) < 0.99),]  #the outcomes have bigger sample size, we are gonna use them as filters for MAF!!
  
  ###################################################################
  #STEP 3.5: let's perform WALD or IVW when we have less than 2 SNPs#
  ###################################################################
  
  if(length(mr_df$SNP) < 3){
    
    mr_df <- mr_df[which(((mr_df$beta.exposure^2)/(mr_df$se.exposure^2)) > 10),] #3/3
    
    #Let's compute MR steiger:
    
    mr_df$r.exposure <- TwoSampleMR::get_r_from_pn(
      p = mr_df$pval.exposure,n = mr_df$samplesize.exposure)
    
    mr_df$r.outcome <- TwoSampleMR::get_r_from_pn(
      p = mr_df$pval.outcome,n = mr_df$samplesize.outcome)
    
    # mr_df$r.outcome <- TwoSampleMR::get_r_from_lor(
    #   lor = mr_df$beta.outcome,
    #   af = mr_df$eaf.outcome,
    #   ncase = mr_df$ncase.outcome,
    #   ncontrol = mr_df$ncontrol.outcome,
    #   prevalence = 0.0038 #whole index population unadjusted period prevalence in % from https://risteys.finngen.fi/endpoints/E4_PCOS
    # )
    
    mr_steiger_df <- TwoSampleMR::steiger_filtering(mr_df)
    mr_steiger_df <- mr_steiger_df[which(mr_steiger_df$steiger_dir == TRUE),] #3/3
    mr_steiger_df$mr_keep <- TRUE
    
    dir.create("output/9_pqtl_results/2_mr/1_mr_res/")
    dir.create("output/9_pqtl_results/2_mr/1_mr_res/t2d")
    dir.create(paste("output/9_pqtl_results/2_mr/1_mr_res/t2d/", protein_, sep = ""))
    
    #save data
    
    saveRDS(mr_steiger_df, paste("output/9_pqtl_results/2_mr/1_mr_res/t2d/", protein_ , "/instruments_after_steiger_df.RDS", sep = ""))
    
    #Compute results:
    
    mr_res <- TwoSampleMR::mr(mr_steiger_df)
    
    print(mr_res)
    
    #And save them:
    
    saveRDS(mr_res, paste("output/9_pqtl_results/2_mr/1_mr_res/t2d/", protein_ , "/mr_res_before_outlier_extraction.RDS", sep = ""))
    
    next()
    
  }
  
  ###########################################################
  #STEP 4: perform F-statistic and steiger filtering section#
  ###########################################################
  
  #Filter for F-statistic
  
  mr_df <- mr_df[which(((mr_df$beta.exposure^2)/(mr_df$se.exposure^2)) > 10),] #3/3

  #Let's compute MR steiger:
  
  mr_df$r.exposure <- TwoSampleMR::get_r_from_pn(
    p = mr_df$pval.exposure,n = mr_df$samplesize.exposure)
  
  #mr_df$r.outcome <- TwoSampleMR::get_r_from_pn(
  #  p = mr_df$pval.outcome,n = mr_df$samplesize.outcome)
  
   mr_df$r.outcome <- TwoSampleMR::get_r_from_lor(
     lor = mr_df$beta.outcome,
     af = mr_df$eaf.outcome,
     ncase = mr_df$ncase.outcome,
     ncontrol = mr_df$ncontrol.outcome,
     prevalence = as.numeric(17.30/100) #whole index population unadjusted period prevalence in % from https://risteys.finngen.fi/endpoints/E4_PCOS
  )
  
  mr_steiger_df <- TwoSampleMR::steiger_filtering(mr_df)
  mr_steiger_df <- mr_steiger_df[which(mr_steiger_df$steiger_dir == TRUE),] #3/3
  mr_steiger_df$mr_keep <- TRUE
  
  dir.create("output/9_pqtl_results/2_mr/1_mr_res/")
  dir.create("output/9_pqtl_results/2_mr/1_mr_res/t2d")
  dir.create(paste("output/9_pqtl_results/2_mr/1_mr_res/t2d/", protein_, sep = ""))
  
  #save data
  
  saveRDS(mr_steiger_df, paste("output/9_pqtl_results/2_mr/1_mr_res/t2d/", protein_ , "/instruments_after_steiger_df.RDS", sep = ""))
  
  ####################################
  #STEP 5: let's compute the analyses#
  ####################################
  
  mr_res <- TwoSampleMR::mr(mr_steiger_df)
  print(mr_res)

  #Let's sensitivity analyses:
  
  rucker <- TwoSampleMR::mr_rucker(mr_steiger_df)
  
  #Here there is heterogeneity and pleitropy
  
  strict_het <- compute_strict_het(mr_steiger_df) #they disagree, but still
  
  #That is quite good:
  
  mr_steiger_df$labels <- NA
  mr_steiger_df$units.exposure <- "SD"
  mr_steiger_df$units.outcome <- "SD"
  
  tiff(filename = paste("output/9_pqtl_results/2_mr/1_mr_res/t2d/", protein_ , "/plots_before_outlier_extraction.tiff", sep = ""), height=5000, width=5000, res=300)
  mr_plots(mr_steiger_df)
  dev.off()
  
  saveRDS(mr_res, paste("output/9_pqtl_results/2_mr/1_mr_res/t2d/", protein_ , "/mr_res_before_outlier_extraction.RDS", sep = ""))
  saveRDS(rucker, paste("output/9_pqtl_results/2_mr/1_mr_res/t2d/", protein_ , "/sensitivity_test_before_outlier_extraction.RDS", sep = ""))
  saveRDS(strict_het, paste("output/9_pqtl_results/2_mr/1_mr_res/t2d/", protein_ , "/strict_het_test_before_outlier_extraction.RDS", sep = ""))
  
  #############################################
  #Let's remove outliers and redo the analyses#
  #############################################
  
  mr_post <- remove_outlier(mr_steiger_df, rucker)
  
  if(length(mr_post) == 1){ #this only happens when it is NA
    
    next()
    
  }
  
  if(length(mr_post$SNP) < 3){ #run WALD or IVW only if less than 3 SNPs
    
    saveRDS(mr_post, paste("output/9_pqtl_results/2_mr/1_mr_res/t2d/", protein_ , "/instruments_after_outlier_removal_df.RDS", sep = ""))
            
    #Compute analyses
    
    mr_res <- TwoSampleMR::mr(mr_post)
    
    print(mr_res)
    
    #And save them:
    
    saveRDS(mr_res, paste("output/9_pqtl_results/2_mr/1_mr_res/t2d/", protein_ , "/mr_res_after_outlier_extraction.RDS", sep = ""))
    
    next()
    
  }
  
  saveRDS(mr_post, paste("output/9_pqtl_results/2_mr/1_mr_res/t2d/", protein_ , "/instruments_after_outlier_removal_df.RDS", sep = ""))
  
  #STEP 6: let's run the results for real now:
  
  mr_res <- TwoSampleMR::mr(mr_post)
  print(mr_res)
  
  rucker <- TwoSampleMR::mr_rucker(mr_post)
  
  strict_het <- compute_strict_het(mr_post) #they disagree, but still
  
  #That is quite good:
  
  mr_post$labels <- NA
  mr_post$units.exposure <- "SD"
  mr_post$units.outcome <- "SD"
  
  tiff(filename = paste("output/9_pqtl_results/2_mr/1_mr_res/t2d/", protein_ , "/plots_after_outlier_extraction.tiff", sep = ""), height=5000, width=5000, res=300)
  mr_plots(mr_post)
  dev.off()
  
  saveRDS(mr_res, paste("output/9_pqtl_results/2_mr/1_mr_res/t2d/", protein_ , "/mr_res_after_outlier_extraction.RDS", sep = ""))
  saveRDS(rucker, paste("output/9_pqtl_results/2_mr/1_mr_res/t2d/", protein_ , "/sensitivity_test_after_outlier_extraction.RDS", sep = ""))
  saveRDS(strict_het, paste("output/9_pqtl_results/2_mr/1_mr_res/t2d/", protein_ , "/strict_het_test_after_outlier_extraction.RDS", sep = ""))
  
}
