##############
#INTRODUCTION#
##############

#This code performs genetic correlations

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)
library(GenomicSEM)

###################
#Loading functions#
###################

cleaning_munge <- function(munged_df, og_df, trait_name){
  
  #STEP 0: make dummy example to run the function:
  
  #munged_df <- whradjbmi_munged
  #og_df <- whradjbmi_og
  
  #STEP 1: let's match the data:
  
  og_match <- og_df[which(og_df$SNP%in%munged_df$SNP),]
  munged_match <- munged_df[which(munged_df$SNP%in%og_df$SNP),]
  
  #STEP 2: let's make sure that the data is OK:
  
  dupl <- og_match$SNP[which(duplicated(og_match$SNP) == TRUE)]
  
  if(length(dupl) != 0){
    
    og_match_non_dupl <- og_match[which(!(og_match$SNP%in%dupl)),]
    og_match_dupl <- og_match[which(og_match$SNP%in%dupl),]
    
    og_match_dupl$id_1 <- paste(og_match_dupl$SNP, og_match_dupl$A1, og_match_dupl$A2, sep= "_")
    og_match_dupl$id_2 <- paste(og_match_dupl$SNP, og_match_dupl$A2, og_match_dupl$A1, sep= "_")
    
    #Now let's make the IDs from the munged data:
    
    munged_match$id <- paste(munged_match$SNP, munged_match$A1, munged_match$A2, sep= "_")
    
    #Finally let's match this properly:
    
    og_match_dupl_solved <- og_match_dupl[which(og_match_dupl$id_1%in%munged_match$id | og_match_dupl$id_2%in%munged_match$id),]
    
    og_match_dupl_solved <- og_match_dupl_solved %>%
      select(-c("id_1", "id_2"))
    
    og_match <- rbind(og_match_non_dupl, og_match_dupl_solved)
    
  } 
  
  #Let's check this out:
  
  og_ordered <- og_match[order(match(og_match$SNP, munged_match$SNP)),]
  
  print(length(which(og_ordered$SNP == munged_match$SNP))) #all
  
  colnames(og_ordered) <- c("SNP",  "A1",   "A2",   "BETA", "SE",   "P",    "MAF",  "N")
  
  new_beta <- ifelse(og_ordered$A1 != munged_match$A1, as.numeric(og_ordered$BETA)*(-1), as.numeric(og_ordered$BETA))
  
  munged_match$BETA <- new_beta
  munged_match$SE <- og_ordered$SE
  munged_match$P <- og_ordered$P
  
  fwrite(munged_match, paste(trait_name, "_clean_sumstats.txt", sep=""), sep=" ", quote=FALSE, col.names = TRUE, row.names = FALSE)
  
}

###############
#Loading paths#
###############

path_2_input <- "/projects/kilpelainen-AUDIT/people/zlc436/IR_GSEM_2023/output/2_models/2_munged_data/fiadjbmi_hdl_tg/"

setwd(path_2_input)

########################################
#STEP 1: munge with the format required#
########################################

files<-c("../../1_data_4_munging/fiadjbmi_hdl_tg/fiadjbmi_4_munging.txt", "../../1_data_4_munging/fiadjbmi_hdl_tg/hdl_4_munging.txt", "../../1_data_4_munging/fiadjbmi_hdl_tg/tg_4_munging.txt")

#define the reference file being used to allign alleles across summary stats
#here we are using hapmap3

hm3<-"../../../../raw_data/ldsc-master/eur_w_ld_chr.tar/eur_w_ld_chr/w_hm3.snplist"

#name the traits 
trait.names<-c("fiadjbmi", "hdl", "tg")

#list the sample sizes. All but PTSD have SNP-specific sum of effective sample sizes so only its
#sample size is listed here
N=c(NA, NA, NA)

#definte the imputation quality filter
info.filter=0

#define the MAF filter
maf.filter=0.01

#run munge
munge(files=files,hm3=hm3,trait.names=trait.names,N=N,info.filter=info.filter,maf.filter=maf.filter)

####################
#Let's run the LDSC#
####################

input_vect <- c("fiadjbmi.sumstats.gz", "hdl.sumstats.gz", "tg.sumstats.gz")

LDSCoutput_SUD <- ldsc(traits = input_vect, sample.prev =  c(NA,NA,NA), population.prev = c(NA,NA,NA), trait.names = c("fiadjbmi","hdl","tg"),
                       ld=("../../../../raw_data/ldsc-master/eur_w_ld_chr.tar/eur_w_ld_chr/"), ## set the directory where the LD scores used the used in the analysis are located
                       wld=("../../../../raw_data/ldsc-master/eur_w_ld_chr.tar/eur_w_ld_chr/"), stand = TRUE)

saveRDS(LDSCoutput_SUD, paste0("../../3_gc/fiadjbmi_hdl_tg/", "fiadjbmi_hdl_tg.rds"))

#################################################################################
#Let's prepare the real data, cuz it needs betas and SE and we do not have them!#
#################################################################################

#First smoking initiation

fiadjbmi_munged <- fread("fiadjbmi.sumstats.gz")
fiadjbmi_og <- fread("../../1_data_4_munging/fiadjbmi_hdl_tg/fiadjbmi_4_munging.txt")

cleaning_munge(fiadjbmi_munged, fiadjbmi_og, "fiadjbmi")

#Now lifetime smoking

hdl_munged <- fread("hdl.sumstats.gz")
hdl_og <- fread("../../1_data_4_munging/fiadjbmi_hdl_tg/hdl_4_munging.txt")

cleaning_munge(hdl_munged, hdl_og, "hdl")

#Now university_degree

tg_munged <- fread("tg.sumstats.gz")
tg_og <- fread("../../1_data_4_munging/fiadjbmi_hdl_tg/tg_4_munging.txt")

cleaning_munge(tg_munged, tg_og, "tg")

#####################################################################
#First, let's define the data and the parameters that we want to use#
#####################################################################

files<-c("fiadjbmi_clean_sumstats.txt", "hdl_clean_sumstats.txt", "tg_clean_sumstats.txt")

ref= "../../../../raw_data/reference.1000G.maf.0.005.txt.gz"

trait.names<-c("fiadjbmi","hdl","tg")

se.logit=c(F,F,F)
linprob=c(F,F,F)
info.filter=0
maf.filter=0.01
OLS=c(T,T,T)
N=NULL
betas=NULL

all_sumstats <-sumstats(files=files,ref=ref,trait.names=trait.names,se.logit=se.logit,OLS=OLS,linprob=linprob,N=N,betas=NULL,info.filter=info.filter,maf.filter=maf.filter,keep.indel=FALSE,parallel=FALSE,cores=NULL)

fwrite(all_sumstats, "../../../3_mv_gwas/fiadjbmi_hdl_tg/fiadjbmi_hdl_tg_4_mvgwas.txt")