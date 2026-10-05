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

path_2_input <- "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/review/review_output/2_mv_gwas/"

setwd(path_2_input)

#################################################################################
#Let's prepare the real data, cuz it needs betas and SE and we do not have them!#
#################################################################################

#First smoking initiation

fiadjbmi_munged <- fread("munged_data/fiadjbmi.sumstats.gz")
fiadjbmi_og <- fread("../1_curated_data/fiadjbmi_4_munging.txt")

cleaning_munge(fiadjbmi_munged, fiadjbmi_og, "fiadjbmi")

#Now lifetime smoking

hdl_munged <- fread("munged_data/hdl.sumstats.gz")
hdl_og <- fread("../1_curated_data//hdl_4_munging.txt")

cleaning_munge(hdl_munged, hdl_og, "hdl")

#Now university_degree

tg_munged <- fread("munged_data/tg.sumstats.gz")
tg_og <- fread("../1_curated_data/tg_4_munging.txt")

cleaning_munge(tg_munged, tg_og, "tg")

#####################################################################
#First, let's define the data and the parameters that we want to use#
#####################################################################

files<-c("fiadjbmi_clean_sumstats.txt", "hdl_clean_sumstats.txt", "tg_clean_sumstats.txt")

ref= "../../../raw_data/reference.1000G.maf.0.005.txt.gz"

trait.names<-c("fiadjbmi","hdl","tg")

se.logit=c(F,F,F)
linprob=c(F,F,F)
info.filter=0
maf.filter=0.01
OLS=c(T,T,T)
N=NULL
betas=NULL

all_sumstats <-sumstats(files=files,ref=ref,trait.names=trait.names,se.logit=se.logit,OLS=OLS,linprob=linprob,N=N,betas=NULL,info.filter=info.filter,maf.filter=maf.filter,keep.indel=FALSE,parallel=FALSE,cores=NULL)

fwrite(all_sumstats, "fiadjbmi_hdl_tg_4_mvgwas.txt")
