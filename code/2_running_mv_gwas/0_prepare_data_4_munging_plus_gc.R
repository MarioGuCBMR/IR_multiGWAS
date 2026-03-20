##############
#INTRODUCTION#
##############

#This code performs genetic correlations and prepares the data for running GSEM

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)
library(GenomicSEM)

###################
#Loading functions#
###################

curated_2_munging <- function(curated_df){

  #STEP 0: let's remove those that have not rsID, cuz LDSC gets confused:

  curated_df <- curated_df[which(str_detect(as.character(curated_df$variant), "rs") == TRUE)]

  curated_df <- curated_df[which(is.na(curated_df$effect_allele_frequency) == FALSE),] #if not, the munging gets very confused!!!

  #STEP 1: get hte MAF:

  curated_df$MAF <- ifelse(as.numeric(curated_df$effect_allele_frequency) > 0.50, 1-as.numeric(curated_df$effect_allele_frequency), as.numeric(curated_df$effect_allele_frequency))

  #STEP 2: get the right columns:

  curated_df_4_ldsc <- curated_df %>%
    select(variant, effect_allele, other_allele, beta, standard_error, p_value, MAF, sample_size)

  colnames(curated_df_4_ldsc) <- c("SNP", "A1", "A2", "BETA", "SE", "P", "MAF", "N")

  #STEP 3: Save the data:

  return(curated_df_4_ldsc)

}

##############
#Loading data#
##############

#Running LDSC:

path_2_input <- "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025"

setwd(path_2_input)

#############################################################################
#STEP 1: let's read the original curated data and prepare it for the munging#
#############################################################################

fiadjbmi <- fread("output/1_curated_gwas/fiadjbmi_curated.txt")
hdl <- fread("output/1_curated_gwas/hdl_curated.txt")
tg <- fread("output/1_curated_gwas/tg_curated.txt")

fiadjbmi_4_munging <- curated_2_munging(fiadjbmi)
hdl_4_munging <- curated_2_munging(hdl)
tg_4_munging <- curated_2_munging(tg)

fwrite(fiadjbmi_4_munging, "output/1_curated_gwas/fiadjbmi_4_munging.txt", col.names = TRUE, row.names = FALSE, sep = " ")
fwrite(hdl_4_munging, "output/1_curated_gwas/hdl_4_munging.txt", col.names = TRUE, row.names = FALSE, sep = " ")
fwrite(tg_4_munging, "output/1_curated_gwas/tg_4_munging.txt", col.names = TRUE, row.names = FALSE, sep = " ")

###############################
#STEP 2: let's run the munging#
###############################

dir.create("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/output/2_mv_gwas/munged_data")

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/output/2_mv_gwas/munged_data") #this works

files<-c("../../1_curated_gwas/fiadjbmi_4_munging.txt", "../../1_curated_gwas/hdl_4_munging.txt", "../../1_curated_gwas/tg_4_munging.txt")

#define the reference file being used to allign alleles across summary stats
#here we are using hapmap3

hm3<-"../../../raw_data/ldsc-master/eur_w_ld_chr.tar/eur_w_ld_chr/w_hm3.snplist"

#name the traits
trait.names<-c("fiadjbmi", "hdl", "tg")

#definte the imputation quality filter
info.filter=0

#define the MAF filter
maf.filter=0.01

#run munge

munge(files=files,hm3=hm3,trait.names=trait.names,info.filter=info.filter,maf.filter=maf.filter)

###########################################
#STEP 3: let's run the genetic correlation#
###########################################

input_vect <- c("fiadjbmi.sumstats.gz", "hdl.sumstats.gz", "tg.sumstats.gz")


LDSCoutput_SUD <- ldsc(traits = input_vect, sample.prev =  c(NA,NA,NA), population.prev = c(NA,NA,NA), trait.names = c("fiadjbmi","hdl","tg"),
                       ld=("../../../raw_data/ldsc-master/eur_w_ld_chr.tar/eur_w_ld_chr"), ## set the directory where the LD scores used the used in the analysis are located
                       wld=("../../../raw_data/ldsc-master/eur_w_ld_chr.tar/eur_w_ld_chr"), stand = TRUE)
# Format decimals
correlationSUD=round(LDSCoutput_SUD$S_Stand, 3)

#fiadjbmi    hdl     tg
#[1,]    1.000 -0.336  0.358
#[2,]   -0.336  1.000 -0.649
#[3,]    0.358 -0.649  1.000

saveRDS(LDSCoutput_SUD, "../fiadjbmi_hdl_tg.rds")

###########################################################################
#Let's study the heritability estimates to ensure we don't have inflation!#
###########################################################################

#First let's make a dummy dataframe!
h2out=data.frame(trait=colnames(correlationSUD),
                 est= rep(NA, length(colnames(correlationSUD))))

#And now fill it up

for (i in 1:length(colnames(correlationSUD))) {
  h2out$est[i]=round(LDSCoutput_SUD$S[i,i],3)
}

h2out

#trait   est
#1 fiadjbmi 0.097
#2      hdl 0.061
#3       tg 0.062

########################################
# Check which GWAS had overinflation!! #
########################################

rownames(correlationSUD)=colnames(correlationSUD)

str(LDSCoutput_SUD)
intercept=LDSCoutput_SUD$I
colnames(intercept)=colnames(LDSCoutput_SUD$S)
rownames(intercept)=colnames(LDSCoutput_SUD$S)
checkGC=data.frame(pheno = colnames(intercept),
                   intercept = diag(intercept)) # get the LDSC intercept per sumstats

checkGC=subset(checkGC, intercept >=1 ) #none had overinflation!! Thus no need to do anything else. I am going to comment the section below.

#######################
#Let's clean the data:#
#######################

# if(is_empty(checkGC$pheno)){
#
#   print("ALL GOOD")
#
# } else {
#
#   for(i in 1:length(checkGC$pheno)) {
#
#     print(paste0("Read in data for ", checkGC$pheno[i]))
#     phenoDirmunge=paste0(checkGC$pheno[i], ".sumstats.gz")
#     phenoDir=paste0("../../1_data_4_munging/fiadjbmi_hdl_tg/", checkGC$pheno[i], "_4_munging.txt")
#     data <- fread(phenoDir,header=T,data.table=F)
#
#     N <- fread(phenoDirmunge,header=T,data.table=F)$N[1]
#     print(paste0("Multiply SE by LDSC intercept of ", checkGC$intercept[i]))
#     data$SE = data$SE * sqrt(checkGC$intercept[i])
#     data$Zscore=as.numeric( data$Beta)/as.numeric( data$SE) # estimate adjusted p-values
#     data$P=2*pnorm(-abs(data$Zscore))
#     data$Zscore=NULL
#     head(data)
#
#     print("Save file on cluster")
#     fwrite(data,
#            file=paste0("../../1_data_4_munging/fiadjbmi_hdl_tg/", checkGC$pheno[i], "_clean"),
#            sep="\t",
#            row.names = FALSE,
#            col.names = TRUE,
#            quote=F)
#     print("Remove non-GC controlled munged file and munge again")
#
#     munge(files = paste0("../../1_data_4_munging/fiadjbmi_hdl_tg/", checkGC$pheno[i]),
#           hm3 = paste0("../../../../raw_data/ldsc-master/eur_w_ld_chr.tar/eur_w_ld_chr/w_hm3.snplist"),
#           trait.names=checkGC$pheno[i],
#           N = N,
#           info.filter = 0.9,
#           maf.filter = 0.01)
#   }
#
#
# }
#
# ##############################################################
# # ===== Re-run LDSC regression including GCed sumstats ===== #
# ##############################################################
#
# if(is_empty(checkGC$pheno)){
#
#   print("SKIPPING REDOING GC")
#
# } else {
#
#   setwd("J:/CBMR/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2023/output/2_models/2_munged_data/fiadjbmi_hdl_tg") #this works
#
#   input_vect <- c("fiadjbmi.sumstats.gz", "hdl.sumstats.gz", "tg.sumstats.gz")
#
#   LDSCoutput_SUD <- ldsc(traits = input_vect, sample.prev =  c(NA,NA,NA), population.prev = c(NA,NA,NA), trait.names = c("fiadjbmi","hdl","tg"),
#                          ld=("../../../../raw_data/ldsc-master/eur_w_ld_chr.tar/eur_w_ld_chr"), ## set the directory where the LD scores used the used in the analysis are located
#                          wld=("../../../../raw_data/ldsc-master/eur_w_ld_chr.tar/eur_w_ld_chr"), stand = TRUE)
#
# }
#
# #Let's save the data now that has run!!
#
# saveRDS(LDSCoutput_SUD, "../../3_gc/fiadjbmi_hdl_tg/fiadjbmi_hdl_tg.rds")



