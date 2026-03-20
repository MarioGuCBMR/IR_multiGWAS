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

curated_2_munging_mv <- function(curated_df){
  
  #STEP 0: let's remove those that have not rsID, cuz LDSC gets confused:
  
  curated_df <- curated_df[which(str_detect(as.character(curated_df$variant), "rs") == TRUE)]
  
  #STEP 1: get hte MAF:
  
  curated_df$MAF <- curated_df$minimum_allele_frequency
  
  #STEP 2: get the right columns:
  
  curated_df_4_ldsc <- curated_df %>%
    select(variant, effect_allele, other_allele, beta, standard_error, p_value, MAF, sample_size)
  
  colnames(curated_df_4_ldsc) <- c("SNP", "A1", "A2", "BETA", "SE", "P", "MAF", "N")
  
  #STEP 3: Save the data:
  
  return(curated_df_4_ldsc)
  
}

curated_2_munging <- function(curated_df){
  
  #STEP 0: let's remove those that have not rsID, cuz LDSC gets confused:
  
  curated_df <- curated_df[which(str_detect(as.character(curated_df$variant), "rs") == TRUE)]
  
  #STEP 1: get hte MAF:
  
  curated_df$MAF <- ifelse(as.numeric(curated_df$effect_allele_frequency) > 0.50, 1-as.numeric(curated_df$effect_allele_frequency), as.numeric(curated_df$effect_allele_frequency))
  
  #STEP 2: get the right columns:
  
  curated_df_4_ldsc <- curated_df %>%
    select(variant, effect_allele, other_allele, beta, standard_error, p_value, MAF, sample_size)
  
  colnames(curated_df_4_ldsc) <- c("SNP", "A1", "A2", "BETA", "SE", "P", "MAF", "N")
  
  #STEP 3: Save the data:
  
  return(curated_df_4_ldsc)
  
}

curated_2_binary <- function(curated_df){
  
  #STEP 0: let's remove those that have not rsID, cuz LDSC gets confused:
  
  curated_df <- curated_df[which(str_detect(as.character(curated_df$variant), "rs") == TRUE)]
  
  #STEP 1: get hte MAF:
  
  curated_df$MAF <- ifelse(as.numeric(curated_df$effect_allele_frequency) > 0.50, 1-as.numeric(curated_df$effect_allele_frequency), as.numeric(curated_df$effect_allele_frequency))
  
  #STEP 2: get the effect sample size:
  
  curated_df$Neff<-4/((2*curated_df$MAF*(1-curated_df$MAF))*curated_df$standard_error^2)
  
  #STEP 2: get the right columns:
  
  curated_df_4_ldsc <- curated_df %>%
    select(variant, effect_allele, other_allele, beta, standard_error, p_value, MAF, Neff)
  
  colnames(curated_df_4_ldsc) <- c("SNP", "A1", "A2", "BETA", "SE", "P", "MAF", "N")
  
  #STEP 3: Save the data:
  
  return(curated_df_4_ldsc)
  
}

#####################################
#Loading data for dichotomous traits#
#####################################

#Running LDSC:

path_2_input <- "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2023"

setwd(path_2_input)

#############################################################################
#STEP 1: let's read the original curated data and prepare it for the munging#
#############################################################################

ir_fiadjbmi <- fread("output/3_mv_gwas/fiadjbmi_hdl_tg/fiadjbmi_hdl_tg_commmon_dwls_curated.txt")

fiadjbmi <- fread("output/1_curated_data/fiadjbmi_curated.txt")
hdl <- fread("output/1_curated_data/hdl_curated.txt")
tg <- fread("output/1_curated_data/tg_curated.txt")
tg_hdl_ratio <- fread("output/1_curated_data/tg_hdl_ratio_curated.txt")

#Anthropometric...

whr <- fread("output/1_curated_data/whr_curated.txt")
bmi <- fread("output/1_curated_data/bmi_curated.txt")
whradjbmi <- fread("output/1_curated_data/whradjbmi_curated.txt")
wc <- fread("output/1_curated_data/wc_curated.txt")
wcadjbmi <- fread("output/1_curated_data/wcadjbmi_giant_curated.txt")
hc <- fread("output/1_curated_data/hc_curated.txt")
hcadjbmi <- fread("output/1_curated_data/hcadjbmi_giant_curated.txt")

#Glycemic...
#fgadjbmi <- fread("../../Team projects/Hermina&Mario&MariaJose/output/1_curated_data/fgadjbmi_curated.txt")
#thgadjbmi <- fread("../../Team projects/Hermina&Mario&MariaJose/output/1_curated_data/thgadjbmi_curated.txt")
isiadjbmi <- fread("../../Team projects/Hermina&Mario&MariaJose/output/1_curated_data/isiadjbmi_curated.txt")
ifcadjbmi <- fread("../../Team projects/Hermina&Mario&MariaJose/output/1_curated_data/ifcadjbmi_curated.txt") #needed to recover variants! - we won't use it here

#Diseases...

t2d <- fread("output/1_curated_data/t2d_curated.txt")
chd <- fread("output/1_curated_data/chd_curated.txt")
nafld <- fread("output/1_curated_data/nafld_curated.txt")
ckd <- fread("output/1_curated_data/ckd_curated.txt")
pcos <- fread("output/1_curated_data/pcos_curated.txt")
hypertension <- fread("output/1_curated_data/hypertension_curated.txt")

##########################################################################
#For the munging I am going to need some help for ISIadjBMI and IFCadjBMI#
##########################################################################

#Let's get RSIDs for ISIadjBMI

fiadjbmi_match <- fiadjbmi[which(fiadjbmi$chr_pos%in%isiadjbmi$chr_pos),]
fiadjbmi_match <- fiadjbmi_match[which(duplicated(fiadjbmi_match$chr_pos) == FALSE),]

isiadjbmi_match <- isiadjbmi[which(isiadjbmi$chr_pos%in%fiadjbmi_match$chr_pos),]
isiadjbmi_match <- isiadjbmi_match[which(duplicated(isiadjbmi_match$chr_pos) == FALSE),]

isiadjbmi_match_ordered <- isiadjbmi_match[order(match(isiadjbmi_match$chr_pos, fiadjbmi_match$chr_pos)),]

length(which(isiadjbmi_match_ordered$chr_pos == fiadjbmi_match$chr_pos)) #perfect

isiadjbmi_match_ordered$variant <- fiadjbmi_match$variant

#Now the same for IFCadjBMI

fiadjbmi_match <- fiadjbmi[which(fiadjbmi$chr_pos%in%ifcadjbmi$chr_pos),]
fiadjbmi_match <- fiadjbmi_match[which(duplicated(fiadjbmi_match$chr_pos) == FALSE),]

ifcadjbmi_match <- ifcadjbmi[which(ifcadjbmi$chr_pos%in%fiadjbmi_match$chr_pos),]
ifcadjbmi_match <- ifcadjbmi_match[which(duplicated(ifcadjbmi_match$chr_pos) == FALSE),]

ifcadjbmi_match_ordered <- ifcadjbmi_match[order(match(ifcadjbmi_match$chr_pos, fiadjbmi_match$chr_pos)),]

length(which(ifcadjbmi_match_ordered$chr_pos == fiadjbmi_match$chr_pos)) #perfect

ifcadjbmi_match_ordered$variant <- fiadjbmi_match$variant

########################################
#One more push for the new tg_hdl ratio#
########################################

colnames(tg_hdl_ratio) <- c("chromosome",  "base_pair_location", "effect_allele", "other_allele", "beta", "standard_error", "effect_allele_frequency","neg_log_10_p_value", "id","variant",  "sample_size",  "p_value",  "chr_pos")

######################################
#Now let's get this ready for munging#
######################################

ir_fiadjbmi_4_munging <- curated_2_munging_mv(ir_fiadjbmi)

fiadjbmi_4_munging <- curated_2_munging(fiadjbmi)
hdl_4_munging <- curated_2_munging(hdl)
tg_4_munging <- curated_2_munging(tg)
tg_hdl_ratio_4_munging <- curated_2_munging(tg_hdl_ratio)

#Anthropometric traits:

whr_4_munging <- curated_2_munging(whr)
bmi_4_munging <- curated_2_munging(bmi)
whradjbmi_4_munging <- curated_2_munging(whradjbmi)
wc_4_munging <- curated_2_munging(wc)
wcadjbmi_4_munging <- curated_2_munging(wcadjbmi)
hc_4_munging <- curated_2_munging(hc)
hcadjbmi_4_munging <- curated_2_munging(hcadjbmi)

#Glycemic traits:

isiadjbmi_4_munging <- curated_2_munging(isiadjbmi_match_ordered)

#Disease traits:

t2d_4_munging <- curated_2_binary(t2d)
nafld_4_munging <- curated_2_binary(nafld)
chd_4_munging <- curated_2_binary(chd)
ckd_4_munging <- curated_2_binary(ckd)
pcos_4_munging <- curated_2_binary(pcos)
hypertension_4_munging <- curated_2_binary(hypertension)

dir.create("output/5_comparing_fiadjbmi_fi/3_gc/")
dir.create("output/5_comparing_fiadjbmi_fi/3_gc/1_data_4_munging/")

#Let's save the data:

fwrite(ir_fiadjbmi_4_munging, "output/5_comparing_fiadjbmi_fi/3_gc/1_data_4_munging/ir_fiadjbmi_4_munging.txt", sep = " ", col.names = TRUE, row.names = FALSE)

fwrite(fiadjbmi_4_munging, "output/5_comparing_fiadjbmi_fi/3_gc/1_data_4_munging/fiadjbmi_4_munging.txt", sep = " ", col.names = TRUE, row.names = FALSE)
fwrite(hdl_4_munging, "output/5_comparing_fiadjbmi_fi/3_gc/1_data_4_munging/hdl_4_munging.txt", sep = " ", col.names = TRUE, row.names = FALSE)
fwrite(tg_4_munging, "output/5_comparing_fiadjbmi_fi/3_gc/1_data_4_munging/tg_4_munging.txt", sep = " ", col.names = TRUE, row.names = FALSE)
fwrite(tg_hdl_ratio_4_munging, "output/5_comparing_fiadjbmi_fi/3_gc/1_data_4_munging/tg_hdl_ratio_4_munging.txt", sep = " ", col.names = TRUE, row.names = FALSE)

#Now anthropoemetric:

fwrite(bmi_4_munging, "output/5_comparing_fiadjbmi_fi/3_gc/1_data_4_munging/bmi_4_munging.txt", sep = " ", col.names = TRUE, row.names = FALSE)
fwrite(whr_4_munging, "output/5_comparing_fiadjbmi_fi/3_gc/1_data_4_munging/whr_4_munging.txt", sep = " ", col.names = TRUE, row.names = FALSE)
fwrite(whradjbmi_4_munging, "output/5_comparing_fiadjbmi_fi/3_gc/1_data_4_munging/whradjbmi_4_munging.txt", sep = " ", col.names = TRUE, row.names = FALSE)
fwrite(wc_4_munging, "output/5_comparing_fiadjbmi_fi/3_gc/1_data_4_munging/wc_4_munging.txt", sep = " ", col.names = TRUE, row.names = FALSE)
fwrite(wcadjbmi_4_munging, "output/5_comparing_fiadjbmi_fi/3_gc/1_data_4_munging/wcadjbmi_4_munging.txt", sep = " ", col.names = TRUE, row.names = FALSE)
fwrite(hc_4_munging, "output/5_comparing_fiadjbmi_fi/3_gc/1_data_4_munging/hc_4_munging.txt", sep = " ", col.names = TRUE, row.names = FALSE)
fwrite(hcadjbmi_4_munging, "output/5_comparing_fiadjbmi_fi/3_gc/1_data_4_munging/hcadjbmi_4_munging.txt", sep = " ", col.names = TRUE, row.names = FALSE)

#Now glycemic:

#fwrite(fgadjbmi_4_munging, "output/5_comparing_fiadjbmi_fi/3_gc/1_data_4_munging/fgadjbmi_4_munging.txt", sep = " ", col.names = TRUE, row.names = FALSE)
#fwrite(thgadjbmi_4_munging, "output/5_comparing_fiadjbmi_fi/3_gc/1_data_4_munging/thgadjbmi_4_munging.txt", sep = " ", col.names = TRUE, row.names = FALSE)
fwrite(isiadjbmi_4_munging, "output/5_comparing_fiadjbmi_fi/3_gc/1_data_4_munging/isiadjbmi_4_munging.txt", sep = " ", col.names = TRUE, row.names = FALSE)
#fwrite(ifcadjbmi_4_munging, "output/5_comparing_fiadjbmi_fi/3_gc/1_data_4_munging/ifcadjbmi_4_munging.txt", sep = " ", col.names = TRUE, row.names = FALSE)

#Now diseases:

fwrite(t2d_4_munging, "output/5_comparing_fiadjbmi_fi/3_gc/1_data_4_munging/t2d_4_munging.txt", sep = " ", col.names = TRUE, row.names = FALSE)
fwrite(nafld_4_munging, "output/5_comparing_fiadjbmi_fi/3_gc/1_data_4_munging/nafld_4_munging.txt", sep = " ", col.names = TRUE, row.names = FALSE)
fwrite(chd_4_munging, "output/5_comparing_fiadjbmi_fi/3_gc/1_data_4_munging/chd_4_munging.txt", sep = " ", col.names = TRUE, row.names = FALSE)
fwrite(ckd_4_munging, "output/5_comparing_fiadjbmi_fi/3_gc/1_data_4_munging/ckd_4_munging.txt", sep = " ", col.names = TRUE, row.names = FALSE)
fwrite(pcos_4_munging, "output/5_comparing_fiadjbmi_fi/3_gc/1_data_4_munging/pcos_4_munging.txt", sep = " ", col.names = TRUE, row.names = FALSE)
fwrite(hypertension_4_munging, "output/5_comparing_fiadjbmi_fi/3_gc/1_data_4_munging/hypertension_4_munging.txt", sep = " ", col.names = TRUE, row.names = FALSE)

###############################
#STEP 2: let's run the munging#
###############################

dir.create("output/5_comparing_fiadjbmi_fi/3_gc/2_munged_data")

setwd("output/5_comparing_fiadjbmi_fi/3_gc/2_munged_data") #this works

files<-c("../1_data_4_munging/ir_fiadjbmi_4_munging.txt",
         "../1_data_4_munging/fiadjbmi_4_munging.txt",
         "../1_data_4_munging/hdl_4_munging.txt",
         "../1_data_4_munging/tg_4_munging.txt",
         "../1_data_4_munging/tg_hdl_ratio_4_munging.txt",
         "../1_data_4_munging/whradjbmi_4_munging.txt",
         "../1_data_4_munging/whr_4_munging.txt",
         "../1_data_4_munging/bmi_4_munging.txt",
         "../1_data_4_munging/wc_4_munging.txt",
         "../1_data_4_munging/wcadjbmi_4_munging.txt",
         "../1_data_4_munging/hc_4_munging.txt",
         "../1_data_4_munging/hcadjbmi_4_munging.txt",
         "../1_data_4_munging/isiadjbmi_4_munging.txt",
         "../1_data_4_munging/t2d_4_munging.txt",
         "../1_data_4_munging/nafld_4_munging.txt",
         "../1_data_4_munging/chd_4_munging.txt",
         "../1_data_4_munging/ckd_4_munging.txt",
         "../1_data_4_munging/pcos_4_munging.txt",
         "../1_data_4_munging/hypertension_4_munging.txt")

#files<-c("../tg_hdl_ratio_4_munging.txt") #added here cuz I found out about this GWAS later...

#define the reference file being used to allign alleles across summary stats
#here we are using hapmap3

hm3<-"../../../../raw_data/ldsc-master/eur_w_ld_chr.tar/eur_w_ld_chr/w_hm3.snplist"

#name the traits 
trait.names<-c("ir_fiadjbmi", "fiadjbmi", "hdl", "tg", "tg_hdl_ratio", "whradjbmi", "whr", "bmi", "wc", "wcadjbmi", "hc", "hcadjbmi", "isiadjbmi", "t2d", "nafld", "chd", "ckd", "pcos", "hypertension")

#trait.names<-c("tg_hdl_ratio") #added here cuz I found out about this GWAS later...

#definte the imputation quality filter
info.filter=0

#define the MAF filter
maf.filter=0.01

#run munge

munge(files=files,hm3=hm3,trait.names=trait.names,info.filter=info.filter,maf.filter=maf.filter)

###########################################
#STEP 3: let's run the genetic correlation#
###########################################

input_vect <- c("ir_fiadjbmi.sumstats.gz", "fiadjbmi.sumstats.gz",  "hdl.sumstats.gz", "tg.sumstats.gz", "tg_hdl_ratio.sumstats.gz", "isiadjbmi.sumstats.gz",
                "bmi.sumstats.gz", "hc.sumstats.gz", "hcadjbmi.sumstats.gz", "wc.sumstats.gz", "wcadjbmi.sumstats.gz","whr.sumstats.gz", "whradjbmi.sumstats.gz",  
                "t2d.sumstats.gz", "nafld.sumstats.gz", "chd.sumstats.gz", "ckd.sumstats.gz", "pcos.sumstats.gz", "hypertension.sumstats.gz")


LDSCoutput_SUD <- ldsc(traits = input_vect, sample.prev =  c(NA, NA, NA,NA, NA, NA,
                                                             NA, NA, NA, NA, NA, NA, NA,
                                                             0.5, 0.5, 0.5, 0.5, 0.5, 0.5), population.prev = c(NA, NA, NA,NA, NA, NA,
                                                                                                                     NA, NA, NA, NA, NA, NA, NA,
                                                                                                                     0.1665, 0.0066, 0.1452, 0.0298, 0.1505, 0.3111), trait.names = c("ir_fiadjbmi",
                                                                                                                                                                                      "fiadjbmi",
                                                                                                                                                                                      "hdl",
                                                                                                                                                                                      "tg",
                                                                                                                                                                                      "tg_hdl_ratio",
                                                                                                                                                                                      "isiadjbmi",
                                                                                                                                                                                      "bmi",
                                                                                                                                                                                      "hc",
                                                                                                                                                                                      "hcadjbmi",
                                                                                                                                                                                      "wc",
                                                                                                                                                                                      "wcadjbmi",
                                                                                                                                                                                      "whr",
                                                                                                                                                                                      "whradjbmi",
                                                                                                                                                                                      "t2d",
                                                                                                                                                                                      "nafld",
                                                                                                                                                                                      "chd",
                                                                                                                                                                                      "ckd",
                                                                                                                                                                                      "pcos",
                                                                                                                                                                                      "hypertension"),
                       ld=("../../../../raw_data/ldsc-master/eur_w_ld_chr.tar/eur_w_ld_chr"), ## set the directory where the LD scores used the used in the analysis are located
                       wld=("../../../../raw_data/ldsc-master/eur_w_ld_chr.tar/eur_w_ld_chr"), stand = TRUE)
# Format decimals
correlationSUD=round(LDSCoutput_SUD$S_Stand, 3) 

#ir_fiadjbmi fiadjbmi    hdl     tg tg_hdl_ratio whradjbmi    whr    bmi     wc     hc fgadjbmi thgadjbmi isiadjbmi ifcadjbmi    t2d  nafld    chd    ckd   pcos hypertension
#[1,]       1.000    0.468 -0.847  0.844        0.916     0.471  0.610  0.396  0.544  0.334    0.159     0.224    -0.353     0.163  0.562  0.657  0.306  0.367  0.235        0.338

saveRDS(LDSCoutput_SUD, "../3_gc/genetic_correlation_res.RDS")

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

checkGC=subset(checkGC, intercept >=1 ) #some had overinflation!!

#######################
#Let's clean the data:#
#######################

if(is_empty(checkGC$pheno)){
  
  print("ALL GOOD")
  
} else {
  
  for(i in 1:length(checkGC$pheno)) {
    
    print(paste0("Read in data for ", checkGC$pheno[i]))
    phenoDirmunge=paste0(checkGC$pheno[i], ".sumstats.gz") 
    phenoDir=paste0("../1_data_4_munging/", checkGC$pheno[i], "_4_munging.txt") 
    data <- fread(phenoDir,header=T,data.table=F)
    
    N <- fread(phenoDirmunge,header=T,data.table=F)$N[1]
    print(paste0("Multiply SE by LDSC intercept of ", checkGC$intercept[i]))
    data$SE = data$SE * sqrt(checkGC$intercept[i])
    data$Zscore=as.numeric( data$BETA)/as.numeric( data$SE) # estimate adjusted p-values
    data$P=2*pnorm(-abs(data$Zscore))
    data$Zscore=NULL
    head(data)
    
    print("Save file on cluster")
    fwrite(data, 
           file=paste0("../1_data_4_munging/", checkGC$pheno[i], "_clean"), 
           sep="\t", 
           row.names = FALSE, 
           col.names = TRUE, 
           quote=F) 
    print("Remove non-GC controlled munged file and munge again")
    
    munge(files = paste0("../1_data_4_munging/", checkGC$pheno[i], "_clean"), 
          hm3 = paste0("../../../../raw_data/ldsc-master/eur_w_ld_chr.tar/eur_w_ld_chr/w_hm3.snplist"),
          trait.names=checkGC$pheno[i],
          N = N,
          info.filter = 0.9, 
          maf.filter = 0.01) 
  }
  
  
}

##############################################################
# ===== Re-run LDSC regression including GCed sumstats ===== #
##############################################################

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2023/output/5_comparing_fiadjbmi_fi/3_gc/2_munged_data") #this works

input_vect <- c("ir_fiadjbmi.sumstats.gz", "fiadjbmi.sumstats.gz",  "hdl.sumstats.gz", "tg.sumstats.gz", "tg_hdl_ratio.sumstats.gz", "isiadjbmi.sumstats.gz",
                "bmi.sumstats.gz", "hc.sumstats.gz", "hcadjbmi.sumstats.gz", "wc.sumstats.gz", "wcadjbmi.sumstats.gz","whr.sumstats.gz", "whradjbmi.sumstats.gz",  
                "t2d.sumstats.gz", "nafld.sumstats.gz", "chd.sumstats.gz", "ckd.sumstats.gz", "pcos.sumstats.gz", "hypertension.sumstats.gz")


LDSCoutput_SUD <- ldsc(traits = input_vect, sample.prev =  c(NA, NA, NA,NA, NA, NA,
                                                             NA, NA, NA, NA, NA, NA, NA,
                                                             0.5, 0.5, 0.5, 0.5, 0.5, 0.5), population.prev = c(NA, NA, NA,NA, NA, NA,
                                                                                                                NA, NA, NA, NA, NA, NA, NA,
                                                                                                                0.1665, 0.0066, 0.1452, 0.0298, 0.1505, 0.3111), trait.names = c("ir_fiadjbmi",
                                                                                                                                                                                 "fiadjbmi",
                                                                                                                                                                                 "hdl",
                                                                                                                                                                                 "tg",
                                                                                                                                                                                 "tg_hdl_ratio",
                                                                                                                                                                                 "isiadjbmi",
                                                                                                                                                                                 "bmi",
                                                                                                                                                                                 "hc",
                                                                                                                                                                                 "hcadjbmi",
                                                                                                                                                                                 "wc",
                                                                                                                                                                                 "wcadjbmi",
                                                                                                                                                                                 "whr",
                                                                                                                                                                                 "whradjbmi",
                                                                                                                                                                                 "t2d",
                                                                                                                                                                                 "nafld",
                                                                                                                                                                                 "chd",
                                                                                                                                                                                 "ckd",
                                                                                                                                                                                 "pcos",
                                                                                                                                                                                 "hypertension"),
                       ld=("../../../../raw_data/ldsc-master/eur_w_ld_chr.tar/eur_w_ld_chr"), ## set the directory where the LD scores used the used in the analysis are located
                       wld=("../../../../raw_data/ldsc-master/eur_w_ld_chr.tar/eur_w_ld_chr"), stand = TRUE)
# Format decimals
correlationSUD=round(LDSCoutput_SUD$S_Stand, 3) 

#Check that the results are sligthly different:

round(LDSCoutput_SUD$S_Stand, 3) #yes!! a lil bit.

saveRDS(LDSCoutput_SUD, "../3_gc/genetic_correlation_clean.RDS")

  

