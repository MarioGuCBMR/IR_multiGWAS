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

get_lower_tri_index <- function(i, j) {
  if (i < j) {
    tmp <- i
    i <- j
    j <- tmp
  }
  return(i * (i - 1) / 2 + j)
}

#####################################
#Loading data for dichotomous traits#
#####################################

#Running LDSC:

path_2_input <- "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025"

setwd(path_2_input)

#############################################################################
#STEP 1: let's read the original curated data and prepare it for the munging#
#############################################################################

ir_fiadjbmi <- fread("output/2_mv_gwas/fiadjbmi_hdl_tg_common_dwls_curated.txt")

fiadjbmi <- fread("output/1_curated_gwas/fiadjbmi_curated.txt")
hdl <- fread("output/1_curated_gwas/hdl_curated.txt")
tg <- fread("output/1_curated_gwas/tg_curated.txt")
tg_hdl_ratio <- fread("output/1_curated_gwas/tg_hdl_ratio_curated.txt")
isiadjbmi <- fread("output/1_curated_gwas/isiadjbmi_curated.txt")
ifcadjbmi <- fread("output/1_curated_gwas/ifcadjbmi_curated.txt") 
fgadjbmi <- fread("output/1_curated_gwas/fgadjbmi_curated.txt")
thgadjbmi <- fread("output/1_curated_gwas/thgadjbmi_curated.txt") 

#Anthropometric...

whr <- fread("output/1_curated_gwas/whr_curated.txt")
bmi <- fread("output/1_curated_gwas/bmi_curated.txt")
whradjbmi <- fread("output/1_curated_gwas/whradjbmi_curated.txt")
wc <- fread("output/1_curated_gwas/wc_curated.txt")
wcadjbmi <- fread("output/1_curated_gwas/wcadjbmi_giant_curated.txt")
hc <- fread("output/1_curated_gwas/hc_curated.txt")
hcadjbmi <- fread("output/1_curated_gwas/hcadjbmi_giant_curated.txt")

#Diseases...

t2d <- fread("output/1_curated_gwas/t2d_curated.txt")
chd <- fread("output/1_curated_gwas/chd_curated.txt")
nafld <- fread("output/1_curated_gwas/nafld_curated.txt")
ckd <- fread("output/1_curated_gwas/ckd_curated.txt")
pcos <- fread("output/1_curated_gwas/pcos_curated.txt")
hypertension <- fread("output/1_curated_gwas/hypertension_curated.txt")

##########################################################################
#For the munging I am going to need some help for ISIadjBMI and IFCadjBMI#
##########################################################################

#Let's get RSIDs for ISIadjBMI

yes_vect <- c("A", "G", "T", "C")
fiadjbmi <- fiadjbmi[which(fiadjbmi$effect_allele%in%yes_vect & fiadjbmi$other_allele%in%yes_vect),]

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

#########################################################################
#For the munging I am going to need some help for FGadjBMI and ISIadjBMI#
#########################################################################

#Let's get RSIDs for FIadjBMI

yes_vect <- c("A", "G", "T", "C")
fiadjbmi <- fiadjbmi[which(fiadjbmi$effect_allele%in%yes_vect & fiadjbmi$other_allele%in%yes_vect),]

fiadjbmi_match <- fiadjbmi[which(fiadjbmi$chr_pos%in%fgadjbmi$chr_pos),]
fiadjbmi_match <- fiadjbmi_match[which(duplicated(fiadjbmi_match$chr_pos) == FALSE),]

fgadjbmi_match <- fgadjbmi[which(fgadjbmi$chr_pos%in%fiadjbmi_match$chr_pos),]
fgadjbmi_match <- fgadjbmi_match[which(duplicated(fgadjbmi_match$chr_pos) == FALSE),]

fgadjbmi_match_ordered <- fgadjbmi_match[order(match(fgadjbmi_match$chr_pos, fiadjbmi_match$chr_pos)),]

length(which(fgadjbmi_match_ordered$chr_pos == fiadjbmi_match$chr_pos)) #perfect

fgadjbmi_match_ordered$variant <- fiadjbmi_match$variant

#Now the same for thgadjbmi

fiadjbmi_match <- fiadjbmi[which(fiadjbmi$chr_pos%in%thgadjbmi$chr_pos),]
fiadjbmi_match <- fiadjbmi_match[which(duplicated(fiadjbmi_match$chr_pos) == FALSE),]

thgadjbmi_match <- thgadjbmi[which(thgadjbmi$chr_pos%in%fiadjbmi_match$chr_pos),]
thgadjbmi_match <- thgadjbmi_match[which(duplicated(thgadjbmi_match$chr_pos) == FALSE),]

thgadjbmi_match_ordered <- thgadjbmi_match[order(match(thgadjbmi_match$chr_pos, fiadjbmi_match$chr_pos)),]

length(which(thgadjbmi_match_ordered$chr_pos == fiadjbmi_match$chr_pos)) #perfect

thgadjbmi_match_ordered$variant <- fiadjbmi_match$variant

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
ifcadjbmi_4_munging <- curated_2_munging(ifcadjbmi_match_ordered)
fgadjbmi_4_munging <- curated_2_munging(fgadjbmi_match_ordered)
thgadjbmi_4_munging <- curated_2_munging(thgadjbmi_match_ordered)

#Disease traits:

t2d_4_munging <- curated_2_binary(t2d)
nafld_4_munging <- curated_2_binary(nafld)
chd_4_munging <- curated_2_binary(chd)
ckd_4_munging <- curated_2_binary(ckd)
pcos_4_munging <- curated_2_binary(pcos)
hypertension_4_munging <- curated_2_binary(hypertension)

dir.create("output/4_ir_loci_discovery/5_gc/")
dir.create("output/4_ir_loci_discovery/5_gc/1_data_4_munging/")

#Let's save the data:

fwrite(ir_fiadjbmi_4_munging, "output/4_ir_loci_discovery/5_gc/1_data_4_munging/ir_fiadjbmi_4_munging.txt", sep = " ", col.names = TRUE, row.names = FALSE)

fwrite(fiadjbmi_4_munging, "output/4_ir_loci_discovery/5_gc/1_data_4_munging/fiadjbmi_4_munging.txt", sep = " ", col.names = TRUE, row.names = FALSE)
fwrite(hdl_4_munging, "output/4_ir_loci_discovery/5_gc/1_data_4_munging/hdl_4_munging.txt", sep = " ", col.names = TRUE, row.names = FALSE)
fwrite(tg_4_munging, "output/4_ir_loci_discovery/5_gc/1_data_4_munging/tg_4_munging.txt", sep = " ", col.names = TRUE, row.names = FALSE)
fwrite(tg_hdl_ratio_4_munging, "output/4_ir_loci_discovery/5_gc/1_data_4_munging/tg_hdl_ratio_4_munging.txt", sep = " ", col.names = TRUE, row.names = FALSE)

#Now anthropoemetric:

fwrite(bmi_4_munging, "output/4_ir_loci_discovery/5_gc/1_data_4_munging/bmi_4_munging.txt", sep = " ", col.names = TRUE, row.names = FALSE)
fwrite(whr_4_munging, "output/4_ir_loci_discovery/5_gc/1_data_4_munging/whr_4_munging.txt", sep = " ", col.names = TRUE, row.names = FALSE)
fwrite(whradjbmi_4_munging, "output/4_ir_loci_discovery/5_gc/1_data_4_munging/whradjbmi_4_munging.txt", sep = " ", col.names = TRUE, row.names = FALSE)
fwrite(wc_4_munging, "output/4_ir_loci_discovery/5_gc/1_data_4_munging/wc_4_munging.txt", sep = " ", col.names = TRUE, row.names = FALSE)
fwrite(wcadjbmi_4_munging, "output/4_ir_loci_discovery/5_gc/1_data_4_munging/wcadjbmi_4_munging.txt", sep = " ", col.names = TRUE, row.names = FALSE)
fwrite(hc_4_munging, "output/4_ir_loci_discovery/5_gc/1_data_4_munging/hc_4_munging.txt", sep = " ", col.names = TRUE, row.names = FALSE)
fwrite(hcadjbmi_4_munging, "output/4_ir_loci_discovery/5_gc/1_data_4_munging/hcadjbmi_4_munging.txt", sep = " ", col.names = TRUE, row.names = FALSE)

#Now glycemic:

fwrite(fgadjbmi_4_munging, "output/4_ir_loci_discovery/5_gc/1_data_4_munging/fgadjbmi_4_munging.txt", sep = " ", col.names = TRUE, row.names = FALSE)
fwrite(thgadjbmi_4_munging, "output/4_ir_loci_discovery/5_gc/1_data_4_munging/thgadjbmi_4_munging.txt", sep = " ", col.names = TRUE, row.names = FALSE)
fwrite(isiadjbmi_4_munging, "output/4_ir_loci_discovery/5_gc/1_data_4_munging/isiadjbmi_4_munging.txt", sep = " ", col.names = TRUE, row.names = FALSE)
fwrite(ifcadjbmi_4_munging, "output/4_ir_loci_discovery/5_gc/1_data_4_munging/ifcadjbmi_4_munging.txt", sep = " ", col.names = TRUE, row.names = FALSE)

#Now diseases:

fwrite(t2d_4_munging, "output/4_ir_loci_discovery/5_gc/1_data_4_munging/t2d_4_munging.txt", sep = " ", col.names = TRUE, row.names = FALSE)
fwrite(nafld_4_munging, "output/4_ir_loci_discovery/5_gc/1_data_4_munging/nafld_4_munging.txt", sep = " ", col.names = TRUE, row.names = FALSE)
fwrite(chd_4_munging, "output/4_ir_loci_discovery/5_gc/1_data_4_munging/chd_4_munging.txt", sep = " ", col.names = TRUE, row.names = FALSE)
fwrite(ckd_4_munging, "output/4_ir_loci_discovery/5_gc/1_data_4_munging/ckd_4_munging.txt", sep = " ", col.names = TRUE, row.names = FALSE)
fwrite(pcos_4_munging, "output/4_ir_loci_discovery/5_gc/1_data_4_munging/pcos_4_munging.txt", sep = " ", col.names = TRUE, row.names = FALSE)
fwrite(hypertension_4_munging, "output/4_ir_loci_discovery/5_gc/1_data_4_munging/hypertension_4_munging.txt", sep = " ", col.names = TRUE, row.names = FALSE)

###############################
#STEP 2: let's run the munging#
###############################

dir.create("output/4_ir_loci_discovery/5_gc/2_munged_data")

setwd("output/4_ir_loci_discovery/5_gc/2_munged_data") #this works

#NOTE: GSEM sometimes dies when the munging is done for a lot of traits. Dunno why, but redoing the analyses helps, hence why I have commented some traits.

files<-c(#"../1_data_4_munging/ir_fiadjbmi_4_munging.txt",
         #"../1_data_4_munging/fiadjbmi_4_munging.txt",
         #"../1_data_4_munging/hdl_4_munging.txt",
         #"../1_data_4_munging/tg_4_munging.txt"
         #,"../1_data_4_munging/tg_hdl_ratio_4_munging.txt",
         #"../1_data_4_munging/whradjbmi_4_munging.txt",
         #"../1_data_4_munging/whr_4_munging.txt",
         #"../1_data_4_munging/bmi_4_munging.txt",
         "../1_data_4_munging/wc_4_munging.txt"
         #,
         #"../1_data_4_munging/wcadjbmi_4_munging.txt",
         #"../1_data_4_munging/hc_4_munging.txt",
         #"../1_data_4_munging/hcadjbmi_4_munging.txt",
         #"../1_data_4_munging/isiadjbmi_4_munging.txt",
         #"../1_data_4_munging/ifcadjbmi_4_munging.txt",
         #"../1_data_4_munging/fgadjbmi_4_munging.txt",
         #"../1_data_4_munging/thgadjbmi_4_munging.txt"
         #,
         #"../1_data_4_munging/t2d_4_munging.txt",
         #"../1_data_4_munging/nafld_4_munging.txt",
         #"../1_data_4_munging/chd_4_munging.txt",
         #"../1_data_4_munging/ckd_4_munging.txt",
         #"../1_data_4_munging/pcos_4_munging.txt",
         #"../1_data_4_munging/hypertension_4_munging.txt"
         )

#define the reference file being used to allign alleles across summary stats
#here we are using hapmap3

hm3<-"../../../../raw_data/ldsc-master/eur_w_ld_chr.tar/eur_w_ld_chr/w_hm3.snplist"

#name the traits 
trait.names<-c(#"ir_fiadjbmi", "fiadjbmi", "hdl", 
               #"tg"
               #, "tg_hdl_ratio", "whradjbmi", "whr", "bmi", 
               "wc"
               #, "wcadjbmi", "hc", "hcadjbmi", "isiadjbmi",   "ifcadjbmi", 
               #"fgadjbmi", "thgadjbmi"
               #, 
               #"t2d", "nafld", "chd", "ckd", "pcos", "hypertension"
               )

#trait.names<-c("tg_hdl_ratio") #added here cuz I found out about this GWAS later...

#definte the imputation quality filter
info.filter=0

#define the MAF filter
maf.filter=0.01

#run munge

#munge(files=files,hm3=hm3,trait.names=trait.names,info.filter=info.filter,maf.filter=maf.filter,log.name = "IR_munging") #run again since it died in ISIadjBMI the first time...
munge(files=files,hm3=hm3,trait.names=trait.names,info.filter=info.filter,maf.filter=maf.filter,log.name = "IR_munging_9") #run again since it died in ISIadjBMI the first time...

###########################################
#STEP 3: let's run the genetic correlation#
###########################################

input_vect <- c("ir_fiadjbmi.sumstats.gz", "fiadjbmi.sumstats.gz",  "hdl.sumstats.gz", "tg.sumstats.gz", "tg_hdl_ratio.sumstats.gz", "isiadjbmi.sumstats.gz", "ifcadjbmi.sumstats.gz", "fgadjbmi.sumstats.gz", "thgadjbmi.sumstats.gz",
                "bmi.sumstats.gz", "hc.sumstats.gz", "hcadjbmi.sumstats.gz", "wc.sumstats.gz", "wcadjbmi.sumstats.gz","whr.sumstats.gz", "whradjbmi.sumstats.gz",  
                "t2d.sumstats.gz", "nafld.sumstats.gz", "chd.sumstats.gz", "ckd.sumstats.gz", "pcos.sumstats.gz", "hypertension.sumstats.gz")


LDSCoutput_SUD <- ldsc(traits = input_vect, sample.prev =  c(NA, NA, NA,NA, NA, NA, NA, NA, NA,
                                                             NA, NA, NA, NA, NA, NA, NA,
                                                             0.5, 0.5, 0.5, 0.5, 0.5, 0.5), population.prev = c(NA, NA, NA,NA, NA, NA, NA, NA, NA,
                                                                                                                     NA, NA, NA, NA, NA, NA, NA,
                                                                                                                     0.1665, 0.0066, 0.1452, 0.0298, 0.1505, 0.3111), trait.names = c("ir_fiadjbmi",
                                                                                                                                                                                      "fiadjbmi",
                                                                                                                                                                                      "hdl",
                                                                                                                                                                                      "tg",
                                                                                                                                                                                      "tg_hdl_ratio",
                                                                                                                                                                                      "isiadjbmi",
                                                                                                                                                                                      "ifcadjbmi",
                                                                                                                                                                                      "fgadjbmi",
                                                                                                                                                                                      "thgadjbmi",
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
                       wld=("../../../../raw_data/ldsc-master/eur_w_ld_chr.tar/eur_w_ld_chr"), stand = TRUE, ldsc.log = "IR_GC")
# Format decimals
correlationSUD=round(LDSCoutput_SUD$S_Stand, 3) 

# ir_fiadjbmi fiadjbmi    hdl     tg tg_hdl_ratio isiadjbmi ifcadjbmi fgadjbmi thgadjbmi    bmi     hc hcadjbmi     wc wcadjbmi    whr whradjbmi    t2d  nafld    chd    ckd   pcos hypertension
# [1,]       1.000    0.468 -0.847  0.844        0.932    -0.353     0.163    0.159     0.224  0.467  0.334   -0.128  0.544    0.238  0.610     0.396  0.562  0.657  0.306  0.367  0.235        0.338

dir.create("../3_gc")

saveRDS(LDSCoutput_SUD, "../3_gc/ldsc_output.RDS")
saveRDS(correlationSUD, "../3_gc/genetic_correlation_res.RDS")

LDSCoutput_SUD = readRDS("../3_gc/ldsc_output.RDS")
correlationSUD = readRDS("../3_gc/genetic_correlation_res.RDS")

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

checkGC=subset(checkGC, intercept >=1 ) #some had overinflation, but very little. From previous analyses we know that it does not affect. Also, all are from publically available data. Hence all is good

#########################################################
#Finally, let's make the table and the plot for the GC!!#
#########################################################

#First, we are going to make a table, but to do that, we need to compute to extract the GC, the standard errors and the p-values

gc_output <- readRDS("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/output/4_ir_loci_discovery/5_gc/3_gc/ldsc_output.RDS")

#Let's get matrices of LDSC and p-values:
cor_mat <- gc_output$S_Stand  # 22 x 22 correlation matrix
V_stand <- gc_output$V_Stand    # 253 x 253 covariance matrix for free parameters

traits <- colnames(cor_mat)
k <- nrow(cor_mat)

# Get indices of free parameters (lower triangle + diagonal)
lt_idx <- lower.tri(cor_mat, diag = TRUE)  # Logical index

#Get the correlations and error for those indeces:

free_params <- cor_mat[lt_idx]          # length 253
SE_free_params <- sqrt(diag(V_stand))   # length 253

# Calculate Z and p-values
z_scores <- free_params / SE_free_params
p_values <- 2 * pnorm(-abs(z_scores))

# Initialize empty p-value matrix
p_mat <- matrix(NA, nrow = k, ncol = k)
rownames(p_mat) <- traits
colnames(p_mat) <- traits

# Fill lower triangle + diagonal with p-values
p_mat[lt_idx] <- p_values

# Fill upper triangle by symmetry
p_mat[upper.tri(p_mat)] <- t(p_mat)[upper.tri(p_mat)]

rownames(cor_mat) <- colnames(cor_mat)

#Let's get only the IR traits

cor_mat_4_plot <- cor_mat[c(1:9, 17),]
cor_mat_4_plot <- as.data.frame(cor_mat_4_plot) %>%
  dplyr::select("ir_fiadjbmi", "fiadjbmi", "hdl", "tg", "tg_hdl_ratio", "isiadjbmi", "ifcadjbmi", "fgadjbmi", "thgadjbmi", "t2d")

#Same for p-values:

p_mat_4_plot <- p_mat[c(1:9, 17),]
p_mat_4_plot <- as.data.frame(p_mat_4_plot) %>%
  dplyr::select("ir_fiadjbmi", "fiadjbmi", "hdl", "tg", "tg_hdl_ratio", "isiadjbmi", "ifcadjbmi", "fgadjbmi", "thgadjbmi", "t2d")

####################################
#Let's make a beautiful figure here#
####################################

rownames(cor_mat_4_plot) <- colnames(cor_mat_4_plot) <- c("Common factor GWAS", "FIadjBMI", "HDL", "TG", "TG/HDL", "ISIadjBMI", "IFCadjBMI", "FGadjBMI", "2hGadjBMI", "T2D")
rownames(p_mat_4_plot) <- colnames(p_mat_4_plot) <- c("Common factor GWAS", "FIadjBMI", "HDL", "TG", "TG/HDL", "ISIadjBMI", "IFCadjBMI", "FGadjBMI", "2hGadjBMI", "T2D")

cor_df_melted <- reshape2::melt(as.matrix(cor_mat_4_plot))
colnames(cor_df_melted) <- c("Trait1", "Trait2", "Correlation")

p_df_melted <- reshape2::melt(as.matrix(p_mat_4_plot))
colnames(p_df_melted) <- c("Trait1", "Trait2", "Pvalue")

# Merge
plot_df <- merge(cor_df_melted, p_df_melted, by = c("Trait1", "Trait2"))

# Size scale: invert p-value for size (lower p = bigger square)
library(ggplot2)
library(dplyr)
library(tidyr)
library(forcats)

# Define trait order
trait_order <- c("Common factor GWAS", "FIadjBMI", "HDL", "TG", "TG/HDL", "ISIadjBMI", "IFCadjBMI", "FGadjBMI", "2hGadjBMI", "T2D")

# Add diagonal entries
traits <- unique(c(plot_df$Trait1, plot_df$Trait2))
diag_df <- data.frame(Trait1 = traits, Trait2 = traits, Correlation = 1, Pvalue = 0)
plot_df <- bind_rows(plot_df, diag_df)

# Create annotation labels
plot_df <- plot_df %>%
  mutate(
    stars = case_when(
      Pvalue < 5e-08 ~ "**",
      Pvalue < 0.05 ~ "*",
      TRUE ~ ""
    ),
    label = paste0(sprintf("%.2f", Correlation), stars)
  )

# Ensure trait order
plot_df$Trait1 <- factor(plot_df$Trait1, levels = trait_order)
plot_df$Trait2 <- factor(plot_df$Trait2, levels = trait_order)

# Plot
plotio <- ggplot(plot_df, aes(x = Trait1, y = fct_rev(Trait2), fill = Correlation)) +
  geom_tile(color = "grey90", linewidth = 0.3) +
  geom_text(aes(label = label), size = 5, family = "Times", color = "black", fontface = "bold") +
  scale_fill_gradient2(
    low = "#4A90E2",    # softer blue
    mid = "white",
    high = "#D64541",   # muted red
    midpoint = 0,
    limits = c(-1, 1),
    name = "Genetic Correlation",
  ) +
  theme_minimal(base_size = 14, base_family = "Times") +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, size = 8, color = "black", face = "bold"),
    axis.text.y = element_text(size = 10, color = "black", face = "bold"),
    panel.grid = element_blank(),
    #panel.border = element_rect(color = "grey85", fill = NA, linewidth = 0.5),
    axis.ticks = element_line(color = "grey85", linewidth = 0.3),
    axis.title = element_blank(),
    legend.title = element_text(face = "bold", size = 12),
    legend.text = element_text(size = 11)
  )

ggsave(
  filename = "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/manuscript/figures/gc_glycemic_traits.svg",
  plot = plotio,
  width = 12,
  height = 10,
  units = "in",  # or "in" depending on your preference, but cm is common
  device = "svg"
)

#####################################
#Let's use this data to make a table#
#####################################

#Let's get the traits ready for this section:

#Clean names

traits <-  c("Common factor GWAS", "FIadjBMI", "HDL", "TG", "TG/HDL", "ISIadjBMI", "IFCadjBMI", "FGadjBMI", "2hGadjBMI", "BMI", "HC", "HCadjBMI", "WC", "WCadjBMI", "WHR", "WHRadjBMI", "T2D", "NAFLD",  "CHD", "CKD", "PCOS", "Hypertension")

#Use all DF because we want all correlations

rownames(cor_mat) <- colnames(cor_mat) <- traits
rownames(p_mat) <- colnames(p_mat) <- traits

#Set up variables for the loop...

k=length(traits)
first_trait <- traits[1]

# including an empty dataframe to fill...
results <- data.frame(
  Trait1 = character(),
  Trait2 = character(),
  GC = numeric(),
  SE = numeric(),
  CI_lower = numeric(),
  CI_upper = numeric(),
  P = numeric(),
  stringsAsFactors = FALSE
)

# Fill results for pairs involving the first trait
for (j in 1:k) {
  trait2 <- traits[j]
  if (trait2 == first_trait) next  # skip self-comparison
  
  # Get matrix index
  i <- which(traits == first_trait)
  idx <- get_lower_tri_index(i, j)
  
  
  # Extract values
  gc <- cor_mat[i, j]
  se <- SE_free_params[idx]
  z  <- gc / se
  p  <- 2 * pnorm(-abs(z))
  ci_lower <- gc - 1.96 * se
  ci_upper <- gc + 1.96 * se
  
  results <- rbind(results, data.frame(
    Trait1 = first_trait,
    Trait2 = trait2,
    GC = gc,
    SE = se,
    CI_lower = ci_lower,
    CI_upper = ci_upper,
    P = p
  ))
}

# Optional: format p-values and CIs for nicer display
results$CI <- sprintf("[%.2f, %.2f]", results$CI_lower, results$CI_upper)
results$P_formatted <- signif(results$P, 2)

# Final table
final_table <- results[, c("Trait1", "Trait2", "GC", "SE", "CI", "P_formatted")]

fwrite(final_table, "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/output/4_ir_loci_discovery/5_gc/3_gc/genetic_correlations_res_clean.csv")
