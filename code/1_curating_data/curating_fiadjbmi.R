##############
#INTRODUCTION#
##############

#This is a code to curate FIadjBMI data from Che et al 2021 (European data)
#Data was downloaded 03/08/2022

#Let's go:

###########
#libraries#
###########

library(data.table)
library(tidyverse)

#######################
#Loading FIadjBMI data#
#######################

#We are gonna load the FIAdjBMI from 2021. 

project_path <- "N:/SUN-CBMR-Kilpelainen-Group/insulin_resistance_variants_common_info/" #change it with your own path.

setwd(project_path)

fiadjbmi <- fread("raw_data/MAGIC1000G_FI_EUR.tsv.gz")

###################
#Curating the data#
###################

#We are gonna follow the following steps:

#1. Create and ID for the variants joining chromosome and position.

#2. Check the dataframe to see anything suspicious.

#3. Check the sample sizes and remove variants with low sample size.

#4. Check the allele frequencies and remove those with MAF < 0.01.

######################
#STEP 1: GENERATE IDs#
######################

#We are gonna set the chr_pos:

fiadjbmi$chr_pos <- paste("chr", fiadjbmi$chromosome, ":", fiadjbmi$base_pair_location, sep = "")

###############################
#STEP 2: CHECK FOR WEIRD STUFF#
###############################

#1. Always order the dataframe with the variant ID. You might find interesting things.

fiadjbmi <- fiadjbmi[order(fiadjbmi$variant),]

head(fiadjbmi) #some of them do not have rsID! Importantly, insertions and deletions are transformed to DI!
tail(fiadjbmi) #we most likely have sexual chromosome data.

#Let's check it out:

summary(as.numeric(fiadjbmi$chromosome)) #exactly: we have NAs.

#Let's see if it is only X and Y:

chromosomes <- fiadjbmi$chromosome[which(duplicated(fiadjbmi$chromosome) == FALSE)] #we have "X" and "XY"!!

fiadjbmi_autosomal <- fiadjbmi[which(!(fiadjbmi$chromosome == "XY" | fiadjbmi$chromosome == "X")),]

#Did this solve the issue?

summary(as.numeric(fiadjbmi_autosomal$chromosome)) #Yes!! No more NAs

#Let's convert the column into a numeric by the end, it is annoying if they are characters.

#####################################
#STEP 3: REMOVE LOW SAMPLE SIZE SNPs#
#####################################

#We are going to remove those that are < 10.000

fiadjbmi_corrected <- fiadjbmi_autosomal[which(fiadjbmi_autosomal$sample_size > 10000),]

###########################
#STEP 4: REMOVE MAF < 0.01#
###########################

summary(as.numeric(fiadjbmi_corrected$effect_allele_frequency)) #We need to remove the NAs, and those with EAF > 0.99 or EAF < 0.01

#We have 4128819 that are NAs.

#################################################################
#We are gonna separate the dataframes to: with NA and without NA#
#################################################################

fiadjbmi_eaf_NA <- fiadjbmi_corrected[which(is.na(fiadjbmi_corrected$effect_allele_frequency) == TRUE),]
fiadjbmi_eaf_OK <- fiadjbmi_corrected[which(is.na(fiadjbmi_corrected$effect_allele_frequency) == FALSE),]

summary(fiadjbmi_eaf_OK$effect_allele_frequency) #No NAs, but we have to correct.
 
fiadjbmi_eaf_OK <- fiadjbmi_eaf_OK[which(fiadjbmi_eaf_OK$effect_allele_frequency > 0.01),]
fiadjbmi_eaf_OK <- fiadjbmi_eaf_OK[which(fiadjbmi_eaf_OK$effect_allele_frequency < 0.99),]

summary(fiadjbmi_eaf_OK$effect_allele_frequency) #worked like a charm.

fiadjbmi_corrected_eaf <- rbind(fiadjbmi_eaf_NA, fiadjbmi_eaf_OK) #13,526,969

##############################
#CURATION OF THE DATA IS DONE#
##############################

#Let's check now if there are any in the MHC region.

fiadjbmi_mhc <- fiadjbmi_corrected_eaf[which(as.numeric(fiadjbmi_corrected_eaf$chromosome) == 6 & as.numeric(fiadjbmi_corrected_eaf$base_pair_location) >= 26000000 & as.numeric(fiadjbmi_corrected_eaf$base_pair_location) <= 34000000),]

summary(as.numeric(fiadjbmi_mhc$chromosome)) #perfect!!
summary(as.numeric(fiadjbmi_mhc$base_pair_location)) #perfect!!

fiadjbmi_end <- fiadjbmi_corrected_eaf[which(!(fiadjbmi_corrected_eaf$variant%in%fiadjbmi_mhc$variant)),]

#########################
#We can save this data!!#
#########################

fwrite(fiadjbmi_end, "output/1_curated_data/fiadjbmi_curated.txt")
