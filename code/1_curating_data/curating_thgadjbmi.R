##############
#INTRODUCTION#
##############

#This is a code to curate thgadjbmi data from Che et al 2021 (European data)
#Data was downloaded 03/08/2022

#Let's go:

###########
#libraries#
###########

library(data.table)
library(tidyverse)

#######################
#Loading thgadjbmi data#
#######################

#We are gonna load the thgadjbmi from 2021. 

project_path <- "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina&Mario&MariaJose//" #change it with your own path.

setwd(project_path)

thgadjbmi <- fread("raw_data/Glycemic_traits_raw_data/MAGIC1000G_2hGlu_EUR.tsv.gz")

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

thgadjbmi$chr_pos <- paste("chr", thgadjbmi$chromosome, ":", thgadjbmi$base_pair_location, sep = "")

###############################
#STEP 2: CHECK FOR WEIRD STUFF#
###############################

#1. Always order the dataframe with the variant ID. You might find interesting things.

thgadjbmi <- thgadjbmi[order(thgadjbmi$variant),]

head(thgadjbmi) #some of them do not have rsID! Importantly, insertions and deletions are transformed to DI!
tail(thgadjbmi) #we most likely have sexual chromosome data.

#Let's check it out:

summary(as.numeric(thgadjbmi$chromosome)) #exactly: we have NAs.

#Let's see if it is only X and Y:

chromosomes <- thgadjbmi$chromosome[which(duplicated(thgadjbmi$chromosome) == FALSE)] #we have "X" and "XY"!!

thgadjbmi_autosomal <- thgadjbmi[which(!(thgadjbmi$chromosome == "XY" | thgadjbmi$chromosome == "X")),]

#Did this solve the issue?

summary(as.numeric(thgadjbmi_autosomal$chromosome)) #Yes!! No more NAs

#Let's convert the column into a numeric by the end, it is annoying if they are characters.

#####################################
#STEP 3: REMOVE LOW SAMPLE SIZE SNPs#
#####################################

#We are going to remove those that are < 10.000

thgadjbmi_corrected <- thgadjbmi_autosomal[which(thgadjbmi_autosomal$sample_size > 10000),]

###########################
#STEP 4: REMOVE MAF < 0.01#
###########################

summary(as.numeric(thgadjbmi_corrected$effect_allele_frequency)) #We need to remove the NAs, and those with EAF > 0.99 or EAF < 0.01

#We have 4128819 that are NAs.

#################################################################
#We are gonna separate the dataframes to: with NA and without NA#
#################################################################

thgadjbmi_eaf_NA <- thgadjbmi_corrected[which(is.na(thgadjbmi_corrected$effect_allele_frequency) == TRUE),]
thgadjbmi_eaf_OK <- thgadjbmi_corrected[which(is.na(thgadjbmi_corrected$effect_allele_frequency) == FALSE),]

summary(thgadjbmi_eaf_OK$effect_allele_frequency) #No NAs, but we have to correct.
 
#thgadjbmi_eaf_OK <- thgadjbmi_eaf_OK[which(thgadjbmi_eaf_OK$effect_allele_frequency > 0.01),]
#thgadjbmi_eaf_OK <- thgadjbmi_eaf_OK[which(thgadjbmi_eaf_OK$effect_allele_frequency < 0.99),]

summary(thgadjbmi_eaf_OK$effect_allele_frequency) #worked like a charm.

thgadjbmi_corrected_eaf <- rbind(thgadjbmi_eaf_NA, thgadjbmi_eaf_OK) #13,526,969

##############################
#CURATION OF THE DATA IS DONE#
##############################

#Let's check now if there are any in the MHC region.

thgadjbmi_mhc <- thgadjbmi_corrected_eaf[which(as.numeric(thgadjbmi_corrected_eaf$chromosome) == 6 & as.numeric(thgadjbmi_corrected_eaf$base_pair_location) >= 26000000 & as.numeric(thgadjbmi_corrected_eaf$base_pair_location) <= 34000000),]

summary(as.numeric(thgadjbmi_mhc$chromosome)) #perfect!!
summary(as.numeric(thgadjbmi_mhc$base_pair_location)) #perfect!!

thgadjbmi_end <- thgadjbmi_corrected_eaf[which(!(thgadjbmi_corrected_eaf$variant%in%thgadjbmi_mhc$variant)),]

#########################
#We can save this data!!#
#########################

fwrite(thgadjbmi_end, "output/1_curated_data/thgadjbmi_curated.txt")
