##############
#INTRODUCTION#
##############

#This is a code to curate fgadjbmi data from Che et al 2021 (European data)
#Data was downloaded 03/08/2022

#Let's go:

###########
#libraries#
###########

library(data.table)
library(tidyverse)

#######################
#Loading fgadjbmi data#
#######################

#We are gonna load the fgadjbmi from 2021. 

project_path <- "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina&Mario&MariaJose//" #change it with your own path.

setwd(project_path)

fgadjbmi <- fread("raw_data/Glycemic_traits_raw_data/MAGIC1000G_FG_EUR.tsv.gz")

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

fgadjbmi$chr_pos <- paste("chr", fgadjbmi$chromosome, ":", fgadjbmi$base_pair_location, sep = "")

###############################
#STEP 2: CHECK FOR WEIRD STUFF#
###############################

#1. Always order the dataframe with the variant ID. You might find interesting things.

fgadjbmi <- fgadjbmi[order(fgadjbmi$variant),]

head(fgadjbmi) #some of them do not have rsID! Importantly, insertions and deletions are transformed to DI!
tail(fgadjbmi) #we most likely have sexual chromosome data.

#Let's check it out:

summary(as.numeric(fgadjbmi$chromosome)) #exactly: we have NAs.

#Let's see if it is only X and Y:

chromosomes <- fgadjbmi$chromosome[which(duplicated(fgadjbmi$chromosome) == FALSE)] #we have "X" and "XY"!!

fgadjbmi_autosomal <- fgadjbmi[which(!(fgadjbmi$chromosome == "XY" | fgadjbmi$chromosome == "X")),]

#Did this solve the issue?

summary(as.numeric(fgadjbmi_autosomal$chromosome)) #Yes!! No more NAs

#Let's convert the column into a numeric by the end, it is annoying if they are characters.

#####################################
#STEP 3: REMOVE LOW SAMPLE SIZE SNPs#
#####################################

#We are going to remove those that are < 10.000

fgadjbmi_corrected <- fgadjbmi_autosomal[which(fgadjbmi_autosomal$sample_size > 10000),]

###########################
#STEP 4: REMOVE MAF < 0.01#
###########################

summary(as.numeric(fgadjbmi_corrected$effect_allele_frequency)) #We need to remove the NAs, and those with EAF > 0.99 or EAF < 0.01

#We have 4128819 that are NAs.

#################################################################
#We are gonna separate the dataframes to: with NA and without NA#
#################################################################

fgadjbmi_eaf_NA <- fgadjbmi_corrected[which(is.na(fgadjbmi_corrected$effect_allele_frequency) == TRUE),]
fgadjbmi_eaf_OK <- fgadjbmi_corrected[which(is.na(fgadjbmi_corrected$effect_allele_frequency) == FALSE),]

#summary(fgadjbmi_eaf_OK$effect_allele_frequency) #No NAs, but we have to correct.
 
#fgadjbmi_eaf_OK <- fgadjbmi_eaf_OK[which(fgadjbmi_eaf_OK$effect_allele_frequency > 0.01),]
#fgadjbmi_eaf_OK <- fgadjbmi_eaf_OK[which(fgadjbmi_eaf_OK$effect_allele_frequency < 0.99),]

#summary(fgadjbmi_eaf_OK$effect_allele_frequency) #worked like a charm.

fgadjbmi_corrected_eaf <- rbind(fgadjbmi_eaf_NA, fgadjbmi_eaf_OK) #13,526,969

##############################
#CURATION OF THE DATA IS DONE#
##############################

#Let's check now if there are any in the MHC region.

fgadjbmi_mhc <- fgadjbmi_corrected_eaf[which(as.numeric(fgadjbmi_corrected_eaf$chromosome) == 6 & as.numeric(fgadjbmi_corrected_eaf$base_pair_location) >= 26000000 & as.numeric(fgadjbmi_corrected_eaf$base_pair_location) <= 34000000),]

summary(as.numeric(fgadjbmi_mhc$chromosome)) #perfect!!
summary(as.numeric(fgadjbmi_mhc$base_pair_location)) #perfect!!

fgadjbmi_end <- fgadjbmi_corrected_eaf[which(!(fgadjbmi_corrected_eaf$variant%in%fgadjbmi_mhc$variant)),]

#########################
#We can save this data!!#
#########################

fwrite(fgadjbmi_end, "output/1_curated_data/fgadjbmi_curated.txt")
