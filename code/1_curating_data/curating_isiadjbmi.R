##############
#INTRODUCTION#
##############

#This is a code to curate isiadjbmi data from Che et al 2021 (European data)
#Data was downloaded 03/08/2022

#Let's go:

###########
#libraries#
###########

library(data.table)
library(tidyverse)

#######################
#Loading isiadjbmi data#
#######################

#We are gonna load the isiadjbmi from 2021. 

project_path <- "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina&Mario&MariaJose//" #change it with your own path.

setwd(project_path)

isiadjbmi <- fread("raw_data/Glycemic_traits_raw_data/MAGIC_postchallengeIR_ISI_adjBMI_EUR.tsv.gz")

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

isiadjbmi$chr_pos <- paste("chr", isiadjbmi$chromosome, ":", isiadjbmi$base_pair_location, sep = "")

###############################
#STEP 2: CHECK FOR WEIRD STUFF#
###############################

#1. Always order the dataframe with the variant ID. You might find interesting things.

isiadjbmi <- isiadjbmi[order(isiadjbmi$variant),]

head(isiadjbmi) #some of them do not have rsID! Importantly, insertions and deletions are transformed to DI!
tail(isiadjbmi) #we most likely have sexual chromosome data.

#Let's check it out:

summary(as.numeric(isiadjbmi$chromosome)) #we don't have NA, but we have chromosome 23! Let's remove it

#Let's see if it is only X and Y:

chromosomes <- isiadjbmi$chromosome[which(duplicated(isiadjbmi$chromosome) == FALSE)] #we have "X" and "XY"!!

isiadjbmi_autosomal <- isiadjbmi[which(!(isiadjbmi$chromosome == 23)),]

#Did this solve the issue?

summary(as.numeric(isiadjbmi_autosomal$chromosome)) #Yes!! No more 23s

#Let's convert the column into a numeric by the end, it is annoying if they are characters.

#####################################
#STEP 3: REMOVE LOW SAMPLE SIZE SNPs#
#####################################

#We are going to remove those that are < 10.000

isiadjbmi_corrected <- isiadjbmi_autosomal[which(isiadjbmi_autosomal$n> 10000),]

###########################
#STEP 4: REMOVE MAF < 0.01#
###########################

summary(as.numeric(isiadjbmi_corrected$effect_allele_frequency)) #We need to remove the NAs, and those with EAF > 0.99 or EAF < 0.01

#We have 4128819 that are NAs.

#################################################################
#We are gonna separate the dataframes to: with NA and without NA#
#################################################################

isiadjbmi_eaf_NA <- isiadjbmi_corrected[which(is.na(isiadjbmi_corrected$effect_allele_frequency) == TRUE),]
isiadjbmi_eaf_OK <- isiadjbmi_corrected[which(is.na(isiadjbmi_corrected$effect_allele_frequency) == FALSE),]

summary(isiadjbmi_eaf_OK$effect_allele_frequency) #No NAs, but we have to correct.
 
#isiadjbmi_eaf_OK <- isiadjbmi_eaf_OK[which(isiadjbmi_eaf_OK$effect_allele_frequency > 0.01),]
#isiadjbmi_eaf_OK <- isiadjbmi_eaf_OK[which(isiadjbmi_eaf_OK$effect_allele_frequency < 0.99),]

summary(isiadjbmi_eaf_OK$effect_allele_frequency) #worked like a charm.

isiadjbmi_corrected_eaf <- rbind(isiadjbmi_eaf_NA, isiadjbmi_eaf_OK) #13,526,969

##############################
#CURATION OF THE DATA IS DONE#
##############################

#Let's check now if there are any in the MHC region.

isiadjbmi_mhc <- isiadjbmi_corrected_eaf[which(as.numeric(isiadjbmi_corrected_eaf$chromosome) == 6 & as.numeric(isiadjbmi_corrected_eaf$base_pair_location) >= 26000000 & as.numeric(isiadjbmi_corrected_eaf$base_pair_location) <= 34000000),]

summary(as.numeric(isiadjbmi_mhc$chromosome)) #perfect!!
summary(as.numeric(isiadjbmi_mhc$base_pair_location)) #perfect!!

isiadjbmi_end <- isiadjbmi_corrected_eaf[which(!(isiadjbmi_corrected_eaf$variant%in%isiadjbmi_mhc$variant)),]

colnames(isiadjbmi_end) <- c("chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", "standard_error", "effect_allele_frequency", "p_value", "sample_size", "het_p_value", "variant", "chr_pos")

isiadjbmi_end <- isiadjbmi_end %>%
  select(variant, chromosome, base_pair_location, effect_allele, other_allele, beta, standard_error, effect_allele_frequency, p_value, sample_size, chr_pos)

#########################
#We can save this data!!#
#########################

fwrite(isiadjbmi_end, "output/1_curated_data/isiadjbmi_curated.txt")
