##############
#INTRODUCTION#
##############

#This is a code to curate ifcadjbmi data from Che et al 2021 (European data)
#Data was downloaded 03/08/2022

#Let's go:

###########
#libraries#
###########

library(data.table)
library(tidyverse)

#######################
#Loading ifcadjbmi data#
#######################

#We are gonna load the ifcadjbmi from 2021. 

project_path <- "J:/CBMR/SUN-CBMR-Kilpelainen-Group/Team projects/Hermina&Mario&MariaJose//" #change it with your own path.

setwd(project_path)

ifcadjbmi <- fread("raw_data/Glycemic_traits_raw_data/MAGIC_postchallengeIR_IFC_adjBMI_EUR.tsv.gz")

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

ifcadjbmi$chr_pos <- paste("chr", ifcadjbmi$chromosome, ":", ifcadjbmi$base_pair_location, sep = "")

###############################
#STEP 2: CHECK FOR WEIRD STUFF#
###############################

#1. Always order the dataframe with the variant ID. You might find interesting things.

ifcadjbmi <- ifcadjbmi[order(ifcadjbmi$variant),]

head(ifcadjbmi) #some of them do not have rsID! Importantly, insertions and deletions are transformed to DI!
tail(ifcadjbmi) #we most likely have sexual chromosome data.

#Let's check it out:

summary(as.numeric(ifcadjbmi$chromosome)) #we don't have NA, but we have chromosome 23! Let's remove it

#Let's see if it is only X and Y:

chromosomes <- ifcadjbmi$chromosome[which(duplicated(ifcadjbmi$chromosome) == FALSE)] #we have "X" and "XY"!!

ifcadjbmi_autosomal <- ifcadjbmi[which(!(ifcadjbmi$chromosome == 23)),]

#Did this solve the issue?

summary(as.numeric(ifcadjbmi_autosomal$chromosome)) #Yes!! No more 23s

#Let's convert the column into a numeric by the end, it is annoying if they are characters.

#####################################
#STEP 3: REMOVE LOW SAMPLE SIZE SNPs#
#####################################

#We are going to remove those that are < 10.000

ifcadjbmi_corrected <- ifcadjbmi_autosomal[which(ifcadjbmi_autosomal$n> 10000),]

###########################
#STEP 4: REMOVE MAF < 0.01#
###########################

summary(as.numeric(ifcadjbmi_corrected$effect_allele_frequency)) #We need to remove the NAs, and those with EAF > 0.99 or EAF < 0.01

#We have 4128819 that are NAs.

#################################################################
#We are gonna separate the dataframes to: with NA and without NA#
#################################################################

ifcadjbmi_eaf_NA <- ifcadjbmi_corrected[which(is.na(ifcadjbmi_corrected$effect_allele_frequency) == TRUE),]
ifcadjbmi_eaf_OK <- ifcadjbmi_corrected[which(is.na(ifcadjbmi_corrected$effect_allele_frequency) == FALSE),]

summary(ifcadjbmi_eaf_OK$effect_allele_frequency) #No NAs, but we have to correct.
 
#ifcadjbmi_eaf_OK <- ifcadjbmi_eaf_OK[which(ifcadjbmi_eaf_OK$effect_allele_frequency > 0.01),]
#ifcadjbmi_eaf_OK <- ifcadjbmi_eaf_OK[which(ifcadjbmi_eaf_OK$effect_allele_frequency < 0.99),]

summary(ifcadjbmi_eaf_OK$effect_allele_frequency) #worked like a charm.

ifcadjbmi_corrected_eaf <- rbind(ifcadjbmi_eaf_NA, ifcadjbmi_eaf_OK) #13,526,969

##############################
#CURATION OF THE DATA IS DONE#
##############################

#Let's check now if there are any in the MHC region.

ifcadjbmi_mhc <- ifcadjbmi_corrected_eaf[which(as.numeric(ifcadjbmi_corrected_eaf$chromosome) == 6 & as.numeric(ifcadjbmi_corrected_eaf$base_pair_location) >= 26000000 & as.numeric(ifcadjbmi_corrected_eaf$base_pair_location) <= 34000000),]

summary(as.numeric(ifcadjbmi_mhc$chromosome)) #perfect!!
summary(as.numeric(ifcadjbmi_mhc$base_pair_location)) #perfect!!

ifcadjbmi_end <- ifcadjbmi_corrected_eaf[which(!(ifcadjbmi_corrected_eaf$variant%in%ifcadjbmi_mhc$variant)),]

colnames(ifcadjbmi_end) <- c("chromosome", "base_pair_location", "effect_allele", "other_allele", "beta", "standard_error", "effect_allele_frequency", "p_value", "sample_size", "het_p_value", "variant", "chr_pos")

ifcadjbmi_end <- ifcadjbmi_end %>%
  select(variant, chromosome, base_pair_location, effect_allele, other_allele, beta, standard_error, effect_allele_frequency, p_value, sample_size, chr_pos)

#########################
#We can save this data!!#
#########################

fwrite(ifcadjbmi_end, "output/1_curated_data/ifcadjbmi_curated.txt")

