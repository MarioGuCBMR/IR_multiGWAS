##############
#INTRODUCTION#
##############

#This code produces proxies for the 282 IR variants! We use Haploreg 4.2 website to obtain them! 

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)

###################
#Loading functions#
###################

get_snp_positions <- function(build = "GRCh38") {
  if (build == "GRCh38") {
    mart <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp")  # Default Ensembl (GRCh38)
  } else if (build == "GRCh37") {
    mart <- useMart("ENSEMBL_MART_SNP", dataset = "hsapiens_snp",
                    host = "grch37.ensembl.org", path = "/biomart/martservice", 
                    port = 80)  # Ensembl archive for GRCh37
  } else {
    stop("Invalid genome build. Use 'GRCh37' or 'GRCh38'.")
  }
  
  # Query rsIDs and genomic positions
  snp_positions <- getBM(attributes = c("refsnp_id", "chr_name", "chrom_start"),
                         filters = "snp_filter",
                         values = rsids,
                         mart = mart)
  
  snp_positions$build <- build  # Add build info
  return(snp_positions)
}


##############
#Loading data#
##############

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/")

ir_variants <- fread("output/4_ir_loci_discovery/2_ir_variants/282_ir_fiadjbmi_hdl_tg_dwls_snps_w_cpassoc.txt")

done_proxies <- fread("output/4_ir_loci_discovery/4_proxies/proxies_build_37_and_38_4_282_variants.txt")

non_ir_variants <- ir_variants[which(!(ir_variants$variant%in%done_proxies$query_snp_rsid)),] #70 proxies! perfect

###########################################################
#Let's save this data so that we can make haploreg read it#
###########################################################

dir.create("output/4_ir_loci_discovery/4_proxies")

fwrite(as.data.frame(non_ir_variants$variant), "output/4_ir_loci_discovery/4_proxies/70_non_ir_4_proxies.txt", quote=FALSE, row.names = FALSE, col.names = FALSE, sep = " ")

#########################################################################################################
#Let's load all the proxies dataframes, see if they have differences and add those that might be missing#
#########################################################################################################

proxies_1 <- fread("output/4_ir_loci_discovery/4_proxies/70_proxies_core_15.txt", fill = TRUE)

#Are all the query SNPs there?

length(which(non_ir_variants$variant%in%proxies_1$query_snp_rsid)) #yes
length(which(non_ir_variants$variant%in%proxies_1$rsID)) #yes

#All of them are found in this version!! 
#Fantastic, let's use this dataframe as a reference:

########################################################################################
#Formatting the data so that we can have map of leads and proxies that we can use later#
########################################################################################

proxies_clean <- proxies_1 %>%
  select(rsID, chr, pos_hg38, query_snp_rsid, r2)

#We have here some mismtaches from how the data is being cleaned:

proxies_clean <- proxies_clean[which(str_detect(proxies_clean$rsID, "rs")),] #We only lose 16

#Let's ensure we have not fucked up anything:

length(which(ir_variants$variant%in%proxies_clean$query_snp_rsid)) #yes
length(which(ir_variants$variant%in%proxies_clean$rsID)) #yes

#Amazing, let's find build 37:

proxies_clean$pos_hg19 <- NA

for(index in seq(1, length(proxies_clean$rsID))){
  
    print(proxies_clean$rsID[index])
    
    variant= proxies_clean$rsID[index]
    query = proxies_clean$query_snp_rsid[index]
      
    skip_to_next <- FALSE
    
    # Let's make sure we retrieve the build37 data...
    
    variant_info <- tryCatch(otargen::variantInfo(variant), error = function(e) { skip_to_next <<- TRUE})
    
    if(skip_to_next){
      
      skip_to_next <- FALSE
      next()
      
    }
    
    chr_37 <- variant_info$chromosomeB37
    bp_37 = variant_info$positionB37
    
    #Then let's update the data:
    
    proxies_clean$pos_hg19[index] <- bp_37
  
}

############################################
#Let's add some more positions just in case#
############################################

library(SNPlocs.Hsapiens.dbSNP155.GRCh37)
library(GenomicRanges)
library(biomaRt)

proxies_missing <- proxies_clean[which(is.na(proxies_clean$pos_hg19)),] #104
proxies_good <- proxies_clean[which(is.na(proxies_clean$pos_hg19)==FALSE),] 

rsids <- proxies_missing$rsID

# Query rsIDs and their genomic positions
snp_positions_grch38 <- get_snp_positions(build = "GRCh38")
snp_positions_grch37 <- get_snp_positions(build = "GRCh37")

#We have almost all of them:

for(index_missing in seq(1, length(proxies_missing$rsID))){
  
  #STEP 1: get variant:
  
  rsid_ <- proxies_missing$rsID[index_missing]
  
  #STEP 2: get the data
  
  snps_38 <- snp_positions_grch38[which(snp_positions_grch38$refsnp_id == rsid_),]
  snps_37 <- snp_positions_grch37[which(snp_positions_grch37$refsnp_id == rsid_),]
  
  proxies_missing$pos_hg38[index_missing] <- ifelse(is_empty(snps_38$chrom_start), NA, snps_38$chrom_start)
  proxies_missing$pos_hg19[index_missing] <- ifelse(is_empty(snps_37$chrom_start), NA, snps_37$chrom_start)
  proxies_missing$chr[index_missing] <- ifelse(is_empty(snps_37$chr_name), NA, snps_37$chr_name)
  
}

proxies_missing_clean <- proxies_missing[which(is.na(proxies_missing$pos_hg19) == FALSE),]

proxies_solved <- rbind(proxies_good, proxies_missing_clean)

#########################################################################################################################
#Let's update the data of proxies_1 with the info that we now have, because we might want to add some more relevant info#
#########################################################################################################################

proxies_end <- proxies_1 %>%
  dplyr::select(rsID, chr, pos_hg38, query_snp_rsid, alt, ref, EUR, r2)

proxies_end <- proxies_end[which(str_detect(proxies_end$rsID, "rs")),] #We only lose 16

proxies_end <- proxies_end[which(proxies_end$rsID%in%proxies_solved$rsID),]

#Let's order the data...

proxies_end <- proxies_end[order(match(proxies_end$rsID, proxies_solved$rsID)),]

length(which(proxies_end$rsID == proxies_solved$rsID))

proxies_end$pos_hg19 <- proxies_solved$pos_hg19

##################
#Last thing to do#
##################

#INDELS

yes_vect <- c("A", "C","G", "T")

proxies_end_ <- proxies_end[which(proxies_end$alt%in%yes_vect & proxies_end$ref%in%yes_vect),]

#Let's clean the array thingie:

chr_ <- unlist(str_split(proxies_end_$chr, "Array"))
chr_ <- chr_[which(chr_ != "")]

proxies_end_$chr <- chr_

length(which(ir_variants$variant%in%proxies_end_$query_snp_rsid)) #yes
length(which(ir_variants$variant%in%proxies_end_$rsID)) #yes

#####################
#LET'S SAVE THE DATA#
#####################

fwrite(proxies_end_, "output/4_ir_loci_discovery/4_proxies/proxies_build_37_and_38_4_70_variants.txt")

check <- data.table::fread("output/4_ir_loci_discovery/4_proxies/proxies_build_37_and_38_4_70_variants.txt")

colnames(check) <- c("rsID", "chr", "pos_hg38", "query_snp_rsid", "alt", "ref", "EUR", "r2", "pos_hg19")

#Careful there is weird ones we have found ad-hoc

check_clean <- check[which(is.na(check$pos_hg19)==FALSE),]
check_weird <- check[which(is.na(check$pos_hg19)),] #0

#We can stop here!!! These work without an issue!

# check_clean <- check_clean[which(!(check_clean$rsID%in%check_weird$rsID)),] #from repeated analyses looking for them we put them twice...
# 
# colnames(check_weird) <- c("rsID", "chr", "pos_hg38", "alt", "ref", "EUR", "r2", "pos_hg19", "query_snp_rsid")
# 
# #those with EUR=0.3 rs13316065
# #EUR=0.04 is rs9813811
# #Chr 12 is rs11171739
# 
# check_weird$query_snp_rsid[which(check_weird$EUR == 0.3)] <- "rs13316065"
# check_weird$query_snp_rsid[which(check_weird$EUR == 0.04)] <- "rs9813811"
# check_weird$query_snp_rsid[which(check_weird$chr == 12)] <- "rs11171739"
# 
# #SOLVED!!
# 
# final_check <- rbind(check_clean, check_weird)
# 
# data.table::fwrite(final_check, "output/4_ir_loci_discovery/4_proxies/proxies_build_37_and_38_4_70_variants.txt")
# 
