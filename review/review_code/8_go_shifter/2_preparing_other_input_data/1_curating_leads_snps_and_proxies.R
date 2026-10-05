##############
#INTRODUCTION#
##############

#This code will prepare the proxies input for Lotta et al 53 IR loci for Go-shifter.

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)

###################
#Loading functions#
###################

chr_parser = function(chr_pos){
  
  #STEP 1: strip away by :
  
  info_vect = unlist(str_split(chr_pos, ":"))
  
  chr = info_vect[1]

  #STEP 3: clean chr
  
  chr_clean = unlist(str_split(chr, "chr"))[2]
  
  return(as.numeric(chr_clean))
    
}

pos_parser = function(chr_pos){
  
  #STEP 1: strip away by :
  
  info_vect = unlist(str_split(chr_pos, ":"))
  
  pos_clean = info_vect[2]
  
  return(as.numeric(pos_clean))
  
}

#####################
#Let's read the data#
#####################

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/")

lotta <- readxl::read_xlsx("raw_data/previous_loci/lotta_et_al.xlsx")
lotta_clean = lotta[which(lotta$SNP != "rs6822892"),] #removing tri-alleic SNPs for which we cannot recover proxies

for(index in seq(1, length(lotta_clean$SNP))){
  
  print(index)
  
  #STEP 1: get the rsID
  
  rsid = lotta_clean$SNP[index]
  
  #STEP 2: get the proxies
  
  proxies_tmp = LDlinkR::LDproxy(rsid, pop = "EUR", token="04cad4ca4374") #the rest of default settings are build 37 and r2, just what we need 
  
  proxies_tmp$lead = rsid
  
  if(!(exists("proxies_df"))){
    
    proxies_df = proxies_tmp
    
  } else {
    
    proxies_df = rbind(proxies_df, proxies_tmp)
    
  }
  
}


#########################################
#Let's format and save the lead variants#
#########################################

leads = proxies_df[which(proxies_df$Distance == 0),] #52
leads$chromosome = as.numeric(as.character(unlist(sapply(leads$Coord, chr_parser))))
leads$base_pair_location = as.numeric(as.character(unlist(sapply(leads$Coord, pos_parser))))

leads <- leads %>%
  select(RS_Number, chromosome, base_pair_location)

colnames(leads) <- c("SNP","Chrom","BP")
leads$Chrom <- paste0("chr",leads$Chrom)

dir.create("review/review_output/5_enrichment_analyses/")
dir.create("review/review_output/5_enrichment_analyses/1_go_shifter")
dir.create("review/review_output/5_enrichment_analyses/1_go_shifter/input")
dir.create("review/review_output/5_enrichment_analyses/1_go_shifter/output")
fwrite(leads, "review/review_output/5_enrichment_analyses/1_go_shifter/input/lotta_input.txt", sep = "\t")

##########################
#Let's get the proxy data#
##########################

check = fread("output/5_enrichment_analyses/3_go_shifter/ld.txt", sep = "\t") #we have all the data in our proxies already! Easy set up

proxies_r2 = proxies_df[which(as.numeric(proxies_df$R2) >= 0.8),]

#Let's add the chr and bp properly:

proxies_r2$chromosome = as.numeric(as.character(unlist(sapply(proxies_r2$Coord, chr_parser))))
proxies_r2$base_pair_location = as.numeric(as.character(unlist(sapply(proxies_r2$Coord, pos_parser))))

proxies_r2$lead_chr = NA
proxies_r2$lead_pos = NA

for(lead_ in unique(proxies_r2$lead)){
  
  #STEP 1: get the lead SNP:
  
  chr_ = proxies_r2$chromosome[which(proxies_r2$RS_Number == lead_)]
  pos_ = proxies_r2$base_pair_location[which(proxies_r2$RS_Number == lead_)]
  
  #STEP 2: add the info on the column
  
  proxies_r2$lead_chr[which(proxies_r2$lead == lead_)] = chr_
  proxies_r2$lead_pos[which(proxies_r2$lead == lead_)] = pos_
  
  
} #this worked!! 

#Now we have all the data, we just need to shape it!

proxies_lead_info_37 <- proxies_r2 %>% 
  dplyr::select("lead_chr","lead_pos","lead","chromosome","base_pair_location","RS_Number","Distance","R2","Dprime")

colnames(proxies_lead_info_37) <- c("ChromA","PosA","RsIdA","ChromB","PosB","RsIdB","Distance","RSquared","Dprime") 

#Only thing to do is add Chr at the beginning of the chromosome and we are done!

proxies_lead_info_37$ChromA = paste("Chr", proxies_lead_info_37$ChromA, sep = "")
proxies_lead_info_37$ChromB = paste("Chr", proxies_lead_info_37$ChromB, sep = "")

fwrite(proxies_lead_info_37, "review/review_output/5_enrichment_analyses/1_go_shifter/input/lotta_ld.txt", sep = "\t")
