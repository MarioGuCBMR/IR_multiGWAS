##############
#INTRODUCTION#
##############

#Let's run hyprcoloc for IR KLK1 hit and FIadjBMI. This one in particular is for FUT2 stop-gain variant. It is associated with TG, especially, but also with FIadjBMI and HDL in the right direction.
#It is also a proxy of our lead SNP, which is a synonymous variant for FUT2.

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)
library(TwoSampleMR)
library(hyprcoloc)
library(GenomicRanges)
library(rtracklayer)

###################
#Loading functions#
###################

protein_parser = function(protein_path){
  
  #STEP 1: take the string and split it:
  
  vect_=unlist(str_split(protein_path, "_"))
  
  id_=vect_[1]
  
  return(id_)
  
}


chr_parser = function(variant_id){
  
  #STEP 1: take the string and split it:
  
  vect_=unlist(str_split(variant_id, "[:]"))
  
  chr_=vect_[1]
  
  return(chr_)
  
}

pos_parser = function(variant_id){
  
  #STEP 1: take the string and split it:
  
  vect_=unlist(str_split(variant_id, "[:]"))
  
  pos_=vect_[2]
  
  return(pos_)
  
}

trait_aligner <- function(query_ss, other_ss){
  
  #We are gonna put the example so that we can know what happens:
  
  #whradjbmi_ss <- exp_df_found_pos
  #other_ss <- hdl_005
  #query_ss <- triangulated_fi_all
  #other_ss <- ebmd_ss
  
  #STEP 0: let's run the matching with TwoSampleMR so we need the data:
  
  #Let's first check if we have the effect_allele_frequency column:
  
  check <- which(colnames(query_ss) == "effect_allele_frequency")
  
  if(is_empty(check)){
    
    query_ss$effect_allele_frequency <- NA
    
  }
  
  #Now we can proceed
  
  exposure <- query_ss %>%
    select(chr_pos, chromosome, base_pair_location, effect_allele, other_allele,
           effect_allele_frequency, beta, standard_error, p_value, sample_size, variant)
  
  colnames(exposure) <- c("SNP", "chr.exposure", "pos.exposure", "effect_allele.exposure",
                          "other_allele.exposure", "eaf.exposure", "beta.exposure", "se.exposure",
                          "pval.exposure", "samplesize.exposure", "rsid.exposure")
  
  exposure$exposure <- "fat_distr"
  exposure$id.exposure <- "fat_distr"
  
  #Now with the outcome:
  
  check <- which(colnames(other_ss) == "effect_allele_frequency")
  
  if(is_empty(check)){
    
    other_ss$effect_allele_frequency <- NA
    
  }
  
  #Now we can proceed
  
  outcome <- other_ss %>%
    select(chr_pos, chromosome, base_pair_location, effect_allele, other_allele, effect_allele_frequency, beta, standard_error, p_value, sample_size, variant)
  
  colnames(outcome) <- c("SNP", "chr.outcome", "pos.outcome", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome", "pval.outcome", "samplesize.outcome", "rsid.outcome")
  
  outcome$outcome <- "Outcome"
  outcome$id.outcome <- "Outcome"
  
  #STEP 1: match the data. This will probably fail with tri-allelic SNPs. Here I think we are gonna be OK...
  
  merged_df <- harmonise_data(exposure, outcome, action=3)
  
  merged_df <- merged_df[which(merged_df$remove == FALSE),] #removing incompatible alleles
  #merged_df <- merged_df[which(merged_df$mr_keep),] #removing incompatible alleles
  
  #I checked that all is working great. whradjbmi betas is positive. The rest is great.
  
  #STEP 3: reorganize the dataframe, cuz we need something clean:
  
  other_ss_aligned <- merged_df %>%
    select("rsid.outcome", "chr.outcome", "pos.outcome", "effect_allele.outcome", "other_allele.outcome", "eaf.outcome", "beta.outcome", "se.outcome", "pval.outcome", "samplesize.outcome", "SNP")
  
  colnames(other_ss_aligned) <- c("variant", "chromosome", "base_pair_location",  "effect_allele", "other_allele", "effect_allele_frequency", "beta", "standard_error", "p_value", "sample_size", "chr_pos")
  
  return(other_ss_aligned)
  
}

organizing_common_data <- function(ss_df, common_snps){
  
  ss_end <- ss_df[which(ss_df$chr_pos%in%common_snps),]
  ss_ordered <- ss_end[order(as.numeric(ss_end$chromosome), as.numeric(ss_end$base_pair_location)),]
  
  return(ss_ordered)
  
}

duplicate_remover <- function(aligned_df){
  
  #Let's remove the duplicates, because with the way that the data is built there is no way to separate those...
  
  duplicated_snps=aligned_df$chr_pos[which(duplicated(aligned_df$chr_pos))]
  
  if(is_empty(duplicated_snps)){
    
    return(aligned_df)
    
  }
  
  aligned_df=aligned_df[which(!(aligned_df$chr_pos%in%duplicated_snps)),]
  
  return(aligned_df)
  
}

#############################
#STEP 0: Let's load the data#
#############################

#This time we are going to use the KLK1 as exposure actually! 
#But we require the lead SNP, so let's get a hold of that:

#Let's get the lead SNPs info just in case:

all_variants=readxl::read_xlsx("/projects/kilpelainen-AUDIT/people/zlc436/IR_GSEM_2025/review/manuscript/03062026/Supplementary Tables_20March2026.xlsx", sheet = 3)
all_variants = as.data.frame(all_variants)
all_variants = all_variants[2:nrow(all_variants),] #remove extra first row that explains stuff
colnames(all_variants) = all_variants[1,]  #change column names
all_variants = all_variants[2:nrow(all_variants),] #remove columns names in the first line now that we have used them

#We do not need to many columns AND we have functions ready for a certain subset of variables.

all_variants = all_variants[,1:12]
colnames(all_variants) = c("variant", "chromosome", "base_pair_location", "minimum_allele_frequency", "effect_allele", "other_allele", "novel_or_reported", "bmi_subgroup", "VEP", "beta", "standard_error", "pvalue")
all_variants$sample_size = 2e06 #an approximation, not needed for the actual analyses, just to make the functions work

################
#Let's get KLK1#
################

list_of_files = list.files("/projects/kilpelainen-AUDIT/people/zlc436/IR_GSEM_2025/review/review_raw_data/pqtl_data", full.names=TRUE)
list_of_genes = unlist(str_split(list_of_files, "/"))
list_of_genes = list_of_genes[str_detect(list_of_genes, "_v")]
prot_list = as.character(unlist(sapply(list_of_genes, protein_parser)))

#Good, across the files in this list, we need KLK1 in the cis-region. Meaning: chr19

lead_="rs601338"
prot_ = prot_list[14]
chr_lead=19

#Good:

folder_for_prot = list_of_files[which(str_detect(list_of_files, paste(prot_, "_", sep= "")))]
list_of_chr_files_for_prot = list.files(folder_for_prot, full.names=TRUE)
chr_2_detect= paste("chr", chr_lead, "_", sep = "") #with this we remove weird stuff!
our_chr_prot_file = list_of_chr_files_for_prot[which(str_detect(list_of_chr_files_for_prot, chr_2_detect))] #worked like a fucking charm
print(our_chr_prot_file)

qtl_tmp = fread(our_chr_prot_file)

#Let's narrow down to the KLK1 locus, the data is in build 38 (thought b37 can be extracted, but it is better to do this later)

pos_lead=48703417
     
start_ = as.numeric(pos_lead)-500000
end_ = as.numeric(pos_lead)+500000

#Let's get that info 

klk1_locus = qtl_tmp[which(as.numeric(qtl_tmp$GENPOS) >= start_ & as.numeric(qtl_tmp$GENPOS) <= end_),]
summary(as.numeric(klk1_locus$GENPOS))

#Let's get the chromosome and the positions in build 37 for FIadjBMI

klk1_locus$pos = as.numeric(as.character(unlist(sapply(klk1_locus$ID, pos_parser)))) 

#We got the data! Let's realign the KLK1-increasing hits, since KLK1 is the exposure

new_a1=ifelse(as.numeric(klk1_locus$BETA) < 0, klk1_locus$ALLELE0, klk1_locus$ALLELE1)
new_a2=ifelse(as.numeric(klk1_locus$BETA) < 0, klk1_locus$ALLELE1, klk1_locus$ALLELE0)
new_beta=ifelse(as.numeric(klk1_locus$BETA) < 0, as.numeric(klk1_locus$BETA)*(-1), as.numeric(klk1_locus$BETA))
new_eaf=ifelse(as.numeric(klk1_locus$BETA) < 0, 1-as.numeric(klk1_locus$A1FREQ), as.numeric(klk1_locus$A1FREQ))

klk1_pos = klk1_locus
klk1_pos$ALLELE1 = new_a1
klk1_pos$ALLELE0 = new_a2
klk1_pos$BETA = new_beta
klk1_pos$A1FREQ = new_eaf

#Great!! The data makes sense; conversions worked. Next we are gonna get the right column for the exposure

colnames(klk1_pos) = c("chromosome", "pos_38", "id", "other_allele", "effect_allele", "effect_allele_frequency", "info", "sample_size", "test", "beta", "standard_error", "chisq", "log10p", "extra", "base_pair_location")
   
klk1_pos$p_value = 10^-as.numeric(klk1_pos$log10p) #checked in OTARGEN - conversion is correct
klk1_pos$chr_pos = paste("chr", klk1_pos$chromosome, ":", klk1_pos$base_pair_location, sep = "")
klk1_pos$variant = klk1_pos$chr_pos

#Maybe we should remove rare alleles here, let's see...

summary(as.numeric(klk1_pos$effect_allele_frequency)) #we do have rare ones.

klk1_pos = klk1_pos[which(as.numeric(klk1_pos$effect_allele_frequency) > 0.01 & as.numeric(klk1_pos$effect_allele_frequency) < 0.99)]

summary(as.numeric(klk1_pos$effect_allele_frequency)) #we do have rare ones.

#####################################################################################
#STEP 1: Let's run liftover for our dataset - the GTEx data is in build 38 after all#
#####################################################################################

fiadjbmi=fread("/maps/projects/kilpelainen-AUDIT/people/zlc436/IR_PLACO_2026/output/1_curated_data/fiadjbmi_eur.txt") #latest release of FIadjBMI data.
fiadjbmi_match = fiadjbmi[which(fiadjbmi$chr_pos%in%klk1_pos$chr_pos),] #quite good match!! 
klk1_match = klk1_pos[which(klk1_pos$chr_pos%in%fiadjbmi_match$chr_pos),] #great, we have triallelic SNPs, nothing we cannot solve with our harmonise_data function. Time to use it.
fiadjbmi_match$variant=fiadjbmi_match$chr_pos

fiadjbmi_aligned = trait_aligner(klk1_match, fiadjbmi_match) #helps removing unnecessary duplicates

klk1_match$id = paste(klk1_match$chr_pos, "_", klk1_match$effect_allele, "_", klk1_match$other_allele, sep = "")
fiadjbmi_aligned$id = paste(fiadjbmi_aligned$chr_pos, "_", fiadjbmi_aligned$effect_allele, "_", fiadjbmi_aligned$other_allele, sep = "")

klk1_aligned=klk1_match[which(klk1_match$id%in%fiadjbmi_aligned$id),] #this worked!! Perfect match

klk1_aligned = klk1_aligned[order(as.numeric(klk1_aligned$chromosome), as.numeric(klk1_aligned$base_pair_location)),]
fiadjbmi_aligned = fiadjbmi_aligned[order(match(fiadjbmi_aligned$chr_pos, klk1_aligned$chr_pos)),]

print(length(which(klk1_aligned$id == fiadjbmi_aligned$id)))
print(length(which(klk1_aligned$effect_allele == fiadjbmi_aligned$effect_allele)))

#Surprisngly, the first variant is missanotated in the original GWAS! The base-pair location across build has a mismatch.
#We will manually remove it

klk1_aligned=klk1_aligned[2:length(klk1_aligned$base_pair_location),]
fiadjbmi_aligned=fiadjbmi_aligned[2:length(fiadjbmi_aligned$base_pair_location),]

#This works! we are ready to run HypRcoloc

betas <- as.matrix(as.data.frame(cbind(as.numeric(klk1_aligned$beta), 
                                           as.numeric(fiadjbmi_aligned$beta))))
    
ses <- as.matrix(as.data.frame(cbind(as.numeric(klk1_aligned$standard_error), 
                                         as.numeric(fiadjbmi_aligned$standard_error))))
    
traits <- c("KLK1", "FIadjBMI")
colnames(betas) <- traits
colnames(ses) <- traits
rownames(betas) <- klk1_aligned$variant
rownames(ses) <- klk1_aligned$variant
    
res <- hyprcoloc::hyprcoloc(effect.est = betas, effect.se = ses, snp.id = row.names(betas), trait.names = traits)
res_df <- res$results
res_df$lead_snp <- lead_
    
#Let's check the plots:

saveRDS(res_df, "/projects/kilpelainen-AUDIT/people/zlc436/IR_GSEM_2025/review/review_output/4_colocalization_and_gtex_checks/klk1_fiadjbmi_rs601338_df.RDS")
