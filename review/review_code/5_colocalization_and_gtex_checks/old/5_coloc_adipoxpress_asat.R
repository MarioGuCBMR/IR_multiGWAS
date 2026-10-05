##############
#INTRODUCTION#
##############

#Let's run hyprcoloc for WHRadjBMI against all fertility, PCOS and other related traits first. In another code we will do this for all of them combined.

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

gene_parser = function(gene_data){
  
  #STEP 1: take the string and split it:
  
  vect_=unlist(str_split(gene_data, "[.]"))
  
  id_=vect_[1]
  
  return(id_)
  
}


chr_parser = function(variant_id){
  
  #STEP 1: take the string and split it:
  
  vect_=unlist(str_split(variant_id, "_"))
  
  chr_=vect_[1]
  
  return(chr_)
  
}

pos_parser = function(variant_id){
  
  #STEP 1: take the string and split it:
  
  vect_=unlist(str_split(variant_id, "_"))
  
  pos_=vect_[2]
  
  return(pos_)
  
}

ref_parser = function(variant_id){
  
  #STEP 1: take the string and split it:
  
  vect_=unlist(str_split(variant_id, "_"))
  
  ref_=vect_[3]
  
  return(ref_)
  
}

alt_parser = function(variant_id){
  
  #STEP 1: take the string and split it:
  
  vect_=unlist(str_split(variant_id, "_"))
  
  alt_=vect_[4]
  
  return(alt_)
  
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

#PERFECT!!
#Now let's move to get the proxies for each SNP:

proxies <- fread("/projects/kilpelainen-AUDIT/people/zlc436/IR_GSEM_2025/output/4_ir_loci_discovery/4_proxies/proxies_build_37_and_38_4_282_variants.txt")

yes_vect <- c("A", "G", "C", "T")

proxies <- proxies[which(proxies$query_snp_rsid%in%all_variants$variant),] #3188
proxies <- proxies[which(proxies$ref%in%yes_vect & proxies$alt%in%yes_vect),]

length(which(duplicated(proxies$query_snp_rsid) == FALSE)) #282

#PERFECT! Now we need the full sumstats, of course:

iradjbmi = fread("/projects/kilpelainen-AUDIT/people/zlc436/IR_GSEM_2025/output/2_mv_gwas/fiadjbmi_hdl_tg_common_dwls_curated.txt")

#Last thing: the list of genes:

qtls=readxl::read_xlsx("/projects/kilpelainen-AUDIT/people/zlc436/IR_GSEM_2025/review/manuscript/03062026/Supplementary Tables_20March2026.xlsx", sheet = 12)

#Let's pass this through the same filters:

qtls = as.data.frame(qtls)
qtls = qtls[2:nrow(qtls),] #remove extra first row that explains stuff
colnames(qtls) = qtls[1,]  #change column names
qtls = qtls[2:nrow(qtls),] #remove columns names in the first line now that we have used them

#The columns here did not work as good as I thought due to the table's structure. 
#No matter, we are gonna try changing some of the columns:

colnames(qtls)[1:11] = c("Lead IR SNP", "Variant", "R2", "Chromosome", "Position (b37)", "Position (b38)", "Effect Allele Frequency", "IR allele", "IS allele", "Gene ID", "Gene Symbol")

#####################################################################################
#STEP 1: Let's run liftover for our dataset - the GTEx data is in build 38 after all#
#####################################################################################

chain <- import.chain("/projects/kilpelainen-AUDIT/people/zlc436/IR_GSEM_2025/review/review_raw_data/hg19ToHg38.over.chain")

gr <- GRanges(
  seqnames = paste("chr", iradjbmi$chromosome, sep = ""),
  ranges = IRanges(
    start = iradjbmi$base_pair_location,
    end   = iradjbmi$base_pair_location
  ),
  strand = "*"
)

# 3) Perform liftOver
lifted <- liftOver(gr, chain)
lifted_unlisted = unlist(lifted)

len <- lengths(lifted)
idx_mapped      <- which(len == 1L)

data_converted <- iradjbmi[idx_mapped,]

converted_df <- as.data.frame(lifted_unlisted)

data_converted$pos_38 = converted_df$start #check - the conversion worked ;) - the order is the same since we are using the index 

iradjbmi=data_converted #this worked fantastically

#Let's add this info to the leads:

data_match = iradjbmi[which(iradjbmi$variant%in%all_variants$variant),]
data_match = data_match[order(match(data_match$variant, all_variants$variant)),]
print(length(which(data_match$variant == all_variants$variant)))
all_variants$pos_38=data_match$pos_38

#########################################################
#STEP 2: let's go and loop over the leads and the traits#
#########################################################

#STEP 2.1: let's loop over chromosomes to only load the parquet files once:

for(chr in seq(1, 22)){

   print("CHROMOSOME")
   print(chr)

   #Let's load the qtl data here:

   qtl_tmp = fread(paste("/projects/kilpelainen-AUDIT/people/zlc436/IR_GSEM_2025/review/review_raw_data/marginal_byChr_EURonly/EURonly_marginal_local_eQTL_meta_chr", chr, ".txt", sep = ""))
   qtl_tmp = as.data.frame(qtl_tmp)

   #Let's process the data

   ir_tmp = all_variants[which(as.numeric(all_variants$chromosome) == chr),]
  
   #STEP 2.2: let's loop over the lead variants!

   for(index in seq(1, length(ir_tmp$variant))){ #we are looping through each index SNP

       print("LEAD INDEX")
       print(index)
  
       #And retrieving iradjbmi data:
    
       #STEP 1: let's get the data for each:
  
       lead_ = ir_tmp$variant[index]
       chr_=as.numeric(ir_tmp$chromosome[index])
       start_ = as.numeric(ir_tmp$base_pair_location[index])-500000
       end_ = as.numeric(ir_tmp$base_pair_location[index])+500000
  
      ############################################
      #STEP 2: let's filter the data for iradjbmi#
      ############################################
  
      iradjbmi_tmp = iradjbmi[which(as.numeric(iradjbmi$chromosome) == chr_ & as.numeric(iradjbmi$base_pair_location) > start_ & as.numeric(iradjbmi$base_pair_location) < end_),]

      #Let's align it to iradjbmi+ allele:
    
      iradjbmi_pos <- iradjbmi_tmp
    
      #Let's align to positive:
    
      new_a1 <- ifelse(iradjbmi_pos$beta < 0, iradjbmi_pos$other_allele, iradjbmi_pos$effect_allele)
      new_a2 <- ifelse(iradjbmi_pos$beta < 0, iradjbmi_pos$effect_allele, iradjbmi_pos$other_allele)
      new_beta <- ifelse(iradjbmi_pos$beta < 0, as.numeric(iradjbmi_pos$beta)*(-1), as.numeric(iradjbmi_pos$beta))
    
      iradjbmi_pos$effect_allele = new_a1
      iradjbmi_pos$other_allele = new_a2
      iradjbmi_pos$beta = new_beta
  
      iradjbmi_pos = iradjbmi_pos[order(iradjbmi_pos$p_value),]
      iradjbmi_pos=iradjbmi_pos[which(duplicated(iradjbmi_pos$variant)==FALSE),]
      iradjbmi_pos=iradjbmi_pos[order(iradjbmi_pos$chromosome, iradjbmi_pos$base_pair_location),] #worked well!! We can move one
  
      #######################################################
      #STEP 3: we need to filter for the genes in the region#
      #######################################################
  
      #Here is where we are loop over ENSG on the list from the qtls:

      proxies_4_lead = proxies$rsID[which(proxies$query_snp_rsid == lead_)]

      #Let's get the genes:

      qtls_info = qtls[which(qtls$Variant%in%proxies_4_lead),] #we found them all!! This strategy works

      genes_for_loop = qtls_info$"Gene ID"

      if(is_empty(genes_for_loop)){

        next()

      }

      for(gene_ in genes_for_loop){ #Let's loop over the data:

        print("GENE")
        print(gene_)

        gene_tmp = qtl_tmp[which(qtl_tmp$ENSG == gene_),]

        if(is_empty(gene_tmp$ENSG)){

           next()
        }

        #Let's adapt the data, this is way easier than gtx ones
        gene_tmp$chr_pos = paste("chr", gene_tmp$chr, ":", gene_tmp$pos, sep = "") #this chromosome already has chr

        #Let's get both datasets ready for matching:

        ir_match = iradjbmi_tmp %>%
          dplyr::select(variant, chromosome, base_pair_location, effect_allele, other_allele, minimum_allele_frequency, beta, standard_error, p_value, sample_size)

        ir_match$chr_pos=paste("chr", ir_match$chromosome, ":", ir_match$base_pair_location, sep = "")

        #This one is ready, let's go for qtl_tmp:

        colnames(gene_tmp) = c("chromosome", "base_pair_location", "other_allele", "effect_allele", "variant_id", "ENSG", "gene", "studies", "beta", "standard_error", "p_value", "sample_size", "chr_pos")
        gene_tmp$variant = gene_tmp$chr_pos

        #Let's get the match data:

        ir_match = ir_match[which(ir_match$chr_pos%in%gene_tmp$chr_pos),]
        gene_aligned = trait_aligner(ir_match, gene_tmp)

        #Let's make sure:

        ir_match = ir_match[which(ir_match$chr_pos%in%gene_aligned$chr_pos),]
        gene_aligned = trait_aligner(ir_match, gene_aligned) #helps removing unnecessary duplicates

        #Let's order by chr_pos:

       ir_ordered = ir_match[order(as.numeric(ir_match$chromosome), as.numeric(ir_match$base_pair_location)),]
       gene_ordered = gene_aligned[order(match(gene_aligned$chr_pos, ir_ordered$chr_pos)),]

       print(length(which(ir_ordered$chr_pos == gene_ordered$chr_pos)))
       print(length(which(ir_ordered$effect_allele == gene_ordered$effect_allele)))
    
       #Let's perform hyprocoloc
    
       betas <- as.matrix(as.data.frame(cbind(as.numeric(ir_ordered$beta), 
                                           as.numeric(gene_ordered$beta))))
    
       ses <- as.matrix(as.data.frame(cbind(as.numeric(ir_ordered$standard_error), 
                                         as.numeric(gene_ordered$standard_error))))
    
      traits <- c("iradjbmi", gene_)
      colnames(betas) <- traits
      colnames(ses) <- traits
      rownames(betas) <- ir_ordered$variant
      rownames(ses) <- ir_ordered$variant
    
      res <- hyprcoloc::hyprcoloc(effect.est = betas, effect.se = ses, snp.id = row.names(betas), trait.names = traits)
      res_df <- res$results
      res_df$lead_snp <- lead_
    
      print(res_df)
    
      if(!(exists("coloc_df"))){
      
        coloc_df <- res_df
      
      } else {
      
        coloc_df <- rbind(coloc_df, res_df)
      
    }
    
  } #for loop for gene

 } #for loop for all variants

} #for loop for chromosomes

###############
#We are done!!#
###############

saveRDS(coloc_df, "/projects/kilpelainen-AUDIT/people/zlc436/IR_GSEM_2025/review/review_output/4_colocalization_and_gtex_checks/raw_asat_adipoxpress_coloc_df.RDS")
