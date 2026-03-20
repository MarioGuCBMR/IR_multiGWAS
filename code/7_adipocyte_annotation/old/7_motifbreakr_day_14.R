##############
#INTRODUCTION#
##############

#This code performs a E2G pipeline that tries to identify enhancer-promoter interactions and links them to genes.

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)
library(MotifDb)

###################
#Loading functions#
###################

# Define the pattern to include your databases
pattern <- "(cisbp_1.02|HOCOMOCOv10|HOCOMOCOv11|hPDI|JASPAR_2014|JASPAR_CORE|jaspar2016|jaspar2018|jolma2013|SwissRegulon|UniPROBE)"
motifs <- query(MotifDb, pattern)

get_tf_family_annotation <- function(tf_name_query, motif_list = motifs) {

  #STEP 1: get metadata:
  
  motif_metadata <- as.data.frame(mcols(motif_list))
  
  # Filter for human motifs (species is typically stored under "organism")
  human_motifs_idx <- motif_metadata[which(motif_metadata$organism == "Hsapiens"),]
  
  #Extract info for the gene:
  
  gene_tmp <- human_motifs_idx[which(human_motifs_idx$geneSymbol == tf_name_query),]
  
  family_vect <- paste(gene_tmp$tfFamily, collapse = ";")
  
}

family_cleaner <- function(family_list){
  
  #STEP 0: call igraph
  
  library(igraph)
  
  #STEP 1: let's get the list of families we want:
  
  family_list <- unique(family_list)
  
  clean_terms <- lapply(family_list, function(family_list) {
    items <- unlist(strsplit(family_list, ";"))
    items <- gsub("\\{.*?\\}", "", items)  # strip {codes}
    items <- trimws(items)                # trim whitespace
    items[items != "NA" & items != ""]    # drop NA/empty
  })
  
  # === Step 2: Make all pairwise edges for each row ===
  edges <- do.call(rbind, lapply(clean_terms, function(terms) {
    if (length(terms) >= 2) {
      t(combn(terms, 2))  # all unique pairs
    } else {
      NULL  # ignore singleton rows
    }
  }))
  
  # === Step 3: Build the graph ===
  g <- graph_from_edgelist(edges, directed = FALSE)
  
  # === Step 4: Find synonym clusters ===
  comps <- components(g)$membership
  groups <- split(names(comps), comps)
  
  # === Step 5: Output collapsed synonyms ===
  collapsed <- sapply(groups, function(g) paste(sort(unique(g)), collapse = "; "))
  
  
  
}

##############
#Loading data#
##############

#First let's load the data with enhancer-promoter evidence:

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025")

enhancer_df <- fread("output/7_functional_annotation/2_enhancer_gene_links/2_enhancer_gene_links_and_hic/proxies_with_stare_and_hic_data.txt")

ir_variants <- fread("output/4_ir_loci_discovery/3_novel_variants/4_colocalization/352_variants_with_all_info_plus_fiadjbmi_and_tg_hdl_coloc.txt")

ir_variants <- ir_variants[which(ir_variants$variant%in%enhancer_df$query_snp_rsid),]

qtl_df <- readRDS("output/6_qtl_annotation/conditional_and_fine_mapped_asat_and_vat_eqtls_sqtls.txt")

#######################################################
#Let's get the data for the lipodystrophy-like effects#
#######################################################

#Let's load the variants in this cluster, specifically,

bmi_ns <- fread("output/5_enrichment_analyses/3_go_shifter/input_bmi_ns.txt")

enhancer_df <- enhancer_df[which(enhancer_df$query_snp_rsid%in%bmi_ns$SNP),]
enhancer_df <- enhancer_df[which(enhancer_df$sgbs_day14 != ""),]

###############################
#Let's focus on all days match#
###############################

snps_4_break <- enhancer_df$rsID

library(motifbreakR)
library(SNPlocs.Hsapiens.dbSNP155.GRCh37)
library(BSgenome.Hsapiens.UCSC.hg19)

#Let's extract the variants:

variants <- motifbreakR::snps.from.rsid(rsid = snps_4_break,
                                        dbSNP = SNPlocs.Hsapiens.dbSNP155.GRCh37,
                                        search.genome = BSgenome.Hsapiens.UCSC.hg19)

#Let's run motif breaker for only the human motifs

motif_db <- query(MotifDb, "Hsapiens")  # Filter for Hsapiens motifs - NOT EXCLUSIVE!

# Exclude "stamlabs" motifs
motif_db <- motif_db[!grepl("stamlab", names(motif_db))]

# Further filter for specific motif databases
pattern <- "(cisbp_1.02|HOCOMOCOv10|HOCOMOCOv11|hPDI|JASPAR_2014|JASPAR_CORE|jaspar2016|jaspar2018|jolma2013|SwissRegulon|UniPROBE)"

# Filter motifs based on the pattern matching the source names in the motif names
selected_motifs <- motif_db[grepl(pattern, names(motif_db))]

#We cannot do serial parametization here so let's change things a bit

library(BiocParallel)
register(SerialParam())

results <- motifbreakR(snpList = variants, filterp = TRUE,
                       pwmList = selected_motifs,
                       threshold = 1e-4,
                       method = "ic")

saveRDS(results, "output/7_functional_annotation/3_motif_analyses/motif_analyses_of_day_14_loci_variants.RDS")