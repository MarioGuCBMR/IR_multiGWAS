##############
#INTRODUCTION#
##############

#This code performs a series of analyses to retrieve the motif data from the proxies overlapping the accessible regions at day 4. 
#It is one hell of a code with a lot of preparation steps due to how motifs work. Let's go.

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)
library(MotifDb)
library(igraph)

###################
#Loading functions#
###################

get_tf_family_annotation <- function(tf_name_query, motif_list = tf_df) {

  #STEP 1: get metadata:
  
  motif_metadata <- tf_df 
  
  #Extract info for the gene:
  
  gene_tmp <- motif_metadata[which(motif_metadata$geneSymbol == tf_name_query),]
  
  family_vect <- unlist(str_remove_all(gene_tmp$tfFamily, " "))
  family_vect <- family_vect[which(family_vect != "NA")]
  family_vect <- family_vect[order(family_vect)]

  family_vect <- paste(family_vect, collapse = ";")
  
  return(family_vect)
  
}

family_cleaner <- function(family_list){
   
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
   collapsed <- as.data.frame(collapsed)
   collapsed$synonym <- collapsed$collapsed
   
   collapsed <- separate_rows(collapsed, synonym, sep = ";")
   
   return(collapsed)
   
}

library(dplyr)
library(purrr)

get_tf_family_network <- function(gene, tf_df=tf_df) {
  
  # Subset tf_df for the input gene
  gene_tmp <- tf_df %>% filter(as.character(geneSymbol) == gene)
  
  # If gene is not in tf_df, return NULL
  if (nrow(gene_tmp) == 0) {
    return(NULL)
  }
  
  # Initialize search with the gene's TF families
  family_search <- unique(as.character(gene_tmp$tfFamily))
  family_search <- family_search[!is.na(family_search)]
  
  # Vector to hold all discovered families
  seen_families <- character(0)
  
  # Iteratively expand the search
  repeat {
    prev_length <- length(seen_families)
    
    # Update the seen families
    seen_families <- unique(c(seen_families, family_search))
    
    # Get all genes from these families
    tf_match <- tf_df$geneSymbol[tf_df$tfFamily %in% family_search]
    
    # Subset tf_df to get updated families
    tf_search <- tf_df %>% filter(geneSymbol %in% tf_match)
    family_search <- unique(as.character(tf_search$tfFamily))
    family_search <- family_search[!is.na(family_search)]
    
    # Break if no new families found
    if (length(seen_families) == prev_length) {
      break
    }
  }
  
  # Return all geneSymbols connected through shared families
  final_genes <- unique(tf_df$geneSymbol[tf_df$tfFamily %in% seen_families])
  
  return(final_genes)
}

##############
#Loading data#
##############

#STEP 1: get the data as in MOTIF DF with the same settings as the ones used in motifbreakr:

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025")

motif_db <- query(MotifDb, "Hsapiens")  # Filter for Hsapiens motifs - NOT EXCLUSIVE!

# Exclude "stamlabs" motifs
motif_db <- motif_db[!grepl("stamlab", names(motif_db))]

# Further filter for specific motif databases
pattern <- "(cisbp_1.02|HOCOMOCOv10|HOCOMOCOv11|hPDI|JASPAR_2014|JASPAR_CORE|jaspar2016|jaspar2018|jolma2013|SwissRegulon|UniPROBE)"

# Filter motifs based on the pattern matching the source names in the motif names
selected_motifs <- motif_db[grepl(pattern, names(motif_db))]

#Now let' get a list of unique TFs:

tf_df <- as.data.frame(mcols(selected_motifs))

#################################
#Let's load the results data now#
#################################

results <- readRDS("output/7_functional_annotation/3_motif_analyses/motif_analyses_of_day_4_plus_cebpb_variants.RDS")
enhancer_df <- fread("output/7_functional_annotation/2_enhancer_gene_links/2_enhancer_gene_links_and_hic/proxies_with_stare_and_hic_data.txt")

#STEP 2: Let's explore the data:

motif_day_4_df <- as.data.frame(results, row.names = NULL)

#STEP 3: let's explore the strong data:

motif_day_4_df <- motif_day_4_df[which(motif_day_4_df$effect == "strong"),]

#Let's get the combination of IDs that truly match:

motif_day_4_df$id_ <- paste(motif_day_4_df$SNP_id, "_", motif_day_4_df$REF, "_", motif_day_4_df$ALT, sep  = "")
enhancer_df$id_1 <- paste(enhancer_df$rsID, "_", enhancer_df$alt, "_", enhancer_df$ref, sep = "")
enhancer_df$id_2 <- paste(enhancer_df$rsID, "_", enhancer_df$ref, "_", enhancer_df$alt, sep = "")

motif_day_4_df <- motif_day_4_df[which(motif_day_4_df$id_%in%enhancer_df$id_1 | motif_day_4_df$id_%in%enhancer_df$id_2),]

#Let's remove duplicates now:

motif_day_4_df$id_prot <- paste(motif_day_4_df$SNP_id, "_", motif_day_4_df$geneSymbol, sep = "")
motif_day_4_df <- motif_day_4_df[order(abs(motif_day_4_df$alleleEffectSize), decreasing = TRUE),]
motif_day_4_df <- motif_day_4_df[which(duplicated(motif_day_4_df$id_prot) == FALSE),] #782

#Let's do a match for those in enhancer_df that are in motif_day_4_df

enhancer_df <- enhancer_df[which(enhancer_df$id_1%in%motif_day_4_df$id_ | enhancer_df$id_2%in%motif_day_4_df$id_),]

######################################
#Which of these are affecting CEBPB??#
######################################

check <- as.data.frame(table(motif_day_4_df$geneSymbol))

#################################################################
#Okay, let's go and try to solve this issue, variant per variant#
#################################################################

for(rsid in unique(enhancer_df$query_snp_rsid)){
  
  print(rsid)
  
  #STEP 0: get all proxies that disrupt a TF for that lead

  proxies_ <- enhancer_df$rsID[which(enhancer_df$query_snp_rsid == rsid)]
  
  #rsid <- unique(enhancer_df$query_snp_rsid)[1]
  
  #STEP 1: get the genes for the proxy:
  
  genes_ <- motif_day_4_df$geneSymbol[which(motif_day_4_df$SNP_id%in%proxies_)]
  
  #STEP 2: let's get the families for these TFs! Data on this is terrible, so we need a while loop to search for all synonyms in motifDB. This is what get_tf_family_network does
  
  gene_families <- lapply(genes_, function(g) get_tf_family_network(g, tf_df))
  names(gene_families) <- genes_
  gene_families <- gene_families[ lengths(gene_families) > 0 ]
  
  if(is_empty(gene_families)){
    
    next()
    
  }
  
  #Careful, let's remove those that do not have any family
  
  genes_ = genes_[which(genes_%in%names(gene_families))]
  
  #Let's make a matrix which tells us the matches between the TFs. Code from chatGPT, but works fantastically... we need to learn how to do this,
  
  # Determine which genes share at least one family network gene
  shared_family_matrix <- outer(genes_, genes_, Vectorize(function(g1, g2) {
    length(intersect(gene_families[[g1]], gene_families[[g2]])) > 0
  }))
  
  rownames(shared_family_matrix) <- genes_
  colnames(shared_family_matrix) <- genes_ 
  
  #The next step I knew how to do! With igraph we can transform the matrix into a list of interconnected data:
  #And now divide by graph:
  
  library(igraph)
  
  g <- graph_from_adjacency_matrix(shared_family_matrix, mode = "undirected", diag = FALSE)
  
  # Find connected components
  components <- components(g)
  
  # Split genes into clusters
  gene_clusters <- split(names(components$membership), components$membership)
  
  #GREAT now we have the list of interconnected data. 
  #Let's loop over them, arrange the one that has the strongest score and only take that one:
  
  for(index_list in seq(1, length(gene_clusters))){
    
    #STEP 1: get the genes:
    
    clust_tmp <- gene_clusters[[index_list]]
    
    #STEP 2: get the motif data:
    
    motif_tmp <- motif_day_4_df[which(motif_day_4_df$geneSymbol%in%clust_tmp & motif_day_4_df$SNP_id%in%proxies_),]
    
    #STEP 3: order data:
    
    motif_tmp <- motif_tmp[order(abs(motif_tmp$alleleEffectSize), decreasing = TRUE),]
    
    #STEP 4: get the best one:
    
    motif_tmp <- motif_tmp[1,]
    
    print(motif_tmp$geneSymbol)
    
    if(!(exists("motif_parsed_df"))){
      
      motif_parsed_df <- motif_tmp
      
    } else {
      
      motif_parsed_df <- rbind(motif_parsed_df, motif_tmp)
      
    }
    
  }

}

check <- as.data.frame(table(motif_parsed_df$geneSymbol))

#####################################
#Let's do an enrichment of these TFs#
#####################################

prot_all <- as.data.frame(table(motif_parsed_df$geneSymbol))
prot_all <- prot_all[order(prot_all$Freq, decreasing = TRUE),]
prot_all_vect <- unique(prot_all$Var1)

#Let's do the enrichment:

library(STRINGdb)

prot_all_df <- data.frame(gene = prot_all_vect, stringsAsFactors = FALSE)

string_db <- STRINGdb$new(version = "11.5", species = 9606, score_threshold = 400, input_directory = "")
mapped_genes <- string_db$map(prot_all_df, "gene", removeUnmappedRows = TRUE)
enrichment_all <- string_db$get_enrichment(mapped_genes$STRING_id)
