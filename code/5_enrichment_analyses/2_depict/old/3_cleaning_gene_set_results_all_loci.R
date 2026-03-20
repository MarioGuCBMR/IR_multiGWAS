##############
#INTRODUCTION#
##############

##Let's check the biology behind the gene-set enrichments

###################
#Loading libraries#
###################

library(data.table)
library(tidyverse)

###################
#Loading functions#
###################

gene_parser <- function(gene_col){
  
  genes_1 <- unlist(str_split(gene_col, "[ (]"))
  genes_1 <- genes_1[which(str_detect(genes_1, "[)]") == FALSE & genes_1 != "")]
  
  return(genes_1)
  
}

##############
#Loading data#
##############

path_2_input <- "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_functional_annotation/"

setwd(path_2_input)

##############################
#First the original gene-sets#
##############################

gene_set_res <- fread("output/3_depict/output/cluster_all/cluster_all_genesetenrichment.txt") 

gene_set_sign <- gene_set_res[which(as.numeric(gene_set_res$`Nominal P value`) < 1e-04),]

table(gene_set_sign$`False discovery rate`)

##########################################
#Let's see if we have any that are shared#
##########################################

network <- fread("output/3_depict/output/cluster_all/network_plot_cluster_results.txt") 

network <- network[which(network$`Nominal P value` < 1e-04),] #198

#Let's remove the lonely guys:

network_group <- as.data.frame(t(table(network$`Cluster ID`)))

network_not_lonely <- network_group[which(network_group$Freq > 2),] #24

#Let's get which are the most relevant:

network_most_relevant <- network[which(network$`Cluster ID`%in%network_not_lonely$Var2),]

#Let's get the centers:

centers <- network_most_relevant[which(network_most_relevant$`Cluster center (boolean)` == 1),] #19 centers

###############################################################
#Let's see which centers have strong correlations between them#
###############################################################

nodes <- fread("output/3_depict/output/cluster_all/network_plot_network_table.txt") 

nodes_of_interest <- nodes[which(nodes$Source%in%centers$`Original gene set ID`),]

centers_matching_centers <- nodes_of_interest[which(nodes_of_interest$Target%in%centers$`Original gene set ID`),]

