##############
#INTRODUCTION#
##############

#This code makes the graph for the gene-set enrichment.

###################
#Loading libraries#
###################

library(igraph)
library(ggraph)
library(tidygraph)
library(ggplot2)
library(RColorBrewer)
library(data.table)
library(tidyverse)

##############
#Loading data#
##############

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/")

# Load results
full_res <- fread("output/5_enrichment_analyses/2_depict/output/cluster_all/cluster_all_genesetenrichment.txt")
full_res <- full_res[which(full_res$`False discovery rate` == "<0.01" | full_res$`False discovery rate` == "<0.05"),]

edges <- read.table("output/5_enrichment_analyses/2_depict/output/cluster_all/network_plot_network_table.txt", header = TRUE, sep = "\t")
nodes <- read.table("output/5_enrichment_analyses/2_depict/output/cluster_all/network_plot_nodeattributes.txt", header = TRUE, sep = "\t")

#############################################################################################################
#We are going to print the gene centers that are correlated with each other! To do so first we need to match#
#############################################################################################################

# Filter and rename
nodes <- nodes[nodes$Original.gene.set.ID %in% unique(c(edges$Source, edges$Target)), ]
edges <- edges[edges$Source %in% nodes$Original.gene.set.ID & edges$Target %in% nodes$Original.gene.set.ID, ]

nodes <- nodes[which(nodes$Original.gene.set.ID %in% edges$Source & nodes$Original.gene.set.ID %in% edges$Target), ] #exclusively target the nodes that are both, sources and targets to make a connected diagram
edges <- edges[edges$Source %in% nodes$Original.gene.set.ID & edges$Target %in% nodes$Original.gene.set.ID, ] #and find the edges

#Next let's add as sources and targets the original description, which is easier to understand

for(i in seq(1, length(edges$Source))) {
  
  edges$Source[i] <- nodes$Original.gene.set.description[which(nodes$Original.gene.set.ID == edges$Source[i])]
  edges$Target[i] <- nodes$Original.gene.set.description[which(nodes$Original.gene.set.ID == edges$Target[i])]
}

#Also apply it in the nodes, but this one is easy:

nodes$Original.gene.set.ID <- nodes$Original.gene.set.description

##################################
#Now we are ready to make a graph#
##################################

# Build graph
graph <- graph_from_data_frame(edges, directed = FALSE, vertices = nodes)

#As attributes we are setting the pearson correlations, particularly the discrete ones, which can be used as thresholds:

if ("Pearson_correlation_discrete" %in% colnames(edges)) {
  E(graph)$weight <- edges$Pearson_correlation_discrete
} else {
  E(graph)$weight <- 1
}

#The nodes sizes will be based on the on the -logP of the p-value of the enrichment

if ("minuslogten_pval_discrete" %in% colnames(nodes)) {
  V(graph)$size <- nodes$minuslogten_pval_discrete * 20
} else {
  V(graph)$size <- degree(graph)
}

# Updated: increased size to better fit labels
V(graph)$size <- V(graph)$size * 1.5
#V(graph)$size <- V(graph)$size * 0.5


# Add node degree
V(graph)$degree <- degree(graph)

# Assign clusters as factor
V(graph)$cluster <- as.factor(nodes$Cluster.ID)

# Label only gene-sets connected by edges with Pearson correlation discrete > 5 - this will make things easier to parse
strong_edges <- E(graph)[weight > 5]
strong_nodes <- unique(c(ends(graph, strong_edges)))
V(graph)$label <- ifelse(V(graph)$name %in% strong_nodes, V(graph)$name, "")


# Detect gene set type
nodes$gene_set_type <- ifelse(str_detect(nodes$Original.gene.set.description, "REACTOME"), "REACTOME",
                              ifelse(str_detect(nodes$Original.gene.set.description, "PPI"), "PPI", "phenotype"))

# Define palette for types
palette_type <- c("REACTOME" = "#7fcdbb",   # teal for REACTOME
                  "PPI" = "#1d91c0",         # strong blue for PPI
                  "phenotype" = "#bcbddc")   # light blue for phenotypes

# Map type to graph nodes
V(graph)$gene_set_type <- nodes$gene_set_type[match(V(graph)$name, nodes$Original.gene.set.description)]


# Plot with labels inside bubbles using geom_node_label
plotio <- ggraph(graph, layout = "kk", k = 0.01) +
  geom_edge_link(aes(width = weight), color = "lightgrey", alpha = 0.6) +
  geom_node_point(aes(size = size, color = gene_set_type), alpha = 0.9) +
  geom_node_label(aes(label = label, fill = gene_set_type), color = "white", fontface = "bold",
                  repel = TRUE, label.size = 0.1, label.r = unit(0.15, "lines")) + 
  scale_edge_width(range = c(0.1, 1)) +
  scale_color_manual(values = palette_type) +
  scale_fill_manual(values = palette_type) +
  theme_void() +
  #ggtitle("DEPICT Gene Set Network")+
          #subtitle = "Coloring by REACTOME, PPI, or Phenotype categories") +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    plot.subtitle = element_text(hjust = 0.5, size = 11)
  )


#########################################
#Let's plot this!! We have it, I believe#
#########################################

ggsave(
  filename = "N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/manuscript/figures/gene_set_enrichment_depict.svg",
  plot = plotio,
  width = 375,   # mm
  height = 200,  # mm
  units = "mm",
  device = "svg"
)
