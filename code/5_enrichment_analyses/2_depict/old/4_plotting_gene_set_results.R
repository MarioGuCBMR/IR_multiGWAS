##############
#INTRODUCTION#
##############

#This code makes networks of gene set results:

###################
#Loading libraries#
###################

library(tidyverse)
library(data.table)
library(igraph)
library(ggraph)
library(ggplot2)

###############
#Load the data#
###############

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/")

full_res <- fread("output/5_enrichment_analyses/2_depict/output/cluster_all/cluster_all_genesetenrichment.txt")
full_res <- full_res[which(full_res$`False discovery rate` == "<0.01" | full_res$`False discovery rate` == "<0.05"),]

edges <- read.table("output/5_enrichment_analyses/2_depict/output/cluster_all/network_plot_network_table.txt", header = TRUE, sep = "\t")
nodes <- read.table("output/5_enrichment_analyses/2_depict/output/cluster_all/network_plot_nodeattributes.txt", header = TRUE, sep = "\t")

#Let's take the best hits:

#nodes_best <- nodes[which(nodes$False.discovery.rate == "<0.01"),]
#nodes_best <- nodes_best[order(as.numeric(nodes_best$Nominal.P.value) < 1e-04),]
# nodes_centers <- nodes_best[which(nodes_best$Cluster.center..boolean. == "True"),] #43 centers
# nodes_centers <- nodes_centers[nodes_centers$Original.gene.set.ID %in% unique(c(edges$Source, edges$Target)), ]
# nodes_centers <- nodes_centers[order(as.numeric(nodes_centers$Nominal.P.value)),]
#nodes_centers <- nodes_centers[1:20,]

#Let's take as nodes only these:

nodes <- nodes[nodes$Original.gene.set.ID %in% unique(c(edges$Source, edges$Target)), ]

# Ensure edges only include those in nodes
edges <- edges[edges$Source %in% nodes$Original.gene.set.ID & edges$Target %in% nodes$Original.gene.set.ID, ]
nodes <- nodes[which(nodes$Original.gene.set.ID%in%edges$Source &  nodes$Original.gene.set.ID%in%edges$Target), ]
edges <- edges[edges$Source %in% nodes$Original.gene.set.ID & edges$Target %in% nodes$Original.gene.set.ID, ]

#Let's change the names:

for(i in seq(1, length(edges$Source))){
  
  #Let's get to it:
  
  def_1 <- nodes$Original.gene.set.description[which(nodes$Original.gene.set.ID == edges$Source[i])]
  def_2 <- nodes$Original.gene.set.description[which(nodes$Original.gene.set.ID == edges$Target[i])]
  
  edges$Source[i] <- def_1
  edges$Target[i] <- def_2
  
}

nodes$Original.gene.set.ID <- nodes$Original.gene.set.description

# Now create igraph object
graph <- graph_from_data_frame(edges, directed = FALSE, vertices = nodes)

# Set edge weights (if available)
if("Pearson_correlation" %in% colnames(edges)) {
  E(graph)$weight <- edges$Pearson_correlation
} else {
  E(graph)$weight <- 1  # Default weight if missing
}

# Set node size based on a relevant attribute (adjust as needed)
if("minuslogten_pval_discrete" %in% colnames(nodes)) {
  V(graph)$size <- nodes$minuslogten_pval_discrete * 50  # Scale for better visualization
} else {
  V(graph)$size <- degree(graph) # Default to degree if no size column
}


# Plot the network using ggraph
 # ggraph(graph, layout = "fr") +
 #   geom_edge_link(aes(width = weight), alpha = 0.5) +
 #   geom_node_point(aes(size = size)) +
 #   geom_node_text(aes(label = name), repel = TRUE, size = 3) +
 #   scale_edge_width(range = c(0.2, 2)) +
 #   theme_void() +
 #   ggtitle("DEPICT Gene Set Network")


plot <- ggraph(graph, layout = "kk") +  # Using Kamada-Kawai layout
  geom_edge_link(aes(width = weight, color = weight), alpha = 0.5) +
  geom_node_point(aes(size = size, color = as.factor(nodes$Cluster.ID))) +
  geom_node_text(aes(label = name), repel = TRUE, size = 3, fontface = "bold") +
  scale_edge_width(range = c(0.2, 2)) +
  scale_edge_color_gradient(low = "blue", high = "red") +
  scale_color_manual(values = rainbow(length(unique(nodes$Cluster.ID)))) +
  theme_void() +
  ggtitle("DEPICT Gene Set Network") +
  theme(legend.position = "right")

tiff("output/5_enrichment_analyses/2_depict/output/gene_set_results.tiff", res=300, height=5000, width=7000)
plot
dev.off()

###################################################
#Let's check the new code that we are working with#
###################################################

library(igraph)
library(ggraph)
library(tidygraph)
library(ggplot2)
library(RColorBrewer)

setwd("N:/SUN-CBMR-Kilpelainen-Group/Mario_Tools/IR_GSEM_2025/")

full_res <- fread("output/5_enrichment_analyses/2_depict/output/cluster_all/cluster_all_genesetenrichment.txt")
full_res <- full_res[which(full_res$`False discovery rate` == "<0.01" | full_res$`False discovery rate` == "<0.05"),]

edges <- read.table("output/5_enrichment_analyses/2_depict/output/cluster_all/network_plot_network_table.txt", header = TRUE, sep = "\t")
nodes <- read.table("output/5_enrichment_analyses/2_depict/output/cluster_all/network_plot_nodeattributes.txt", header = TRUE, sep = "\t")

# Filter nodes to keep only those used in edges
nodes <- nodes[nodes$Original.gene.set.ID %in% unique(c(edges$Source, edges$Target)), ]
edges <- edges[edges$Source %in% nodes$Original.gene.set.ID & edges$Target %in% nodes$Original.gene.set.ID, ]
nodes <- nodes[which(nodes$Original.gene.set.ID %in% edges$Source & nodes$Original.gene.set.ID %in% edges$Target), ]
edges <- edges[edges$Source %in% nodes$Original.gene.set.ID & edges$Target %in% nodes$Original.gene.set.ID, ]

# Replace node IDs in edges with readable names
for(i in seq_len(nrow(edges))) {
  edges$Source[i] <- nodes$Original.gene.set.description[which(nodes$Original.gene.set.ID == edges$Source[i])]
  edges$Target[i] <- nodes$Original.gene.set.description[which(nodes$Original.gene.set.ID == edges$Target[i])]
}

# Replace IDs in node table as well
nodes$Original.gene.set.ID <- nodes$Original.gene.set.description

# Build graph
graph <- graph_from_data_frame(edges, directed = FALSE, vertices = nodes)

# Set edge weights
E(graph)$weight <- if ("Pearson_correlation" %in% colnames(edges)) {
  edges$Pearson_correlation
} else {
  1
}

# Set node size based on enrichment
V(graph)$size <- if ("minuslogten_pval_discrete" %in% colnames(nodes)) {
  nodes$minuslogten_pval_discrete * 50
} else {
  degree(graph)
}

# Add node attributes: degree + cluster
V(graph)$degree <- degree(graph)
V(graph)$cluster <- as.factor(nodes$Cluster.ID)

# Label filtering: only label highly connected or enriched nodes
#V(graph)$label_flag <- V(graph)$degree >= 6 | V(graph)$size >= quantile(V(graph)$size, 0.9)

# Filter weak edges (optional, adjust threshold)
graph_filtered <-graph 

# Define cluster colors
cluster_ids <- sort(unique(nodes$Cluster.ID))
cluster_colors <- setNames(brewer.pal(n = max(3, min(length(cluster_ids), 8)), "Set2")[1:length(cluster_ids)],
                           cluster_ids)

# Final network plot
plot <- ggraph(graph_filtered, layout = "kk") +
  geom_edge_link(aes(width = weight, color = weight), alpha = 0.6) +
  geom_node_point(aes(size = size, color = cluster, shape = highlight)) +
  geom_node_text(aes(label = ifelse(label_flag, name, "")), repel = TRUE, size = 3, fontface = "bold") +
  scale_edge_width(range = c(0.5, 2)) +
  scale_edge_color_gradient2(low = "lightblue", mid = "orchid", high = "firebrick", midpoint = 0.5) +
  scale_color_manual(values = cluster_colors) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17)) +
  theme_minimal(base_size = 13) +
  ggtitle("DEPICT Gene Set Network",
          subtitle = "Clusters of gene sets and pathways based on DEPICT; edge color = correlation") +
  labs(size = "-log10(P)", color = "Cluster", edge_color = "Pearson correlation") +
  theme(legend.position = "right")

# Print
print(plot)

#############################

# Filter hub nodes: label only the top X by degree
#V(graph)$degree <- degree(graph)
top_n <- 50  # Change this for more/less labels
V(graph)$label <- ifelse(V(graph)$degree >= sort(V(graph)$degree, decreasing = TRUE)[top_n], V(graph)$name, "")

# Define the plot
plot <- ggraph(graph, layout = "fr") +
  geom_edge_link(color = "grey80", alpha = 0.3, width = 0.5) +
  geom_node_point(color = "black", size = 3) +
  geom_node_text(aes(label = label), repel = TRUE, size = 3, fontface = "plain") +
  theme_void() +
  ggtitle("Minimalist Gene Set Network") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  )
