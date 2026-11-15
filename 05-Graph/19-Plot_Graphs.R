# =============================================================================
# 19-Plot_Graphs.R
#
# Plot the calibration of the multi-block networks and the networks before
# and after swimming, using a common layout and colouring edges by whether
# they are found before, after, or both before and after swimming.
#
# Input:  Results/Network/stab.cluster.21.prot.wv{1,2}.rds
#         Results/Annotation/cluster.annotation.rds
# Output: Results/Network/*.pdf
# =============================================================================

rm(list = ls())
setwd(Sys.getenv("PISCINA_DIR", unset = ".")) # project root (see README)

library(dplyr)
library(sharp)
library(igraph)


# Read data ---------------------------------------------------------------
# Read network results
stab.cluster.21.prot.wv1 <- readRDS("Results/Network/stab.cluster.21.prot.wv1.rds")
stab.cluster.21.prot.wv2 <- readRDS("Results/Network/stab.cluster.21.prot.wv2.rds")

# Read annotation
cluster.annotation <- readRDS("Results/Annotation/cluster.annotation.rds")

# Modify annotation displayed in networks
cluster.annotation <- ifelse(cluster.annotation == "Unidentified chlorinated compounds",
  "Unidentified\nchlorinated compounds",
  cluster.annotation
)


# Extract graph objects
graph.cluster.21.prot.wv1 <- Graph(stab.cluster.21.prot.wv1, satellites = TRUE)
graph.cluster.21.prot.wv2 <- Graph(stab.cluster.21.prot.wv2, satellites = TRUE)


# Calibration plots -------------------------------------------------------
pdf("Results/Network/calibration_multiomics_sPLS_selected_centroids_21_prepost.pdf",
  height = 8, width = 12
)
par(mfrow = c(2, 3), mar = c(7, 6, 10, 6))

CalibrationPlot(stab.cluster.21.prot.wv1, block_id = 1)
mtext("Metabolomics", line = 6, cex = 1.5)
mtext("A", line = 6, cex = 2, adj = -0.32)
CalibrationPlot(stab.cluster.21.prot.wv1, block_id = 2)
mtext("Metabolomics-proteomics", line = 6, cex = 1.5)
CalibrationPlot(stab.cluster.21.prot.wv1, block_id = 3)
mtext("Proteomics", line = 6, cex = 1.5)


CalibrationPlot(stab.cluster.21.prot.wv2, block_id = 1)
mtext("Metabolomics", line = 6, cex = 1.5)
mtext("B", line = 6, cex = 2, adj = -0.32)
CalibrationPlot(stab.cluster.21.prot.wv2, block_id = 2)
mtext("Metabolomics-proteomics", line = 6, cex = 1.5)
CalibrationPlot(stab.cluster.21.prot.wv2, block_id = 3)
mtext("Proteomics", line = 6, cex = 1.5)
dev.off()

# Characterise networks ---------------------------------------------------
edges.wv1.df <- as.data.frame(do.call(rbind, lapply(1:length(E(graph.cluster.21.prot.wv1)), function(x) ends(graph.cluster.21.prot.wv1, x))))
edges.wv2.df <- as.data.frame(do.call(rbind, lapply(1:length(E(graph.cluster.21.prot.wv2)), function(x) ends(graph.cluster.21.prot.wv2, x))))

edges.wv1.annotated.df <- as.data.frame(apply(edges.wv1.df, 2, function(x) ifelse(is.na(cluster.annotation[x]), x, cluster.annotation[x])))
edges.wv2.annotated.df <- as.data.frame(apply(edges.wv2.df, 2, function(x) ifelse(is.na(cluster.annotation[x]), x, cluster.annotation[x])))

colnames(edges.wv1.df) <- colnames(edges.wv2.df) <- colnames(edges.wv1.annotated.df) <- colnames(edges.wv2.annotated.df) <- c("Start", "End")

edges.wv1 <- paste(edges.wv1.df$Start, edges.wv1.df$End, sep = " - ")
edges.wv2 <- paste(edges.wv2.df$Start, edges.wv2.df$End, sep = " - ")

edges.wv1.annotated <- paste(edges.wv1.annotated.df$Start, edges.wv1.annotated.df$End, sep = " - ")
edges.wv2.annotated <- paste(edges.wv2.annotated.df$Start, edges.wv2.annotated.df$End, sep = " - ")

edges.df <- data.frame(Edge = unique(c(edges.wv1, edges.wv2)))
edges.df <- edges.df %>%
  mutate(
    Found = case_when(
      Edge %in% intersect(edges.wv1, edges.wv2) ~ "Both",
      Edge %in% edges.wv1 ~ "Before",
      TRUE ~ "After"
    ),
    Colour = case_when(
      Found == "Both" ~ "#CCCCCC",
      Found == "Before" ~ "#F4CAE4",
      TRUE ~ "#FDCDAC"
    )
  )

edges.colour <- structure(edges.df$Colour, names = edges.df$Edge)



# Graph parameters: Plot network with edge colours ------------------------
metab.col <- "#B3CDE3"
prot.col <- "#FFF2AE"

V(graph.cluster.21.prot.wv1)$label.cex <- rep(1.15, length(V(graph.cluster.21.prot.wv1)))
V(graph.cluster.21.prot.wv2)$label.cex <- rep(1.15, length(V(graph.cluster.21.prot.wv2)))
V(graph.cluster.21.prot.wv1)$color <- c(rep(metab.col, 21), rep(prot.col, 8))
V(graph.cluster.21.prot.wv2)$color <- c(rep(metab.col, 21), rep(prot.col, 8))
V(graph.cluster.21.prot.wv1)$frame.color <- V(graph.cluster.21.prot.wv1)$color
V(graph.cluster.21.prot.wv2)$frame.color <- V(graph.cluster.21.prot.wv2)$color
E(graph.cluster.21.prot.wv1)$color <- edges.colour[edges.wv1]
E(graph.cluster.21.prot.wv2)$color <- edges.colour[edges.wv2]
E(graph.cluster.21.prot.wv1)$width <- 3
E(graph.cluster.21.prot.wv2)$width <- 3


# Plotting ----------------------------------------------------------------
nodes.remov <- which(igraph::degree(graph.cluster.21.prot.wv1) == 0 & igraph::degree(graph.cluster.21.prot.wv2) == 0)
pdf("Results/Network/multiomics_sPLS_selected_centroids_21_prepost_col_edges.pdf",
  width = 13, height = 8
)

# Define multiple plotting grid to allow legends at the bottom
mat <- matrix(c(1, 2, 3, 3), byrow = TRUE, nrow = 2)
layout(mat, heights = c(1, 0.3))

set.seed(474)
# Create a fixed layout to be applied for both networks
layout <- layout_with_fr(induced_subgraph(graph.cluster.21.prot.wv1, V(graph.cluster.21.prot.wv1)[-nodes.remov]))
# Only keep nodes with at least one edge
grp.sub <- induced_subgraph(graph.cluster.21.prot.wv1, V(graph.cluster.21.prot.wv1)[-nodes.remov])
# Plot network
par(mar = c(0, 4, 5, 1))
plot(grp.sub,
  layout = layout,
  vertex.label = ifelse(is.na(cluster.annotation[V(grp.sub)$label]), V(grp.sub)$label, cluster.annotation[V(grp.sub)$label])
)
mtext("A", cex = 2.8, adj = -0.05, line = 0.5)
mtext("Before swimming", cex = 1.5, outer = FALSE, cex.main = 2, line = 0.5)

grp.sub <- induced_subgraph(graph.cluster.21.prot.wv2, V(graph.cluster.21.prot.wv2)[-nodes.remov])
plot(grp.sub,
  layout = layout,
  vertex.label = ifelse(is.na(cluster.annotation[V(grp.sub)$label]), V(grp.sub)$label, cluster.annotation[V(grp.sub)$label])
)
mtext("B", cex = 2.8, adj = -0.05, line = 0.5)
mtext("After swimming", cex = 1.5, outer = FALSE, cex.main = 2, line = 0.5)

par(mar = c(0, 0, 0, 0))
plot(1, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")

legend(
  x = 0.8, y = 1.4, pch = 19,
  col = c(metab.col, prot.col),
  legend = c("Clusters of metabolic features", "Immune marker"),
  horiz = FALSE,
  bty = "n",
  cex = 1.8,
  pt.cex = 2.5,
  inset = 0.3,
  title = "Node",
  xjust = 0.5
)

legend(
  x = 1.2, y = 1.4, lty = 1,
  col = c("#CCCCCC", "#F4CAE4", "#FDCDAC"),
  legend = c("Both before and after", "Only before swimming", "Only after swimming"),
  horiz = FALSE,
  bty = "n",
  cex = 1.8,
  lwd = 3.5,
  inset = 0.5,
  title = "Appearance of edge",
  xjust = 0.5
)

dev.off()
