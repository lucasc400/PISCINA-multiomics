# =============================================================================
# 11-Clustering_results.R
#
# Summarise each stable cluster of metabolic features by its first principal
# component (centroid), describe the clusters (explained variance, ranges of
# retention time and mass) and compare feature correlations before and after
# clustering.
#
# Input:  Data/Cleaned_Data/metab.piscina.denoised.rds
#         Results/Clustering/consensus_clustering_broad*.rds
# Output: Results/Clustering/metabolites_summarised_centroid.rds,
#         Results/Clustering/Description_summary.rds (if write.res)
#         Results/Figures/* (if write.fig)
# =============================================================================

rm(list = ls())
setwd(Sys.getenv("PISCINA_DIR", unset = ".")) # project root (see README)

library(dplyr)
library(stringr)
library(ggplot2)
library(sharp)
library(pheatmap)
library(patchwork)
library(RColorBrewer)
library(gridExtra)

source("Scripts/functions.R")

metab <- readRDS("Data/Cleaned_Data/metab.piscina.denoised.rds")
clust_stab <- readRDS("Results/Clustering/consensus_clustering_broad.rds")
calib <- readRDS("Results/Clustering/consensus_clustering_broad_calibration.rds")

write.res <- FALSE
write.fig <- FALSE

# Summarise clusters ------------------------------------------------------
# Extracting stable clusters
myclusters <- Clusters(clust_stab)

# Summarising the data by cluster (centroids)
summarised <- GroupingSummary(
  x = metab,
  group = myclusters,
  summary = "centroid"
)
metab_summarised_centroids <- summarised$xsummarised
description_summary <- summarised$description

# Saving the outputs
if (write.res) {
  write.csv(description_summary, "Results/Tables/Description_summary.csv")
  saveRDS(description_summary, "Results/Clustering/Description_summary.rds")
  saveRDS(metab_summarised_centroids, "Results/Clustering/metabolites_summarised_centroid.rds")
}

# Plot calibration --------------------------------------------------------
# Stability score by number of clusters (candidates up to nc = 1,600); the
# dashed line marks the selected number of clusters (nc = 550)
if (write.fig) {
  pdf("Results/Figures/Clustering_calibration_plot_rm_last.pdf", width = 10, height = 6)
  par(mar = c(5, 6, 4, 4))
  plot(as.numeric(names(calib)[1:281]), as.numeric(calib[1:281]),
    xlab = "Number of clusters", ylab = "Stability score", pch = 19,
    col = "#0019a8", cex = 0.3, cex.axis = 1, cex.lab = 1
  )
  abline(v = 550, col = "indianred1", lty = 2)
  dev.off()
}

# Ranges by cluster -------------------------------------------------------
# Extracting the multivariate clusters (>1 members)
multiv_clusters <- sort(unique(myclusters[duplicated(myclusters)]))
mymultiv_clusters <- myclusters[myclusters %in% multiv_clusters]

# Extracting mass and retention time
mass <- as.numeric(str_split_fixed(names(mymultiv_clusters), "@", 2)[, 1])
time <- as.numeric(str_split_fixed(str_split_fixed(names(mymultiv_clusters), ":", 2)[, 1], "@", 2)[, 2])

# Save in a dataframe
cluster.stats <- data.frame(
  group = multiv_clusters,
  mass.range = rep(NA, length(multiv_clusters)),
  mass.IQR = rep(NA, length(multiv_clusters)),
  rt.range = rep(NA, length(multiv_clusters)),
  rt.IQR = rep(NA, length(multiv_clusters))
)

mass.list <- split(mass, f = mymultiv_clusters)
time.list <- split(time, f = mymultiv_clusters)

cluster.stats$mass.range <- sapply(mass.list, function(x) {
  diff(range(x))
})
cluster.stats$mass.IQR <- sapply(mass.list, IQR)
cluster.stats$rt.range <- sapply(time.list, function(x) {
  diff(range(x))
})
cluster.stats$rt.IQR <- sapply(time.list, IQR)


# Get description summary for clusters with multiple metabolites
description_summary.multi <- description_summary[, c("group", "min_cor", "max_cor", "mean_cor", "ev")]
description_summary.multi <- description_summary.multi[!duplicated(description_summary.multi), ]

cluster.stats <- merge(cluster.stats, description_summary.multi, all.x = TRUE, by = "group")
cluster.stats$ev.pct <- cluster.stats$ev * 100

# Panel plot
p.ev <- cluster.stats %>%
  mutate(range = cut(ev.pct, breaks = seq(40, 100, by = 10), include.lowest = TRUE)) %>%
  ggplot(aes(ev.pct, fill = range)) +
  geom_histogram(breaks = seq(40, 100, by = 2.5), col = "white") +
  theme_bw() +
  scale_fill_manual(values = rev(RColorBrewer::brewer.pal(9, "YlGnBu")[3:8])) +
  xlab("Explained Variance (%) by PC1") +
  ylab("Count") +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

p.rt.range <- cluster.stats %>%
  ggplot(aes(rt.range)) +
  geom_histogram(fill = "grey75", col = "white", bins = 20) +
  theme_bw() +
  xlab("Range of retention time (min)") +
  ylab("Count") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

p.mass.range <- cluster.stats %>%
  ggplot(aes(mass.range)) +
  geom_histogram(fill = "grey75", col = "white", bins = 20) +
  theme_bw() +
  xlab("Range of monoisotopic mass") +
  ylab("Count") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

p.scatter <- cluster.stats %>%
  ggplot(aes(rt.range, mass.range, col = cut(ev.pct, breaks = seq(40, 100, 10)))) +
  geom_point() +
  theme_bw() +
  scale_color_manual(values = rev(RColorBrewer::brewer.pal(9, "YlGnBu")[3:8]), name = "Explained\nVariance (%)\nby PC1") +
  xlab("Range of retention time (min)") +
  ylab("Range of monoisotopic mass") +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

panel <- (p.rt.range + p.mass.range) / (p.ev + p.scatter) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 20))

if (write.fig) {
  ggsave("Results/Figures/PISCINA_cluster_statistics.pdf", panel, width = 10, height = 7)
}


# Heatmaps for before and after clustering --------------------------------
cor.metab <- pheatmap(cor(metab),
  breaks = seq(-1, 1, length.out = 100),
  main = "Before clustering (N=1,696)", fontsize = 6,
  show_rownames = FALSE, show_colnames = FALSE,
  treeheight_row = FALSE, treeheight_col = FALSE
)

cor.metab.summarised <- pheatmap(cor(metab_summarised_centroids),
  breaks = seq(-1, 1, length.out = 100),
  main = "After clustering (N=550)", fontsize = 6,
  show_rownames = FALSE, show_colnames = FALSE,
  treeheight_row = FALSE, treeheight_col = FALSE
)

plot_list_metab <- list()
plot_list_metab[["p1"]] <- cor.metab[[4]]
plot_list_metab[["p2"]] <- cor.metab.summarised[[4]]
g <- arrangeGrob(grobs = plot_list_metab, ncol = 2, nrow = 1)
if (write.fig) {
  ggsave(filename = "Results/Figures/cor_heatmap_metab_cluster.pdf", g, width = 7, height = 4)
}
