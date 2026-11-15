# =============================================================================
# 18-Stability_Graphs.R
#
# Multi-block conditional independence networks (stability-calibrated
# graphical LASSO, sharp::GraphicalModel) integrating the 21 sPLS-selected
# metabolic cluster centroids and 8 immune markers, estimated separately
# before and after swimming.
#
# Input:  Data/Cleaned_Data/*.rds
#         Results/Clustering/metabolites_summarised_centroid.rds
#         Results/sPLS/metab.centroids.selected.clusters.21.rds
# Output: Results/Network/stab.cluster.21.prot.wv{1,2}.rds
# =============================================================================

rm(list = ls())
setwd(Sys.getenv("PISCINA_DIR", unset = ".")) # project root (see README)

library(dplyr)
library(sharp)

metab.centroid <- readRDS("Results/Clustering/metabolites_summarised_centroid.rds")
covars <- readRDS("Data/Cleaned_Data/covars.piscina.rds")
prot <- readRDS("Data/Cleaned_Data/protein.piscina.rds")

metab.centroids.selected.clusters.21 <- readRDS("Results/sPLS/metab.centroids.selected.clusters.21.rds")

dir.create("Results/Network", showWarnings = FALSE, recursive = TRUE)

# Prepare data ------------------------------------------------------------
# Immune markers previously associated with THM exposure
prot.selected <- prot %>%
  dplyr::select(`IL-8`, VEGF, CCL22, CCL11, CRP, CXCL10, Periostin, IL1ra)

metab.wv1.clusters.21 <- metab.centroid[covars$wave == 1, metab.centroids.selected.clusters.21]
metab.wv2.clusters.21 <- metab.centroid[covars$wave == 2, metab.centroids.selected.clusters.21]

# Samples with both metabolomics and proteomics data
metab.prot.wv1.id <- intersect(rownames(prot), rownames(metab.wv1.clusters.21))
metab.prot.wv2.id <- intersect(rownames(prot), rownames(metab.wv2.clusters.21))

metab.wv1.clusters.21.112 <- metab.wv1.clusters.21[metab.prot.wv1.id, ]
metab.wv2.clusters.21.112 <- metab.wv2.clusters.21[metab.prot.wv2.id, ]

prot.wv1.112 <- prot.selected[metab.prot.wv1.id, ]
prot.wv2.112 <- prot.selected[metab.prot.wv2.id, ]

metab.cluster.21.prot.wv1 <- cbind(metab.wv1.clusters.21.112, prot.wv1.112)
metab.cluster.21.prot.wv2 <- cbind(metab.wv2.clusters.21.112, prot.wv2.112)


# Networks: 21 metabolic clusters + 8 immune markers ----------------------
# Same penalty grid for the three blocks (metab-metab, metab-prot, prot-prot)
grid.customised <- matrix(rep(LambdaSequence(lmin = 0.08, lmax = 0.9, cardinal = 100), 3), ncol = 3)

stab.cluster.21.prot.wv1 <- GraphicalModel(metab.cluster.21.prot.wv1,
  pk = c(21, 8), K = 100,
  n_cat = 3, pi_list = seq(0.5, 1, by = 0.01),
  PFER_thr = 50, Lambda = grid.customised
)
stab.cluster.21.prot.wv2 <- GraphicalModel(metab.cluster.21.prot.wv2,
  pk = c(21, 8), K = 100,
  n_cat = 3, pi_list = seq(0.5, 1, by = 0.01),
  PFER_thr = 50, Lambda = grid.customised
)

saveRDS(stab.cluster.21.prot.wv1, "Results/Network/stab.cluster.21.prot.wv1.rds")
saveRDS(stab.cluster.21.prot.wv2, "Results/Network/stab.cluster.21.prot.wv2.rds")
