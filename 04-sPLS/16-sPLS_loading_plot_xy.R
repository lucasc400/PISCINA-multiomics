# =============================================================================
# 16-sPLS_loading_plot_xy.R
#
# Characterise the final sPLS models: loadings and variable importance of the
# exposures, variance explained by the components, loadings of the selected
# metabolic clusters, comparison of sPLS selection with the MWAS results and
# a circular plot of loadings vs. MWAS significance.
#
# Input:  Data/Cleaned_Data/*.rds, Results/Clustering/*, Results/MVN/*,
#         Results/sPLS/*, Results/Annotation/*
# Output: Results/Figures/sPLS_XY_ev_loading_panel.pdf
#         Results/Figures/circular_sPLS_loading_MVN_pval.pdf
#         Results/Figures/* (other figures, if write.fig)
# =============================================================================

rm(list = ls())
setwd(Sys.getenv("PISCINA_DIR", unset = ".")) # project root (see README)

library(dplyr)
library(tidyr)
library(ggplot2)
library(mixOmics)
library(stringr)
library(patchwork)
library(circlize)
library(ComplexHeatmap)

source("Scripts/functions.R")


# Read data ---------------------------------------------------------------
# Read and format data
expo <- readRDS("Data/Cleaned_Data/expo.piscina.rds")
metab.centroid <- readRDS("Results/Clustering/metabolites_summarised_centroid.rds")
covars <- readRDS("Data/Cleaned_Data/covars.piscina.rds")

# Within-subject variation (multilevel decomposition)
designID <- as.matrix(as.numeric(as.factor(covars$ID)))
withinVarX <- withinVariation(expo, design = designID)
metab.centroid.wv <- withinVariation(metab.centroid, design = designID)

# Read sPLS models and results
spls.centroid.2.21 <- readRDS("Results/sPLS/spls.centroid.2.21.rds")
spls.centroid.4.26 <- readRDS("Results/sPLS/spls.centroid.4.26.rds")
metab.centroids.selected.clusters.21 <- readRDS("Results/sPLS/metab.centroids.selected.clusters.21.rds")
metab.centroids.selected.clusters.26 <- readRDS("Results/sPLS/metab.centroids.selected.clusters.26.rds")
metab.features.selected.centroids.21 <- readRDS("Results/sPLS/metab.features.selected.centroids.21.rds")

# Clustering results
cluster.description <- readRDS("Results/Clustering/Description_summary.rds")

# Annotation
annotation <- readRDS("Results/Annotation/annotation.sPLS.selected.metab.rds")
cluster.annotation <- readRDS("Results/Annotation/cluster.annotation.rds")

write.fig <- FALSE

# Compute metrics ---------------------------------------------------------
# Exposure
vip.spls.centroid.2.21 <- vip(spls.centroid.2.21)
vip.spls.centroid.4.26 <- vip(spls.centroid.4.26)

vip.loading.X <- data.frame(
  expo = colnames(withinVarX),
  loading_2 = as.numeric(spls.centroid.2.21$loadings$X),
  loading_4 = as.numeric(spls.centroid.4.26$loadings$X),
  vip_2 = as.numeric(vip.spls.centroid.2.21),
  vip_4 = as.numeric(vip.spls.centroid.4.26)
)
vip.loading.X$expo <- factor(vip.loading.X$expo, levels = c("CHCl3", "BDCM", "DBCM", "CHBr3"))

# Explained variance
ev.centroid.2.21 <- quality.var(spls.centroid.2.21, withinVarX, metab.centroid.wv)
ev.centroid.4.26 <- quality.var(spls.centroid.4.26, withinVarX, metab.centroid.wv)
ev.X.on.X.2.21 <- explained_variance(spls.centroid.2.21$X, spls.centroid.2.21$variates$X, ncomp = 1)
ev.X.on.X.4.26 <- explained_variance(spls.centroid.4.26$X, spls.centroid.4.26$variates$X, ncomp = 1)
ev.Y.on.Y.2.21 <- explained_variance(spls.centroid.2.21$Y, spls.centroid.2.21$variates$Y, ncomp = 1)
ev.Y.on.Y.4.26 <- explained_variance(spls.centroid.4.26$Y, spls.centroid.4.26$variates$Y, ncomp = 1)

# Rows: model 2.21 (sPLS on X and Y) then model 4.26 (sPLS on Y), for each
# of: selected Y explained by C_X, X explained by C_X, Y explained by C_Y
ev.df <- data.frame(
  var.type = rep(c("Y~by~C[X]", "X~by~C[X]", "Y~by~C[Y]"), each = 2),
  model = rep(c("sPLS on X and Y", "sPLS on Y"), 3),
  value = c(
    ev.centroid.2.21$PerExplain, ev.centroid.4.26$PerExplain,
    ev.X.on.X.2.21, ev.X.on.X.4.26,
    ev.Y.on.Y.2.21, ev.Y.on.Y.4.26
  )
)

# Prepare for plotting ----------------------------------------------------
ev.loading <- data.frame(
  cluster = metab.centroids.selected.clusters.26,
  ev_21 = ev.centroid.2.21$PerExplainY[metab.centroids.selected.clusters.26],
  ev_26 = ev.centroid.4.26$PerExplainY[metab.centroids.selected.clusters.26],
  loading_21 = spls.centroid.2.21$loadings$Y[metab.centroids.selected.clusters.26, ],
  loading_26 = spls.centroid.4.26$loadings$Y[metab.centroids.selected.clusters.26, ],
  cluster_annotation = ifelse(is.na(cluster.annotation[metab.centroids.selected.clusters.26]),
    metab.centroids.selected.clusters.26,
    cluster.annotation[metab.centroids.selected.clusters.26]
  )
)

# Compare sPLS with univariate --------------------------------------------
metab.piscina.MVN.res <- readRDS("Results/MVN/metab.piscina.cleaned.res.rds")
PISCINA_cor_centroid <- readRDS("Results/sPLS/PISCINA_sPLS_tune_cor_metab_centroid.rds")

# Selection proportions and mean loadings of each cluster across the 500
# resampled models, for every (nX, nY) combination

MySelPro.centroid <- as.data.frame(do.call(rbind, lapply(PISCINA_cor_centroid, function(x) x$y.selection.proportion)))
MySelPro.centroid <- cbind(as.data.frame(expand.grid(nX = c(1:4), nY = seq(1, 550, 1))), MySelPro.centroid)
colnames(MySelPro.centroid) <- c("nX", "nY", colnames(metab.centroid))
MyLoading.centroid <- as.data.frame(do.call(rbind, lapply(PISCINA_cor_centroid, function(x) x$yloadings.mean)))
MyLoading.centroid <- cbind(as.data.frame(expand.grid(nX = c(1:4), nY = seq(1, 550, 1))), MyLoading.centroid)
colnames(MyLoading.centroid) <- c("nX", "nY", colnames(metab.centroid))

# nX=2, nY=21
spls.mvn.res.21 <- cluster.description %>%
  dplyr::select(feature, group) %>%
  rename(metab = feature) %>%
  mutate(group = paste0("C", group)) %>%
  left_join(., MySelPro.centroid %>% filter(nX == 2, nY == 21) %>% pivot_longer(., everything(), names_to = "group", values_to = "SelectionProportion")) %>%
  left_join(., MyLoading.centroid %>% filter(nX == 2, nY == 21) %>% pivot_longer(., everything(), names_to = "group", values_to = "Loading")) %>%
  left_join(., metab.piscina.MVN.res %>%
    filter(Expo %in% c("Chloroform", "BDCM", "DBCM", "Bromoform")) %>%
    pivot_wider(id_cols = OMIC_ID, values_from = pval, names_from = Expo),
  by = c("metab" = "OMIC_ID")
  ) %>%
  mutate(spls_selected = ifelse(group %in% metab.centroids.selected.clusters.21, "sPLS selected", "Not selected"))

# nX=4, nY=26
spls.mvn.res.26 <- cluster.description %>%
  dplyr::select(feature, group) %>%
  rename(metab = feature) %>%
  mutate(group = paste0("C", group)) %>%
  left_join(., MySelPro.centroid %>% filter(nX == 4, nY == 26) %>% pivot_longer(., everything(), names_to = "group", values_to = "SelectionProportion")) %>%
  left_join(., MyLoading.centroid %>% filter(nX == 4, nY == 26) %>% pivot_longer(., everything(), names_to = "group", values_to = "Loading")) %>%
  left_join(., metab.piscina.MVN.res %>%
    filter(Expo %in% c("Chloroform", "BDCM", "DBCM", "Bromoform")) %>%
    pivot_wider(id_cols = OMIC_ID, values_from = pval, names_from = Expo),
  by = c("metab" = "OMIC_ID")
  ) %>%
  mutate(spls_selected = ifelse(group %in% metab.centroids.selected.clusters.26, "Selected", "Not selected"))

# Create a combined df ----------------------------------------------------
sel.pro.df <- spls.mvn.res.21 %>%
  rename(sel.pro_21 = SelectionProportion) %>%
  left_join(
    spls.mvn.res.26 %>%
      rename(sel.pro_26 = SelectionProportion) %>%
      dplyr::select(metab, sel.pro_26),
    by = "metab"
  ) %>%
  dplyr::select(group, sel.pro_21, sel.pro_26) %>%
  rename(cluster = group) %>%
  filter(!duplicated(cluster))

# Arrange clusters
cluster.annot.levels <- sel.pro.df %>%
  filter(cluster %in% metab.centroids.selected.clusters.26) %>%
  arrange(desc(sel.pro_26)) %>%
  mutate(cluster_annotation = as.character(ifelse(is.na(cluster.annotation[cluster]),
    cluster, cluster.annotation[cluster]
  ))) %>%
  mutate(cluster_annotation = ifelse(cluster_annotation != cluster, paste0(cluster, " (", cluster_annotation, ")"), cluster_annotation)) %>%
  pull(cluster_annotation)

ev.loading.selpro.long <- ev.loading %>%
  left_join(., sel.pro.df) %>%
  pivot_longer(c("ev_21", "ev_26", "loading_21", "loading_26", "sel.pro_21", "sel.pro_26"),
    names_sep = "_", names_to = c("metric", "nY")
  ) %>%
  mutate(model = ifelse(nY == 21, "sPLS on X and Y", "sPLS on Y")) %>%
  mutate(cluster_annotation = ifelse(cluster_annotation != cluster, paste0(cluster, " (", cluster_annotation, ")"), cluster_annotation)) %>%
  mutate(cluster_annotation = factor(cluster_annotation, levels = cluster.annot.levels))


# Explained variance, feature importance ----------------------------------
# Prepare data
ev.loading.selpro.wide <- ev.loading.selpro.long %>%
  dplyr::select(-nY) %>%
  pivot_wider(names_from = c("metric", "model"), values_from = "value", id_cols = c("cluster", "cluster_annotation")) %>%
  mutate(cluster_annotation = gsub("\\bF", "C", cluster_annotation)) %>%
  as.data.frame()

vip.X <- t(as.matrix(vip.loading.X[c("vip_2", "vip_4")]))
colnames(vip.X) <- vip.loading.X$expo
loading.X <- t(as.matrix(vip.loading.X[c("loading_2", "loading_4")]))
colnames(loading.X) <- vip.loading.X$expo


ev.feature.matrix <- as.matrix(t(ev.loading.selpro.wide[, c("ev_sPLS on X and Y", "ev_sPLS on Y")]))
loading.Y.matrix <- as.matrix(t(ev.loading.selpro.wide[, c("loading_sPLS on X and Y", "loading_sPLS on Y")]))
colnames(ev.feature.matrix) <- colnames(loading.Y.matrix) <- ev.loading.selpro.wide$cluster_annotation

# Rank clusters by absolute loading in the sPLS on X and Y model
metab.order <- rev(order(abs(loading.Y.matrix[1, ])))
ev.feature.matrix <- ev.feature.matrix[, metab.order]
loading.Y.matrix <- loading.Y.matrix[, metab.order]


# Panel figure for the sPLS on X and Y model ------------------------------
pdf("Results/Figures/sPLS_XY_ev_loading_panel.pdf", width = 9, height = 10)

dbp.col1 <- "#D2F9F0"
dbp.col2 <- "#FEE1CD"
col1 <- "#B3CDE3"

mat <- matrix(c(1, 2, 3, 3, 4, 4, 5, 5, 6, 6),
  byrow = TRUE,
  nrow = 5
)
layout(mat, heights = c(1.1, 1.5, 1.8, 0.9, 0.3))

# Loading, VIP
par(mar = c(1, 10, 2, 2))
barplot(loading.X[1, ],
  beside = TRUE, ylab = "Loading",
  col = c(dbp.col1, dbp.col2, dbp.col1, dbp.col2),
  names.arg = c("Chloroform", "BDCM", "DBCM", "Bromoform"),
  xaxs = "i", yaxt = "n",
  cex.axis = 1.2,
  cex.lab = 1.2,
  cex.names = 1.2
)
axis(2, at = seq(-0.6, 0, by = 0.2), labels = seq(-0.6, 0, by = 0.2), cex.axis = 1.2)

barplot(vip.X[1, ],
  beside = TRUE, ylab = "Variable importance",
  col = c(dbp.col1, dbp.col2, dbp.col1, dbp.col2),
  names.arg = c("Chloroform", "BDCM", "DBCM", "Bromoform"),
  xaxs = "i",
  cex.axis = 1.2,
  cex.lab = 1.2,
  cex.names = 1.2
)
abline(h = 1, lty = 2, col = "grey50")

# EV, loading
par(mar = c(2, 10, 2, 2))
plot(seq(1, 42, by = 2) + 0.5, ev.feature.matrix[1, 1:21] * 100,
  pch = 19, col = col1,
  ylab = expression("Variance explained by" ~ C[X] ~ "(%)"),
  xlab = "",
  xaxt = "n",
  bty = "n",
  xaxs = "i",
  xlim = c(0, 42),
  ylim = c(35, 60),
  cex = 1.5,
  cex.axis = 1.2,
  cex.lab = 1.2
)
abline(v = seq(1.5, 42.5, by = 2), lty = 2, col = "grey80")

par(mar = c(14, 10, 0, 2))
barplot(loading.Y.matrix[1, 1:21],
  ylab = "Loading",
  width = 1,
  space = 1,
  col = col1,
  xaxt = "n",
  xaxs = "i",
  xlim = c(0, 42),
  cex.axis = 1.2,
  cex.lab = 1.2
)
abline(v = seq(1.5, 42.5, by = 2), lty = 2, col = "grey80")
text(
  x = seq(1.5, 42.5, by = 2),
  y = -0.5,
  labels = colnames(ev.feature.matrix)[1:21],
  srt = 45,
  adj = 1,
  xpd = TRUE,
  cex = 1.2
)

# EV by components
par(mar = c(2, 10, 3, 2))
plot(1,
  type = "n",
  yaxt = "n", xaxt = "n",
  ylab = "", xlab = "",
  xlim = c(0, 100),
  ylim = c(0.5, 8.5),
  xaxs = "i"
)
abline(h = c(3, 6), lty = 1, col = "grey50")
abline(h = c(1.5, 4.5, 7.5), lty = 2, col = "grey90")
points(ev.df$value[c(2, 4, 6)] * 100, c(1.5, 4.5, 7.5),
  col = col1,
  cex = 1.5,
  pch = 19
)
axis(
  side = 2, at = c(1.5, 4.5, 7.5),
  labels = c(expression(Selected ~ Y ~ by ~ C[X]), expression(X ~ by ~ C[X]), expression(Y ~ by ~ C[Y])),
  las = 1,
  tick = FALSE,
  cex.axis = 1.2
)
axis(side = 3, at = seq(0, 100, 10), cex.axis = 1.2)
mtext("Explained variance (%)", side = 3, line = 2.5, cex = 0.9)

# Legend
par(mar = c(0, 10, 1, 2))
plot(1, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")

dev.off()


if (write.fig) {
  pdf("Results/Figures/sPLS_ev_loading_panel.pdf", width = 9, height = 10)
  {
    mat <- matrix(c(1, 2, 3, 3, 4, 4, 5, 5, 6, 6),
      byrow = TRUE,
      nrow = 5
    )
    layout(mat, heights = c(1.1, 1.5, 1.8, 0.9, 0.3))

    # Loading, VIP
    par(mar = c(1, 10, 2, 2))
    barplot(loading.X,
      beside = TRUE, ylab = "Loading",
      col = rep(c(col1, col2), 3),
      names.arg = c("Chloroform", "BDCM", "DBCM", "Bromoform"),
      xaxs = "i", xlim = c(0, 12), yaxt = "n",
      cex.axis = 1.2,
      cex.lab = 1.2,
      cex.names = 1.2
    )
    axis(2, at = seq(-0.6, 0, by = 0.2), labels = seq(-0.6, 0, by = 0.2), cex.axis = 1.2)

    barplot(vip.X,
      beside = TRUE, ylab = "Variable importance",
      col = rep(c(col1, col2), 3),
      names.arg = c("Chloroform", "BDCM", "DBCM", "Bromoform"),
      xaxs = "i", xlim = c(0, 12),
      cex.axis = 1.2,
      cex.lab = 1.2,
      cex.names = 1.2
    )
    abline(h = 1, lty = 2, col = "grey50")

    # EV, loading
    par(mar = c(2, 10, 2, 2))
    plot(seq(1, 76, by = 3) + 0.5, ev.feature.matrix[1, ] * 100,
      pch = 19, col = col1,
      ylab = expression("Variance explained by" ~ C[X] ~ "(%)"),
      xlab = "",
      xaxt = "n",
      bty = "n",
      xaxs = "i",
      xlim = c(0, 78),
      cex = 1.2,
      cex.axis = 1.2,
      cex.lab = 1.2
    )
    points(seq(2, 77, by = 3) + 0.5, ev.feature.matrix[2, ] * 100,
      cex = 1.2,
      pch = 19, col = col2
    )
    abline(v = seq(3.5, 78.5, by = 3), lty = 2, col = "grey80")

    par(mar = c(14, 10, 0, 2))
    barplot(loading.Y.matrix,
      beside = TRUE, ylab = "Loading",
      col = rep(c(col1, col2), 3),
      xaxt = "n",
      xaxs = "i",
      xlim = c(0, 78),
      cex.axis = 1.2,
      cex.lab = 1.2
    )
    abline(v = seq(3.5, 78.5, by = 3), lty = 2, col = "grey80")
    text(
      x = seq(2, 77, by = 3),
      y = -0.5,
      labels = colnames(ev.feature.matrix),
      srt = 45,
      adj = 1,
      xpd = TRUE,
      cex = 1.2
    )

    # EV by components
    par(mar = c(2, 10, 3, 2))
    plot(1,
      type = "n",
      yaxt = "n", xaxt = "n",
      ylab = "", xlab = "",
      xlim = c(0, 100),
      ylim = c(0.5, 8.5),
      xaxs = "i"
    )
    abline(h = c(3, 6), lty = 1, col = "grey50")
    abline(h = c(1, 2, 4, 5, 7, 8), lty = 2, col = "grey90")
    points(ev.df$value * 100, c(1, 2, 4, 5, 7, 8),
      col = rep(c(col1, col2), 3),
      cex = 1.5,
      pch = 19
    )
    axis(
      side = 2, at = c(1.5, 4.5, 7.5),
      labels = c(expression(Selected ~ Y ~ by ~ C[X]), expression(X ~ by ~ C[X]), expression(Y ~ by ~ C[Y])),
      las = 1,
      tick = FALSE,
      cex.axis = 1.2
    )
    axis(side = 3, at = seq(0, 100, 10), cex.axis = 1.2)
    mtext("Explained variance (%)", side = 3, line = 2.5, cex = 0.9)

    # Legend
    par(mar = c(0, 10, 1, 2))
    plot(1, type = "n", bty = "n", xaxt = "n", yaxt = "n", xlab = "", ylab = "")

    legend("top",
      fill = c(col2, col1),
      legend = c("sPLS on Y", "sPLS on X and Y"),
      horiz = TRUE,
      inset = -0.1,
      bty = "n",
      cex = 1.3
    )
  }
  dev.off()
}


# Circlize ----------------------------------------------------------------
spls.mvn.res.selected <- spls.mvn.res.21 %>%
  filter(spls_selected == "sPLS selected") %>%
  arrange(desc(SelectionProportion))

# Sort out annotation
spls.mvn.res.selected.annotation <- left_join(spls.mvn.res.selected, annotation[, c("feature", "Annotation")],
  by = c("metab" = "feature")
)
spls.mvn.res.selected.annotation$Annotation <- ifelse(grepl("Unknown", spls.mvn.res.selected.annotation$Annotation),
  spls.mvn.res.selected.annotation$metab,
  spls.mvn.res.selected.annotation$Annotation
)


# Clusters
split.cluster <- factor(spls.mvn.res.selected$group, levels = unique(spls.mvn.res.selected$group))

# Loading
loading <- as.matrix(spls.mvn.res.selected$Loading)

loading.named <- loading
rownames(loading.named) <- spls.mvn.res.selected.annotation$Annotation

loading.cat <- ifelse(loading.named > 0.2, ">0.2",
  ifelse(loading.named > 0, ">0",
    ifelse(loading.named < -0.2, "<-0.2", "<0")
  )
)

# -log10 transformed p-values
spls.mvn.res.selected.pval <- apply(spls.mvn.res.selected[, c("Chloroform", "BDCM", "DBCM", "Bromoform")], 2, function(x) -log10(x))
rownames(spls.mvn.res.selected.pval) <- spls.mvn.res.selected.annotation$Annotation
colnames(spls.mvn.res.selected.pval) <- c("Chloroform", "BDCM", "DBCM", "Bromoform")

# Categorical
# Bonferroni significance
spls.mvn.res.selected.sign <- apply(spls.mvn.res.selected[, c("Chloroform", "BDCM", "DBCM", "Bromoform")], 2, function(x) ifelse(x < 0.05 / 1696, "Sign", "Not sign"))
rownames(spls.mvn.res.selected.sign) <- spls.mvn.res.selected.annotation$Annotation
colnames(spls.mvn.res.selected.sign) <- c("Chloroform", "BDCM", "DBCM", "Bromoform")

# Directionalities of association
spls.mvn.res.selected.dir <- metab.piscina.MVN.res %>%
  filter(
    OMIC_ID %in% metab.features.selected.centroids.21,
    Expo %in% c("Chloroform", "BDCM", "DBCM", "Bromoform")
  ) %>%
  mutate(sign = ifelse(coef > 0, 1, -1)) %>%
  pivot_wider(id_cols = OMIC_ID, names_from = Expo, values_from = sign) %>%
  as.data.frame()

rownames(spls.mvn.res.selected.dir) <- spls.mvn.res.selected.dir$OMIC_ID
spls.mvn.res.selected.dir <- as.matrix(spls.mvn.res.selected.dir[spls.mvn.res.selected$metab, c("Chloroform", "BDCM", "DBCM", "Bromoform")])

# Categorise p-values
spls.mvn.res.selected.pval.category <- apply(
  spls.mvn.res.selected[, c("Chloroform", "BDCM", "DBCM", "Bromoform")], 2,
  function(x) {
    ifelse(x < 1e-10, 2,
      ifelse(x < 0.05 / 1696, 1, 0)
    )
  }
)
spls.mvn.res.selected.cat.dir <- spls.mvn.res.selected.pval.category * spls.mvn.res.selected.dir
spls.mvn.res.selected.cat.dir <- as.matrix(spls.mvn.res.selected.cat.dir)
rownames(spls.mvn.res.selected.cat.dir) <- spls.mvn.res.selected.annotation$Annotation


{
  pdf("Results/Figures/circular_sPLS_loading_MVN_pval.pdf", width = 7, height = 9)
  plot.new()
  circle_size <- unit(2, "snpc")
  # Grid parameters
  pushViewport(viewport(
    x = 0, y = 0, width = circle_size, height = circle_size,
    just = c("left", "center")
  ))
  par(omi = c(1.8, 0, 0, 0), new = TRUE)

  # Plotting
  circos.clear()

  col_fun1 <- structure(c("#2E45B8", "#D6DBF5", "#D2F9F0", "#1DC9A4"),
    names = c("<-0.2", "<0", ">0", ">0.2")
  )
  circos.par(gap.after = c(rep(2, 20), 20))

  # Heatmap: Loading
  circos.heatmap(loading.cat,
    split = split.cluster, col = col_fun1, cluster = FALSE,
    track.height = 0.05, bg.border = TRUE, bg.lwd = 0.25
  )
  # Add column name
  circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
    if (CELL_META$sector.numeric.index == 21) { # the last sector
      cn <- "Loading"
      n <- length(cn)
      circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"),
        1:n - 0.5, cn,
        cex = 0.7, adj = c(0, 0.5), facing = "inside"
      )
    }
  }, bg.border = NA)

  # Add cluster name
  circos.track(track.index = get.current.track.index(), panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$cell.ylim[2] + convert_y(3, "mm"),
      CELL_META$sector.index,
      facing = "bending.inside", cex = 0.5,
      adj = c(0.5, 0), niceFacing = TRUE
    )
  }, bg.border = NA)

  # Heatmap: significance
  col_fun3 <- colorRamp2(c(-2, -1, 0, 1, 2), c("#1DC9A4", "#D2F9F0", "grey90", "#D6DBF5", "#2E45B8"))
  circos.heatmap(spls.mvn.res.selected.cat.dir,
    split = split.cluster, col_fun3,
    cluster = FALSE, cell.border = TRUE, cell.lwd = 0.25, rownames.side = "inside"
  )

  # Add column name
  circos.track(track.index = get.current.track.index() - 1, panel.fun = function(x, y) {
    if (CELL_META$sector.numeric.index == 21) { # the last sector
      cn <- c("Chloroform", "BDCM", "DBCM", "Bromoform")
      n <- length(cn)
      circos.text(rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"),
        1:n - 0.5, cn,
        cex = 0.7, adj = c(0, 0.5), facing = "inside"
      )
    }
  }, bg.border = NA)

  # upViewport
  upViewport()

  # Add legend
  lgd.loading <- Legend(
    title = "Loading in sPLS on X and Y model",
    labels = c("< -0.2", "< 0", "> 0", "> 0.2"),
    legend_gp = gpar(fill = c("#2E45B8", "#D6DBF5", "#D2F9F0", "#1DC9A4")),
    ncol = 4
  )
  lgd.sign <- Legend(
    title = "Significance in MVN model (Direction of association)",
    labels = c(
      expression("p-value <" ~ 10^-10 ~ "(Positive)"),
      "Bonferroni significant (Positive)",
      "Not significant",
      "Bonferroni significant (Negative)",
      expression("p-value <" ~ 10^-10 ~ "(Negative)")
    ),
    legend_gp = gpar(fill = c("#2E45B8", "#D6DBF5", "grey90", "#D2F9F0", "#1DC9A4")),
    ncol = 3, nrow = 2, by_row = TRUE
  )

  lgd.list <- packLegend(lgd.loading, lgd.sign, direction = "vertical")
  draw(lgd.list, y = 0.13 * circle_size, just = "top")
  dev.off()
}


# Compare sPLS selection proportion with MWAS p-value ---------------------
sPLS.2.21.pval.vs.sel.pro <- spls.mvn.res.21 %>%
  mutate(spls_selected = ifelse(spls_selected == "sPLS selected", "sPLS selected (Model 2, N=21 clusters of metabolic features)", spls_selected)) %>%
  pivot_longer(Chloroform:Bromoform, names_to = "Expo", values_to = "p-value") %>%
  mutate(SelectionProportion = 100 * SelectionProportion) %>%
  ggplot(aes(-log10(`p-value`), SelectionProportion)) +
  geom_line(aes(group = group), linewidth = 0.15, col = "grey70") +
  geom_point(aes(col = spls_selected), size = 0.5) +
  facet_wrap(~Expo, labeller = labeller(Expo = label_parsed)) +
  theme_bw() +
  scale_colour_manual(values = c("black", "#CC3333"), name = "") +
  theme(
    legend.position = "bottom",
    legend.text = element_text(size = 12),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_y_continuous(breaks = seq(0, 100, 20)) +
  ylab("Selection Proportion (% of 500 models)") +
  xlab(expression(paste("-log"[10] * "(", italic("p"), "-value) of individual feature in MVN model"))) +
  geom_vline(aes(xintercept = -log10(0.05 / 6040)), linetype = "dashed", linewidth = 0.25)

if (write.fig) {
  ggsave("Results/Figures/sPLS_2_21_selection_proportion_vs_pvalue.pdf", sPLS.2.21.pval.vs.sel.pro)
}
