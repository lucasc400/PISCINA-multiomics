# =============================================================================
# 14-sPLS_Results.R
#
# Results of the sPLS calibration: select the best models, fit the two final
# multilevel sPLS models (sPLS on X and Y: nX = 2, nY = 21; sPLS on Y:
# nX = 4, nY = 26), extract the selected clusters/features, evaluate their
# performance (Q2, MSEP, out-of-sample correlation) and plot the calibration.
#
# Input:  Data/Cleaned_Data/*.rds, Results/Clustering/*, Results/MVN/*,
#         Results/sPLS/PISCINA_sPLS_*.rds (from scripts 12 and 13)
# Output: Results/sPLS/*.rds (if write.res / run.perf)
#         Results/Figures/PISCINA_sPLS_calibration_*.pdf (if write.fig)
# =============================================================================

rm(list = ls())
setwd(Sys.getenv("PISCINA_DIR", unset = ".")) # project root (see README)

library(dplyr)
library(tidyr)
library(ggplot2)
library(mixOmics)
library(stringr)
library(patchwork)
library(sharp)

source("Scripts/functions.R")

expo <- readRDS("Data/Cleaned_Data/expo.piscina.rds")
metab <- readRDS("Data/Cleaned_Data/metab.piscina.denoised.rds")
covars <- readRDS("Data/Cleaned_Data/covars.piscina.rds")

# Within-subject variation (multilevel decomposition)
designID <- as.matrix(as.numeric(as.factor(covars$ID)))
withinVarX <- withinVariation(expo, design = designID)

description_summary <- readRDS("Results/Clustering/Description_summary.rds") # Clustering results
metab.sign.list <- readRDS("Results/MVN/sign.metab.piscina.188.rds") # Univariate results

# Reshape the output of tune_spls_cor() into a long data frame
clean_cor <- function(dat) {
  data.frame(
    nX = do.call(rbind, lapply(dat, function(x) x$nX)),
    nY = do.call(rbind, lapply(dat, function(x) x$nY)),
    cor.tu_all_mean = do.call(rbind, lapply(dat, function(x) x$cor.tu)),
    cor.tu_train_mean = do.call(rbind, lapply(dat, function(x) x$cor.tu.train.mean)),
    cor.tu_test_mean = do.call(rbind, lapply(dat, function(x) x$cor.tu.test.mean)),
    cor.tu_train_moe = do.call(rbind, lapply(dat, function(x) x$cor.tu.train.moe)),
    cor.tu_test_moe = do.call(rbind, lapply(dat, function(x) x$cor.tu.test.moe))
  ) %>%
    pivot_longer(
      cols = cor.tu_all_mean:cor.tu_test_moe,
      names_to = c(".value", "set", "stat"),
      names_sep = "_",
      values_drop_na = TRUE
    ) %>%
    pivot_wider(names_from = stat, values_from = cor.tu) %>%
    mutate(
      set = factor(
        case_when(
          set == "all" ~ "All data",
          set == "train" ~ "Training",
          TRUE ~ "Test"
        ),
        levels = c("All data", "Training", "Test")
      ),
      nX = as.factor(nX)
    )
}

write.res <- FALSE
write.tab <- FALSE
write.fig <- FALSE
run.perf <- FALSE # Re-run cross-validation on the best models

# Clustering + sPLS -------------------------------------------------------
# Read data
metab.centroid <- readRDS("Results/Clustering/metabolites_summarised_centroid.rds")
PISCINA_cor_centroid <- readRDS("Results/sPLS/PISCINA_sPLS_tune_cor_metab_centroid.rds")
PISCINA_perf_centroid <- readRDS("Results/sPLS/PISCINA_sPLS_centroid_perf.rds")
ClusterStab <- readRDS("Results/Clustering/consensus_clustering_broad.rds")

summarised <- GroupingSummary(x = metab, group = Clusters(ClusterStab))
cluster.description <- summarised$description

# WithinVariation
metab.centroid.wv <- withinVariation(metab.centroid, design = designID)

# Results
PISCINA_cor_centroid_cleaned <- clean_cor(PISCINA_cor_centroid)

# Best nY for each nX and data split
sPLS.cor.res <- PISCINA_cor_centroid_cleaned %>%
  group_by(nX, set) %>%
  filter(mean == max(mean)) %>%
  arrange(set)

# Final models and selected clusters -------------------------------------
spls.centroid.2.21 <- spls(withinVarX, metab.centroid.wv, keepX = 2, keepY = 21, ncomp = 1, multilevel = designID, mode = "regression")
spls.centroid.4.26 <- spls(withinVarX, metab.centroid.wv, keepX = 4, keepY = 26, ncomp = 1, multilevel = designID, mode = "regression")

# Selected centroids
metab.centroids.selected.clusters.21 <- rownames(spls.centroid.2.21$loadings$Y)[spls.centroid.2.21$loadings$Y != 0]
metab.centroids.selected.clusters.26 <- rownames(spls.centroid.4.26$loadings$Y)[spls.centroid.4.26$loadings$Y != 0]

metab.features.selected.centroids.21 <- cluster.description[cluster.description$group %in% as.numeric(gsub("C", "", metab.centroids.selected.clusters.21)), "feature"]
metab.features.selected.centroids.26 <- cluster.description[cluster.description$group %in% as.numeric(gsub("C", "", metab.centroids.selected.clusters.26)), "feature"]


description_summary.21 <- description_summary %>%
  filter(group %in% as.numeric(gsub("C", "", metab.centroids.selected.clusters.21))) %>%
  mutate(
    mass = as.numeric(str_split_fixed(feature, "@", 2)[, 1]),
    rt = as.numeric(str_split_fixed(feature, "@", 2)[, 2]),
    mz = mass + 1.0078,
    MVN.significance = ifelse(feature %in% metab.sign.list, 1, 0)
  ) %>%
  arrange(group)


if (write.res) {
  saveRDS(PISCINA_cor_centroid_cleaned, "Results/sPLS/PISCINA_cor_centroid_cleaned.rds")
  saveRDS(spls.centroid.2.21, "Results/sPLS/spls.centroid.2.21.rds")
  saveRDS(spls.centroid.4.26, "Results/sPLS/spls.centroid.4.26.rds")
  saveRDS(metab.centroids.selected.clusters.21, "Results/sPLS/metab.centroids.selected.clusters.21.rds")
  saveRDS(metab.centroids.selected.clusters.26, "Results/sPLS/metab.centroids.selected.clusters.26.rds")
  saveRDS(metab.features.selected.centroids.21, "Results/sPLS/metab.features.selected.centroids.21.rds")
  saveRDS(metab.features.selected.centroids.26, "Results/sPLS/metab.features.selected.centroids.26.rds")
}

if (write.tab) {
  write.csv(description_summary.21, "Results/Tables/Description_summary.21.csv")
}

# Perf for best models ----------------------------------------------------
if (run.perf) {
  spls.centroid.2.21.perf <- perf(spls.centroid.2.21, validation = "Mfold", folds = 5, nrepeat = 500, criterion = "all")
  spls.centroid.4.26.perf <- perf(spls.centroid.4.26, validation = "Mfold", folds = 5, nrepeat = 500, criterion = "all")

  # MSEP and Q2
  spls.centroid.2.21.MSEP <- spls.centroid.2.21.perf$measures$MSEP$values %>%
    group_by(nrep) %>%
    summarise(mean_val = mean(value)) %>%
    pull(mean_val)

  spls.centroid.4.26.MSEP <- spls.centroid.4.26.perf$measures$MSEP$values %>%
    group_by(nrep) %>%
    summarise(mean_val = mean(value)) %>%
    pull(mean_val)

  spls.centroid.2.21.Q2 <- spls.centroid.2.21.perf$measures$Q2$values %>%
    group_by(nrep) %>%
    summarise(mean_val = mean(value)) %>%
    pull(mean_val)

  spls.centroid.4.26.Q2 <- spls.centroid.4.26.perf$measures$Q2$values %>%
    group_by(nrep) %>%
    summarise(mean_val = mean(value)) %>%
    pull(mean_val)

  spls.centroid.2.21.perf.cor <- perf_cor(spls.centroid.2.21, designID = designID, id.vec = 1:58, nrepeat = 500, keep_value = TRUE)
  spls.centroid.4.26.perf.cor <- perf_cor(spls.centroid.4.26, designID = designID, id.vec = 1:58, nrepeat = 500, keep_value = TRUE)

  spls.centroid.2.21.cor.test <- spls.centroid.2.21.perf.cor$cor.test
  spls.centroid.4.26.cor.test <- spls.centroid.4.26.perf.cor$cor.test

  spls.centroid.models.perf <- data.frame(
    model = rep(c("sPLS on X and Y", "sPLS on Y"), each = 500),
    Q2 = c(spls.centroid.2.21.Q2, spls.centroid.4.26.Q2),
    MSEP = c(spls.centroid.2.21.MSEP, spls.centroid.4.26.MSEP),
    cor_test = c(spls.centroid.2.21.cor.test, spls.centroid.4.26.cor.test)
  )

  saveRDS(spls.centroid.models.perf, "Results/sPLS/spls.centroid.models.perf.rds")
}


# Plot --------------------------------------------------------------------
# Perf
PISCINA_perf_centroid_df <- as.data.frame(do.call(rbind, PISCINA_perf_centroid))

# Read results for the calibration performances of the selected models
spls.centroid.models.perf <- readRDS("Results/sPLS/spls.centroid.models.perf.rds")


# Q2 and correlation
Q2.calibration.plot <- PISCINA_perf_centroid_df %>%
  pivot_longer(Q2_comp1:Q2_comp4, names_to = "comp") %>%
  mutate(comp = str_replace(comp, "Q2_comp", "Component ")) %>%
  ggplot(aes(nY, value, col = as.factor(nX))) +
  geom_point(size = 0.08, alpha = 0.5) +
  facet_wrap(~ factor(comp), scales = "free_y") +
  scale_colour_manual(values = c("#E6A6C7", "#56B4E9", "#69C9B9", "#FD8D3C")) +
  theme_bw() +
  theme(
    legend.position = "top",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  guides(colour = guide_legend(
    override.aes = list(size = 2),
    title = "Number of DBPs (nX)"
  )) +
  scale_x_continuous(breaks = seq(0, 550, 100)) +
  xlab("") +
  ylab("Mean Q2")

cor.calibration.plot <- PISCINA_cor_centroid_cleaned %>%
  filter(set == "Test") %>%
  ggplot(aes(as.numeric(nY), mean, col = as.factor(nX))) +
  geom_point(size = 0.1) +
  theme_bw() +
  scale_colour_manual(values = c("#E6A6C7", "#56B4E9", "#69C9B9", "#FD8D3C")) +
  theme(
    legend.position = "none",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(size = 3), title = "Number of DBPs (nX)")) +
  scale_x_continuous(breaks = seq(0, 550, 100), limits = c(-5, 560), expand = c(0, 0)) +
  scale_y_continuous(breaks = seq(0.7, 0.86, 0.05), limits = c(0.69, 0.87), expand = c(0, 0)) +
  xlab("Number of clusters of metabolic features (nY)") +
  ylab("Mean correlation between\ncomponent 1 of X and Y (test data)") +
  annotate("segment", x = 21, xend = 21, y = 0.69, yend = 0.851, linetype = "dashed", col = "#56B4E9", linewidth = 0.25) +
  annotate("segment", x = 26, xend = 26, y = 0.69, yend = 0.817, linetype = "dashed", col = "#FD8D3C", linewidth = 0.25) +
  annotate("segment", x = -5, xend = 21, y = 0.851, yend = 0.851, linetype = "dashed", col = "#56B4E9", linewidth = 0.25) +
  annotate("segment", x = -5, xend = 26, y = 0.817, yend = 0.817, linetype = "dashed", col = "#FD8D3C", linewidth = 0.25) +
  annotate("text", x = 70, y = 0.75, label = "nX=4, nY=26", col = "#FD8D3C", size = 4) +
  annotate("text", x = 70, y = 0.72, label = "nX=2, nY=21", col = "#56B4E9", size = 4)

# Add in performance metrics
p3.1 <- spls.centroid.models.perf %>%
  ggplot(aes(model, cor_test, fill = model)) +
  geom_boxplot() +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  ) +
  scale_fill_manual(values = c("#56B4E9", "#FD8D3C"), name = "") +
  xlab("") +
  ylab("Correlation in component\nX and component Y") +
  scale_x_discrete(labels = c("sPLS on X\n and Y", "sPLS on Y"))

p3.2 <- spls.centroid.models.perf %>%
  ggplot(aes(model, MSEP, fill = model)) +
  geom_boxplot() +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  ) +
  scale_fill_manual(values = c("#56B4E9", "#FD8D3C"), name = "") +
  xlab("") +
  scale_x_discrete(labels = c("sPLS on X\n and Y", "sPLS on Y"))


p3.3 <- spls.centroid.models.perf %>%
  ggplot(aes(model, Q2, fill = model)) +
  geom_boxplot() +
  theme_bw() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
  ) +
  scale_fill_manual(values = c("#56B4E9", "#FD8D3C"), name = "") +
  xlab("") +
  scale_x_discrete(labels = c("sPLS on X\n and Y", "sPLS on Y"))

p3 <- p3.1 + p3.2 + p3.3

sPLS.calibration <- Q2.calibration.plot / cor.calibration.plot / p3 + plot_layout(heights = c(1, 1, 1)) +
  plot_annotation(tag_levels = "A") & theme(plot.tag = element_text(size = 20))

if (write.fig) {
  ggsave("Results/Figures/PISCINA_sPLS_calibration_Q2_cor_calibrated_models.pdf", sPLS.calibration, height = 9, width = 8)
}
