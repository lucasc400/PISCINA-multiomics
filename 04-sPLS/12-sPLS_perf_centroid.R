# =============================================================================
# 12-sPLS_perf_centroid.R
#
# Calibration of multilevel sPLS models (X = four THMs, Y = metabolic cluster
# centroids). For every combination of the number of selected exposures (nX)
# and clusters (nY), Q2 and MSEP are estimated by 5-fold cross-validation
# (500 repeats). Runs in parallel on an HPC node (see spls.sh).
#
# Input:  Data/Cleaned_Data/expo.piscina.rds
#         Data/Cleaned_Data/covars.piscina.rds
#         Results/Clustering/metabolites_summarised_centroid.rds
# Output: Results/sPLS/PISCINA_sPLS_centroid_perf.rds
# =============================================================================

rm(list = ls())
setwd(Sys.getenv("PISCINA_DIR", unset = ".")) # project root (see README)

library(mixOmics)
library(parallel)
library(dplyr)

source("Scripts/functions.R")

expo <- readRDS("Data/Cleaned_Data/expo.piscina.rds")
metab <- readRDS("Results/Clustering/metabolites_summarised_centroid.rds")
covars <- readRDS("Data/Cleaned_Data/covars.piscina.rds")

# Within-subject variation (multilevel decomposition)
designID <- as.matrix(as.numeric(as.factor(covars$ID)))
withinVarX <- withinVariation(expo, design = designID)
withinVarY <- withinVariation(metab, design = designID)


# Tuning ------------------------------------------------------------------
param.grid <- as.matrix(expand.grid(nX = c(1:4), nY = seq(1, 550, 1), comp = 4))

t0 <- Sys.time()
no_cores <- 128
cl <- makeCluster(no_cores)
clusterExport(cl, c("pls.perf.par", "param.grid", "designID", "withinVarX", "withinVarY"))
clusterEvalQ(cl, library(mixOmics))
clusterEvalQ(cl, library(dplyr))
res <- parLapply(
  cl = cl, X = 1:nrow(param.grid), fun = pls.perf.par, param = param.grid,
  Xdata = withinVarX, Ydata = withinVarY, ID = designID
)
stopCluster(cl)
t1 <- Sys.time()
print(t1 - t0)

dir.create("Results/sPLS", showWarnings = FALSE, recursive = TRUE)
saveRDS(res, "Results/sPLS/PISCINA_sPLS_centroid_perf.rds")
