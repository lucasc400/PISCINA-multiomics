# =============================================================================
# 13-sPLS_cor_centroid.R
#
# Calibration of multilevel sPLS models based on out-of-sample correlation.
# For every combination of nX and nY, the model is refitted on 80% of the
# subjects and the correlation between the first X and Y components is
# computed on the remaining 20% (500 repeats). Runs in parallel on an HPC node
# (see spls.sh).
#
# Input:  Data/Cleaned_Data/expo.piscina.rds
#         Data/Cleaned_Data/covars.piscina.rds
#         Results/Clustering/metabolites_summarised_centroid.rds
# Output: Results/sPLS/PISCINA_sPLS_tune_cor_metab_centroid.rds
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
param.grid <- as.matrix(expand.grid(nX = c(1:4), nY = seq(1, 550, 1)))

t0 <- Sys.time()
no_cores <- 64
cl <- makeCluster(no_cores)
clusterExport(cl, c("tune_spls_cor", "perf_cor", "moe", "param.grid", "designID", "withinVarX", "withinVarY"))
clusterEvalQ(cl, library(mixOmics))
clusterEvalQ(cl, library(dplyr))
res <- parLapply(
  cl = cl, X = 1:nrow(param.grid), fun = tune_spls_cor, param = param.grid,
  Xdata = withinVarX, Ydata = withinVarY, designID = designID, n_subject = 58,
  nrepeat = 500
)
stopCluster(cl)
t1 <- Sys.time()
print(t1 - t0)

dir.create("Results/sPLS", showWarnings = FALSE, recursive = TRUE)
saveRDS(res, "Results/sPLS/PISCINA_sPLS_tune_cor_metab_centroid.rds")
