# =============================================================================
# 10-Merge_clustering_results.R
#
# Combine the stability scores of the 300 consensus clustering jobs and keep
# the number of clusters with the highest stability score.
#
# Input:  Results/Clustering/consensus_clustering_broad_<1..300>.rds
# Output: Results/Clustering/consensus_clustering_broad_calibration.rds
#         Results/Clustering/consensus_clustering_broad.rds
# =============================================================================

rm(list = ls())
setwd(Sys.getenv("PISCINA_DIR", unset = ".")) # project root (see README)

clust_dir <- "Results/Clustering"
n_jobs <- 300 # number of HPC array jobs in 09-Clustering_metab_hpc.R

# Calibration of hyper-parameters
nc <- sc <- rep(NA, n_jobs)
for (i in 1:n_jobs) {
  out <- readRDS(file.path(clust_dir, paste0("consensus_clustering_broad_", i, ".rds")))
  nc[i] <- out$nc[1]
  sc[i] <- out$Sc[1]
}
calib <- sc
names(calib) <- nc
saveRDS(calib, file.path(clust_dir, "consensus_clustering_broad_calibration.rds"))


# Extracting the calibrated results (the last candidate, nc = 1695, is
# excluded)
hat_id <- which.max(sc[1:(n_jobs - 1)])
print(hat_id)
print(nc[hat_id])
print(sc[hat_id])
out <- readRDS(file.path(clust_dir, paste0("consensus_clustering_broad_", hat_id, ".rds")))
saveRDS(out, file.path(clust_dir, "consensus_clustering_broad.rds"))
