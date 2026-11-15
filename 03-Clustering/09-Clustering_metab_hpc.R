# =============================================================================
# 09-Clustering_metab_hpc.R
#
# Stability-calibrated consensus clustering of the PISCINA-II metabolic
# features (samples collected before swimming), using sharp::Clustering().
# Each HPC array job evaluates one candidate number of clusters
# (nc = 200, 205, ..., 1695; 300 jobs). See run_consensus_clustering.sh.
#
# Usage:  Rscript Scripts/03-Clustering/09-Clustering_metab_hpc.R <job_id>
# Input:  Data/Cleaned_Data/metab.piscina.denoised.rds
#         Data/Cleaned_Data/covars.piscina.rds
# Output: Results/Clustering/consensus_clustering_broad_<job_id>.rds
# =============================================================================

rm(list = ls())
setwd(Sys.getenv("PISCINA_DIR", unset = ".")) # project root (see README)

library(sharp)

filepath <- "Results/Clustering"
dir.create(filepath, showWarnings = FALSE, recursive = TRUE)

# Parameters
args <- commandArgs(trailingOnly = TRUE)
id <- as.numeric(args[1])
nc_min <- 200
nc_max <- 1695
step <- 5

# Loading data
metab <- readRDS("Data/Cleaned_Data/metab.piscina.denoised.rds")
covars <- readRDS("Data/Cleaned_Data/covars.piscina.rds")

metab.wave1 <- metab[covars$wave == 1, ]

# Consensus clustering of the features (columns) for one number of clusters
print(system.time({
  stab <- Clustering(
    xdata = metab.wave1,
    nc = seq(nc_min, nc_max, by = step)[id],
    row = FALSE
  )
}))
print(stab$nc)
print(stab$Sc)

# Saving the output
saveRDS(stab, paste0(filepath, "/consensus_clustering_broad_", id, ".rds"))
