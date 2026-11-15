# =============================================================================
# 01-Imputation.R
#
# PISCINA-II metabolomics: remove features missing in more than 60% of samples,
# log2-transform intensities and impute remaining missing values with QRILC.
#
# Input:  Data/New dataset/metabolites_raw.rds
# Output: Data/Cleaned_Data/metabolites.piscina.log.impute.rds
# =============================================================================

rm(list = ls(all = TRUE))
setwd(Sys.getenv("PISCINA_DIR", unset = ".")) # project root (see README)

library(imputeLCMD)


# Read data ---------------------------------------------------------------
metab.piscina <- readRDS("Data/New dataset/metabolites_raw.rds")

# Keep features with at most 60% missing values
metab.piscina.60 <- metab.piscina[, apply(metab.piscina, 2, function(x) sum(is.na(x)) / length(x) <= 0.6)]

# Imputation of missing values (QRILC on log2 intensities) -----------------
metab.piscina.60.log <- log2(metab.piscina.60)
metab.piscina.60.log.t <- t(metab.piscina.60.log)

set.seed(400)
impute1 <- impute.QRILC(metab.piscina.60.log.t)
metab.piscina.imp <- t(impute1[[1]])

# Save --------------------------------------------------------------------
dir.create("Data/Cleaned_Data", showWarnings = FALSE)
saveRDS(metab.piscina.imp, "Data/Cleaned_Data/metabolites.piscina.log.impute.rds")
