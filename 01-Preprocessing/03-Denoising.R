# =============================================================================
# 03-Denoising.R
#
# Remove plate (batch) effects from the PISCINA-II metabolomics data. For each
# feature, a linear mixed model with plate as a random intercept is fitted and
# the residuals (after removing the plate random intercept) are kept as the
# "denoised" feature levels.
#
# Input:  Data/Cleaned_Data/metab.piscina.rds
#         Data/Cleaned_Data/metab.piscina.tech.rds
# Output: Data/Cleaned_Data/metab.piscina.denoised.rds
# =============================================================================

rm(list = ls(all = TRUE))
setwd(Sys.getenv("PISCINA_DIR", unset = ".")) # project root (see README)

library(lme4)
source("Scripts/functions.R")

metab.piscina <- readRDS("Data/Cleaned_Data/metab.piscina.rds")
metab.tech <- readRDS("Data/Cleaned_Data/metab.piscina.tech.rds")

# Metabolomics ------------------------------------------------------------
lmer.metab.PISCINA <- Multilmer(ydata = metab.piscina, xdata = metab.tech, random = "Plate", residuals = TRUE)
metab.piscina.denoised <- lmer.metab.PISCINA$residuals
rownames(metab.piscina.denoised) <- metab.tech$XP_ID

saveRDS(metab.piscina.denoised, "Data/Cleaned_Data/metab.piscina.denoised.rds")
