# =============================================================================
# 06-MWAS_PISCINA.R
#
# Metabolome-wide association study (MWAS) in PISCINA-II. Each denoised
# metabolic feature is regressed on each THM exposure (four individual THMs,
# BrTHM and TotalTHM) using a multivariate normal model (generalised least
# squares with an unstructured within-subject correlation), adjusting for the
# subject's mean exposure, age, sex and BMI. The same models are fitted on the
# original metabolomics dataset for comparison.
#
# Input:  Data/Cleaned_Data/*.rds
# Output: Results/MVN/metab.piscina.MVN.res.rds
#         Results/MVN/metab.piscina.MVN.res.old.rds
# =============================================================================

rm(list = ls())
setwd(Sys.getenv("PISCINA_DIR", unset = ".")) # project root (see README)

library(doBy)
library(nlme)

source("Scripts/functions.R")

expo.composite <- readRDS("Data/Cleaned_Data/expo.piscina.composite.rds")
metab <- readRDS("Data/Cleaned_Data/metab.piscina.denoised.rds")
metab.old <- readRDS("Data/Cleaned_Data/metab.piscina.old.rds")
covars <- readRDS("Data/Cleaned_Data/covars.piscina.rds")

data_gls <- covars[, c(
  "ID", "XP_ID", "Per", "wave", "age", "sex", "BMI", "data_month", "ph", "tempw",
  "tempa", "perc_HRmax_69", "kcal"
)]

dir.create("Results/MVN", showWarnings = FALSE, recursive = TRUE)

# Fit the MVN model for every feature in `ydata` and every exposure
run_mwas <- function(ydata) {
  res <- list()
  for (CurrExpo in colnames(expo.composite)) {
    print(CurrExpo)
    data_gls$Expo <- expo.composite[, CurrExpo]
    res[[CurrExpo]] <- multigls(
      ydata = ydata, xdata = data_gls, covars = c("age", "sex", "BMI"),
      progressbar = FALSE, standardise_outcome = FALSE
    )
  }
  res
}

# MWAS --------------------------------------------------------------------
coef_denoised <- run_mwas(metab)
saveRDS(coef_denoised, "Results/MVN/metab.piscina.MVN.res.rds")

# MWAS on the original dataset --------------------------------------------
coef_denoised_old <- run_mwas(metab.old)
saveRDS(coef_denoised_old, "Results/MVN/metab.piscina.MVN.res.old.rds")
