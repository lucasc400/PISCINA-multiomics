# =============================================================================
# 07-MWAS_MCC.R
#
# MWAS in MCC-Spain (validation data). Each metabolic feature is regressed on
# each residential THM exposure (recent and long-term) with linear models
# adjusted for sex, age, BMI, smoking and education. Sensitivity analyses
# additionally adjust for the WCRF score or residential area.
#
# Input:  Data/Cleaned_Data/covars.expo.MCC.rds
#         Data/Cleaned_Data/metab.MCC.imputed.rds
# Output: Results/MCC_linear/*.rds (if write.res)
# =============================================================================

rm(list = ls())
setwd(Sys.getenv("PISCINA_DIR", unset = ".")) # project root (see README)

library(dplyr)

covars_expo <- readRDS("Data/Cleaned_Data/covars.expo.MCC.rds")
metab <- readRDS("Data/Cleaned_Data/metab.MCC.imputed.rds")

source("Scripts/functions.R")

dir.create("Results/MCC_linear", showWarnings = FALSE, recursive = TRUE)
write.res <- FALSE

# Prepare data ------------------------------------------------------------
# Columns 10-21: recent (no suffix) and long-term ("_18") THM exposures
expo.list <- colnames(covars_expo)[10:21]
covars_expo$Sex <- as.factor(ifelse(covars_expo$Sex == "Hombre", "Male", "Female"))
covars_expo$Ever_smoking <- as.factor(covars_expo$Ever_smoking)
covars_expo$Education <- as.factor(covars_expo$Education)

# Main model --------------------------------------------------------------
MWAS.res.list <- list()
for (expo in expo.list) {
  myMWAS <- RunLinear(covars_expo, expo, metab, conf = c("Sex", "Age", "BMI", "Ever_smoking", "Education"))
  myMWAS <- cbind(
    expo = expo,
    expo_type = ifelse(grepl("_18", expo), "Long-term", "Recent"),
    myMWAS
  )
  MWAS.res.list[[expo]] <- myMWAS
}
MWAS.res.df <- do.call(rbind, MWAS.res.list)
colnames(MWAS.res.df)[which(colnames(MWAS.res.df) == "outcome")] <- "metab"
rownames(MWAS.res.df) <- c()

# Save results
if (write.res) {
  saveRDS(MWAS.res.df, "Results/MCC_linear/MCC.MWAS.res.rds")
}

# Significant features (BH-adjusted p < 0.05, recent exposures) -----------
sign.metab.FDR.recent <- MWAS.res.df %>%
  filter(expo_type == "Recent") %>%
  group_by(expo) %>%
  mutate(pval.adjusted = p.adjust(`p-value`, method = "BH")) %>%
  filter(pval.adjusted < 0.05) %>%
  ungroup() %>%
  pull(metab) %>%
  unique()

if (write.res) {
  saveRDS(sign.metab.FDR.recent, "Results/MCC_linear/sign.fdr.metab.recent.rds")
}

# Sensitivity analysis: adjusting for WCRF score --------------------------
MWAS.WCRF.res.list <- list()

for (expo in expo.list) {
  myMWAS <- RunLinear(covars_expo[complete.cases(covars_expo$WCRF), ],
    expo,
    metab[complete.cases(covars_expo$WCRF), ],
    conf = c("Sex", "Age", "BMI", "Ever_smoking", "Education", "WCRF")
  )
  myMWAS <- cbind(
    expo = expo,
    expo_type = ifelse(grepl("_18", expo), "Recent", "Long-term"),
    myMWAS
  )
  MWAS.WCRF.res.list[[expo]] <- myMWAS
}

MWAS.WCRF.res.df <- do.call(rbind, MWAS.WCRF.res.list)
colnames(MWAS.WCRF.res.df)[which(colnames(MWAS.WCRF.res.df) == "outcome")] <- "metab"
rownames(MWAS.WCRF.res.df) <- c()

# Save results
if (write.res) {
  saveRDS(MWAS.WCRF.res.df, "Results/MCC_linear/MCC.MWAS.WCRF.res.rds")
}

# Sensitivity analysis: adjusting for residential area -------------------
MWAS.area.res.list <- list()

for (expo in expo.list) {
  myMWAS <- RunLinear(covars_expo,
    expo,
    metab,
    conf = c("Sex", "Age", "BMI", "Ever_smoking", "Education", "Area")
  )
  myMWAS <- cbind(
    expo = expo,
    expo_type = ifelse(grepl("_18", expo), "Recent", "Long-term"),
    myMWAS
  )
  MWAS.area.res.list[[expo]] <- myMWAS
}

MWAS.area.res.df <- do.call(rbind, MWAS.area.res.list)
colnames(MWAS.area.res.df)[which(colnames(MWAS.area.res.df) == "outcome")] <- "metab"
rownames(MWAS.area.res.df) <- c()

# Save results
if (write.res) {
  saveRDS(MWAS.area.res.df, "Results/MCC_linear/MCC.MWAS.area.res.rds")
}
