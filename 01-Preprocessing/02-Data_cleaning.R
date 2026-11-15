# =============================================================================
# 02-Data_cleaning.R
#
# Clean and align all datasets used in the analysis.
#   PISCINA-II: metabolomics (with plate/well information), covariates,
#               exhaled-breath THM exposures, transcriptomics, proteomics and
#               the original (previously processed) metabolomics dataset.
#   MCC-Spain:  covariates and residential THM exposures for cancer-free
#               controls, and metabolomics (missingness filter, QRILC
#               imputation).
#
# Input:  raw files under Data/ (see paths below)
# Output: Data/Cleaned_Data/*.rds
# =============================================================================

rm(list = ls(all = TRUE))
setwd(Sys.getenv("PISCINA_DIR", unset = ".")) # project root (see README)

library(dplyr)
library(stringr)
library(readxl)
library(foreign)
library(imputeLCMD)


# PISCINA-II --------------------------------------------------------------
# Metabolomics sample worklist and sample ID mapping ----------------------
metab.worklist <- read_excel("Data/New dataset/AnalsisWorklist.xlsx", sheet = 1, col_names = FALSE)
metab.mapping <- read.csv("Data/New dataset/new_feature_table.csv", nrows = 3, header = FALSE)

# Clean mapping
metab.mapping <- as.data.frame(t(metab.mapping[, 2:121]))
colnames(metab.mapping) <- c("Metx_data_file_name_POS", "ID_PREPOST", "XP_ID_original")
metab.mapping$Per <- sapply(metab.mapping$ID_PREPOST, function(x) str_split(x, "_")[[1]][3])
metab.mapping$ID <- as.integer(sapply(metab.mapping$ID_PREPOST, function(x) str_split(x, "_")[[1]][2]))
metab.mapping$XP_ID <- paste0(metab.mapping$ID, ifelse(metab.mapping$Per == "PRE", 1, 2))

metab.worklist <- metab.worklist[, 4:5]
colnames(metab.worklist) <- c("Plate_row_column", "Metx_data_file_name_POS")
metab.mapping <- metab.mapping[, !(colnames(metab.mapping) %in% "XP_ID_original")]
rownames(metab.mapping) <- metab.mapping$XP_ID

# Clean worklist: remove calibration (CC), blank and QC (LQC) injections and
# extract plate, row and column positions
metab.worklist.clean <- metab.worklist[-grep("CC|LQC|Blank", metab.worklist$Metx_data_file_name_POS), ]
metab.worklist.clean$Plate <- as.character(sapply(metab.worklist.clean$Plate_row_column, function(x) substr(str_split(x, "-")[[1]][1], 2, 2)))
metab.worklist.clean$Row <- as.character(sapply(metab.worklist.clean$Plate_row_column, function(x) substr(str_split(x, "-")[[1]][2], 1, 1)))
metab.worklist.clean$Column <- as.character(sapply(metab.worklist.clean$Plate_row_column, function(x) substr(str_split(x, "-")[[1]][2], 2, 2)))

# Join worklist and mapping
metab.worklist.clean.mapped <- merge(metab.worklist.clean, metab.mapping, all.x = TRUE, by = "Metx_data_file_name_POS")
rownames(metab.worklist.clean.mapped) <- metab.worklist.clean.mapped$XP_ID

# Metabolomics, covariates and exposures ----------------------------------
metab.piscina <- readRDS("Data/Cleaned_Data/metabolites.piscina.log.impute.rds")
covars <- readRDS("Data/New dataset/PISCINA_covars_imputed.rds")
expo_ID <- readRDS("Data/New dataset/Full_expo.RDS")

rownames(covars) <- covars$Row.names
drop.cols <- c("Row.names", "ID.y", "XP_ID.y", "Per.y")
covars <- covars %>%
  dplyr::select(-all_of(drop.cols)) %>%
  rename("XP_ID" = "XP_ID.x", "Per" = "Per.x", "ID" = "ID.x") %>%
  arrange(XP_ID)

# Align metabolites and covariates
metab.piscina <- metab.piscina[as.character(rownames(covars)), ]

# Check all rownames are the same
all(rownames(covars) == rownames(expo_ID))

# Keep the four individual THMs (drop ID columns)
expo <- expo_ID[, 1:4]
colnames(expo) <- c("Chloroform", "BDCM", "DBCM", "Bromoform")

# Create composite exposure variables
expo.composite <- expo %>%
  mutate(
    BrTHM = BDCM + DBCM + Chloroform,
    TotalTHM = Chloroform + BDCM + DBCM + Bromoform
  )

# Add a "wave" for covariates where 1 = "Pre", 2 = "Post"
covars$wave <- as.numeric(covars$Per == "post") + 1

# Set energy expenditure (kcal) to 0 before swimming
covars$kcal <- ifelse(covars$Per == "pre", 0, covars$kcal)

# Align plate with metabolites
metab.worklist.clean.mapped <- metab.worklist.clean.mapped[rownames(metab.piscina), ]

# Original (previously processed) metabolomics dataset --------------------
load("Data/Old dataset/PISCINA_data_metab_imputed_60%.RData")

# Columns 6472-6475 are not metabolic features
metab.piscina.old <- Full_peaks_imputed[rownames(covars), -c(6472:6475)]

# Transcriptomics ---------------------------------------------------------
mRNA.raw <- read.table("Data/Transcripts/PISC2_Tx_mRNA_procd_background_corrected.txt", header = TRUE)
rownames(mRNA.raw) <- mRNA.raw[, 1]
mRNA.raw <- mRNA.raw[, 2:ncol(mRNA.raw)]
colnames(mRNA.raw) <- gsub("X", "", colnames(mRNA.raw))

mRNA.id <- as.data.frame(t(read.csv("Data/Transcripts/PISCINA_mRNA_IDs.csv")))
colnames(mRNA.id) <- mRNA.id[1, ]
mRNA.id <- mRNA.id[2:nrow(mRNA.id), ]
rownames(mRNA.id) <- mRNA.id[, 2]

mRNA.id <- mRNA.id[colnames(mRNA.raw), ]
all(rownames(mRNA.id) == colnames(mRNA.raw))
mRNA.id$ID <- as.integer(str_split_fixed(mRNA.id$UniqueSampleID, "_", n = 3)[, 2])
mRNA.id$status <- str_split_fixed(mRNA.id$UniqueSampleID, "_", n = 3)[, 3]
mRNA.id$XID <- paste0(mRNA.id$ID, ifelse(mRNA.id$status == "PRE", 1, 2))

colnames(mRNA.raw) <- mRNA.id$XID
mRNA.clean <- t(mRNA.raw)

# Proteomics (immune markers) ---------------------------------------------
prot.raw <- read.csv("Data/Prot/Data/PISC2_PrePost_cov_expos_prot13.csv", header = TRUE)
prot.clean <- prot.raw %>%
  arrange(UniqueSubjectID, desc(timepoint)) %>%
  mutate(XP_ID = paste0(as.numeric(gsub("PISC2_", "", UniqueSubjectID)), ifelse(timepoint == "PRE", 1, 2))) %>%
  relocate(XP_ID) %>%
  dplyr::select(XP_ID, MCP.1_original:IL.1ra_imputed) %>%
  rename_with(.col = MCP.1_original:IL.1ra_imputed, .fn = ~ str_remove(., "_.*")) %>%
  tibble::column_to_rownames(var = "XP_ID") %>%
  rename(
    CCL2 = "MCP.1",
    `IL-8` = "IL.8.C",
    `EGF` = "EGF.2",
    MPO = "MPO.5",
    `IL-17` = "IL.17",
    CCL22 = "MDC.CC",
    `G-CSF` = "G.CSF",
    CCL11 = "Eotaxi",
    CXCL10 = "IP.10",
    Periostin = "Perios",
    IL1ra = "IL.1ra"
  ) %>%
  as.data.frame()

prot.clean <- prot.clean[as.character(rownames(covars)), ]
prot.clean <- prot.clean[complete.cases(prot.clean), ]

# Save PISCINA-II data ----------------------------------------------------
dir.create("Data/Cleaned_Data", showWarnings = FALSE)

saveRDS(expo, "Data/Cleaned_Data/expo.piscina.rds")
saveRDS(expo.composite, "Data/Cleaned_Data/expo.piscina.composite.rds")
saveRDS(covars, "Data/Cleaned_Data/covars.piscina.rds")
saveRDS(metab.piscina, "Data/Cleaned_Data/metab.piscina.rds")
saveRDS(metab.worklist.clean.mapped, "Data/Cleaned_Data/metab.piscina.tech.rds")
saveRDS(mRNA.clean, "Data/Cleaned_Data/mRNA.piscina.rds")
saveRDS(prot.clean, "Data/Cleaned_Data/protein.piscina.rds")
saveRDS(metab.piscina.old, "Data/Cleaned_Data/metab.piscina.old.rds")


# MCC-Spain ---------------------------------------------------------------
covars.MCC.raw <- read.dta("Data/MCC/BD_exposomics_feb2017_v2.dta")
metab.MCC <- read.csv("Data/MCC/MCC_Spain feature table 020223.csv", header = TRUE)

metab.MCC <- t(metab.MCC)
metab.ID <- rownames(metab.MCC)
colnames(metab.MCC) <- metab.MCC[1, ]
metab.MCC <- as.data.frame(apply(metab.MCC[-1, ], 2, as.numeric))
rownames(metab.MCC) <- str_replace_all(str_remove(str_split_fixed(metab.ID[-1], "_", 2)[, 1], "SP."), "\\.", "")

# Covariates and residential THM exposures: recent (no suffix) and
# long-term ("_18" suffix) exposures
covars.MCC.clean <- with(covars.MCC.raw, data.frame(
  ID = newid,
  Cancer = tipo_cancer,
  Area = area,
  Age = edad,
  Sex = sexo,
  Education = education_basic_final,
  BMI = BMIcomplement,
  Ever_smoking = fumador_ever,
  WCRF = score_WCRF_simple,
  Chloroform = chcl3_ex2ye_ly,
  DBCM = chclbr2_ex2ye_ly,
  BDCM = chcl2br_ex2ye_ly,
  Bromoform = chbr3_ex2ye_ly,
  BrTHM = chclbr2_ex2ye_ly + chcl2br_ex2ye_ly + chbr3_ex2ye_ly,
  TotalTHM = chcl3_ex2ye_ly + chclbr2_ex2ye_ly + chcl2br_ex2ye_ly + chbr3_ex2ye_ly,
  Chloroform_18 = chcl3_18act2ye,
  DBCM_18 = chclbr2_18act2ye,
  BDCM_18 = chcl2br_18act2ye,
  Bromoform_18 = chbr3_18act2ye,
  BrTHM_18 = chclbr2_18act2ye + chcl2br_18act2ye + chbr3_18act2ye,
  TotalTHM_18 = chcl3_18act2ye + chclbr2_18act2ye + chcl2br_18act2ye + chbr3_18act2ye
))

rownames(covars.MCC.clean) <- covars.MCC.clean$ID

# Only use cancer-free controls
covars.MCC.clean <- covars.MCC.clean[covars.MCC.clean$Cancer == 0, ]

# Exclude subjects with missing data in any variable other than WCRF
# (1 subject with missing exposure). The 26 subjects with missing WCRF are kept
# here; WCRF is only used in a sensitivity analysis.
covars.MCC.clean <- covars.MCC.clean[complete.cases(covars.MCC.clean[, -which(colnames(covars.MCC.clean) == "WCRF")]), ]

# Align ID
ID.selected <- intersect(rownames(metab.MCC), covars.MCC.clean$ID)
metab.MCC <- metab.MCC[ID.selected, ]
covars.MCC.clean <- covars.MCC.clean[ID.selected, ]


# MCC-Spain imputation ----------------------------------------------------
# Number of features missing in more than 60% of samples
table(apply(metab.MCC, 2, function(x) sum(is.na(x)) / length(x)) < 0.6)

metab.MCC.60 <- metab.MCC[, apply(metab.MCC, 2, function(x) sum(is.na(x)) / length(x)) < 0.6]

# Set negative intensities to missing
which(apply(metab.MCC.60, 2, function(x) any(x < 0, na.rm = TRUE)))

metab.MCC.60 <- apply(metab.MCC.60, 2, function(x) {
  ifelse(x < 0, NA, x)
})

metab.MCC.60.log2 <- log2(metab.MCC.60)
metab.MCC.60.log2.t <- t(metab.MCC.60.log2)

set.seed(125)
impute1 <- impute.QRILC(metab.MCC.60.log2.t)
metab.MCC.imputed <- t(impute1[[1]])

# Save MCC-Spain data -----------------------------------------------------
saveRDS(covars.MCC.clean, "Data/Cleaned_Data/covars.expo.MCC.rds")
saveRDS(metab.MCC.imputed, "Data/Cleaned_Data/metab.MCC.imputed.rds")
