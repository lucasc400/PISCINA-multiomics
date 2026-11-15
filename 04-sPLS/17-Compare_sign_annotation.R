# =============================================================================
# 17-Compare_sign_annotation.R
#
# Validation in MCC-Spain. The metabolic features of the 21 sPLS-selected
# clusters are matched to MCC-Spain features (mass difference <= 7 ppm and
# retention time difference <= 0.1 min), combined with their annotation and
# detection rates, and the MCC-Spain MWAS results are tabulated for the
# matched annotated compounds (main model and sensitivity analyses).
#
# Input:  Data/Cleaned_Data/*.rds, Data/New dataset/metabolites_raw.rds,
#         Results/MVN/*, Results/MCC_linear/*, Results/sPLS/*,
#         Results/Clustering/*, Results/Annotation/*
# Output: Results/Validation/*.csv, Results/Tables/*.csv (if write.tab)
# =============================================================================

rm(list = ls(all = TRUE))
setwd(Sys.getenv("PISCINA_DIR", unset = ".")) # project root (see README)

library(dplyr)
library(tidyr)
library(ggplot2)
library(stringr)

source("Scripts/functions.R")

# Read data ---------------------------------------------------------------
expo.piscina <- readRDS("Data/Cleaned_Data/expo.piscina.rds")
metab.mcc <- readRDS("Data/Cleaned_Data/metab.MCC.imputed.rds")

sign.metab.piscina.188 <- readRDS("Results/MVN/sign.metab.piscina.188.rds") # PISCINA-II MWAS (Bonferroni)
sign.fdr.metab.recent <- readRDS("Results/MCC_linear/sign.fdr.metab.recent.rds") # MCC-Spain MWAS (FDR)
sign.fdr.metab.longterm <- readRDS("Results/MCC_linear/sign.fdr.metab.longterm.rds")

spls.selected.metab.centroids.21 <- readRDS("Results/sPLS/metab.features.selected.centroids.21.rds")
cluster.description <- readRDS("Results/Clustering/Description_summary.rds")
annotation <- readRDS("Results/Annotation/annotation.sPLS.selected.metab.rds")

dir.create("Results/Validation", showWarnings = FALSE, recursive = TRUE)

write.tab <- FALSE


# Match sPLS-selected features with MCC-Spain -----------------------------
spls.matched <- data.frame(compare_metab(spls.selected.metab.centroids.21, colnames(metab.mcc)))
colnames(spls.matched) <- c("PISCINA", "MCC")
spls.matched <- spls.matched %>%
  separate(PISCINA, sep = "@", into = c("PISCINA.mass", "PISCINA.rt"), remove = FALSE, convert = TRUE) %>%
  separate(MCC, sep = "@", into = c("MCC.mass", "MCC.rt"), remove = FALSE, convert = TRUE) %>%
  mutate(
    mass.delta = PISCINA.mass - MCC.mass,
    mass.delta.ppm = (PISCINA.mass - MCC.mass) / PISCINA.mass * 1e6,
    rt.delta = PISCINA.rt - MCC.rt
  )

spls_cluster_matched_df <- spls.matched %>%
  mutate(
    MCC.FDR.sign.recent = ifelse(MCC %in% sign.fdr.metab.recent, "Significant", "Not significant"),
    MCC.FDR.sign.longterm = ifelse(MCC %in% sign.fdr.metab.longterm, "Significant", "Not significant")
  ) %>%
  relocate(PISCINA, MCC) %>%
  arrange(PISCINA.rt)

spls_cluster_matched_MCC_df <- spls_cluster_matched_df %>%
  arrange(desc(MCC.FDR.sign.longterm)) %>%
  dplyr::select(PISCINA, MCC, mass.delta, mass.delta.ppm, rt.delta, MCC.FDR.sign.longterm) %>%
  left_join(., cluster.description %>% dplyr::select(feature, group),
    by = c("PISCINA" = "feature")
  )

spls_selected_metab_validation <- cluster.description %>%
  filter(feature %in% spls.selected.metab.centroids.21) %>%
  dplyr::select(feature, group) %>%
  left_join(., spls_cluster_matched_MCC_df %>% dplyr::select(PISCINA, MCC), by = c("feature" = "PISCINA")) %>%
  mutate(
    `Significant in PISCINA` = ifelse(feature %in% sign.metab.piscina.188, "Significant in PISCINA", "Not significant in PISCINA"),
    `Matched in MCC` = ifelse(feature %in% spls_cluster_matched_MCC_df$PISCINA, "Matched", "Not matched"),
    `Significant in MCC` = ifelse(feature %in% (spls_cluster_matched_MCC_df %>% filter(MCC.FDR.sign.longterm == "Significant") %>% pull(PISCINA)),
      "Significant in MCC", "Not significant in MCC"
    )
  ) %>%
  rename(`Feature in PISCINA` = feature, `Feature in MCC` = MCC) %>%
  arrange(group)

if (write.tab) {
  write.csv(spls_cluster_matched_MCC_df, "Results/Validation/sPLS_cluster_matched_MCC_df.csv", row.names = FALSE)
  write.csv(spls_selected_metab_validation, "Results/Validation/spls_selected_metab_validation.csv", row.names = FALSE)
}

# Matched features significant in either study
spls.matched.centroids.any.df <- spls.matched %>%
  mutate(
    MCC.significant = ifelse(MCC %in% sign.fdr.metab.recent | MCC %in% sign.fdr.metab.longterm, "Significant", "Not significant"),
    PISCINA.significant = ifelse(PISCINA %in% sign.metab.piscina.188, "Significant", "Not significant")
  ) %>%
  relocate(PISCINA, MCC) %>%
  arrange(desc(PISCINA.significant), desc(MCC.significant), PISCINA.rt)


# Tables of selected features ---------------------------------------------
spls.selected.metab.centroids.21.df <- data.frame(metab = spls.selected.metab.centroids.21) %>%
  mutate(
    `Retention time` = as.numeric(str_split_fixed(metab, "@", 2)[, 2]),
    Mass = as.numeric(str_split_fixed(metab, "@", 2)[, 1])
  ) %>%
  arrange(`Retention time`)

sign.metab.piscina.188.df <- data.frame(metab = sign.metab.piscina.188) %>%
  mutate(
    `Retention time` = as.numeric(str_split_fixed(metab, "@", 2)[, 2]),
    Mass = as.numeric(str_split_fixed(metab, "@", 2)[, 1])
  ) %>%
  arrange(`Retention time`)

sign.metab.piscina.spls.df <- data.frame(
  metab = unique(c(spls.selected.metab.centroids.21))
) %>%
  mutate(
    `Retention time` = as.numeric(str_split_fixed(metab, "@", 2)[, 2]),
    Mass = as.numeric(str_split_fixed(metab, "@", 2)[, 1]),
    sPLS.clusters = ifelse(metab %in% spls.selected.metab.centroids.21, "Selected", "Not selected"),
    `Significance in univariate` = ifelse(metab %in% sign.metab.piscina.188, "Significant", "Not significant")
  ) %>%
  arrange(desc(`Significance in univariate`), desc(sPLS.clusters), `Retention time`)

if (write.tab) {
  write.csv(sign.metab.piscina.188.df, "Results/Tables/PISCINA_univariate_selected_metab.csv")
  write.csv(sign.metab.piscina.spls.df, "Results/Tables/PISCINA_sPLS_selected_metab.csv")
  write.csv(spls.selected.metab.centroids.21.df, "Results/Tables/PISCINA_sPLS_selected_centroids_metab.csv")
  write.csv(spls.matched.centroids.any.df, "Results/Tables/PISCINA_matched_metab_any.csv")
}


# Annotation of the sPLS-selected features --------------------------------
spls_selected_metab_validation_annotation <- spls_selected_metab_validation %>%
  left_join(., annotation %>%
    dplyr::select(-c(Annotation_group, group)),
  by = c(`Feature in PISCINA` = "feature")
  ) %>%
  rename(
    Cluster = group,
    `MSI level` = MSI.confidence
  )

spls_selected_metab_validation_annotation_formatted <- spls_selected_metab_validation_annotation %>%
  mutate(
    `Significant in MCC` = ifelse(`Significant in MCC` == "Significant in MCC", "Yes", ""),
    `Significant in PISCINA` = ifelse(`Significant in PISCINA` == "Significant in PISCINA", "Yes", ""),
    `Matched in MCC` = ifelse(`Matched in MCC` == "Matched", "Yes", ""),
    Annotation = ifelse(grepl("Unknown", Annotation), "", Annotation),
    `MSI level` = ifelse(is.na(`MSI level`), "", `MSI level`),
    `Feature in MCC` = ifelse(is.na(`Feature in MCC`), "", `Feature in MCC`)
  )

formatted.colnames <- c(
  "Feature in PISCINA", "Cluster", "Matched in MCC", "Feature in MCC",
  "Significant in PISCINA", "Significant in MCC", "Annotation", "MSI level"
)
spls_selected_metab_validation_annotation_formatted <- spls_selected_metab_validation_annotation_formatted[, formatted.colnames]

# Detection rate (before imputation) before and after swimming
metab.piscina.raw <- readRDS("Data/New dataset/metabolites_raw.rds")

# Align ID
metab.piscina.raw.116 <- metab.piscina.raw[rownames(expo.piscina), ]
metab.raw.72 <- metab.piscina.raw.116[, spls.selected.metab.centroids.21]

n.subjects <- 58
metab.nmiss.before <- apply(metab.raw.72[seq(1, 118, 2), ], 2, function(x) sum(is.na(x)))
metab.nmiss.after <- apply(metab.raw.72[seq(2, 118, 2), ], 2, function(x) sum(is.na(x)))

metab.det.df <- data.frame(
  metab = names(metab.nmiss.before),
  before = round((n.subjects - metab.nmiss.before) / n.subjects * 100, 1),
  after = round((n.subjects - metab.nmiss.after) / n.subjects * 100, 1)
)

spls_selected_metab_validation_annotation_formatted_detect <-
  merge(spls_selected_metab_validation_annotation_formatted,
    metab.det.df,
    by.x = "Feature in PISCINA",
    by.y = "metab"
  ) %>%
  arrange(Cluster)

if (write.tab) {
  write.csv(spls_selected_metab_validation_annotation, "Results/Tables/spls_selected_metab_validation_annotation.csv", row.names = FALSE)
  write.csv(spls_selected_metab_validation_annotation_formatted, "Results/Tables/spls_selected_metab_validation_annotation_formatted.csv", row.names = FALSE)
  write.csv(spls_selected_metab_validation_annotation_formatted_detect, "Results/Tables/spls_selected_metab_validation_annotation_detection_formatted.csv", row.names = FALSE)
}


# MCC-Spain results for the matched annotated compounds -------------------
MCC.MWAS.res <- readRDS("Results/MCC_linear/MCC.MWAS.res.rds")
MCC.MWAS.WCRF.res <- readRDS("Results/MCC_linear/MCC.MWAS.WCRF.res.rds")
MCC.MWAS.area.res <- readRDS("Results/MCC_linear/MCC.MWAS.area.res.rds")
MCC.MWAS.noed.res <- readRDS("Results/MCC_linear/MCC.MWAS.noed.res.rds")
MCC.MWAS.noba.res <- readRDS("Results/MCC_linear/MCC.MWAS.noba.res.rds")

# Annotated PISCINA-II compounds matched to MCC-Spain features
matched.annotation <- c(
  "208.0851@1.8680028" = "L-Kynurenine",
  "501.2873@6.9139724" = "LysoPE (20:4)",
  "285.1375@5.988003" = "Piperine",
  "360.1941@5.0279937" = "Cortisone",
  "217.1312@1.3240031" = "Propionylcarnitine (C3:0)",
  "161.1051@0.6729988" = "L-Carnitine",
  "168.0284@1.095001" = "Uric acid"
)
metab.MCC.matched <- names(matched.annotation)

# Keep the matched compounds and the recent exposures of an MCC-Spain MWAS
select_matched <- function(res) {
  res %>%
    mutate(expo = str_remove(expo, "_18")) %>%
    filter(
      expo_type == "Recent",
      metab %in% metab.MCC.matched
    ) %>%
    dplyr::select(-expo_type) %>%
    mutate(Annotation = unname(matched.annotation[as.character(metab)])) %>%
    arrange(metab)
}

# Number formatting for the tables
format_pvalue <- function(p) {
  case_when(
    p >= 0.1 ~ as.character(round(p, 2)),
    p >= 0.01 ~ as.character(round(p, 3)),
    p >= 0.001 ~ as.character(round(p, 4)),
    TRUE ~ format(p, digits = 2)
  )
}

format_estimate_main <- function(x) {
  case_when(
    abs(x) <= 0.001 ~ format(x, digits = 1),
    abs(x) <= 0.01 ~ format(x, digits = 2),
    TRUE ~ format(x, digits = 3)
  )
}

format_estimate_sensitivity <- function(x) {
  case_when(
    abs(x) <= 0.001 ~ as.character(format(x, digits = 1)),
    TRUE ~ as.character(round(x, digits = 3))
  )
}

# Wide table: one row per compound and metric, one column per exposure
tabulate_matched <- function(res, format_estimate) {
  res %>%
    dplyr::select(-c(nobs, SE, lci, uci)) %>%
    mutate(
      estimate = format_estimate(estimate),
      `p-value` = format_pvalue(`p-value`)
    ) %>%
    pivot_longer(cols = c(estimate, `p-value`), names_to = "Metric") %>%
    pivot_wider(names_from = expo, values_from = value, id_cols = c(Annotation, Metric)) %>%
    relocate(Annotation)
}

# Main model and sensitivity analyses (WCRF, area, no education or smoking,
# excluding Barcelona)
metab.MCC.matched.MCC.res <- select_matched(MCC.MWAS.res)
metab.MCC.matched.MCC.WCRF.res <- select_matched(MCC.MWAS.WCRF.res)
metab.MCC.matched.MCC.area.res <- select_matched(MCC.MWAS.area.res)
metab.MCC.matched.MCC.noed.res <- select_matched(MCC.MWAS.noed.res)
metab.MCC.matched.MCC.noba.res <- select_matched(MCC.MWAS.noba.res)

MCC.res.tab <- tabulate_matched(metab.MCC.matched.MCC.res, format_estimate_main)
MCC.res.WCRF.tab <- tabulate_matched(metab.MCC.matched.MCC.WCRF.res, format_estimate_sensitivity)
MCC.res.area.tab <- tabulate_matched(metab.MCC.matched.MCC.area.res, format_estimate_sensitivity)
MCC.res.noed.tab <- tabulate_matched(metab.MCC.matched.MCC.noed.res, format_estimate_sensitivity)
MCC.res.noba.tab <- tabulate_matched(metab.MCC.matched.MCC.noba.res, format_estimate_sensitivity)

if (write.tab) {
  write.csv(MCC.res.tab, "Results/Tables/MCC.MWAS.res.csv", row.names = FALSE)
  write.csv(MCC.res.WCRF.tab, "Results/Tables/MCC.MWAS.res.WCRF.csv", row.names = FALSE)
  write.csv(MCC.res.area.tab, "Results/Tables/MCC.MWAS.res.area.csv", row.names = FALSE)
  write.csv(MCC.res.noed.tab, "Results/Tables/MCC.MWAS.res.noed.csv", row.names = FALSE)
}

# Exploratory plots: models with vs. without education and smoking --------
res.noed.vs.main <- metab.MCC.matched.MCC.noed.res %>%
  mutate(Model = "No education or smoking") %>%
  bind_rows(metab.MCC.matched.MCC.res %>% mutate(Model = "Adjusted for education and smoking"))

pval.noed.plot <- res.noed.vs.main %>%
  pivot_wider(names_from = Model, values_from = `p-value`, id_cols = c("metab", "expo", "Annotation")) %>%
  ggplot(aes(-log10(`No education or smoking`), -log10(`Adjusted for education and smoking`), col = Annotation)) +
  geom_point(size = 1.5) +
  theme_bw() +
  facet_wrap(~expo) +
  theme(legend.position = "top") +
  scale_colour_brewer(palette = "Dark2") +
  xlab("-log10(p-value, no education or smoking)") +
  ylab("-log10(p-value, adjusted for education and smoking)")

coef.noed.plot <- res.noed.vs.main %>%
  pivot_wider(names_from = Model, values_from = estimate, id_cols = c("metab", "expo", "Annotation")) %>%
  ggplot(aes(`No education or smoking`, `Adjusted for education and smoking`, col = Annotation)) +
  geom_point(size = 1.5) +
  theme_bw() +
  facet_wrap(~expo) +
  theme(legend.position = "top") +
  scale_colour_brewer(palette = "Dark2") +
  xlab("Estimate (no education or smoking)") +
  ylab("Estimate (adjusted for education and smoking)")

forest.noed.plot <- res.noed.vs.main %>%
  ggplot(aes(Annotation, estimate, col = Model)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymax = uci, ymin = lci), position = position_dodge(width = 0.5)) +
  facet_wrap(~expo, scales = "free_x") +
  theme_bw() +
  coord_flip() +
  theme(legend.position = "top")
