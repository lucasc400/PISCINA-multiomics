# =============================================================================
# 08-Analysis_univariate.R
#
# Summarise the PISCINA-II MWAS results: features significant after
# Bonferroni correction, Manhattan, volcano and UpSet plots, and the
# equivalent summary for the original metabolomics dataset.
#
# Input:  Data/Cleaned_Data/metab.piscina.denoised.rds
#         Results/MVN/metab.piscina.MVN.res.rds
#         Results/MVN/metab.piscina.MVN.res.old.rds
# Output: Results/MVN/*, Results/Figures/* (if write.res / write.tab / write.fig)
# =============================================================================

rm(list = ls())
setwd(Sys.getenv("PISCINA_DIR", unset = ".")) # project root (see README)

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(ComplexUpset)
library(RColorBrewer)
source("Scripts/functions.R")

metab.PISCINA <- readRDS("Data/Cleaned_Data/metab.piscina.denoised.rds")
metab.PISCINA.res.list <- readRDS("Results/MVN/metab.piscina.MVN.res.rds")

write.res <- FALSE
write.tab <- FALSE
write.fig <- FALSE

# Bonferroni threshold
n.features <- 1696 # number of metabolic features tested
bonferroni.thr <- 0.05 / n.features

# Stack the per-exposure results into a single data frame
stack_results <- function(res.list) {
  res <- data.frame()
  for (dbp in names(res.list)) {
    tmp <- as.data.frame(res.list[[dbp]])
    tmp$Expo <- rep(dbp, nrow(tmp))
    res <- rbind(res, tmp)
  }
  res
}

# Exposure labels used in the facets (parsed as plotmath)
label_expo <- function(expo) {
  factor(
    case_when(
      expo == "BrTHM" ~ "Brominated~THM",
      expo == "TotalTHM" ~ "Total~THM",
      TRUE ~ expo
    ),
    levels = c("Chloroform", "BDCM", "DBCM", "Bromoform", "Brominated~THM", "Total~THM")
  )
}


# Results on the reprocessed dataset --------------------------------------
metab.PISCINA.res <- stack_results(metab.PISCINA.res.list)

metab.PISCINA.sign.list <- unique(metab.PISCINA.res$OMIC_ID[metab.PISCINA.res$pval < bonferroni.thr])

# Number of significant features per exposure
metab.PISCINA.res.df <- metab.PISCINA.res %>%
  dplyr::select(OMIC_ID, pval, Expo) %>%
  pivot_wider(id_cols = OMIC_ID, values_from = pval, names_from = Expo) %>%
  mutate(Significant = ifelse(OMIC_ID %in% metab.PISCINA.sign.list, 1, 0)) %>%
  arrange(desc(Significant))

if (write.res) {
  saveRDS(metab.PISCINA.res, "Results/MVN/metab.piscina.cleaned.res.rds")
  saveRDS(metab.PISCINA.sign.list, "Results/MVN/sign.metab.piscina.188.rds")
}

if (write.tab) {
  write.csv(metab.PISCINA.sign.list, "Results/MVN/PISCINA_sign_metab_list.csv")
  write.csv(metab.PISCINA.res.df, "Results/MVN/PISCINA_metab_1696_pvalues.csv")
}

# Figures -----------------------------------------------------------------
# Number of exposures each feature is significantly associated with
metab.PISCINA.res.plot <- metab.PISCINA.res %>%
  mutate(bf.sign = ifelse(pval < bonferroni.thr, "Bonferroni significant", "Insignificant")) %>%
  group_by(OMIC_ID) %>%
  mutate(n.significant = sum(bf.sign == "Bonferroni significant")) %>%
  ungroup() %>%
  mutate(Expo = label_expo(Expo))

# Manhattan plot
Manhattan.plot <- metab.PISCINA.res.plot %>%
  mutate(retention.time = as.numeric(str_split_fixed(str_split_fixed(OMIC_ID, ":", 2)[, 1], "@", 2)[, 2])) %>%
  ggplot(aes(retention.time, -log10(pval) * sign(coef), col = as.factor(n.significant))) +
  geom_point(size = 0.5) +
  facet_wrap(~Expo, ncol = 2, labeller = labeller(Expo = label_parsed)) +
  scale_colour_manual(
    values = c("grey90", brewer.pal(6, "Paired")[c(1, 3, 5, 2, 4, 6)]),
    name = "Number of significantly associated DBP"
  ) +
  scale_y_continuous(n.breaks = 9) +
  scale_x_continuous(n.breaks = 7, limits = c(0.4, 10)) +
  geom_hline(aes(yintercept = 0), linetype = "dashed", linewidth = 0.3, col = "black") +
  theme_bw() +
  theme(
    legend.position = "top",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  xlab("Retention time (min)") +
  ylab(expression(-log[10] * (italic(p) * "-value") * " * sign" * (italic(beta)))) +
  guides(colour = guide_legend(override.aes = list(size = 2)))

# Volcano plot
Volcano.plot <- metab.PISCINA.res.plot %>%
  ggplot(aes(coef, -log10(pval), col = as.factor(n.significant))) +
  geom_point(size = 0.5) +
  facet_wrap(~Expo, ncol = 2, labeller = labeller(Expo = label_parsed), scales = "free_x") +
  scale_colour_manual(
    values = c("grey90", brewer.pal(6, "Paired")[c(1, 3, 5, 2, 4, 6)]),
    name = "Number of significantly associated DBP"
  ) +
  theme_bw() +
  theme(legend.position = "top") +
  xlab("Coefficient") +
  ylab(expression(-log[10] * (italic(p) * "-value"))) +
  guides(colour = guide_legend(override.aes = list(size = 2))) +
  geom_hline(yintercept = -log10(bonferroni.thr), linetype = "dashed", col = "grey")

# UpSet plot
metab.PISCINA.upset.dat <- metab.PISCINA.res %>%
  dplyr::select(-c(coef, coef.se, tval)) %>%
  pivot_wider(names_from = Expo, values_from = pval, id_cols = OMIC_ID) %>%
  filter(OMIC_ID %in% metab.PISCINA.sign.list) %>%
  dplyr::select(-OMIC_ID) %>%
  mutate(across(everything(), ~ .x < bonferroni.thr))
colnames(metab.PISCINA.upset.dat) <- c("Chloroform", "BDCM", "DBCM", "Bromoform", "TTHM", "BrTHM")

Upset.plot <- upset(metab.PISCINA.upset.dat, colnames(metab.PISCINA.upset.dat),
  base_annotations = list(
    "Intersection size" = intersection_size(
      text = list(size = 3)
    )
  ),
  set_sizes = (
    upset_set_size() + geom_text(aes(label = after_stat(count)), hjust = 0, stat = "count", col = "white", size = 4)),
  themes = upset_default_themes(text = element_text(size = 15))
)

if (write.fig) {
  ggsave("Results/Figures/PISCINA_Manhattan.pdf", Manhattan.plot, width = 8, height = 8)
  ggsave("Results/Figures/PISCINA_Volcano.pdf", Volcano.plot, width = 8, height = 6)
  ggsave("Results/Figures/PISCINA_metab_upset.pdf", Upset.plot, width = 12, height = 6)
}

# Results on the original dataset -----------------------------------------
metab.PISCINA.old.list <- readRDS("Results/MVN/metab.piscina.MVN.res.old.rds")
metab.PISCINA.old.res <- stack_results(metab.PISCINA.old.list)

n.features.old <- 6471 # number of features in the original dataset
metab.PISCINA.sign.old.list <- unique(metab.PISCINA.old.res$OMIC_ID[metab.PISCINA.old.res$pval < 0.05 / n.features.old])

if (write.res) {
  saveRDS(metab.PISCINA.old.res, "Results/MVN/metab.piscina.cleaned.res.old.rds")
  saveRDS(metab.PISCINA.sign.old.list, "Results/MVN/sign.metab.piscina.old.308.rds")
}
