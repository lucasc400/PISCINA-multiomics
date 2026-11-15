# =============================================================================
# 05-Reprocessing.R
#
# Compare the original (n = 6,471 features) and reprocessed (n = 1,696
# features) PISCINA-II metabolomics datasets: distribution of features by
# retention time and monoisotopic mass.
#
# Input:  Data/Cleaned_Data/metab.piscina.rds
#         Data/Cleaned_Data/metab.piscina.old.rds
# Output: Results/Figures/Reprocessing.pdf
# =============================================================================

rm(list = ls(all = TRUE))
setwd(Sys.getenv("PISCINA_DIR", unset = ".")) # project root (see README)

library(dplyr)
library(ggplot2)
library(stringr)
library(patchwork)

metab.piscina <- readRDS("Data/Cleaned_Data/metab.piscina.rds")
metab.piscina.old <- readRDS("Data/Cleaned_Data/metab.piscina.old.rds")

# Feature IDs are formatted as "mass@retention_time"
mass.old <- as.numeric(str_split_fixed(colnames(metab.piscina.old), "@", 2)[, 1])
rt.old <- as.numeric(str_split_fixed(str_split_fixed(colnames(metab.piscina.old), ":", 2)[, 1], "@", 2)[, 2])
mass.new <- as.numeric(str_split_fixed(colnames(metab.piscina), "@", 2)[, 1])
rt.new <- as.numeric(str_split_fixed(str_split_fixed(colnames(metab.piscina), ":", 2)[, 1], "@", 2)[, 2])

mass.rt.old <- data.frame(omic = colnames(metab.piscina.old), mass = mass.old, rt = rt.old, dataset = "Original")
mass.rt.new <- data.frame(omic = colnames(metab.piscina), mass = mass.new, rt = rt.new, dataset = "New")

mass.rt.df <- rbind(mass.rt.old, mass.rt.new)
mass.rt.df$dataset <- factor(mass.rt.df$dataset, levels = c("Original", "New"))

# Plots -------------------------------------------------------------------
scatter <- mass.rt.df %>%
  mutate(dataset = ifelse(dataset == "Original", "Original (N=6,471)", "New (N=1,696)")) %>%
  mutate(dataset = factor(dataset, levels = c("Original (N=6,471)", "New (N=1,696)"))) %>%
  ggplot(aes(rt, mass, col = dataset)) +
  geom_point(size = 0.06, alpha = 0.6) +
  facet_wrap(~dataset) +
  theme_bw() +
  scale_colour_brewer(palette = "Pastel1") +
  theme(legend.position = "none") +
  xlab("Retention time (min)") +
  ylab("Mass")

histo.mass <- ggplot(aes(mass, fill = dataset), data = mass.rt.df) +
  geom_histogram(col = "grey50") +
  facet_wrap(~dataset) +
  theme_bw() +
  theme(
    strip.text = element_blank(),
    strip.background = element_blank(),
    legend.position = "none"
  ) +
  scale_fill_brewer(palette = "Pastel1") +
  xlab("Mass") +
  ylab("Count")

histo.rt <- ggplot(aes(rt, fill = dataset), data = mass.rt.df) +
  geom_histogram(col = "grey50") +
  facet_wrap(~dataset) +
  theme_bw() +
  theme(
    strip.text = element_blank(),
    strip.background = element_blank(),
    legend.position = "none"
  ) +
  scale_fill_brewer(palette = "Pastel1") +
  xlab("Retention time (min)") +
  ylab("Count")

reprocessing.plot <- scatter / histo.mass / histo.rt + plot_layout(heights = c(2, 1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(
    plot.tag = element_text(size = 20),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

dir.create("Results/Figures", showWarnings = FALSE, recursive = TRUE)
ggsave("Results/Figures/Reprocessing.pdf", reprocessing.plot, height = 6, width = 8)
