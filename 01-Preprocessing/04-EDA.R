# =============================================================================
# 04-EDA.R
#
# Descriptive analyses: baseline characteristics and THM exposure tables for
# PISCINA-II and MCC-Spain, exposure plots and correlation heatmaps of the
# exposures and denoised metabolomics data.
#
# Input:  Data/Cleaned_Data/*.rds
# Output: Results/Tables/*, Results/Figures/* (if write.tab / write.fig)
# =============================================================================

rm(list = ls())
setwd(Sys.getenv("PISCINA_DIR", unset = ".")) # project root (see README)

library(gt)
library(gtsummary)
library(dplyr)
library(ggplot2)
library(tidyr)
library(patchwork)
library(stringr)
library(pheatmap)
library(RColorBrewer)
library(grid)
library(gridExtra)


# Create folders ----------------------------------------------------------
dir.create("Results", showWarnings = FALSE)
dir.create("Results/Figures", showWarnings = FALSE)
dir.create("Results/Tables", showWarnings = FALSE)

# Read data ---------------------------------------------------------------
covars <- readRDS("Data/Cleaned_Data/covars.piscina.rds")
expo <- readRDS("Data/Cleaned_Data/expo.piscina.rds")
expo.composite <- readRDS("Data/Cleaned_Data/expo.piscina.composite.rds")
covars.expo.mcc <- readRDS("Data/Cleaned_Data/covars.expo.MCC.rds")

covars.mcc <- covars.expo.mcc %>% dplyr::select(ID:Ever_smoking)
expo.mcc <- covars.expo.mcc %>% dplyr::select(Chloroform:TotalTHM_18)


# Write tables/figures ----------------------------------------------------
write.tab <- FALSE
write.fig <- FALSE

# Table 1 -----------------------------------------------------------------
### PISCINA
PISCINA.Table1 <- covars %>%
  filter(Per == "post") %>%
  select(age, Weight, Height, BMI, sex, distance, ph, kcal) %>%
  rename(Age = age, Sex = sex, Distance = distance) %>%
  mutate(Sex = ifelse(Sex == "male", "Male", "Female")) %>%
  relocate(Sex) %>%
  tbl_summary(
    statistic = list(all_continuous() ~ "{mean} ({sd})"),
    missing = "no"
  ) %>%
  as_gt() %>%
  tab_style(style = list(cell_borders(
    color = "white",
    sides = c("bottom", "top", "left", "right"),
    weight = px(1)
  )), locations = cells_column_labels(everything())) %>%
  cols_align(align = "left", columns = everything())


PISCINA.expo.Table <- covars %>%
  select(Per, Cl3CHeb:Br3CHeb, ID) %>%
  mutate(
    Per = factor(ifelse(Per == "pre", "Before swimming", "After swimming"),
      levels = c("Before swimming", "After swimming")
    ),
    BrTHM = BDCMeb + DBCMeb + Br3CHeb,
    TotalTHM = Cl3CHeb + BrTHM
  ) %>%
  rename(
    `Chloroform (CHCl3)` = Cl3CHeb,
    `Bromodichloromethane (BDCM)` = BDCMeb,
    `Dibromochloromethane (DBCM)` = DBCMeb,
    `Bromoform (CHBr3)` = Br3CHeb,
    `Brominated THMs (BrTHM)` = BrTHM,
    `Total THMs (TTHM)` = TotalTHM
  ) %>%
  mutate(ID = as.factor(ID)) %>%
  relocate(ID) %>%
  tbl_summary(
    by = Per,
    type = all_continuous() ~ "continuous2",
    statistic = list(all_continuous() ~ c("{mean} ({sd})", "{min}, {max}")),
    include = `Chloroform (CHCl3)`:`Total THMs (TTHM)`,
    digits = c(`Total THMs (TTHM)`) ~ 1
  ) %>%
  add_p(
    group = ID,
    test = list(`Chloroform (CHCl3)`:`Total THMs (TTHM)` ~ "paired.t.test"),
    pvalue_fun = function(x) format(x, scientific = TRUE, digits = 3)
  ) %>%
  modify_header(label = "Disinfection By-Products (DBPs)") %>%
  as_gt() %>%
  tab_style(style = list(cell_borders(
    color = "white",
    sides = c("bottom", "top", "left", "right"),
    weight = px(1)
  )), locations = cells_column_labels(everything())) %>%
  cols_align(align = "left", columns = everything())


# Write outputs
if (write.tab) {
  gtsave(PISCINA.Table1, filename = "Results/Tables/PISCINA_baseline_characteristics.html")
  gtsave(PISCINA.Table1, filename = "Results/Tables/PISCINA_baseline_characteristics.tex")
  gtsave(PISCINA.expo.Table, filename = "Results/Tables/PISCINA_DBP_levels.tex")
}

### MCC-Spain
# baseline characteristics
MCC.Table1 <- covars.mcc %>%
  mutate(
    Sex = ifelse(Sex == "Hombre", "Male", "Female"),
    Smoking = ifelse(Ever_smoking == 1, "Ever smoker", "Never smoker")
  ) %>%
  select(Sex, Age, Area, Education, BMI, Smoking) %>%
  tbl_summary(
    statistic = list(all_continuous() ~ "{mean} ({sd})"),
    missing = "no"
  ) %>%
  as_gt() %>%
  tab_style(style = list(cell_borders(
    color = "white",
    sides = c("bottom", "top", "left", "right"),
    weight = px(1)
  )), locations = cells_column_labels(everything())) %>%
  cols_align(align = "left", columns = everything())

# Exposure
MCC.expo.Table <- expo.mcc %>%
  pivot_longer(
    cols = everything(),
    names_to = c(".value", "Exposure_type"),
    names_sep = "_"
  ) %>%
  mutate(Exposure_type = ifelse(is.na(Exposure_type), "Recent", "Long-term")) %>%
  rename(
    `Chloroform (CHCl3)` = Chloroform,
    `Bromodichloromethane (BDCM)` = BDCM,
    `Dibromochloromethane (DBCM)` = DBCM,
    `Bromoform (CHBr3)` = Bromoform,
    `Brominated THMs (BrTHM)` = BrTHM,
    `Total THMs (TTHM)` = TotalTHM
  ) %>%
  tbl_summary(
    by = Exposure_type,
    type = all_continuous() ~ "continuous2",
    statistic = list(all_continuous() ~ c("{mean} ({sd})", "{min}, {max}")),
    include = `Chloroform (CHCl3)`:`Total THMs (TTHM)`,
    digits = c(`Total THMs (TTHM)`) ~ 1
  ) %>%
  add_p(
    test = list(`Chloroform (CHCl3)`:`Total THMs (TTHM)` ~ "t.test"),
    pvalue_fun = function(x) format(x, scientific = TRUE, digits = 3)
  ) %>%
  modify_header(label = "Disinfection By-Products (DBPs)") %>%
  as_gt() %>%
  tab_style(style = list(cell_borders(
    color = "white",
    sides = c("bottom", "top", "left", "right"),
    weight = px(1)
  )), locations = cells_column_labels(everything())) %>%
  cols_align(align = "left", columns = everything())

# Write outputs
if (write.tab) {
  gtsave(MCC.Table1, filename = "Results/Tables/MCC_baseline_characteristics.html")
  gtsave(MCC.Table1, filename = "Results/Tables/MCC_baseline_characteristics.tex")
  gtsave(MCC.expo.Table, filename = "Results/Tables/MCC_exposure.tex")
}


# Visualise exposure levels -----------------------------------------------
expo.line <- expo.composite %>%
  mutate(
    Per = covars$Per,
    ID = covars$ID
  ) %>%
  pivot_longer(Chloroform:TotalTHM, names_to = "DBP") %>%
  mutate(
    DBP = factor(DBP, levels = c("Chloroform", "BDCM", "DBCM", "Bromoform", "BrTHM", "TotalTHM")),
    Per = factor(ifelse(Per == "pre", "Before", "After"),
      levels = c("Before", "After")
    )
  ) %>%
  ggplot(aes(Per, value)) +
  geom_point(size = 0.5) +
  geom_line(aes(group = ID)) +
  facet_wrap(~DBP, scales = "free") +
  theme_bw() +
  scale_colour_brewer(palette = "Set1") +
  xlab("") +
  ylab(expression(paste("Expression level (", mu, "g/L)")))

# MCC exposure
covars.mcc.long <- covars.expo.mcc %>%
  pivot_longer(Chloroform:TotalTHM_18, names_to = "DBP") %>%
  mutate(
    Exposure_type = ifelse(str_detect(DBP, "_18"), "Long-term", "Recent"),
    DBP = factor(str_remove(DBP, "_18"), levels = c("Chloroform", "BDCM", "DBCM", "Bromoform", "BrTHM", "TotalTHM"))
  )

p1 <- covars.mcc.long %>%
  filter(!(DBP %in% c("BrTHM", "TotalTHM"))) %>%
  filter(Exposure_type == "Recent") %>%
  ggplot(aes(DBP, value)) +
  geom_boxplot() +
  theme_bw() +
  scale_fill_brewer(palette = "Set1") +
  xlab("") +
  scale_y_continuous(limits = c(0, 60), breaks = seq(0, 60, 20)) +
  ylab(expression(paste("Concentration (", mu, "g/L)"))) +
  theme(text = element_text(size = 15), legend.position = "top")

p2 <- covars.mcc.long %>%
  filter(DBP == "BrTHM") %>%
  filter(Exposure_type == "Recent") %>%
  ggplot(aes(DBP, value)) +
  geom_boxplot() +
  theme_bw() +
  scale_fill_brewer(palette = "Set1") +
  xlab("") +
  ylab("") +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 25)) +
  theme(text = element_text(size = 15), legend.position = "none")

p3 <- covars.mcc.long %>%
  filter(DBP == "TotalTHM") %>%
  filter(Exposure_type == "Recent") %>%
  ggplot(aes(DBP, value)) +
  geom_boxplot() +
  theme_bw() +
  scale_fill_brewer(palette = "Set1", name = "Exposure type") +
  xlab("") +
  ylab("") +
  scale_y_continuous(limits = c(0, 120), breaks = seq(0, 120, 30)) +
  theme(text = element_text(size = 15), legend.position = "none")

MCC.expo.recent.plot <- p1 + p2 + p3 + plot_layout(widths = c(3, 1, 1))

# MCC exposure stratified by location
MCC.expo.by.area <- covars.mcc.long %>%
  ggplot(aes(DBP, value, fill = Area)) +
  geom_boxplot() +
  facet_grid(rows = vars(Exposure_type)) +
  theme_bw() +
  xlab("") +
  ylab(expression(paste("DBP levels (", mu, "g/L)")))

if (write.fig) {
  ggsave("Results/Figures/PISCINA_expo_lineplot.pdf", expo.line)
  ggsave("Results/Figures/MCC_exposure_recent_boxplot.pdf", MCC.expo.recent.plot, width = 8, height = 4)
  ggsave("Results/Figures/MCC_exposure_boxplot_by_area.pdf", MCC.expo.by.area)
}


# Correlation heatmaps ----------------------------------------------------
metab.denoised <- readRDS("Data/Cleaned_Data/metab.piscina.denoised.rds")

# Samples are ordered by subject, alternating before (odd rows) and after
# (even rows) swimming
expo.pre <- expo[seq(1, 115, 2), c("Chloroform", "BDCM", "DBCM", "Bromoform")]
expo.post <- expo[seq(2, 116, 2), c("Chloroform", "BDCM", "DBCM", "Bromoform")]

cor.expo.pre <- pheatmap(cor(expo.pre, method = "spearman"),
  cluster_cols = FALSE, cluster_rows = FALSE,
  col = colorRampPalette(brewer.pal(n = 9, name = "RdYlBu"))(100),
  breaks = seq(-1, 1, length.out = 100),
  main = "Pre swimming", fontsize = 4,
  display_numbers = TRUE,
  border_color = "white",
  cellwidth = 21,
  cellheight = 21,
  fontsize_number = 4
)

cor.expo.post <- pheatmap(cor(expo.post, method = "spearman"),
  cluster_cols = FALSE, cluster_rows = FALSE,
  col = colorRampPalette(brewer.pal(n = 9, name = "RdYlBu"))(100),
  breaks = seq(-1, 1, length.out = 100),
  main = "Post swimming", fontsize = 4,
  display_numbers = TRUE,
  border_color = "white",
  cellwidth = 21,
  cellheight = 21,
  fontsize_number = 4
)


cor.metab <- pheatmap(cor(metab.denoised),
  treeheight_row = FALSE, treeheight_col = FALSE, show_rownames = FALSE, show_colnames = FALSE,
  col = colorRampPalette(brewer.pal(n = 9, name = "RdYlBu"))(100),
  breaks = seq(-1, 1, length.out = 100),
  width = 8,
  main = "Metabolomics", fontsize = 4,
  height = 10,
  border_color = "white",
  cellwidth = 0.15,
  cellheight = 0.15
)

plot_list_expo <- list()
plot_list_expo[["p1"]] <- cor.expo.pre[[4]]
plot_list_expo[["p2"]] <- cor.expo.post[[4]]
g <- arrangeGrob(grobs = plot_list_expo, ncol = 2, nrow = 2)
if (write.fig) {
  ggsave(filename = "Results/Figures/cor_heatmap_expo.pdf", g, width = 4, height = 3)
}


plot_list_all <- list()
plot_list_all[["p1"]] <- cor.expo.pre[[4]]
plot_list_all[["p2"]] <- cor.expo.post[[4]]
plot_list_all[["p3"]] <- cor.metab[[4]]
layout <- matrix(c(1, 2, 3, 3), byrow = TRUE, nrow = 2)
g <- arrangeGrob(grobs = plot_list_all, ncol = 2, nrow = 2, heights = c(1, 2), layout_matrix = layout)

if (write.fig) {
  ggsave(filename = "Results/Figures/cor_heatmap_combined.pdf", g, width = 4, height = 6)
  ggsave(filename = "Results/Figures/cor_heatmap_combined.jpeg", g, width = 4, height = 6)
  ggsave(filename = "Results/Figures/cor_heatmap_combined.png", g, width = 4, height = 6)
}
