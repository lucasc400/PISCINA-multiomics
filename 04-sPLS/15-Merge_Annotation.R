# =============================================================================
# 15-Merge_Annotation.R
#
# Merge the annotation of the sPLS-selected metabolic features with their
# cluster membership, and label each cluster by an annotated compound (or by
# its cluster ID when no feature in the cluster could be annotated).
#
# Input:  Annotation/Annotation_011223.csv
#         Results/Clustering/Description_summary.rds
# Output: Results/Annotation/annotation.sPLS.selected.metab.rds
#         Results/Annotation/cluster.annotation.rds
# =============================================================================

rm(list = ls())
setwd(Sys.getenv("PISCINA_DIR", unset = ".")) # project root (see README)

library(dplyr)

dir.create("Results/Annotation", showWarnings = FALSE, recursive = TRUE)

annotation <- read.csv("Annotation/Annotation_011223.csv")
cluster.description <- readRDS("Results/Clustering/Description_summary.rds")

annotation <- annotation %>%
  left_join(., cluster.description %>% dplyr::select(feature, group)) %>%
  arrange(group) %>%
  mutate(group = paste0("C", group)) %>%
  group_by(group) %>%
  mutate(Cluster_annotated = ifelse(all(grepl("Unknown", Annotation)), 0, 1)) %>%
  ungroup() %>%
  mutate(Annotation_group = ifelse(grepl("Unknown", Annotation), group, Annotation))

# Represent each cluster with a single compound where available, or cluster ID
cluster.annotation.df <- annotation %>%
  filter(Cluster_annotated == 0 | (Cluster_annotated == 1 & !grepl("Unknown", Annotation))) %>%
  filter(!duplicated(group))

cluster.annotation <- structure(cluster.annotation.df$Annotation_group, names = cluster.annotation.df$group)

# Clean annotation names
cluster.annotation <- gsub("\\*", "", cluster.annotation)
cluster.annotation <- sapply(cluster.annotation, function(x) {
  gsub(" \\(.*|\\*.*", "", x)
})
cluster.annotation["C148"] <- "Chlorinated compound"

saveRDS(annotation, "Results/Annotation/annotation.sPLS.selected.metab.rds")
saveRDS(cluster.annotation, "Results/Annotation/cluster.annotation.rds")
