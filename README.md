# Multi-omic signatures of exposure to water disinfection by-products

Analysis code for:

> Cheng SL, Keski-Rahkonen P, Chatziioannou AC, *et al.* (2026). **Exploring molecular signatures of exposure to water disinfection by-products: A multi-omic analysis across multiple exposure routes.** *Environment International*, 215, 110442. [doi:10.1016/j.envint.2026.110442](https://doi.org/10.1016/j.envint.2026.110442)

## Overview

Trihalomethanes (THMs) are disinfection by-products (DBPs) formed when chlorine reacts with organic matter in water, and they have been linked to bladder and colorectal cancer. This project looks for molecular signatures of THM exposure in two studies with very different exposure routes and levels:

- **PISCINA-II (discovery):** 58 healthy adults swam for 40 minutes in a chlorinated indoor pool. Serum metabolomics (UHPLC-QTOF-MS), serum immune markers and THMs in exhaled breath were measured before and after swimming.
- **MCC-Spain (validation):** 293 cancer-free controls from a multicentre case-control study, with residential tap-water THM exposure estimates and serum metabolomics.

The analysis combines univariate and multivariate methods:

1. **Pre-processing:** missingness filtering, QRILC imputation of left-censored values, and removal of plate effects with linear mixed models.
2. **Metabolome-wide association study (MWAS):** multivariate normal models for repeated measures (generalised least squares with an unstructured within-subject correlation), one model per feature and THM. Features are significant after Bonferroni correction.
3. **Stability-calibrated consensus clustering** ([`sharp`](https://cran.r-project.org/package=sharp)) of correlated features, which typically come from adducts and fragments of the same compound. Each cluster is summarised by its first principal component.
4. **Multilevel sparse PLS** ([`mixOmics`](https://mixomics.org)) with the four THMs as predictors and the metabolic clusters as outcomes. The model is calibrated on cross-validated Q² and on out-of-sample correlation between components.
5. **Stability-calibrated graphical LASSO networks** linking the selected metabolic clusters and immune markers, estimated before and after swimming.
6. **Validation in MCC-Spain:** features are matched across studies (≤ 7 ppm, ± 0.1 min) and the annotated compounds are tested with linear regression.

## Requirements

The analyses were run with R 4.1.3. Required packages:

- **CRAN:** `dplyr`, `tidyr`, `stringr`, `tibble`, `ggplot2`, `patchwork`, `RColorBrewer`, `pheatmap`, `gridExtra`, `gt`, `gtsummary`, `ComplexUpset`, `readxl`, `foreign`, `imputeLCMD`, `lme4`, `nlme`, `doBy`, `FactoMineR`, `sharp`, `igraph`, `circlize`
- **Bioconductor:** `mixOmics`, `ComplexHeatmap`; `pcaMethods` and `impute` are needed by `imputeLCMD`

```r
install.packages("BiocManager")
BiocManager::install(c("mixOmics", "ComplexHeatmap", "pcaMethods", "impute"))
install.packages(c(
  "dplyr", "tidyr", "stringr", "tibble", "ggplot2", "patchwork", "RColorBrewer",
  "pheatmap", "gridExtra", "gt", "gtsummary", "ComplexUpset", "readxl", "foreign",
  "imputeLCMD", "lme4", "nlme", "doBy", "FactoMineR", "sharp", "igraph", "circlize"
))
```

## Citation

If you use this code, please cite the paper:

```bibtex
@article{cheng2026dbp,
  title   = {Exploring molecular signatures of exposure to water disinfection by-products:
             A multi-omic analysis across multiple exposure routes},
  author  = {Cheng, Sibo Lucas and Keski-Rahkonen, Pekka and Chatziioannou, Anastasia Chrysovalantou
             and Vuckovic, Dragana and Barron, Leon P. and Dehghan, Abbas and Vlaanderen, Jelle
             and Portengen, L{\"u}tzen and Jim{\'e}nez Zabala, Ana and Guevara, Marcela
             and Gomez Acebo, Ines and Molina de la Torre, Antonio Jos{\'e} and Scalbert, Augustin
             and Vineis, Paolo and Elliott, Paul and Vermeulen, Roel and Kogevinas, Manolis
             and Villanueva, Cristina M. and Dagnino, Sonia and Chadeau-Hyam, Marc},
  journal = {Environment International},
  volume  = {215},
  pages   = {110442},
  year    = {2026},
  doi     = {10.1016/j.envint.2026.110442}
}
```

## Author

Sibo Lucas Cheng ([ORCID 0009-0009-8930-6148](https://orcid.org/0009-0009-8930-6148)), Department of Epidemiology and Biostatistics, School of Public Health, Imperial College London.
