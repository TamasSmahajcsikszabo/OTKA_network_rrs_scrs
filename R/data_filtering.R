---
title: "Data Filtering: Removing redundant items"
output: html_document
---

```{r}
library(tidyverse)
library(networktools)
```
# sourcing the data import script to get the dataset as tibble
source("~/repos/rumination_lasso/R/data_import.R")

# 1. searching for highly correlated pairs of items [>.5]
var_names <- names(dataset)

#### RRS #####
RRS_items <- grep("RRS_[1,2,3,4,5,6,7,8,9,10]", var_names)
RRS_subset <- dataset[, RRS_items]

RRS_goldbricker <- goldbricker(RRS_subset, p = 0.05, method = "hittner2003", threshold = 0.25, corMin = 0.5, progressbar = TRUE)
RRS_goldbricker$suggested_reductions

# SCRS
SCRS_items <- grep("SCRS_[1,2,3,4,5,6,7,8,9,10]", var_names)
SCRS_subset <- dataset[, SCRS_items]

SCRS_goldbricker <- goldbricker(SCRS_subset, p = 0.05, method = "hittner2003", threshold = 0.25, corMin = 0.5, progressbar = TRUE)
SCRS_goldbricker$proportion_matrix
suggested_SCRS <- SCRS_goldbricker$suggested_reductions
SCRS_subset_reduced <- net_reduce(SCRS_subset, SCRS_goldbricker, method = c("PCA", "best_goldbricker"))

# MHC
MHC_items <- grep("MHC_[1,2,3,4,5,6,7,8,9,10]", var_names)
MHC_subset <- dataset[, MHC_items]

MHC_goldbricker <- goldbricker(MHC_subset, p = 0.05, method = "hittner2003", threshold = 0.25, corMin = 0.5, progressbar = TRUE)
MHC_goldbricker$proportion_matrix
MHC_goldbricker$suggested_reductions
MHC_subset_reduced <- net_reduce(MHC_subset, MHC_goldbricker, method = c("PCA", "best_goldbricker"))
