library(tidyverse)
library(networktools)
library(qgraph)
library(glasso)
library(ggraph)
library(igraph)
library(bootnet)
library(modelr)
library(scales)
library(knitr)
library(ggpubr)
library(broom)
library(ergm)
library(fdrtool)
library(eigenmodel)
library(ROCR)
library(bookdown)
library(foreign)
library(tidyverse)
library(haven)
library(tibble)
library(ergm)
library(fdrtool)
library(eigenmodel)
library(ROCR)
library(bookdown)

if (!"captioner" %in% installed.packages()) {
  devtools::install_github("adletaw/captioner")
} else {
  library(captioner)
}

if (!"networkAnalysisTools" %in% installed.packages()) {
  devtools::install_github("TamasSmahajcsikszabo/networkAnalysisTools")
} else {
  library(networkAnalysisTools)
}