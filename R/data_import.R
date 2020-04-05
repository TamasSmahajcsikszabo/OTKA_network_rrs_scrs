# libraries
library(foreign)
library(tidyverse)

# data importing
dataset <- read.spss("~/repos/rumination_lasso/data/scrs.sav", to.data.frame = TRUE)
dataset <- as_tibble(dataset) # transform raw dataframe into tibble dataframe
