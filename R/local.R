library(tidyverse)

sample_data <- readRDS("../output/sample_data.RDS")
bootstrap_aggr <- readRDS("../output/bootstrap_data_aggr.RDS")

# greater than 0:
data <- bootstrap_aggr %>%
  filter(lower > 0)

write_csv(data, "../output/accuracy.csv")
