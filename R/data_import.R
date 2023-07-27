# libraries
library(foreign)
library(tidyverse)
library(haven)
library(igraph)
library(sand)
library(qgraph)
library(tibble)
library(tidyverse)

# data importing
dataset <- read.spss("data/scrs888.sav", to.data.frame = TRUE, use.value.labels = TRUE)
dataset <- as_tibble(dataset) # transform raw dataframe into tibble dataframe

# data labels
data_labels <- tibble(var = names(dataset), labels = attributes(dataset)$variable.labels) %>%
  filter(labels != "")


# node_labels <- readRDS("data/node_labels.RDS")
# node_labels <- node_labels[node_labels$var %in% c(paste0('RRS_',seq(1:10)), paste0('SCRS_',seq(1:10))),]
# rumi <- dataset[,c(13:32)]
# rumi <- na.omit(rumi)
# names(rumi)<-c(paste0('RRS_',seq(1:10)), paste0('SCRS_',seq(1:10)))
# cormatrix <-cor_auto(rumi, detectOrdinal=TRUE)
# set.seed(42)
# rumi_data <- data.frame(qgraph::EBICglasso(cormatrix, n = nrow(rumi), threshold =TRUE,nlambda=1000))
# mycorr <- as.matrix(rumi_data)
# rumi_data[lower.tri(rumi_data, diag=TRUE)] <- NA
# saveRDS(rumi_data, "output/mycorr.RDS")
# rumi_data <- tibble(data.frame(rumi_data))
# rumi_data <- bind_cols(name = colnames(rumi_data), rumi_data)
# rumi_graph <- graph_from_data_frame(rumi_data[,1:2], directed=FALSE)
