---
title: "Network Analysis on Rumination"
output: html_document
---

```{r include = FALSE}
# libraries declared
library(tidyverse)
library(networktools)
library(qgraph)
library(glasso)
library(ggraph)
library(igraph)

```
## 2. Methods

### 2.1 Network analysis 

#### 2.1.1 Removing redundant items

As a preparation step, we applied the *goldbricker* function from the R package *networktools* to search for pairs of items in the RRS, SCRS and MHS items which were highly correlated (r > 0.50) and showed highly the same pattern with other items. When such "bad pairs" were identified, reduction techniques were put into use (with the *net_reduce* function in R): instead of the pair, the first principal component of the two variables were kept in the data as a new, merged variable.


```{r include = FALSE}

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
MHC_items <- grep("MHC_[1,2,3,4,5,6,7,8,9,10,11,12,14]", var_names)
MHC_subset <- dataset[, MHC_items]

MHC_goldbricker <- goldbricker(MHC_subset, p = 0.05, method = "hittner2003", threshold = 0.25, corMin = 0.5, progressbar = TRUE)
MHC_goldbricker$proportion_matrix
MHC_goldbricker$suggested_reductions
MHC_subset_reduced <- net_reduce(MHC_subset, MHC_goldbricker, method = c("PCA", "best_goldbricker"))
```

The RRS items required no reduction steps as the highly correlating items did not extend their pattern to more nodes (their correlations to other nodes did not violate the criteria that > 75% of these correlations did not differ significantly).
Among the SCRS scales the `r names(suggested_SCRS)` scales required reductions steps using PCA. Among the MHC nodes `r names(MHC_goldbricker$suggested_reductions)` items were subject to reduction.

#### 2.1.2 Network Analysis with LASSO

For the network analysis, the tuning parameter of *gamma* in the grahical LASSO alogrithm has been set to 0.5 value in accordance with Bernstein et al in order to achieve sparse graphs. The final networks were chosen by the lowest values of the *extended Bayesian Information Criterion* (EBIC).

#### 2.1.3. Community detection
In order to get an understanding of community of nodes and to see if nodes can be grouped into neighbourhoods, we used the *igraph* R package. When detecting possible communities of nodes, the detection happened with the following parameters in the *spinglass algorithm*:
TODO: try the walktrap algorithm from igraph!

## 3. Results

### 3.1 Network plots
RSS network

```{r echo = FALSE, message = FALSE, warning = FALSE}

# RRS data
# orevious attempt
##RSScov <- cov(RRS_subset, use = "pairwise.complete") 
#RRSgLASSO <- glasso(RSScov, rho = 2^-2) 
#RSSnetwork <-  qgraph(RRSgLASSO, theme="TeamFortress")
#RSSnetwork$graphAttributes$Nodes$labels <- names(dataset[,RRS_items]) 

RSScov2 <- qgraph::cor_auto(RRS_subset) 
rss_network <- qgraph(RSScov2, graph="glasso", tuning=0.5, layout="spring", sampleSize=530, theme="TeamFortress", details = TRUE, threshold = TRUE)
```


SCRS network
```{r echo = FALSE, message = FALSE, warning = FALSE}


SCRScov <- qgraph::cor_auto(SCRS_subset_reduced) 
qgraph(SCRScov, graph="glasso", tuning=0.5, layout="spring", sampleSize=530, theme="TeamFortress", details = TRUE, threshold = TRUE)
```

MHC network
```{r echo = FALSE, message = FALSE, warning = FALSE}

# MHC data

MHCcov <- qgraph::cor_auto(MHC_subset_reduced) 
qgraph(MHCcov, graph="glasso", tuning=0.5, layout="spring", sampleSize=530, theme="TeamFortress", details = TRUE, threshold = TRUE)
```
 
### 3.2 Community detection 

```{r echo = FALSE, message = FALSE, warning = FALSE}
# community detection for RRS
rss_graph <-  as.igraph(rss_network, attributes=TRUE)
#spinglass_rss <- spinglass.community(rss_graph) 
#rss_communities <- spinglass_rss$membership
rss_groups <-  cluster_spinglass(
                  rss_graph,
                  weights = NULL,
                  vertex = NULL,
                  spins = 25,
                  parupdate = FALSE,
                  start.temp = 1,
                  stop.temp = 0.01,
                  cool.fact = 0.99,
                  update.rule = c("config", "random", "simple"),
                  gamma = 0.5,
                  implementation = c("orig", "neg"),
                  gamma.minus = 1
)

rss_graph %>%
    geom_node_circle(aes(fill = rss_groups))

```
