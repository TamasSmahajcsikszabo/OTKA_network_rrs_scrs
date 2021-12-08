library(tibble)
library(bnlearn)
library(stringr)
library(igraph)
library(sand)

source("~/repos/rumination_lasso/R/data_import.R")

selected_data <- dataset[,str_detect(names(dataset), "RRS") | str_detect(names(dataset), "SCRS")]
selected_data <- dataset[,str_detect(names(dataset), "RRS")]
selected_data <- selected_data[, 1:(ncol(selected_data)-4)]
selected_data <- na.omit(selected_data)

dagfit <- gs(selected_data)
dagfit2 <- iamb(selected_data)
dagfit3 <- hc(selected_data)

graphviz.plot(dagfit)
graphviz.plot(dagfit2)
graphviz.plot(dagfit3)

bnlearnObject <- dagfit3

create_vertex_dataframe <- function(bnlearnObject, ...) {
    nodes <- bnlearnObject$nodes
    network <- tibble()

    for (n in nodes) {
        subnetwork <- tibble()
        node  <- names(n)
        obj <- n[[1]]
        parents <- obj$parents
        children <- obj$children
        for (parent in parents) {
            actual_network <- tibble(V1 = parent, V2 = node)
            subnetwork <- bind_rows(subnetwork, actual_network)
        }
        for (child in children) {
            actual_network <- tibble(V1 = node, V2 = child)
            subnetwork <- bind_rows(subnetwork, actual_network)
        }
    }
    
}

?igraph::as_adjacency_matrix()
data(lazega)
names(elist.lazega)
