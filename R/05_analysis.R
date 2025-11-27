#Additional computations and model building

# ergm

library(statnet)
library(tidyverse)
library(ergm)

rumi<-readRDS("data/rumi_graph_full.RDS")
A <-igraph::as_adjacency_matrix(rumi)
v.attrs <- igraph::as_data_frame(rumi, what="vertices")
rumi.s <- network::as.network(as.matrix(A), directed=FALSE)
network::set.vertex.attribute(rumi.s, "Tool", ifelse(v.attrs$tool=="RRS",1,2))
network::set.vertex.attribute(rumi.s, "Scale", v.attrs$subscale_enum)
rumi.ergm <- formula(rumi.s ~ edges + gwdegree(3, fixed=TRUE) + gwesp(c(1,2,3), fixed=TRUE) + match("Tool") + match("Scale"))
ergm.fit<-ergm(rumi.ergm,estimate='MLE', control=control.ergm(MCMLE.maxit = 10))
saveRDS(ergm.fit, "output/ergmfit.RDS")
