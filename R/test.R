library(qgraph)

AM <- matrix(0, 10, 10)
AM[1, 2] <- AM[2, 1] <- AM[2, 3] <- AM[3, 2] <- AM[2, 4] <- AM[4, 2] <- AM[3, 4] <- AM[4, 3] <- AM[3, 6] <- AM[6, 3] <- AM[3, 8] <- AM[8, 3] <- AM[3, 9] <- AM[9, 3] <- AM[3, 10] <- AM[10, 3] <- AM[4, 7] <- AM[7, 4] <- AM[5, 7] <- AM[7, 5] <- AM[5, 10] <- AM[10, 5] <- AM[9, 10] <- AM[10, 9] <- 1
gr <- list(c(1, 2, 4:10), 3)
names <- c("1", "3", "6", "3", "2", "1", "2", "1", "2", "3")
N <- qgraph(AM,
  groups = gr, color = c("#cccccc", "#3CB371"), labels = names,
  border.width = 3, edge.width = 2, vsize = 9,
  border.color = "#555555", edge.color = "#555555", label.color = "#555555"
)
