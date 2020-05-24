library(bootnet)
library(tidyverse)

source("~/repos/rumination_lasso/R/data_import.R")

var_names <- names(dataset)

#### RRS #####
RRS_items <- grep("RRS_[1,2,3,4,5,6,7,8,9,10]", var_names)
RRS_subset <- dataset[, RRS_items]

RRS_network <- estimateNetwork(RRS_subset, default = "EBICglasso")
RRS_edge_weights <- bootnet(
  data = RRS_network,
  type = "nonparametric",
  nCore = 4
)

plot(RRS_edge_weights,
  decreasing = TRUE,
  sampleColor = "cornflowerblue"
)


# bootstrap difference test
# id1 <- RRS_edge_weights$bootTable$node1
# id2 <- RRS_edge_weights$bootTable$node2
ids <- RRS_edge_weights$sampleTable$id
ids <- ids[grep("--", ids)]

matrix_values <- tibble(expand.grid(ids, ids))

# setup session for parallel computing
library(foreach)
library(doParallel)
cores <- detectCores()
cl <- makeCluster(cores[1] - 1)
registerDoParallel(cl)


test <- tibble()

for (i in seq(1, nrow(matrix_values))) {
  diff_i <- differenceTest(RRS_edge_weights,
    x = matrix_values[i, 1][[1]],
    y = matrix_values[i, 2][[1]],
    "edge",
    verbose = FALSE
  )
  test[i, "id"] <- i
  nested_result <- nest(diff_i)
  test[i, "model"] <- nested_result
  cat(paste0("\r", "Progress: ", round((i / nrow(matrix_values)) * 100, 0), "%"))
}
stopCluster()

saveRDS(test, "test.RDS")

unnested <- test %>%
  unnest() %>%
  mutate(fill_flag = if_else(significant & !lower == 0, "sig", if_else(!significant & !lower == 0, "nonsig", if_else(lower == 0 & upper == 0 & !significant, "itself", "nonsig"))))


unnested %>%
  filter(fill_flag == "valid") %>%
  filter(significant)
# order by number of significant
level_order <- unnested %>%
  group_by(id1) %>%
  summarize(sig = sum(significant)) %>%
  arrange(desc(sig))

level_order <- as.character(level_order$id1)

ggplot(unnested) +
  geom_tile(aes(x = factor(id1, level = level_order), y = factor(id2, level = level_order), fill = fill_flag, linejoin = "round", width = 0.9, height = 0.9)) +
  geom_text(aes(id1, id2, color = fill_flag), label = "*", vjust = 0.75, size = 10, fontface = "bold") +
  scale_fill_manual(values = c("white", "gray70", "grey70")) +
  scale_color_manual(values = c("white", "gray70", "grey20")) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(x = "Edges", y = "Edge pairs", title = "Difference test significances")
