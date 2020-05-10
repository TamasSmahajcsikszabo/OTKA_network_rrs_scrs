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


unnested <- test %>%
  unnest()  %>%
  mutate(fill_flag = if_else(significant, "sig", if_else(!significant, "nonsig", if_else(lower == 0 & upper == 0, "itself"))))

ggplot(unnested) +
  geom_tile(aes(x = id1, y = id2, fill = significant)) +
  scale_fill_manual(values = c("gray70", "coral")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
