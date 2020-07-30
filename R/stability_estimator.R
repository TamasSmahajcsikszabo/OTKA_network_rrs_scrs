
library(tidyverse)
library(networktools)
library(qgraph)
library(glasso)
#library(ggraph)
library(igraph)
library(bootnet)
library(modelr)
library(scales)

setwd("/home/tamas/repos/rumination_lasso")
RRS_stability <- readRDS("output/RRS_stability.RDS")
SCRS_stability <- readRDS("output/SCRS_stability.RDS")
full_stability <- readRDS("output/full_stability.RDS")

network_stability_estimator <- function(network, name) {
  B_levels = seq(1, 2500)
  sample_data  <- t(matrix(network$sampleTable))
  names(sample_data) <- names(network$sampleTable)
  boot <- network$bootTable  %>%
    mutate(sample_rate = nPerson / mean(sample_data$nPerson)) %>%
    mutate(
      row = row_number(),
      bootstrap_indicator = (row %/% 75) + 1,
      bootstrap_indicator = factor(bootstrap_indicator, levels = B_levels))
  
  boot_data <- t(matrix(boot))
  names(boot_data) <- names(boot)
  
  drop_levels  <- unique(boot_data$nPerson)
  measures<- unique(sample_data$type)
  total_operations  <- 2500 * length(measures) * length(drop_levels)
  plot_data  <- matrix(ncol = 4, nrow = 0)
  colnames(plot_data) <- c("b", "measure", "level", "correlation")
  i = 0
  
  for (b in B_levels) {
    for (level in drop_levels) {
      for (m in measures) {
        i = i + 1
        
        boot_sample_bootstrap_mask  <- unlist(boot_data['bootstrap_indicator'], use.names=FALSE) == b
        boot_sample_level_mask  <- unlist(boot_data['nPerson'], use.names=FALSE) == level
        boot_sample_measure_mask  <- unlist(boot_data['type'], use.names=FALSE) == m
        boot_mask  <- boot_sample_bootstrap_mask & boot_sample_level_mask & boot_sample_measure_mask
        
        
        boot_sample <- unlist(boot_data['value'], use.names=FALSE)[boot_mask]
        
        
        original_sample_mask <- unlist(sample_data['type'], use.names=FALSE) == m
        original_sample <- unlist(sample_data['value'], use.names = FALSE)[original_sample_mask]
        
        test  <-all(!is.null(boot_sample), !is.null(original_sample), length(original_sample) == length(boot_sample)) 
        if (test) {
          correlation = cor(boot_sample, original_sample)
          new_record <- c(
            "b" = b,
            "measure" = m,
            "level" = level,
            "correlation" = correlation
          )
          plot_data <- rbind(plot_data, new_record)
          cat(paste0("\rProgress: ", round((i/total_operations)*100,2), "%"))
        }
      }
    }
  }
  plot_df  <- tibble(as.data.frame(plot_data))
  plot_df  <- plot_df %>%
    mutate(level = as.numeric(level),
           correlation = as.numeric(correlation)) %>%
    group_by(b, measure) %>%
    summarise(r = mean(correlation, na.rm = TRUE),
              sample_rate = level / max(unlist(sample_data['nPerson'], use.names = FALSE)))
  
  saveRDS(plot_df, paste0("output/plot_", name, ".RDS"))
  
}

#network_stability_estimator(RRS_stability, "RRS_stability")
#network_stability_estimator(SCRS_stability, "SCRS_stability")
network_stability_estimator(full_stability, "combined_stability")
