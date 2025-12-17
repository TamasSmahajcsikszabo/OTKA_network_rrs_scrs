library(knitr)
library(ltm)
library(tibble)
library(WRS)
library(dplyr)
library(tidyr)
library(rogme)
library(magrittr)
library(stringr)
library(igraph)
library(qgraph)
library(tidyverse)

# import main data
# original data with demographics
dataset <- readRDS("data/original_data.RDS")
data_labels <- readRDS("data/data_labels.RDS")

# data for analysis
dataset_only_scales <- dataset[, names(dataset) %in% c(paste0("RRS_", 1:10), paste0("SCRS_", 1:10))]


# prepared data for analysis
var_names <- names(dataset)

### create networks ###

# 1. RRS Graph
source_data_rrs <- dataset_only_scales[, paste0("RRS_", 1:10)]
RRS_graph_scale_id <- c(1, 2, 1, 2, 2, 1, 1, 1, 2, 2)
RRS_graph <- construct_graph(source_data_rrs, "RRS_graph_2025.RDS", RRS_graph_scale_id, data_labels)

# 2. SCRS Graph
source_data_scrs <- dataset_only_scales[, paste0("SCRS_", 1:10)]
SCRS_graph_scale_id <- rep(3, 10)
SCRS_graph <- construct_graph(source_data_scrs, "SCRS_graph_2025.RDS", SCRS_graph_scale_id, data_labels)

# 3. Full Combined Graph of RRS + SCRS
full_graph_scale_id <- c(1, 2, 1, 2, 2, 1, 1, 1, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 2, 3)
rumi_graph <- construct_graph(dataset_only_scales, "rumi_graph_full_2025.RDS", full_graph_scale_id, data_labels)

# load computed measures
RRS_simulated_communities <- readRDS("output/RRS_communities.RDS")
SCRS_simulated_communities <- readRDS("output/SCRS_communities.RDS")
combined_simulated_communities <- readRDS("output/combined_communities.RDS")
RRS_item_TDR <- readRDS("output/RRS_item_TPR.RDS")
SCRS_item_TDR <- readRDS("output/SCRS_item_TPR.RDS")
combined_item_TDR <- readRDS("output/combined_item_TPR.RDS")
accuracy_rrs <- readRDS("output/accuracy_rrs.RDS")
accuracy_scrs <- readRDS("output/accuracy_scrs.RDS")
accuracy_combined <- readRDS("output/accuracy_combined.RDS")
accuracy_rrs <- network_accuracy_data(accuracy_rrs)
accuracy_scrs <- network_accuracy_data(accuracy_scrs)
accuracy_combined <- network_accuracy_data(accuracy_combined)

# add accuracy to networks
rumi <- add_accuracy_to_graph(rumi_graph, accuracy_combined)
RRS_graph <- add_accuracy_to_graph(RRS_graph, accuracy_rrs)
SCRS_graph <- add_accuracy_to_graph(SCRS_graph, accuracy_scrs)

# network plots
overalltextsize <- 18
textsize <- 7
overallnodesize <- 8
RRS_plot <- beautify(RRS_graph, RRS_simulated_communities, title = "RRS undirected graph", no_caption = TRUE, textsize = textsize, overalltextsize = overalltextsize, item_TDR = RRS_item_TDR, overallnodesize = overallnodesize)
SCRS_plot <- beautify(SCRS_graph, SCRS_simulated_communities, title = "SCRS undirected graph", no_caption = TRUE, textsize = textsize, overalltextsize = overalltextsize, item_TDR = SCRS_item_TDR, overallnodesize = overallnodesize)
combined_plot <- beautify(rumi, combined_simulated_communities, title = "RRS & SCRS undirected graph", force_caption = TRUE, textsize = textsize, overalltextsize = overalltextsize, item_TDR = combined_item_TDR, overallnodesize = overallnodesize)
network_plot <- ggpubr::ggarrange(ggpubr::ggarrange(RRS_plot, SCRS_plot, nrow = 2, heights = c(1, 1)), combined_plot, widths = c(0.95, 1.05))
ggsave("output/networkplot.png", network_plot, dpi = 400, device = "png", width = 35, height = 18)
save_figure(network_plot, "Figure_1_Tableau_View_Graphs", height=18, width=35)

# community detection methods
i <- 1000
RRS_tpr <- suppressMessages(get_TPR(RRS_simulated_communities, RRS_graph))
RRS_tpr_s <- suppressMessages(RRS_tpr %>%
    group_by(Subscale) %>%
    summarise(TPR_m = mean(TPR)))
SCRS_tpr <- suppressMessages(get_TPR(SCRS_simulated_communities, SCRS_graph))
SCRS_tpr_s <- suppressMessages(SCRS_tpr %>%
    group_by(Subscale, Method, Weights) %>%
    summarise(TPR_m = mean(TPR)))
combined_tpr <- suppressMessages(get_TPR(combined_simulated_communities, rumi))
combined_tpr_s <- suppressMessages(combined_tpr %>%
    group_by(Subscale) %>%
    summarise(TPR_m = mean(TPR)))
combinend_tpr_avg <- mean(combined_tpr_s$TPR_m)
combined_tpr_top <- suppressMessages(combined_tpr %>%
    filter(!Method %in% c("walktrap", "optimal", "fast_greedy")) %>%
    group_by(Subscale) %>%
    summarise(TPR_m = mean(TPR)))

# summary table
seed <- 647
weights <- TRUE
method <- "optimal"
RRS_summary <- network_summary(RRS_graph, name = "RRS", seed = seed, method = method, weights = weights, item_TPR = RRS_item_TDR)
SCRS_summary <- network_summary(SCRS_graph, name = "SCRS", seed = seed, method = method, weights = weights, item_TPR = SCRS_item_TDR)
combined_graph_summary <- network_summary(rumi, name = "Combined", seed = seed, method = method, weights = weights, item_TPR = combined_item_TDR)
networksummary_table <- bind_rows(RRS_summary, SCRS_summary, combined_graph_summary)
saveRDS(networksummary_table, "./output/networksummary_table.RDS")


# create a questionnaire scales and add them to data
scales <- data.frame(
    scale = c(rep("reflection", 5), rep("brooding", 5), rep("SCRS", 10)),
    item = c(paste0("RRS_", c(2, 4, 5, 9, 10)), paste0("RRS_", c(1, 3, 6, 7, 8)), paste0("SCRS_", seq(1, 10)))
)


# reliability results
RRS_total <- dataset[paste0("RRS_", seq(1:10))]
crA_RRS_total <- ltm::cronbach.alpha(RRS_total)
RRS_reflection <- dataset[paste0("RRS_", c(2, 4, 5, 9, 10))]
crA_RRS_reflection <- ltm::cronbach.alpha(RRS_reflection)
RRS_brooding <- dataset[paste0("RRS_", c(1, 3, 6, 7, 8, 10))]
crA_RRS_brooding <- ltm::cronbach.alpha(RRS_brooding)

SCRS_total <- dataset[paste0("SCRS_", seq(1:10))]
crA_SCRS_total <- ltm::cronbach.alpha(SCRS_total)
crA_title <- paste0("Cronbach's ", Alpha())

summary_table <- tribble(
    ~Scale, ~crA,
    "RRS brooding", Round(crA_RRS_brooding$alpha[[1]]),
    "RRS reflection", Round(crA_RRS_reflection$alpha[[1]]),
    "RRS", Round(crA_RRS_total$alpha[[1]]),
    "SCRS", Round(crA_SCRS_total$alpha[[1]])
)
colnames(summary_table) <- c("Scale", crA_title)

dataset_long <- dataset %>%
    mutate(ID = row_number()) %>%
    pivot_longer(RRS_1:SCRS_10, names_to = "item", values_to = "value") %>%
    dplyr::select(-nem) %>%
    left_join(scales)


scores <- dataset_long %>%
    group_by(gender, ID, scale) %>%
    summarise(score = sum(value))
scores <- scores %>%
    spread(scale, score) %>%
    mutate(RRS = brooding + reflection)
scores <- scores %>% pivot_longer(3:6)
colnames(scores)[3:4] <- c("scale", "score")
check_homoscedasticity <- function(dataset, group_index_vector, measured) {
    res <- list()
    for (g in unique(group_index_vector)) {
        print(g)
        estimated_var <- var(dataset[c(group_index_vector == g), measured][[1]], na.rm = TRUE)
        res[g] <- estimated_var
    }
    res
}

check_normality <- function(dataset, mapping, method = "histogram_overlay") {
    require("ggplot2")
    g <- ggplot(dataset, mapping = mapping)
    g +
        geom_histogram(color = "black", position = "dodge") +
        scale_fill_manual(values = c("orange", "cornflowerblue"))
}


### Cohen's d
### SCRS
SCRSf <- scores %>%
    filter(gender == "female", scale == "SCRS") %>%
    ungroup() %>%
    dplyr::select(-ID, -gender, -scale, score) %>%
    unlist() %>%
    unname()
SCRSm <- scores %>%
    filter(gender == "male", scale == "SCRS") %>%
    ungroup() %>%
    dplyr::select(-ID, -gender, -scale, score) %>%
    unlist() %>%
    unname()
SCRScohenD <- akp.effect(SCRSf, SCRSm, tr = 0, EQVAR = TRUE)
SCRScohenDtr <- akp.effect(SCRSf, SCRSm, tr = 0.2, EQVAR = TRUE)
SCRScohenDtrboth <- akp.effect(SCRSf, SCRSm, tr = 0.2, EQVAR = FALSE)
SCRScohenDexpl <- yuenv2(SCRSf, SCRSm, tr = 0.2)
SCRSCliff <- cidv2(SCRSf, SCRSm)[c(5, 8)]

### brooding
broodingf <- scores %>%
    filter(gender == "female", scale == "brooding") %>%
    ungroup() %>%
    dplyr::select(-ID, -gender, -scale, score) %>%
    unlist() %>%
    unname()
broodingm <- scores %>%
    filter(gender == "male", scale == "brooding") %>%
    ungroup() %>%
    dplyr::select(-ID, -gender, -scale, score) %>%
    unlist() %>%
    unname()
broodingcohenD <- akp.effect(broodingf, broodingm, tr = 0, EQVAR = TRUE)
brooingcohenDtr <- akp.effect(broodingf, broodingm, tr = 0.2, EQVAR = TRUE)
broodingcohenDtrboth <- akp.effect(broodingf, broodingm, tr = 0.2, EQVAR = FALSE)
broodingcohenDexpl <- yuenv2(broodingf, broodingm, tr = 0.2)
broodingCliff <- cidv2(broodingf, broodingm)[c(5, 8)]

### reflection
reflectionf <- scores %>%
    filter(gender == "female", scale == "reflection") %>%
    ungroup() %>%
    dplyr::select(-ID, -gender, -scale, score) %>%
    unlist() %>%
    unname()
reflectionm <- scores %>%
    filter(gender == "male", scale == "reflection") %>%
    ungroup() %>%
    dplyr::select(-ID, -gender, -scale, score) %>%
    unlist() %>%
    unname()
reflectioncohenD <- akp.effect(reflectionf, reflectionm, tr = 0, EQVAR = TRUE)
reflectioncohenDtr <- akp.effect(reflectionf, reflectionm, tr = 0.2, EQVAR = TRUE)
reflectioncohenDtrboth <- akp.effect(reflectionf, reflectionm, tr = 0.2, EQVAR = FALSE)
reflectioncohenDexpl <- yuenv2(reflectionf, reflectionm, tr = 0.2)
reflectionCliff <- cidv2(reflectionf, reflectionm)[c(5, 8)]


### RRS total
RRSscoresf <- scores %>%
    filter(gender == "female", scale == "RRS") %>%
    ungroup() %>%
    dplyr::select(-ID, -gender, -scale, score) %>%
    unlist() %>%
    unname()
RRSscoresm <- scores %>%
    filter(gender == "male", scale == "RRS") %>%
    ungroup() %>%
    dplyr::select(-ID, -gender, -scale, score) %>%
    unlist() %>%
    unname()
RRSTOTALcohenD <- akp.effect(RRSscoresf, RRSscoresm, tr = 0, EQVAR = TRUE)
RRSTOTALcohenDtr <- akp.effect(RRSscoresf, RRSscoresm, tr = 0.2, EQVAR = TRUE)
RRSTOTALcohenDtrboth <- akp.effect(RRSscoresf, RRSscoresm, tr = 0.2, EQVAR = FALSE)
RRSTOTALcohenDexpl <- yuenv2(RRSscoresf, RRSscoresm, tr = 0.2)
RRSTOTALCliff <- cidv2(RRSscoresf, RRSscoresm)[c(5, 8)]

cohen_summary <- tribble(
    ~`Cohen's d (20% trimmed)`, ~`Cliff's delta`, ~`Explanatory measure of ES`,
    Round(brooingcohenDtr), paste0("P(f>m)=", round(broodingCliff$summary.dvals[1, 3], 2)), Round(broodingcohenDexpl$Effect.Size),
    Round(reflectioncohenDtr), paste0("P(f>m)=", round(reflectionCliff$summary.dvals[1, 3], 2)), Round(reflectioncohenDexpl$Effect.Size),
    Round(SCRScohenDtr), paste0("P(f>m)=", round(SCRSCliff$summary.dvals[1, 3], 2)), Round(SCRScohenDexpl$Effect.Size),
    Round(RRSTOTALcohenDtr), paste0("P(f>m)=", round(RRSTOTALCliff$summary.dvals[1, 3], 2)), Round(RRSTOTALcohenDexpl$Effect.Size),
)

summary_table <- bind_cols(summary_table, cohen_summary)


# score summaries
score_summary <- scores %>%
    rowwise() %>%
    mutate(tool = if_else(str_detect("SCRS", scale), "SCRS", "RRS")) %>%
    mutate(Scale = if_else(str_detect("SCRS", tool), tool, paste0("RRS ", scale))) %>%
    group_by(Scale, gender) %>%
    summarise(
        M = mean(score),
        SD = sd(score)
    )
score_summary <- score_summary %>%
    ungroup() %>%
    mutate(Scale = if_else(Scale == "RRS RRS", "RRS", Scale))

score_summary <- bind_rows(score_summary) %>% arrange(Scale)
score_summary <- score_summary %>%
    mutate(MSD = paste0(Round(M), " (", Round(SD), ")")) %>%
    dplyr::select(-M, -SD) %>%
    pivot_wider(names_from = "gender", values_from = "MSD") %>%
    ungroup()


summary_table <- summary_table %>% arrange(Scale)
summary_table <- summary_table %>% left_join(score_summary)


# corr matrix
cormat <- scores %>%
    pivot_wider(names_from = "scale", values_from = "score") %>%
    ungroup() %>%
    dplyr::select(-gender, -ID) %>%
    mutate(RRS = brooding + reflection) %>%
    dplyr::select(RRS, brooding, reflection, SCRS)
corr <- data.frame(cor(cormat))

corr[upper.tri(cor(cormat), diag = TRUE)] <- ""
corr[lower.tri(corr, diag = FALSE)] <- Round(as.double(corr[lower.tri(corr, diag = FALSE)]))
colnames(corr) <- 1:4

cormat <- pcor(cormat)

cormat$p.value[cormat$p.value == 0] <- 0.001

lower_table <- summary_table %>%
    pivot_longer(crA_title:male, names_to = "col", values_to = "val", values_transform = as.character) %>%
    mutate(val = if_else(!is.na(as.numeric(val)), Round(val), val)) %>%
    pivot_wider(names_from = "Scale", values_from = "val") %>%
    dplyr::select(col, `RRS`, `RRS brooding`, `RRS reflection`, `SCRS`)


lower_table <- lower_table[c(5, 6, 1, 2, 3), ]
colnames(lower_table)[2:5] <- c(1:4)
corr <- bind_cols(tibble(col = c("1. RRS", "2. RRS - brooding", "3. RRS - reflection", "4. SCRS")), corr)
summary_table <- bind_rows(corr, lower_table)
colnames(summary_table) <- c("", paste0(1:4, "."))

summary_table[5, 1] <- "Female Mean (SD)"
summary_table[6, 1] <- "Male Mean (SD)"
opening_line <- tibble_row("Correlations:", `1.` = "", `2.` = "", `3.` = "", `4.` = "")
colnames(opening_line)[1] <- ""
summary_table <- bind_rows(opening_line, summary_table)
colnames(summary_table)[1] <- ""



#### glossary tibble
 glossary_tb <- tribble(
     ~Term, ~Definition,
     "Articulation point", "A measure of vertex influence: V(i) vertex of G graph is an Articulation Point, if removed from the graph, would segment, cut or disintegrate the graph; therefore is considered locally important.",
     "","",
     "Association network","A network where neighbouring vertices when connected are assumed to being associated; in other words, edges represent some measure of correlation or other kind of association; this makes them different from energy or information flow networks, for instance.",
     "","",
      "Community Detection", "A group of algorithmic methods to find clusters of nodes  whose community can be described by systematic variance encoded in data represented in edges and node characteristics, and therefore are non-trivial. In the current paper, we used four selected algorithms to estimate how well they identify the original subscale-structure of the RRS and SCRS questionnaires.",
     "","",
      "Cut vertex", "see 'Articulation point'",
     "","",
     "Edge", "A link (E, arc) between nodes of a network. Undirected if direction of the connection (e.g. flow of energy or information) is not essential to the network; directed (uni- or bidirectional) otherwise. Can be weighted by some measure of association, frequency or other attribute.",
     "","",
     "Exponential Random Graph Models", "ERGM; a family of models to the analogy of generalized linear models adopted to network modeling. They allow the construction, fitting and comparison of network models. Their nature is structural, incorporating different configurations of a graph (sets of possible edges). ERGMs have numerous extensions and variations.",
     "","",
     "Graph", "see 'Network'",
     "","",
     "Markov Chain Monte Carlo","Monte Carlo methods are a collection of techniques which approximate a quantity (mean, variance, some other measure of location, for instance) by generating random values of a random variable. In the present paper, MCMC was used to get the ERGM model parameters; and also in the fitting of the eigenmodel for the latent network modeling.",
     "","",
     "Network", "Graph G = (V,E), a construct of V vertices and E edges.",
     "","",
     "Network, Adjancency matrix", "Given v number of vertices, the adjacency matrix is a v X v matrix, where any [i,j] cell's value reflects with 0 or 1 whether a selected pair of vertices are connected or not",
     "","",
     "Network, Clustering coefficient", "Transitivity, a descriptive measure to reflect the degree of connectedness [or process of clustering] in a graph by giving the ratio of connected triangles (triads of vertices) versus the total number of triangles in the graph.",
     "","",
     "Network, Diameter", "The value of the longest distance (traceable route along edges and incident nodes) in a graph.",
     "","",
     "Network, Order", "Descriptive; order of a network is the number of vertices.",
     "","",
     "Network, Size","Descriptive; size of a network is the total number of edges.",
     "","",
     "Network, Vertex, Betweenness","Descriptive of vertex centrality; quantifies how much a given vertex is between other pairs of vertices, thus has a connecting power.",
     "","",
     "Network, Vertex, Closeness","Descriptive of vertex centrality; assumes that important vertices are closer to others; its measure is inverse of total distances to other vertices.",
     "","",
     "Network, Vertex, Degree","Descriptive; the degree of a vertex is the number of edges incident on the given vertex. Vertex degrees are summarized as Degree Distribution of a graph.",
     "","",
     "Network, Vertex, Strength", "Descriptive of vertex centrality; equals the sum of edge weights incident upon the given vertex, reflecting vertex centrality (importance). Vertex Strength Distribution is a histogram of vertex strength estimates of an entire graph.",
     "","",
     "Node", "Nodes (V,vertices) are representation of constructs, latent variables, subjects or any level of attributes in a network. Like edges, they can inherently bear any attributes.",
     "","",
       "Penalization", "In regression, penalization introduces additional penalty to the sum of squared errors the model is trying to minimize. This way, by allowing a small increase in bias, the model is expected to exhibit  a lower level of variance. The graphical LASSO method, used in the current paper, uses an absolute fraction of coefficietnts as penalties. With the applied penalty, insignificant or weak associations are expected to shrink towards 0."
 )
