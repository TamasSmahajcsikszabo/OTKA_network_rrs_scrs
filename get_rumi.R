setwd('data/')
rumi<-readRDS("rumi_graph_full.RDS")
rumid<-readRDS("rumi_graph_full_din.RDS")
prumi<-readRDS("personified_rumi.RDS")
rumi_data<-readRDS("rumi_data.RDS")
temporal_summary <- readRDS("temporal_summary.RDS")
dynamic <- readRDS("dynamic.RDS")
ids <- readRDS("ids.RDS")
mycorr <- readRDS("mycorr.RDS")
RRS_graph <- readRDS("RRS.RDS")
SCRS_graph <- readRDS("SCRS.RDS")
setwd('../')
RRS_simulated_communities <- readRDS("output/RRS_communities.RDS")
SCRS_simulated_communities <- readRDS("output/SCRS_communities.RDS")
combined_simulated_communities <- readRDS("output/combined_communities.RDS")
RRS_item_TDR <- readRDS('output/RRS_item_TPR.RDS')
SCRS_item_TDR <- readRDS('output/SCRS_item_TPR.RDS')
combined_item_TDR <- readRDS('output/combined_item_TPR.RDS')
accuracy_rrs  <- readRDS("output/accuracy_rrs.RDS")
accuracy_scrs  <- readRDS("output/accuracy_scrs.RDS")
accuracy_combined <- readRDS("output/accuracy_combined.RDS")
accuracy_rrs <- network_accuracy_data(accuracy_rrs)
accuracy_scrs <- network_accuracy_data(accuracy_scrs)
accuracy_combined <- network_accuracy_data(accuracy_combined)
get_itemnames <- function(namevector){
    gridN <- tibble(expand.grid(namevector, namevector))
    gridN <- gridN %>% filter(!Var1==Var2)
    gridN
}
# add accuracy to networks
rumi <- add_accuracy_to_graph(rumi, accuracy_combined)
RRS_graph <- add_accuracy_to_graph(RRS_graph, accuracy_rrs)
SCRS_graph <- add_accuracy_to_graph(SCRS_graph, accuracy_scrs)

overalltextsize <- 24
textsize <- 7
overallnodesize  <- 10
RRS_plot <- beautify(RRS_graph, RRS_simulated_communities, title="RRS undirected graph", no_caption=TRUE, textsize=textsize, overalltextsize=overalltextsize,item_TDR=RRS_item_TDR, overallnodesize=overallnodesize)
SCRS_plot <- beautify(SCRS_graph, SCRS_simulated_communities, title="SCRS undirected graph", no_caption=TRUE, textsize=textsize, overalltextsize=overalltextsize,item_TDR=SCRS_item_TDR, overallnodesize=overallnodesize)
combined_plot <- beautify(rumi, combined_simulated_communities, title="RRS & SCRS undirected graph",force_caption=TRUE, textsize=textsize, overalltextsize=overalltextsize, item_TDR=combined_item_TDR, overallnodesize=overallnodesize)
#network_plot <- ggarrange(ggarrange(RRS_plot, SCRS_plot, nrow=2,heights=c(1,1)), combined_plot, widths=c(0.95,1.05))
#ggsave("output/networkplot.png",network_plot, dpi=400, device='png',width=35, height=18)
RRS_tpr <- suppressMessages(get_TPR(RRS_simulated_communities,RRS_graph))

RRS_tpr_s <- suppressMessages(RRS_tpr %>% 
    group_by(Subscale)  %>%
    summarise(TPR_m = mean(TPR)))

SCRS_tpr <- suppressMessages(get_TPR(SCRS_simulated_communities,SCRS_graph))
SCRS_tpr_s <- suppressMessages(SCRS_tpr %>% 
    group_by(Subscale, Method, Weights)  %>%
    summarise(TPR_m = mean(TPR)))

combined_tpr <- suppressMessages(get_TPR(combined_simulated_communities,rumi))
combined_tpr_s <- suppressMessages(combined_tpr %>% 
    group_by(Subscale)  %>%
    summarise(TPR_m = mean(TPR)))
combinend_tpr_avg <- mean(combined_tpr_s$TPR_m)

combined_tpr_top <- suppressMessages(combined_tpr %>% 
    filter(!Method %in% c("walktrap", "optimal","fast_greedy")) %>%
    group_by(Subscale)  %>%
    summarise(TPR_m = mean(TPR)))



# summary table 
seed <- 647
weights <- TRUE
method <- 'optimal'
networksummary_table <- bind_rows(network_summary(RRS_graph, name="RRS",seed=seed,method=method,weights=weights, item_TPR=RRS_item_TDR),
                            network_summary(SCRS_graph, name="SCRS", seed=seed,method=method,weights=weights, item_TPR=SCRS_item_TDR),
                            network_summary(rumi, name="Combined",seed=seed,method=method,weights=weights, item_TPR=combined_item_TDR))

