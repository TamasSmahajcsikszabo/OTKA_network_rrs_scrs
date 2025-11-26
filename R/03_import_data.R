library(knitr)
library(ltm)
library(tibble)
library(WRS)
library(dplyr)
library(tidyr)
library(rogme)
library(magrittr)
library(stringr)

#import main data
dataset_raw <- read.spss("data/scrs888.sav", to.data.frame = TRUE, use.value.labels = TRUE)
dataset_raw <- as_tibble(dataset_raw) # transform raw dataframe into tibble dataframe

# extract data labels
data_labels <- tibble(var = names(dataset_raw), labels = attributes(dataset_raw)$variable.labels) %>%
  filter(labels != "")

rumidata <- dataset_raw %>% mutate(ID = row_number())
origdata <- rumidata
rumidata <- na.omit(rumidata[,c(4, 13:32, 56)])
filterd_idx <- rumidata$ID
origdata <- origdata %>% filter(ID %in% filterd_idx)
rumidata  <- rumidata[,1:21]
rumidata$gender <- ""
rumidata$gender[rumidata$nem=='lany'] <- 'female'
rumidata$gender[rumidata$nem=='fiu'] <- 'male'
rumidata$gender <- factor(rumidata$gender)

# prepared dataset
dataset <- rumidata

# load computed measures
rumi<-readRDS("data/rumi_graph_full.RDS")
rumid<-readRDS("data/rumi_graph_full_din.RDS")
prumi<-readRDS("data/personified_rumi.RDS")
rumi_data<-readRDS("data/rumi_data.RDS")
temporal_summary <- readRDS("data/temporal_summary.RDS")
dynamic <- readRDS("data/dynamic.RDS")
ids <- readRDS("data/ids.RDS")
mycorr <- readRDS("data/mycorr.RDS")
RRS_graph <- readRDS("data/RRS.RDS")
SCRS_graph <- readRDS("data/SCRS.RDS")
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

# add accuracy to networks
rumi <- add_accuracy_to_graph(rumi, accuracy_combined)
RRS_graph <- add_accuracy_to_graph(RRS_graph, accuracy_rrs)
SCRS_graph <- add_accuracy_to_graph(SCRS_graph, accuracy_scrs)

# network plots
overalltextsize <- 24
textsize <- 7
overallnodesize  <- 10
RRS_plot <- beautify(RRS_graph, RRS_simulated_communities, title="RRS undirected graph", no_caption=TRUE, textsize=textsize, overalltextsize=overalltextsize,item_TDR=RRS_item_TDR, overallnodesize=overallnodesize)
SCRS_plot <- beautify(SCRS_graph, SCRS_simulated_communities, title="SCRS undirected graph", no_caption=TRUE, textsize=textsize, overalltextsize=overalltextsize,item_TDR=SCRS_item_TDR, overallnodesize=overallnodesize)
combined_plot <- beautify(rumi, combined_simulated_communities, title="RRS & SCRS undirected graph",force_caption=TRUE, textsize=textsize, overalltextsize=overalltextsize, item_TDR=combined_item_TDR, overallnodesize=overallnodesize)
network_plot <- ggarrange(ggarrange(RRS_plot, SCRS_plot, nrow=2,heights=c(1,1)), combined_plot, widths=c(0.95,1.05))
ggsave("output/networkplot.png",network_plot, dpi=400, device='png',width=35, height=18)

# community detection methods
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


# create a questionnaire scales and add them to data
scales <- data.frame(scale=c(rep('reflection', 5), rep('brooding',5), rep('SCRS',10)), 
    item=c(paste0('RRS_',c(2,4,5,9,10)), paste0("RRS_",c(1,3,6,7,8)), paste0("SCRS_",seq(1,10))))


# reliability results
RRS_total <- dataset[paste0("RRS_", seq(1:10))]
crA_RRS_total <- ltm::cronbach.alpha(RRS_total)
RRS_reflection <- dataset[paste0("RRS_", c(2,4,5,9,10))]
crA_RRS_reflection <- ltm::cronbach.alpha(RRS_reflection)
RRS_brooding <- dataset[paste0("RRS_", c(1,3,6,7,8,10))]
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
colnames(summary_table)  <-  c('Scale', crA_title)

dataset_long <- dataset %>% 
    mutate(ID = row_number())  %>% 
    pivot_longer(RRS_1:SCRS_10, names_to='item', values_to='value') %>% 
    dplyr::select(-nem) %>% 
    left_join(scales)


scores <- dataset_long%>% group_by(gender, ID, scale) %>% summarise(score=sum(value))
scores  <- scores %>% spread(scale, score) %>% mutate(RRS = brooding + reflection)
scores <- scores %>% pivot_longer(3:6)
colnames(scores)[3:4] <- c("scale", "score")
check_homoscedasticity <- function(dataset, group_index_vector, measured) {
    res <- list()
    for (g in unique(group_index_vector)) {
        print(g)
        estimated_var <- var(dataset[c(group_index_vector == g),measured][[1]],na.rm=TRUE)
        res[g] <- estimated_var
    }
    res
}

check_normality <- function(dataset, mapping, method='histogram_overlay') {
    require("ggplot2")
    g <- ggplot(dataset, mapping=mapping)
    g + 
        geom_histogram(color='black',position='dodge') +
        scale_fill_manual(values=c("orange", "cornflowerblue"))
}



### Cohen's d
### SCRS
SCRSf <- scores %>% filter(gender=='female', scale=='SCRS') %>% ungroup() %>% dplyr::select(-ID,-gender,-scale, score) %>% unlist() %>% unname()
SCRSm <- scores %>% filter(gender=='male', scale=='SCRS') %>% ungroup() %>% dplyr::select(-ID,-gender,-scale, score) %>% unlist() %>% unname()
SCRScohenD <- akp.effect(SCRSf, SCRSm, tr=0, EQVAR=TRUE)
SCRScohenDtr <- akp.effect(SCRSf, SCRSm, tr=0.2, EQVAR=TRUE)
SCRScohenDtrboth <- akp.effect(SCRSf, SCRSm, tr=0.2, EQVAR=FALSE)
SCRScohenDexpl <- yuenv2(SCRSf, SCRSm, tr=0.2)
SCRSCliff <- cidv2(SCRSf, SCRSm)[c(5,8)]

### brooding
broodingf <- scores %>% filter(gender=='female', scale=='brooding') %>% ungroup() %>% dplyr::select(-ID,-gender,-scale, score) %>% unlist() %>% unname()
broodingm <- scores %>% filter(gender=='male', scale=='brooding') %>% ungroup() %>% dplyr::select(-ID,-gender,-scale, score) %>% unlist() %>% unname()
broodingcohenD <- akp.effect(broodingf, broodingm, tr=0, EQVAR=TRUE)
brooingcohenDtr <- akp.effect(broodingf,  broodingm, tr=0.2, EQVAR=TRUE)
broodingcohenDtrboth <- akp.effect(broodingf,  broodingm, tr=0.2, EQVAR=FALSE)
broodingcohenDexpl <- yuenv2(broodingf,  broodingm, tr=0.2)
broodingCliff <- cidv2(broodingf,  broodingm)[c(5,8)]

### reflection
reflectionf <- scores %>% filter(gender=='female', scale=='reflection') %>% ungroup() %>% dplyr::select(-ID,-gender,-scale, score) %>% unlist() %>% unname()
reflectionm <- scores %>% filter(gender=='male', scale=='reflection') %>% ungroup() %>% dplyr::select(-ID,-gender,-scale, score) %>% unlist() %>% unname()
reflectioncohenD <- akp.effect(reflectionf, reflectionm, tr=0, EQVAR=TRUE)
reflectioncohenDtr <- akp.effect(reflectionf,  reflectionm, tr=0.2, EQVAR=TRUE)
reflectioncohenDtrboth <- akp.effect(reflectionf,  reflectionm, tr=0.2, EQVAR=FALSE)
reflectioncohenDexpl <- yuenv2(reflectionf,  reflectionm, tr=0.2)
reflectionCliff <- cidv2(reflectionf,  reflectionm)[c(5,8)]


### RRS total
RRSscoresf <- scores %>% filter(gender=='female', scale=='RRS') %>% ungroup() %>% dplyr::select(-ID,-gender,-scale, score) %>% unlist() %>% unname()
RRSscoresm <- scores %>% filter(gender=='male', scale=='RRS') %>% ungroup() %>% dplyr::select(-ID,-gender,-scale, score) %>% unlist() %>% unname()
RRSTOTALcohenD <- akp.effect(RRSscoresf, RRSscoresm, tr=0, EQVAR=TRUE)
RRSTOTALcohenDtr <- akp.effect(RRSscoresf, RRSscoresm, tr=0.2, EQVAR=TRUE)
RRSTOTALcohenDtrboth <- akp.effect(RRSscoresf,  RRSscoresm, tr=0.2, EQVAR=FALSE)
RRSTOTALcohenDexpl <- yuenv2(RRSscoresf,  RRSscoresm, tr=0.2)
RRSTOTALCliff <- cidv2(RRSscoresf,  RRSscoresm)[c(5,8)]

cohen_summary <- tribble(
    ~`Cohen's d (20% trimmed)`, ~`Cliff's delta`, ~`Explanatory measure of ES`,
    Round(brooingcohenDtr),   paste0("P(f>m)=",round(broodingCliff$summary.dvals[1,3],2)), Round(broodingcohenDexpl$Effect.Size),
    Round(reflectioncohenDtr),   paste0("P(f>m)=",round(reflectionCliff$summary.dvals[1,3],2)), Round(reflectioncohenDexpl$Effect.Size),
    Round(SCRScohenDtr),   paste0("P(f>m)=",round(SCRSCliff$summary.dvals[1,3],2)), Round(SCRScohenDexpl$Effect.Size),
    Round(RRSTOTALcohenDtr),   paste0("P(f>m)=",round(RRSTOTALCliff$summary.dvals[1,3],2)), Round(RRSTOTALcohenDexpl$Effect.Size),
)

summary_table <- bind_cols(summary_table, cohen_summary)


# score summaries
score_summary <- scores %>% 
    rowwise() %>% 
    mutate(tool = if_else(str_detect("SCRS", scale), "SCRS", "RRS")) %>% 
    mutate(Scale = if_else(str_detect("SCRS", tool), tool, paste0("RRS ", scale))) %>%
    group_by(Scale, gender) %>% 
    summarise(M=mean(score),
              SD=sd(score))
score_summary <- score_summary %>%  ungroup() %>% mutate(Scale = if_else(Scale == "RRS RRS", "RRS", Scale))

score_summary <- bind_rows(score_summary) %>% arrange(Scale)
score_summary <- score_summary %>% 
    mutate(MSD = paste0(Round(M), ' (',Round(SD), ')')) %>% 
    dplyr::select(-M, -SD) %>% 
    pivot_wider(names_from='gender', values_from="MSD") %>% 
    ungroup()


summary_table <- summary_table %>% arrange(Scale)
summary_table <- summary_table %>% left_join(score_summary)


# corr matrix
cormat <- scores %>% pivot_wider(names_from='scale', values_from='score') %>% 
    ungroup() %>% 
    dplyr::select(-gender, -ID) %>% 
    mutate(RRS = brooding+reflection) %>% 
    dplyr::select(RRS, brooding, reflection, SCRS)
corr <- data.frame(cor(cormat))

corr[upper.tri(cor(cormat), diag=TRUE)] <- ""
corr[lower.tri(corr, diag=FALSE)] <- Round(as.double(corr[lower.tri(corr, diag=FALSE)] ))
colnames(corr) <- 1:4

cormat <- pcor(cormat)

cormat$p.value[cormat$p.value==0] <- 0.001

lower_table <- summary_table %>% 
    pivot_longer(crA_title:male, names_to='col', values_to='val', values_transform=as.character) %>%
    mutate(val=if_else(!is.na(as.numeric(val)), Round(val), val)) %>% 
pivot_wider(names_from='Scale', values_from='val') %>% 
dplyr::select(col, `RRS`, `RRS brooding`, `RRS reflection`, `SCRS`)


lower_table <- lower_table[c(5,6,1,2,3),]
colnames(lower_table)[2:5] <- c(1:4)
corr <- bind_cols(tibble(col=c('1. RRS','2. RRS - brooding','3. RRS - reflection','4. SCRS')),corr)
summary_table <-bind_rows(corr, lower_table)
colnames(summary_table) <- c("", paste0(1:4, "."))

summary_table[5,1] <- 'Female Mean (SD)'
summary_table[6,1] <- 'Male Mean (SD)'
opening_line <- tibble_row("Correlations:", `1.`="", `2.`="", `3.`="", `4.`="")
colnames(opening_line)[1] <- ''
summary_table <- bind_rows(opening_line, summary_table)
colnames(summary_table)[1] <- ''

