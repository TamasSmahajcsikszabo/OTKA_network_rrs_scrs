library(knitr)
library(ltm)
library(tibble)
library(WRS)
library(dplyr)
library(tidyr)
library(rogme)

rumi <- readRDS("/home/tamas/repos/networks_with_r/full_data.RDS")
rumi <- na.omit(rumi[,c(4, 13:32)])
rumi$gender <- ""
rumi$gender[rumi$nem=='lany'] <- 'female'
rumi$gender[rumi$nem=='fiu'] <- 'male'
rumi$gender <- factor(rumi$gender)
dataset <- rumi

scales <- data.frame(scale=c(rep('brooding', 5), rep('reflection',5), rep('SCRS',10)), 
    item=c(paste0('RRS_',c(2,4,5,9,10)), paste0("RRS_",c(1,3,6,7,8)), paste0("SCRS_",seq(1,10))))

rumi <- rumi  %>% 
    mutate(ID = row_number())  %>% 
    pivot_longer(RRS_1:SCRS_10, names_to='item', values_to='value') %>% 
    dplyr::select(-nem) %>% 
    left_join(scales)

RRS_total <- dataset[paste0("RRS_", seq(1:10))]
crA_RRS_total <- ltm::cronbach.alpha(RRS_total)



RRS_reflection <- dataset[paste0("RRS_", c(2,4,5,9,10))]
crA_RRS_reflection <- ltm::cronbach.alpha(RRS_reflection)


RRS_brooding <- dataset[paste0("RRS_", c(1,3,6,7,8,10))]
crA_RRS_brooding <- ltm::cronbach.alpha(RRS_brooding)

SCRS_total <- dataset[paste0("SCRS_", seq(1:10))]
crA_SCRS_total <- ltm::cronbach.alpha(SCRS_total)

summary_table <- tribble(
    ~Scale, ~`Cronbach's Alpha`,
    "RRS brooding", crA_RRS_brooding$alpha[[1]],
    "RRS reflection", crA_RRS_reflection$alpha[[1]],
    "RRS total", crA_RRS_total$alpha[[1]],
    "SCRS total", crA_SCRS_total$alpha[[1]]
)

scores <- rumi %>% group_by(gender, ID, scale) %>% summarise(score=sum(value))

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


# akp.effect(dataset[!dataset$gender=='male',]$reflection,dataset[dataset$gender=='male',]$reflection,tr=0.2)

check_homoscedasticity(scores[scores$scale=='brooding',], scores[scores$scale=='brooding',]$gender,"score")

check_normality(scores, aes(x=score, group=gender, fill=gender)) + facet_wrap(~gender)

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
RRSscoresf <- rumi %>% 
    group_by(gender, ID, scale) %>% 
    summarise(score=sum(value)) %>% 
    ungroup() %>% 
    filter(!scale == "SCRS") %>% 
    group_by(gender,ID) %>% summarise(score=sum(score)) %>% 
    ungroup() %>% 
    filter(gender=='female') %>% 
    dplyr::select(-gender, -ID, score) %>% unlist()  %>%  unname()
RRSscoresm <- rumi %>% 
    group_by(gender, ID, scale) %>% 
    summarise(score=sum(value)) %>% 
    ungroup() %>% 
    filter(!scale == "SCRS") %>% 
    group_by(gender,ID) %>% summarise(score=sum(score)) %>% 
    ungroup() %>% 
    filter(gender=='male') %>% 
    dplyr::select(-gender, -ID, score) %>% unlist()  %>%  unname()

RRSTOTALcohenD <- akp.effect(RRSscoresf, RRSscoresm, tr=0, EQVAR=TRUE)
RRSTOTALcohenDtr <- akp.effect(RRSscoresf, RRSscoresm, tr=0.2, EQVAR=TRUE)
RRSTOTALcohenDtrboth <- akp.effect(RRSscoresf,  RRSscoresm, tr=0.2, EQVAR=FALSE)
RRSTOTALcohenDexpl <- yuenv2(RRSscoresf,  RRSscoresm, tr=0.2)
RRSTOTALCliff <- cidv2(RRSscoresf,  RRSscoresm)[c(5,8)]

cohen_summary <- tribble(
    ~`Cohen's d (20% trimmed)`, ~`Cliff's delta`, ~`Explanatory measure of ES`,
    brooingcohenDtr,   paste0("P(f>m)=",round(broodingCliff$summary.dvals[1,3],3)), broodingcohenDexpl$Effect.Size,
    reflectioncohenDtr,   paste0("P(f>m)=",round(reflectionCliff$summary.dvals[1,3],3)), reflectioncohenDexpl$Effect.Size,
    SCRScohenDtr,   paste0("P(f>m)=",round(SCRSCliff$summary.dvals[1,3],3)), SCRScohenDexpl$Effect.Size,
    RRSTOTALcohenDtr,   paste0("P(f>m)=",round(RRSTOTALCliff$summary.dvals[1,3],3)), RRSTOTALcohenDexpl$Effect.Size,
)

summary_table <- bind_cols(summary_table, cohen_summary)
