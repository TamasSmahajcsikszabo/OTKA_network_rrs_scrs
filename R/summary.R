library(knitr)
library(ltm)
library(tibble)
library(WRS)
library(dplyr)
library(tidyr)
library(rogme)
library(magrittr)
library(stringr)

Alpha <- function(){
    '$\\alpha$'
}

Round <- function(value) {
        return(format(round(as.numeric(value), 2), nsmall=2))
}

rumidata <- readRDS("/home/tamas/repos/networks_with_r/full_data.RDS")
rumidata <- rumidata %>% mutate(ID = row_number())
origdata <- rumidata
rumidata <- na.omit(rumidata[,c(4, 13:32, 56)])
filterd_idx <- rumidata$ID
origdata <- origdata %>% filter(ID %in% filterd_idx)
rumidata  <- rumidata[,1:21]
rumidata$gender <- ""
rumidata$gender[rumidata$nem=='lany'] <- 'female'
rumidata$gender[rumidata$nem=='fiu'] <- 'male'
rumidata$gender <- factor(rumidata$gender)
dataset <- rumidata

# sanity check
# dataset %>% dplyr::select(nem, c(paste0('RRS_',c(2,4,5,9,10)), paste0("RRS_",c(1,3,6,7,8)))) %>%
#     mutate(RRS = RRS_1 + RRS_2 + RRS_3 + RRS_4 + RRS_5 + RRS_6 + RRS_7 + RRS_8 + RRS_9 + RRS_10) %>%
#     group_by(nem) %>%
#     summarise (m = mean(RRS), SD = sd(RRS))

scales <- data.frame(scale=c(rep('reflection', 5), rep('brooding',5), rep('SCRS',10)), 
    item=c(paste0('RRS_',c(2,4,5,9,10)), paste0("RRS_",c(1,3,6,7,8)), paste0("SCRS_",seq(1,10))))

rumidata <- rumidata  %>% 
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
crA_title <- paste0("Cronbach's ", Alpha())

summary_table <- tribble(
    ~Scale, ~crA,
    "RRS brooding", Round(crA_RRS_brooding$alpha[[1]]),
    "RRS reflection", Round(crA_RRS_reflection$alpha[[1]]),
    "RRS", Round(crA_RRS_total$alpha[[1]]),
    "SCRS", Round(crA_SCRS_total$alpha[[1]])
)
colnames(summary_table)  <-  c('Scale', crA_title)

scores <- rumidata %>% group_by(gender, ID, scale) %>% summarise(score=sum(value))
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


# akp.effect(dataset[!dataset$gender=='male',]$reflection,dataset[dataset$gender=='male',]$reflection,tr=0.2)

# check_homoscedasticity(scores[scores$scale=='brooding',], scores[scores$scale=='brooding',]$gender,"score")

# check_normality(scores, aes(x=score, group=gender, fill=gender)) + facet_wrap(~gender)

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
# RRSscoresf <- rumidata %>% 
#     group_by(gender, ID, scale) %>% 
#     summarise(score=sum(value)) %>% 
#     ungroup() %>% 
#     filter(!scale == "SCRS") %>% 
#     group_by(gender,ID) %>% summarise(score=sum(score)) %>% 
#     ungroup() %>% 
#     filter(gender=='female') %>% 
#     dplyr::select(-gender, -ID, score) %>% unlist()  %>%  unname()
# RRSscoresm <- rumidata %>% 
#     group_by(gender, ID, scale) %>% 
#     summarise(score=sum(value)) %>% 
#     ungroup() %>% 
#     filter(!scale == "SCRS") %>% 
#     group_by(gender,ID) %>% summarise(score=sum(score)) %>% 
#     ungroup() %>% 
#     filter(gender=='male') %>% 
#     dplyr::select(-gender, -ID, score) %>% unlist()  %>%  unname()

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
score_summary = score_summary %>%  ungroup() %>% mutate(Scale = if_else(Scale == "RRS RRS", "RRS", Scale))

# score_summary_RRS <- scores %>% 
#     rowwise() %>% 
#     mutate(tool = if_else(str_detect("SCRS", scale), "SCRS", "RRS")) %>% 
#     filter(tool=='RRS') %>% 
#     group_by(tool, gender) %>% 
#     summarise(M=mean(score),
#               SD=sd(score)) %>% 
#     ungroup() %>% 
#     dplyr::select(-tool) %>% 
#     mutate(Scale='RRS')

score_summary <- bind_rows(score_summary) %>% arrange(Scale)
score_summary <- score_summary %>% 
    mutate(MSD = paste0(Round(M), ' (',Round(SD), ')')) %>% 
    dplyr::select(-M, -SD) %>% 
    pivot_wider(names_from='gender', values_from="MSD") %>% 
    ungroup()
    # dplyr::select(-Scale)


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
# corr <- as.vector(unlist(corr))
# corr <- lapply(corr, function(c){Round(c)})
# corr <- data.frame(matrix(corr, ncol=4))


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
print(summary_table)
