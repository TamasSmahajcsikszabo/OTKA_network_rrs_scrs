library(dplyr)
library(ggplot2)
source('R/functions.R')
source('R/items.R')

RRS <- readRDS('output/RRSsummary.RDS')
SCRS <- readRDS('output/SCRSsummary.RDS')
combined <- readRDS('output/combinedSummary.RDS')
scales <- get_factor_levels()

RRS$network <- 'RRS'
SCRS$network <- 'SCRS'
combined$network <- 'RRS + SCRS'

alldata  <- bind_rows(RRS, SCRS, combined)
colnames(alldata)[1:2] <- c('node1', 'node2')
alldata <- alldata %>% left_join(items_table[,c(1,4)], by=c('node1'='Questionnaire'))
alldata <- alldata %>% left_join(items_table[,c(1,4)], by=c('node2'='Questionnaire'))
colnames(alldata)[13:14]  <- c('node1_label', 'node2_label')

alldata <- alldata %>% mutate(edge= paste0(node1_label, ' - ', node2_label))
alldata <- alldata %>%dplyr::select(1,2,3,4,5,6,7,12,8,9,10,11,13,14,15)

plotdata <- alldata %>% pivot_longer(3:11)
plotdata <- plotdata %>% filter(!is.na(value))
levels <- unique(plotdata$name)
levels <- levels[c(1,2,3,4,5,9,6,7,8)]
plotdata$name <- factor(plotdata$name, levels=levels)

ggplot(plotdata) +
    geom_point(aes(name, value, group=edge, color=edge)) +
    geom_path(aes(name, value, group=edge, color=edge)) +
    facet_wrap(~network)


combinedplotdata %>%
    group_by(network, edge, name) %>% summarize(v = mean(value))


# combined  <- alldata[alldata$network == 'RRS + SCRS',]
# RRS <-  alldata[alldata$network == 'RRS',]
# RRS <-  alldata[alldata$network == 'SCRS',]

# colnames(combined) <- paste0('combined_', colnames(combined))
# colnames(RRS) <- paste0('RRS_', colnames(RRS))
# colnames(SCRS) <- paste0('SCRS_', colnames(SCRS))

comparative <- alldata %>% pivot_longer(3:11) %>% pivot_wider(names_from = network, values_from=value)

comparative <- comparative %>%
    ungroup() %>%
    group_by(edge) %>%
    mutate(diff = if_else(is.na(SCRS), `RRS + SCRS` - RRS, `RRS + SCRS` - SCRS)) %>%
    dplyr::select(edge, name, RRS, SCRS, `RRS + SCRS`, diff)
comparative <- comparative %>% pivot_longer(3:5, names_to='network')
comparative <- comparative %>% filter(!is.na(value))

RRS_comparative <- comparative %>% filter(network == 'RRS')
RRS_comparative <- RRS_comparative %>% filter(str_detect(name, 'Acc'))
RRS_comparative$name <- factor(RRS_comparative$name, levels=c('L.B.Acc.', 'Avg. Acc.', 'U.B.Acc.'))
RRS_comparative <- RRS_comparative %>% filter(!is.na(diff))

ggplot(RRS_comparative) +
    geom_point(aes(name, value), color="skyblue") +
    geom_line(aes(name, value, group=edge), color="skyblue") +
    geom_point(aes(name, diff)) +
    geom_line(aes(name, diff, group=edge)) +
facet_wrap(~edge) + theme_light()


SCRS_comparative <- comparative %>% filter(network == 'SCRS')
SCRS_comparative <- SCRS_comparative %>% filter(str_detect(name, 'Acc'))
SCRS_comparative$name <- factor(SCRS_comparative$name, levels=c('L.B.Acc.', 'Avg. Acc.', 'U.B.Acc.'))
SCRS_comparative <- SCRS_comparative %>% filter(!is.na(diff))

ggplot(SCRS_comparative) +
    geom_point(aes(name, value), color="skyblue") +
    geom_line(aes(name, value, group=edge), color="skyblue") +
    geom_point(aes(name, diff)) +
    geom_line(aes(name, diff, group=edge)) +
facet_wrap(~edge) + theme_light()
