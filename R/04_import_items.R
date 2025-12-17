library(tibble)
library(stringr)
library(igraph)
RRS_graph <- readRDS("data/RRS.RDS")
SCRS_graph <- readRDS("data/SCRS.RDS")

items_table <- tribble(
    ~Questionnaire, ~Scale, ~Item,
    "RRS (1)","Brooding", "Think 'What am I doing to deserve this?'",
    "RRS (2)","Reflection","Analyze recent events to try to understand why you are depressed",
    "RRS (3)","Brooding", "Think 'Why do I always react this way?'",
     "RRS (4)","Reflection", "Go away by yourself and think about why you feel this way",
    "RRS (5)","Reflection", "Write down what you are thinking and analyze it",
    "RRS (6)","Brooding", "Think about a recent situation, wishing it had gone better",
    "RRS (7)","Brooding", "Think 'Why do I have problems other people don’t have?'",
    "RRS (8)","Brooding", "Think 'Why can’t I handle things better?'",
    "RRS (9)","Reflection", "Analyze your personality to try to understand why you are depressed",
    "RRS (10)","Reflection", "Go someplace alone to think about your feelings",
    "SCRS (1)","Self-critical", "My attention is often focused on aspects of myself that I’m ashamed of.",
 "SCRS (2)","Self-critical", "I always seem to be rehashing in my mind stupid things that I’ve said or done.",
"SCRS (3)","Self-critical", "Sometimes it is hard for me to shut off critical thoughts about myself.",
 "SCRS (4)","Self-critical", "I can’t stop thinking about how I should have acted differently in certain situations.",
 "SCRS (5)","Self-critical", "I spend a lot of time thinking about how ashamed I am of some of my personal habits.",
 "SCRS (6)","Self-critical", "I criticize myself a lot for how I act around other people.",
 "SCRS (7)","Self-critical", "I wish I spent less time criticizing myself.",
 "SCRS (8)","Self-critical", "I often worry about all of the mistakes I have made.",
 "SCRS (9)","Self-critical", "I spend a lot of time wishing I were different.",
 "SCRS (10)","Self-critical", "I often berate myself for not being as productive as I should be.")

item_codes <- c(V(RRS_graph)$label)
for (i in 1:10) {
    node_label <- paste0("SCRS_", i)
    item_codes <- c(item_codes, V(SCRS_graph)[node_label]$label)
}
items_table$Label  <- item_codes

lookUpItemName <- function(itemName) {
    items_table[items_table$Label == itemName,]$Item
}

makeItemStats <- function(itemName, summaryTable, stats=c('EI1', 'EI2','Bet.'), names=c('Exp. Inf.', 'Two-Step Exp. Inf.', 'Betweenness'), fullExplain=FALSE, onlyExplain=TRUE) {
    itemLabel = lookUpItemName((itemName))
    colMask = colnames(summaryTable) %in% stats
    colIndex = seq(1,length(colMask))
    colIndex = colIndex[colMask]
    colTable <- data.frame(index = stats, label = names)

    textSummary <- paste0("'",itemLabel,"'-"," [")

    statsRow = summaryTable[summaryTable$Label == itemName,]
    for (j in seq_along(colIndex)) {
        colName = colTable[colTable$index == colnames(summaryTable)[colIndex[j]],]
        if (j == length(colIndex)){
            if (fullExplain){
            textSummary <- paste0(textSummary, colName$index, "(", colName$label, ") = ", format(round(statsRow[1,colIndex[j]][[1]],2),nsmall=2), "")
            } else {
                if (onlyExplain) {
                    textSummary  <- paste0(textSummary, colName$label)
                } else {
                    textSummary  <- paste0(textSummary, colName$index)
                }
            textSummary <- paste0(textSummary, "=", format(round(statsRow[1,colIndex[j]][[1]],2),nsmall=2), "")
            }
        } else {
            if (fullExplain){
                textSummary <- paste0(textSummary, colName$index, "(", colName$label,") = ", format(round(statsRow[1,colIndex[j]][[1]],2),nsmall=2), "; ")
            } else {

                if (onlyExplain) {
                    textSummary  <- paste0(textSummary, colName$label)
                } else {
                    textSummary  <- paste0(textSummary, colName$index)
                }
                textSummary <- paste0(textSummary, "=", format(round(statsRow[1,colIndex[j]][[1]],2),nsmall=2), "; ")

            }
        }
    }
    textSummary <- paste0(textSummary, "]")
    textSummary
}

