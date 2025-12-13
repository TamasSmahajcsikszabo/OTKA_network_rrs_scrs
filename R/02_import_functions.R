library(ggraph)
library(networkAnalysisTools)
library(tidyverse)
library(ggrepel)


#' Constructor function for graphs
construct_graph <- function(data, savename, subscale_id, data_labels) {
    source_data <- data
    cormatrix <- cor_auto(source_data)
    ebic_data <- data.frame(qgraph::EBICglasso(cormatrix, n = nrow(source_data), threshold = TRUE))
    ebic_data[lower.tri(ebic_data, diag=TRUE)] <- NA
    ebic_data <- tibble(data.frame(ebic_data))
    ebic_data <- bind_cols(name = colnames(ebic_data), ebic_data)
    ebic_data <- ebic_data %>% pivot_longer(cols=-name, names_to="pair", values_to="r") %>% unique() %>% filter(!is.na(r)) %>% filter(r > 0 )
    graph <- graph_from_data_frame(ebic_data[,1:2], directed=FALSE)
    E(graph)$strength <-ebic_data[,3][[1]]
    E(graph)$weight <-ebic_data[,3][[1]]
    tool <- unname(unlist(sapply(names(V(graph)), function(x){ strsplit(x,split="_")[[1]][1]})))
    toolcolor <- tool
    toolcolor[toolcolor == "RRS"] <- "cornflowerblue"
    toolcolor[toolcolor == "SCRS"] <- "coral"
    V(graph)$tool <- tool
    tools <- unique(tool)
    tool_enum <- c()
    for (item in tool) {
        for (i in 1:length(tools)) {
            if (item == tools[i]) {
                tool_enum <- c(tool_enum, i)
            }
        }
    }
    V(graph)$tool_enum <- tool_enum
    V(graph)$toolcolor <- toolcolor
    labels <- c()
    for (node_id in V(graph)$name){
        label <- data_labels[data_labels$var == node_id, ]['labels'][[1]]
        labels <- c(labels, label)
    }
    V(graph)$label <- labels
    subscales <- c("brooding", "reflection", "self-critical")
    graph_subscale <- c()
    for (id in subscale_id) {
        graph_subscale <- c(graph_subscale, subscales[id])
    }
    V(graph)$subscale <- graph_subscale
    V(graph)$subscale_enum <- subscale_id

    saveRDS(graph, paste0("data/", savename))
    graph

    }

#' Pastes Alpha character
Alpha <- function() {
    '$\\alpha$'
}

#' Custom formatting round function
Round <- function(value) {
        return(format(round(as.numeric(value), 2), nsmall=2))
}


#' Standardize function
#' @param a numeric vector
#' @returns standardized numeric vector of a
standardize <- function(a) {
    est_mean <- mean(a, na.rm = TRUE)
    est_SD <- sd(a, na.rm = TRUE)
    (a - est_mean) / est_SD
}


#' Helper function to find item name of questionnaire
find_tool <- function(itemname, data_labels) {
    if (str_locate(itemname, " .+ ")[1, 1] > 0 & !is.na(str_locate(itemname, " .+ "))) {
        itemname <- str_sub(itemname, start = 1, end = str_locate(itemname, " .+ ")[1, 1] - 1)
    } else {
        itemname <- itemname
    }
    toolname <- data_labels[itemname == data_labels$labels, ]$var
    toolname <- str_sub(toolname, start = 1, end = str_locate(toolname, "_")[[1, 1]] - 1)
    toolname
}




get_factor_levels <- function() {
    itemLevels <- c(paste0('RRS_', seq(1,10)), paste0('SCRS_', seq(1,10)))
    itemLabels <- c(paste0('RRS (', seq(1,10), ')'), paste0('SCRS (', seq(1,10), ')'))
    scaleLevels <- c('brooding', 'reflection', 'self-critical')
    scaleLabels <- c('Brooding', 'Reflection', 'Self-critism')
    list(
        "scaleLevels" = scaleLevels,
        "scaleLabels" = scaleLabels,
        "itemLevels" = itemLevels,
        "itemLabels" = itemLabels
    )
}

network_summary <- function(graph, dec = 2, name = "Graph", single_scale = FALSE, seed = 1234, method = "optimal", weights = TRUE, item_TPR) {
    summary <- matrix(nrow = vcount(graph))
    summary <- data.frame(summary)
    summary[, 1] <- tibble("Item" = V(graph)$name)
    colnames(summary) <- "Item"
    summary["Scale"] <- tibble("subscale" = V(graph)$subscale)
    summary["Label"] <- tibble("label" = V(graph)$label)
    summary["Deg."] <- tibble("degree" = igraph::degree(graph))
    summary["Str."] <- tibble("strength" = round(igraph::strength(graph), dec))
    summary["Bet."] <- tibble("betweenness" = round(igraph::betweenness(graph), dec))
    summary["Clo."] <- tibble("closeness" = round(igraph::closeness(graph), dec))
    influence_df <- networktools::expectedInf(graph, step = c("both"), directed = FALSE)
    summary["EI1"] <- tibble("EI1" = round(standardize(influence_df$step1), dec))
    summary["EI2"] <- tibble("EI2" = round(standardize(influence_df$step2), dec))
    bridge <- get_bridge_estimate(graph, seed = seed, method = method, weights = weights) %>% as.data.frame()
    names(bridge) <- c("Br.Str.", "Br.Bet.", "Br.Cl.", "Br.EI1", "Br.EI2", "Community")
    summary <- bind_cols(summary, bridge)
    summary["Title"] <- c(name, rep("", vcount(graph) - 1))
    summary["TPR"] <- round(V(add_community_certainty(graph, item_TPR))$TPR,3)



    summary  <-  summary %>% dplyr::select("Graph"="Title", "Item", "Scale", "Label", everything(), -"Community")


    labels <- get_factor_levels()
    summary$Item  <-  factor(summary$Item, levels=labels['itemLevels'][[1]], labels=labels['itemLabels'][[1]])
    summary$Scale  <-  factor(summary$Scale, levels=labels['scaleLevels'][[1]], labels=labels['scaleLabels'][[1]])
    summary  <-  summary %>% arrange(Scale)
    summary

}


get_original_communities <- function(graph) {
    true_table <- tibble(
        Item = V(graph)$name,
        `Community (True)` = V(graph)$subscale_enum
    )
    original_communities <- list()
    for (community in unique(true_table$`Community (True)`)) {
        original_community <- list(true_table[true_table$`Community (True)` == community, ]$Item)
        original_communities[community] <- original_community
    }
    original_communities
}

get_TPR <- function(simulatedCommunities, graph) {
    true_table <- tibble(
        Item = V(graph)$name,
        `Community (True)` = V(graph)$subscale_enum
    )
    original_communities <- list()
    for (community in unique(true_table$`Community (True)`)) {
        original_community <- list(true_table[true_table$`Community (True)` == community, ]$Item)
        original_communities[community] <- original_community
    }
    simulatedCommunities <- simulatedCommunities %>%
        group_by(Item, Method, Weights, Community) %>%
        summarise(N = n()) %>%
        ungroup()
    # group_by(Item, Method, Weights) %>%
    # filter(N == max(N))

    best_found_communities <- list()
    for (method in c("optimal", "spinglass", "walktrap", "fast_greedy")) {
        for (weight in c(TRUE, FALSE)) {
            slice <- simulatedCommunities[simulatedCommunities$Method == method & simulatedCommunities$Weights == weight, ]
            for (community in as.numeric(unique(slice$Community))) {
                found_community <- list(unique(slice[slice$Community == community, ]$Item))

                label <- paste0(method, ", ", weight)
                names(found_community) <- label
                best_found_communities <- append(best_found_communities, found_community)
            }
        }
    }

    results <- tibble()
    subscales <- data.frame(item = V(graph)$name, subscale = V(graph)$subscale)

    for (original_community in original_communities) {
        if (!is.null(original_community)) {
        for (i in seq_along(best_found_communities)) {
            match <- length(intersect(original_community, unlist(best_found_communities[i])))
            label <- strsplit(names(best_found_communities)[i], ", ")[[1]]
            method <- label[1]
            weight <- label[2]
            subscale <- suppressMessages(unique(data.frame(item = original_community) %>% left_join(subscales) %>% dplyr::select(subscale))[[1]])
            if (match > 0) {
                result <- tibble(
                    Method = method,
                    Weights = weight,
                    Subscale = subscale,
                    `# of True community membership` = match,
                    `TPR` = match / length(unlist(original_community)),
                    `Original Subscale` = paste0(original_community, collapse = ", "),
                    `Found Community` = paste0(unlist(best_found_communities[i]), collapse = ", ")
                )
                results <- bind_rows(results, result)
            }
        }
    }
    }
    results <- unique(results)
    results
    results %>%
        group_by(Method, `Original Subscale`) %>%
        filter(TPR == max(TPR))
}



get_match <- function(input) {
    columns <- strsplit(input, ",")
    tibble("Item" = unlist(columns))
}


get_rate <- function(input) {
    values <- strsplit(input, ",")[[1]]
    tibble(
        Item = values,
        Value = rep(1, length(values))
    )
}

item_level_TPR <- function(simulatedCommunities, graph, threshold = 0.85) {
    true_table <- tibble(
        Item = V(graph)$name,
        `Community (True)` = V(graph)$subscale_enum
    )
    original_communities <- list()
    for (community in unique(true_table$`Community (True)`)) {
        original_community <- list(true_table[true_table$`Community (True)` == community, ]$Item)
        original_communities[community] <- original_community
    }
    matching <- tibble()
    rates <- tibble(Item = V(graph)$name, "TP" = 0, "FP" = 0)
    for (method in c("spinglass", "optimal", "fast_greedy", "walktrap")) {
        for (weight in c(TRUE, FALSE)) {
            for (iteration in unique(simulatedCommunities$Iteration)) {
                simulated_community <- simulatedCommunities %>%
                    filter(
                        Method == method,
                        Weights == weight,
                        Iteration == iteration
                    )

                communities <- get_community(simulated_community)
                match <- found_pattern(original_communities, communities)
                matching <- bind_rows(matching, match)
            }
        }
    }
    filtered_matching <- matching %>% filter(match_rate >= threshold)
    total <- nrow(filtered_matching)
    TP <- map_dfr(filtered_matching["TP"][[1]], function(x) {
        get_rate(x)
    })
    TP <- as.data.frame(TP)
    names(TP) <- c("Item", "TP")
    TP <- as_tibble(TP)
    TP <- TP %>%
        group_by(Item) %>%
        summarise(TP = sum(TP, na.rm = TRUE)) %>%
        ungroup()
    FP <- map_dfr(filtered_matching["FP"][[1]], function(x) {
        get_rate(x)
    })
    FP <- as.data.frame(FP)
    names(FP) <- c("Item", "FP")
    FP <- as_tibble(FP)
    FP <- FP %>%
        group_by(Item) %>%
        summarise(FP = sum(FP, na.rm = TRUE)) %>%
        ungroup()
    rates <- bind_rows(rates, FP)
    rates <- bind_rows(rates, TP)
    per_item <- rates %>%
        group_by(Item) %>%
        summarise(
            TP = sum(TP, na.rm = TRUE),
            FP = sum(FP, na.rm = TRUE)
        ) %>%
        rowwise() %>%
        mutate(Total = sum(TP, FP)) %>%
        ungroup() %>%
        mutate(
            TP = TP / Total,
            FP = FP / Total
        )
    per_item
}


add_label <- function(graphPlot, rowVector, styleVector) {
    vjust <- 1.4
    for (i in 1:length(rowVector)) {
        graphPlot <- graphPlot + geom_text_node(aes(label=rowVector[i]), fontface = styleVector[i], vjust=vjust+0.2)

    }
    graphPlot
}




table_nums <- captioner::captioner(prefix = "Tab.")

f.ref <- function(x) {
  stringr::str_extract(table_nums(x), "[^:]*")
}



plot_eigen <- function(evecs, postmeans, graph, namevecs, subscale, labels,textsize=14) {
    plot_data <- data.frame(evecs)
    colnames(plot_data) <- c("dim1", "dim2")
    colors <- data.frame(tool=namevecs)
    plot_data <- bind_cols(plot_data, colors)
    item_name <- namevecs
    plot_data['name'] <- item_name
    plot_data['Sub-scale'] <- subscale
    names(plot_data)[3] <- 'Questionnaire'
    territory <- plot_data %>%
        group_by(Questionnaire) %>%
        slice(chull(dim1, dim2))

    custom_colors <- tibble('Sub-scale'=c('Brooding', 'Reflection', 'Self-critical'), color=c('white', '#CA382A', '#0C38A0'))
    my_color_scale <- plot_data %>%  left_join(custom_colors)
    my_color_scale <- as.character(my_color_scale$color)
    names(my_color_scale) <- plot_data$`Sub-scale`

    ggplot() +
        geom_polygon(data=territory, aes(dim1, dim2, fill=Questionnaire), alpha=1/3,color='black', show.legend=TRUE) +
        geom_point(data=plot_data, aes(dim1, dim2, shape=`Sub-scale`),size=6,alpha=1/2,color='black',fill='white') +
        geom_point(data=plot_data, aes(dim1, dim2, color=`Sub-scale`, shape=`Sub-scale`),size=5) +
        geom_text_repel(data=plot_data,aes(dim1, dim2, label=labels), size=textsize/3)+
        theme_bw() +
        scale_color_manual(values=my_color_scale) +
        scale_fill_manual(values=c("#CA382A", "#0C38A0")) +
        labs(x=paste0('1. Eigen Vector [m=', round(postmeans[1],3), "]"), y=paste0('2. Eigen Vector [m=', round(postmeans[2],3), ']')) +
        labs(caption='* posterior means of eigenvalues are in [] after axis labels')+
        theme(legend.position = "bottom",
            text=element_text(size=textsize))
}

