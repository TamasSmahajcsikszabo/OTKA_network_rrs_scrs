standardize <- function(a) {
    est_mean <- mean(a, na.rm = TRUE)
    est_SD <- sd(a, na.rm = TRUE)

    (a - est_mean) / est_SD
}

tau_estimate <- function(x, y = NULL) {
    if (is.list(x)) {
        input <- x
        x <- input[[1]]
        y <- input[[2]]
    }

    concordants <- c()
    discordants <- c()
    total <- c(0)
    max_op <- length(x) * (length(y) - 1)
    for (index in seq(1, length(x))) {
        for (pair_index in seq(1, length(x))[!seq(1, length(x)) == index]) {
            first <- list(x[index], y[index])
            second <- list(x[pair_index], y[pair_index])
            if ((first[[1]] > second[[1]] & first[[2]] > second[[2]]) |
                (first[[1]] < second[[1]] & first[[2]] < second[[2]])) {
                concordants <- c(concordants, 1)
            } else {
                discordants <- c(discordants, -1)
            }
            total <- total + 1
            cat(paste0("\r", "Progress: ", round(total / max_op, 2) * 100, "%"))
        }
    }
    tau <- (sum(concordants) + sum(discordants)) / total
    if (is.nan(tau) & (length(x) == 1 | length(y) == 1)) {
        warning("Estimation not possible; input values are scalars")
    } else {
        cat(paste0("\n", "Tau estimate = ", tau))
    }
}

get_item_names <- function(covMatrix, labels) {
    new_labels <- c()
    original_names <- unlist(lapply(rownames(covMatrix), function(x) {
        str_replace(x, "PCA.", "")
    }))
    for (i in seq_along(original_names)) {
        if (str_detect(original_names[i], "\\.")) {
            split_labels <- unlist(strsplit(original_names[i], "\\."))
            lookup_names <- paste0(unlist(lapply(split_labels, function(x) {
                labels[names(labels) == x]
            })), collapse = " + ")
            split_labels <- paste0(split_labels, collapse = " + ")
            new_labels[i] <- paste0(split_labels, ": \n", lookup_names)
        } else {
            new_labels[i] <- paste0(original_names[i], ": \n", labels[names(labels) == original_names[i]])
        }
    }
    new_labels
}

make_community_graph <- function(graph, covMatrix, labels) {
    igraph_converted <- as.igraph(graph, attributes = TRUE)
    group_estimation <- cluster_spinglass(
        igraph_converted,
        weights = NULL,
        vertex = NULL,
        spins = 25,
        parupdate = FALSE,
        start.temp = 1,
        stop.temp = 0.01,
        cool.fact = 0.99,
        update.rule = c("config", "random", "simple"),
        gamma = 0.5,
        implementation = c("orig", "neg"),
        gamma.minus = 1
    )


    grouping_order <- data.frame(id = group_estimation$membership) %>%
        left_join(global_colors) %>%
        select(color) %>%
        unlist()
    labels <- get_item_names(covMatrix, labels)
    qgraph(covMatrix, graph = "glasso", tuning = 0.5, layout = "spring", sampleSize = 888, theme = "TeamFortress", details = TRUE, threshold = FALSE, color = grouping_order, labels = labels)
}

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
expected_inf_plot <- function(network, toolname) {
    influence_df <- networktools::expectedInf(network, step = c("both"), directed = FALSE)
    step1 <- unlist(influence_df["step1"])
    step2 <- unlist(influence_df["step2"])

    # z-scores
    step1_mean <- mean(step1, na.rm = TRUE)
    step1_sd <- sd(step1, na.rm = TRUE)
    step1 <- (step1 - step1_mean) / step1_sd
    step2_mean <- mean(step2, na.rm = TRUE)
    step2_sd <- sd(step2, na.rm = TRUE)
    step2 <- (step2 - step2_mean) / step2_sd

    varnames <- names(step1)
    varnames <- stringr::str_replace(varnames, "step1.", "")
    for (i in seq_along(varnames)) {
        tool <- find_tool(varnames[i], data_labels)
        varnames[i] <- paste0(tool, ": ", varnames[i])
    }
    influence_tibble <- tibble(
        OneStep = step1,
        TwoStep = step2,
        item = factor(varnames, levels = varnames)
    ) %>%
        gather(1:2, key = Method, value = inf)

    influence_tibble %>%
        ggplot(aes(item, inf, group = Method, color = Method)) +
        geom_point(shape = 0, size = 2) +
        geom_path() +
        coord_flip() +
        scale_color_manual(values = global_colors$color[1:2]) +
        theme_light() +
        scale_x_discrete(limits = rev(levels(influence_tibble$item))) +
        labs(
            title = paste0("One- and Two-Step Expected Influence estimates for ", toolname),
            x = "Items",
            y = "Expected Influence (z-scores)"
        ) +
        theme(
            plot.title = element_text(size = 11)
        )
}
make_bridge_plot <- function(network, toolname) {
    igraph_converted <- as.igraph(network, attributes = TRUE)
    group_estimation <- cluster_spinglass(
        igraph_converted,
        weights = NULL,
        vertex = NULL,
        spins = 25,
        parupdate = FALSE,
        start.temp = 1,
        stop.temp = 0.01,
        cool.fact = 0.99,
        update.rule = c("config", "random", "simple"),
        gamma = 0.5,
        implementation = c("orig", "neg"),
        gamma.minus = 1
    )

    #    community_vector <- c()
    #    for (comm in seq_along(group_estimation[1:2])) {
    #        comm_items <- c()
    #        for (item in unlist(group_estimation[comm])) {
    #            item <- c(item)
    #            names(item) <- comm
    #            comm_items <- c(comm_items, item)
    #        }
    #        community_vector <- c(community_vector, comm_items)
    #    }
    if (length(unique(group_estimation$membership)) > 1) {
        bridge_estimate <- networktools::bridge(network, communities = group_estimation$membership)
    } else {
        bridge_estimate <- networktools::bridge(network)
    }

    # standardize scores into a tibble
    output_tibble <- tibble(.rows = length(bridge_estimate[[1]]))

    for (i in seq(1, 5)) {
        target_vector <- tibble(standardize(bridge_estimate[[i]]))
        colnames(target_vector) <- names(bridge_estimate)[i]
        output_tibble <- bind_cols(output_tibble, target_vector)
    }
    # add communities
    output_tibble <- bind_cols(output_tibble, bridge_estimate[6])
    output_tibble <- output_tibble[, names(bridge_estimate)]
    output_tibble <- output_tibble %>%
        mutate(item = names(bridge_estimate[[1]]))
    output_tibble <- output_tibble %>%
        mutate(Community = as_factor(communities))

    # reshape
    reshaped_output <- output_tibble %>%
        gather(1:5, key = measure, value = value)

    for (r in seq(1, nrow(reshaped_output))) {
        reshaped_output[r, "var"] <- find_tool(reshaped_output[r, "item"][[1]], data_labels)
        reshaped_output[r, "longitem"] <- paste0(reshaped_output[r, "var"], ": ", reshaped_output[r, "item"])
    }

    names(bridge_estimate[[1]]) <- unique(reshaped_output$longitem)



    # create plot

    reshaped_output %>%
        ggplot(aes(longitem, value)) +
        geom_path(aes(group = measure), color = "cornflowerblue") +
        geom_point(size = 4, alpha = 1 / 2) +
        coord_flip() +
        scale_color_manual(values = global_colors$color) +
        theme_light() +
        scale_x_discrete(limits = rev(levels(as_factor(names(bridge_estimate[[1]]))))) +
        facet_wrap(~measure, ncol = 5) +
        theme(
            legend.position = "bottom"
        ) +
        labs(
            x = "Item",
            y = "Z-scores",
            title = paste0("Bridge Expected Influence Estimates for ", toolname)
        ) +
        theme(
            plot.title = element_text(size = 11)
        )
}


network_accuracy_plot <- function(accuracy, tool = "", multiple = FALSE, save_data = FALSE) {
    sample_data <- accuracy$sampleTable %>%
        filter(nchar(id) > 7) %>%
        arrange(desc(value))

    sample_data <- sample_data %>%
        left_join(data_labels, by = c("node1" = "var")) %>%
        left_join(data_labels, by = c("node2" = "var")) %>%
        rename(
            new_node_1 = labels.x,
            new_node_2 = labels.y
        ) %>%
        filter(!is.na(new_node_2)) %>%
        mutate(
            tool1 = str_sub(node1, 1, str_locate(node1, "_")[, 1] - 1),
            tool2 = str_sub(node2, 1, str_locate(node2, "_")[, 1] - 1)
        ) %>%
        mutate(id = if_else(tool1 == tool2,
            paste0(tool1, ": ", new_node_1, " - ", new_node_2),
            paste0(tool1, ": ", new_node_1, " - ", tool2, ": ", new_node_2)
        ))

    bootstrap_data_aggr <- accuracy$bootTable %>%
        left_join(data_labels, by = c("node1" = "var")) %>%
        left_join(data_labels, by = c("node2" = "var")) %>%
        rename(
            new_node_1 = labels.x,
            new_node_2 = labels.y
        ) %>%
        filter(!is.na(new_node_2)) %>%
        mutate(
            tool1 = str_sub(node1, 1, str_locate(node1, "_")[, 1] - 1),
            tool2 = str_sub(node2, 1, str_locate(node2, "_")[, 1] - 1)
        ) %>%
        mutate(id = if_else(tool1 == tool2,
            paste0(tool1, ": ", new_node_1, " - ", new_node_2),
            paste0(tool1, ": ", new_node_1, " - ", tool2, ": ", new_node_2)
        )) %>%
        filter(nchar(id) > 7) %>%
        arrange(desc(value)) %>%
        group_by(id) %>%
        summarize(
            m = mean(value, na.rm = TRUE),
            lower = m - 1.96 * sd(value, na.rm = TRUE),
            upper = m + 1.96 * sd(value, na.rm = TRUE)
        )



    bootstrap_data <- accuracy$bootTable %>%
        filter(nchar(id) > 7) %>%
        arrange(desc(value))

    if (save_data) {
        saveRDS(sample_data, "../output/sample_data.RDS")
        saveRDS(bootstrap_data, "../output/bootstrap_data.RDS")
        saveRDS(bootstrap_data_aggr, "../output/bootstrap_data_aggr.RDS")
    }

    ggplot() +
        geom_vline(aes(xintercept = 0), linetype = "dotted", size = 0.75) +
        geom_point(data = bootstrap_data_aggr, aes(m, reorder(id, m)), shape = 22, fill = "cornflowerblue", size = 2) +
        geom_errorbar(data = bootstrap_data_aggr, aes(m, id, xmin = lower, xmax = upper), color = "cornflowerblue") +
        geom_point(data = sample_data, aes(value, reorder(id, value)), shape = 21, fill = "coral", size = 2) +
        theme_light() +
        labs(
            title = paste0("Edge Weight Accuracy Estimates with Bootstrap for ", tool),
            x = "Edge Weight",
            y = "Edge",
            caption = "DOT - the sample edge weight, SQUARE and confidence interval - the bootstrap estimates for edge weights (mean and 95% CI)"
        ) +
        scale_x_continuous(breaks = seq(-0.2, 0.8, 0.05))
}

network_stability_estimator <- function(network) {
    B_levels <- seq(1, 2500)
    sample_data <- t(matrix(network$sampleTable))
    names(sample_data) <- names(network$sampleTable)
    boot <- network$bootTable %>%
        mutate(sample_rate = nPerson / mean(sample_data$nPerson)) %>%
        mutate(
            row = row_number(),
            bootstrap_indicator = (row %/% 75) + 1,
            bootstrap_indicator = factor(bootstrap_indicator, levels = B_levels)
        )

    boot_data <- t(matrix(boot))
    names(boot_data) <- names(boot)

    drop_levels <- unique(boot_data$nPerson)
    measures <- unique(sample_data$type)
    total_operations <- 2500 * length(measures) * length(drop_levels)
    plot_data <- matrix(ncol = 4, nrow = 0)
    colnames(plot_data) <- c("b", "measure", "level", "correlation")
    i <- 0

    for (b in B_levels) {
        for (level in drop_levels) {
            for (m in measures) {
                i <- i + 1

                boot_sample_bootstrap_mask <- unlist(boot_data["bootstrap_indicator"], use.names = FALSE) == b
                boot_sample_level_mask <- unlist(boot_data["nPerson"], use.names = FALSE) == level
                boot_sample_measure_mask <- unlist(boot_data["type"], use.names = FALSE) == m
                boot_mask <- boot_sample_bootstrap_mask & boot_sample_level_mask & boot_sample_measure_mask

                boot_sample <- unlist(boot_data["value"], use.names = FALSE)[boot_mask]


                original_sample_mask <- unlist(sample_data["type"], use.names = FALSE) == m
                original_sample <- unlist(sample_data["value"], use.names = FALSE)[original_sample_mask]

                test <- all(!is.null(boot_sample), !is.null(original_sample), length(original_sample) == length(boot_sample))
                if (test) {
                    correlation <- cor(boot_sample, original_sample)
                    new_record <- c(
                        "b" = b,
                        "measure" = m,
                        "level" = level,
                        "correlation" = correlation
                    )
                    plot_data <- rbind(plot_data, new_record)
                    cat(paste0("\rProgress: ", round((i / total_operations) * 100, 2), "%"))
                }
            }
        }
    }
    plot_df <- tibble(as.data.frame(plot_data))
    plot_df <- plot_df %>%
        mutate(
            level = as.numeric(level),
            correlation = as.numeric(correlation)
        ) %>%
        group_by(b, measure) %>%
        summarise(
            r = mean(correlation, na.rm = TRUE),
            sample_rate = level / max(unlist(sample_data["nPerson"], use.names = FALSE))
        )

    saveRDS(plot_df, paste0("../output/plot_", deparse(quote(network)), ".RDS"))
}

network_stability_plot <- function(stability_estimate) {
    stability_estimate %>%
        mutate(b = as.numeric(b)) %>%
        group_by(sample_rate, measure) %>%
        mutate(mean_r = mean(r, na.rm = TRUE)) %>%
        filter(!is.na(r)) %>%
        filter(!is.na(sample_rate)) %>%
        ggplot(aes(reorder(percent(sample_rate), desc(sample_rate)), r, group = measure)) +
        geom_point(aes(y = mean_r, color = measure), size = 3) +
        geom_line(aes(y = mean_r, color = measure), size = 1) +
        geom_jitter(aes(color = measure), alpha = 1 / 5, size = 0.5) +
        geom_hline(aes(yintercept = 0)) +
        ylim(c(-1, 1)) +
        theme_light() +
        scale_color_manual(values = c("coral4", "coral", "cornflowerblue", "darkturquoise")) +
        scale_fill_manual(values = c("coral4", "coral", "cornflowerblue", "darkturquoise")) +
        labs(
            title = "Stability estimates of centrality measures with bootstrap",
            x = "Portion of original sample",
            y = "Average correlation with the original sample",
            fill = "Centrality measure",
            color = "Centrality measure",
            caption = "Large points are average correlation estimates, jittered points are actual bootstrap estimates"
        ) +
        # scale_x_continuous(breaks = seq(1.0, 0.0, -0.1)) +
        scale_y_continuous(breaks = seq(-1.0, 1.0, 0.1))
}


to_vector <- function(input) {
    if (typeof(input) == "list") {
        output <- unlist(input)
    } else {
        output <- input
    }
    output
}

find_subscale <- function(original, predicted) {
    original <- to_vector(original)
    predicted <- to_vector(predicted)
    common <- intersect(original, predicted)
    match_rate <- length(common) / length(original)
    if (identical(common, original)) {
        output <- list(
            "match" = TRUE, "FP" = predicted[!predicted %in% common],
            "TP" = common, "original" = original, "predicted" = predicted, "match_rate" = match_rate
        )
    } else {
        output <- list(
            "match" = FALSE, "FP" = predicted[!predicted %in% common],
            "TP" = common, "original" = original, "predicted" = predicted, "match_rate" = match_rate
        )
    }

    output
}

get_bridge_estimate <- function(network, dec = 2, seed = 1234, method = "fast_greedy", weights = TRUE) {
    if (!is.null(seed)) {
        set.seed(seed)
    }
    if (weights) {
        W <- E(network)$weight
    } else {
        W <- NULL
    }
    if (method == "spinglass") {
        group_estimation <- cluster_spinglass(
            network,
            weights = W,
            vertex = NULL,
            spins = 25,
            parupdate = FALSE,
            start.temp = 1,
            stop.temp = 0.01,
            cool.fact = 0.99,
            update.rule = c("config", "random", "simple"),
            gamma = 0.5,
            implementation = c("orig", "neg"),
            gamma.minus = 1
        )
    } else if (method == "fast_greedy") {
        group_estimation <- cluster_fast_greedy(
            network,
            weights = W
        )
    } else if (method == "walktrap") {
        group_estimation <- cluster_walktrap(network,
            weights = W
        )
    } else if (method == "optimal") {
        group_estimation <- cluster_optimal(network,
            weights = W
        )
    }
    if (length(unique(group_estimation$membership)) > 1) {
        bridge_estimate <- networktools::bridge(network, communities = group_estimation$membership)
    } else {
        bridge_estimate <- networktools::bridge(network)
    }

    itemrownames <- names(bridge_estimate[1][[1]])

    # standardize scores into a tibble
    output_tibble <- tibble(.rows = length(bridge_estimate[[1]]))

    for (i in seq(1, 5)) {
        target_vector <- tibble(round(standardize(bridge_estimate[[i]]), dec))
        colnames(target_vector) <- names(bridge_estimate)[i]
        output_tibble <- bind_cols(output_tibble, target_vector)
    }
    # add communities
    output_tibble <- bind_cols(output_tibble, bridge_estimate[6])
    output_tibble <- output_tibble[, names(bridge_estimate)]
    output_tibble <- output_tibble %>%
        mutate(Community = as_factor(communities)) %>%
        dplyr::select(-communities)
    # mutate(Item = itemrownames)
    # true_table <- tibble(
    #     Item = V(network)$name,
    #     `Community (True)` = V(network)$subscale_enum
    # )
    # output_tibble <- output_tibble %>%
    #     left_join(true_table) %>%
    #     mutate(TP = Community == `Community (True)`)

    # as.data.frame(output_tibble)
    output_tibble
}
network_summary <- function(graph, dec = 2, name = "Graph", single_scale = FALSE, seed = 1234, method = "optimal", weights = TRUE, item_TPR) {
    summary <- matrix(nrow = vcount(graph))
    summary <- data.frame(summary)
    summary[, 1] <- tibble("Item" = V(graph)$name)
    colnames(summary) <- "Item"
    summary["Scale"] <- tibble("subscale" = V(graph)$subscale)
    summary["Label"] <- tibble("label" = V(graph)$label)
    summary["Deg."] <- tibble("degree" = degree(graph))
    summary["Str."] <- tibble("strength" = round(strength(graph), dec))
    summary["Bet."] <- tibble("betweenness" = round(betweenness(graph), dec))
    summary["Clo."] <- tibble("closeness" = round(closeness(graph), dec))
    influence_df <- networktools::expectedInf(graph, step = c("both"), directed = FALSE)
    summary["EI1"] <- tibble("EI1" = round(standardize(influence_df$step1), dec))
    summary["EI2"] <- tibble("EI2" = round(standardize(influence_df$step2), dec))
    bridge <- get_bridge_estimate(graph, seed = seed, method = method, weights = weights) %>% as.data.frame()
    names(bridge) <- c("Br.Str.", "Br.Bet.", "Br.Cl.", "Br.EI1", "Br.EI2", "Community")
    summary <- bind_cols(summary, bridge)
    summary["Title"] <- c(name, rep("", vcount(graph) - 1))
    summary["TPR"] <- round(V(add_community_certainty(graph, item_TPR))$TPR,3)

    summary %>% select("Title", "Item", "Scale", "Label", everything(), -"Community")
}

degree_distribution_summary <- function(graph) {
    tibble(
        "Degree" = 0:max(degree(graph), na.rm = TRUE),
        "Degree dist." = degree_distribution(graph, mode = "all")
    )
}

simulate_communities <- function(graph, i = 1000, path = "simulate.RDS") {
    results <- tibble()
    totalN <- 4 * 2 * i

    for (method in c("optimal", "spinglass", "walktrap", "fast_greedy")) {
        for (w in c(TRUE, FALSE)) {
              for (iteration in 1:i) {
                  estimate <- get_bridge_estimate(graph, seed = NULL, method = method, weights = w)
                  estimate["Item"] <- V(graph)$name
                  estimate <- estimate[c("Item", "Community")]
                  estimate["Iteration"] <- iteration
                  estimate["Method"] <- method
                  estimate["Weights"] <- w
                  results <- bind_rows(results, estimate)
                  progress <- paste0("\r", round(((nrow(results) / vcount(graph)) / totalN) * 100), "%")
                  cat(progress)
              }
          }
    }
    saveRDS(results, path)
}

# simulatedCommunities <- RRS_simulated_communities
# graph <- rumi
# weight = TRUE
# method='fast_greedy'
# community <- 1
# a <- c(1,2,3,4,5,6)
# b <- c(2,3,4,5,6,7)
# length(intersect(a,b))/length(a)
# i <- 1
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
    results <- unique(results)
    results
    results %>%
        group_by(Method, `Original Subscale`) %>%
        filter(TPR == max(TPR))
}

add_articulation_point <- function(graph) {
    V(graph)[igraph::articulation_points(graph)]$articulation_point <- "*"
    V(graph)$articulation_point[is.na(V(graph)$articulation_point)] <- ""
    graph
}

add_network_descriptives <- function(graph) {
    o <- vcount(graph)
    s <- ecount(graph)
    t <- transitivity(graph)
    d <- diameter(graph)
    md <- mean_distance(graph)
    paste0("Order=", o, "; ", "Size=", s, "; ", "Clustering Coefficient=", round(t, 3), "; ", "Diameter=", round(d, 3), "; ", "Mean distance=", round(md, 3))
}

add_toolname <- function(graph) {
    V(graph)$tool <- unlist(lapply(strsplit(V(graph)$name, "_"), function(x) {
        x[1]
    }))
    V(graph)$item_number <- unlist(lapply(strsplit(V(graph)$name, "_"), function(x) {
        x[2]
    }))
    graph
}

# item_analysis <- function(simulated_community, graph) {
#     true_table <- tibble(
#             Item = V(graph)$name,
#             `Community (True)` = V(graph)$subscale_enum
#         )

#     simulated_community %>%
#         filter(Method %in% method) %>%
#         group_by(Item, Community, Method) %>%
#         summarise(N = n()) %>%
#         ungroup() %>%
#         group_by(Item) %>%
#         mutate(Total = sum(N)) %>%
#         ungroup() %>%
#         rowwise() %>%
#         mutate(Certainty = round(N / Total,3)) %>%
#         group_by(Item, Method) %>%
#         left_join(true_table) %>%
#         as.data.frame()
# }

# item_analysis_TPR <- function(simulatedCommunities, graph) {
#    true_table <- tibble(
#        Item = V(graph)$name,
#        `Community (True)` = V(graph)$subscale_enum
#    )
#    original_communities <- list()
#    for (community in unique(true_table$`Community (True)`)) {
#        original_community <- list(true_table[true_table$`Community (True)` == community, ]$Item)
#        original_communities[community] <- original_community
#    }

#    #TODO: special case if original scale is just one community!

#    results <- tibble()
#    for (iteration in unique(simulatedCommunities$Iteration)) {
#        iteration_data <- simulatedCommunities[simulatedCommunities$Iteration==i,]

#        for (weight in c(TRUE, FALSE)) {
#            for (method in unique(simulatedCommunities$Method)){
#                iteration_data_subset <- iteration_data |>
#                    filter(Method==method, Weights==weight)
#                for (community in unique(iteration_data_subset$Community)){
#                    subscale <- iteration_data_subset[iteration_data_subset$Community==community,]
#                    subscale <- subscale$Item
#                    for (original in original_communities){
#                        match <- find_subscale(original, subscale)
#                        if (match['match'][[1]]) {
#                            for (item in match['TP']){
#                                res <- tibble('Method'=method,
#                                            'Weights'=weight,
#                                            'Item'=item,
#                                            'Community'=community,
#                                            'Iteration'=iteration,
#                                            'TP'=1,
#                                            'FP'=0)
#                                results <- bind_rows(results, res)
#                            }
#                            for (item in match['FP']) {
#                                res <- tibble('Method'=method,
#                                            'Weights'=weight,
#                                            'Item'=item,
#                                            'Community'=community,
#                                            'Iteration'=iteration,
#                                            'TP'=0,
#                                            'FP'=1)
#                                results <- bind_rows(results, res)
#                            }
#                        }
#                    }
#                }
#            }
#        }
#    }
#    results
get_community <- function(dataset) {
    communities <- unique(dataset$Community)
    lapply(communities, function(x) {
        dataset[dataset$Community == x, ]$Item
    })
}

found_pattern <- function(original_communities, communities) {
    results <- tibble()
    for (original in original_communities) {
        for (community in communities) {
            res <- find_subscale(original, community)
            result <- tibble(
                match = res["match"][[1]],
                FP = paste0(res["FP"][[1]], collapse = ","),
                TP = paste0(res["TP"][[1]], collapse = ","),
                original = paste0(res["original"][[1]], collapse = ","),
                predicted = paste0(res["predicted"][[1]], collapse = ","),
                match_rate = res["match_rate"][[1]]
            )
            results <- bind_rows(results, result)
        }
    }
    results
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

add_community_certainty <- function(graph, item_TPR) {
    for (i in 1:vcount(graph)) {
        V(graph)[i]$TPR <- item_TPR %>%
            filter(Item == V(graph)[i]$name) %>%
            ungroup() %>%
            dplyr::select(TP) %>%
            unlist()
    }
    graph
}
beautify <- function(graph, simulated_community, title = "Graph", no_caption = FALSE, force_caption = FALSE, textsize = 6, overalltextsize = 18, item_TDR) {
    require("ggraph")
    S <- nrow(simulated_community) / vcount(graph) / 8
    set.seed(42)
    caption <- paste0("TPR is average True Positive Rate with ", S, " times reruns of community detection")
    if (length(articulation_points(graph)) > 0) {
        caption <- paste0("* - articulation point (cut vertex; when removed disconnects the graph)", "\n", caption)
    }
    if (no_caption) {
        caption <- ""
    }
    if (force_caption) {
        caption <- paste0("TPR is average True Positive Rate with ", S, " times reruns of community detection")
        caption <- paste0("* marks Articulation Points (cut vertices; when such vertices are  removed disconnect the graph)", "\n", caption)
        caption <- paste0(caption, "\n Order is # of vertices; Size is # of edges")
        caption <- paste0(caption, "\n Edge width reflect edge weight (penalized part. corr.), while edge shade reflects lower bound of 95% CI of bootstrap accuracy estimate")
        caption <- paste0(caption, "\n Dashed edge line indicates the 95% CI of accuracy estimate ranges below 0.0")
    }
    custom_colors <- tibble("subscale" = c("brooding", "reflection", "self-critical"), color = c("white", "#CA382A", "#0C38A0"))
    my_color_scale <- tibble("subscale" = V(graph)$subscale) %>% left_join(custom_colors)
    my_color_scale <- as.character(my_color_scale$color)
    names(my_color_scale) <- V(graph)$subscale
    add_community_certainty(graph, item_TDR) %>%
        add_articulation_point() %>%
        add_toolname() %>%
        ggraph(layout = "fr") +
        geom_edge_density(edge_fill = "grey100") +
        # scale_edge_color_manual(values=c("grey40", "grey0"))+
        geom_edge_fan(aes(alpha = accuracy, width = weight, linetype = accuracy < 0), show.legend = FALSE) +
        scale_color_manual(values = my_color_scale, name = "Sub-scale") +
        geom_node_point(color = "black", size = 12) +
        geom_node_point(aes(color = subscale), size = 10) +
        geom_node_point(color = "white", size = 5) +
        geom_node_point(aes(alpha = TPR), size = 5) +
        # geom_node_label(aes(label = label), alpha=1/5, color='grey70', size = 5, vjust=-0.6) +
        geom_node_text(aes(label = label), size = textsize, vjust = -1.4, fontface = "bold") +
        # geom_node_label(aes(label = paste0(name, "-", subscale)),alpha=1/5, color='grey70', size = 5,vjust=-1.1) +
        geom_node_text(aes(label = paste0(tool, "(", item_number, ")", "-", subscale)), size = textsize, vjust = -2.6) +
        geom_node_text(aes(label = articulation_point), size = 13, hjust = -2.9, vjust = -1.0) +
        # scale_color_manual(values = c("white", "grey30", "grey60"), name = "Sub-scale") +
        labs(
            caption = caption,
            title = title,
            subtitle = add_network_descriptives(graph)
        ) +
        theme(
            legend.position = "right",
            text = element_text(size = overalltextsize),
            panel.background = element_rect(color = "black", fill = "white")
        )
}


network_accuracy_data <- function(accuracy, tool = "", multiple = FALSE, save_data = FALSE) {
    sample_data <- accuracy$sampleTable %>%
        filter(nchar(id) > 7) %>%
        arrange(desc(value))

    sample_data <- sample_data %>%
        left_join(data_labels, by = c("node1" = "var")) %>%
        left_join(data_labels, by = c("node2" = "var")) %>%
        rename(
            new_node_1 = labels.x,
            new_node_2 = labels.y
        ) %>%
        filter(!is.na(new_node_2)) %>%
        mutate(
            tool1 = str_sub(node1, 1, str_locate(node1, "_")[, 1] - 1),
            tool2 = str_sub(node2, 1, str_locate(node2, "_")[, 1] - 1)
        ) %>%
        mutate(id = if_else(tool1 == tool2,
            paste0(tool1, ": ", new_node_1, " - ", new_node_2),
            paste0(tool1, ": ", new_node_1, " - ", tool2, ": ", new_node_2)
        ))

    bootstrap_data_aggr <- accuracy$bootTable %>%
        # left_join(data_labels, by = c("node1" = "var")) %>%
        # left_join(data_labels, by = c("node2" = "var")) %>%
        # rename(
        #     new_node_1 = labels.x,
        #     new_node_2 = labels.y
        # ) %>%
        filter(!node2 == "") %>%
        # mutate(
        #     tool1 = str_sub(node1, 1, str_locate(node1, "_")[, 1] - 1),
        #     tool2 = str_sub(node2, 1, str_locate(node2, "_")[, 1] - 1)
        # ) %>%
        # mutate(id = if_else(tool1 == tool2,
        #     paste0(tool1, ": ", new_node_1, " - ", new_node_2),
        #     paste0(tool1, ": ", new_node_1, " - ", tool2, ": ", new_node_2)
        # )) %>%
        # filter(nchar(id) > 7) %>%
        # arrange(desc(value)) %>%
        group_by(node1, node2) %>%
        summarize(
            m = mean(value, na.rm = TRUE),
            lower = m - 1.96 * sd(value, na.rm = TRUE),
            upper = m + 1.96 * sd(value, na.rm = TRUE)
        )



    bootstrap_data <- accuracy$bootTable %>%
        filter(nchar(id) > 7) %>%
        arrange(desc(value))

    list(sample_data, bootstrap_data_aggr, bootstrap_data)
}

add_accuracy_to_graph <- function(graph, accuracy_data, metric = "lower") {
    namevector <- V(graph)$name
    for (i in 1:nrow(get_itemnames(namevector))) {
        pair <- get_itemnames(namevector)[i, ] %>%
            unlist() %>%
            as.character()
        accuracy <- accuracy_data[2][[1]] %>%
            ungroup() %>%
            filter((node1 == pair[1] & node2 == pair[2]) | (node1 == pair[2] & node2 == pair[1]))
        node1 <- accuracy$node1
        node2 <- accuracy$node2
        mask <- paste0(node1, "|", node2)
        test <- attributes(E(graph))$vnames == mask
        if (sum(test) > 0){
            E(graph)[test]$accuracy <- accuracy$lower
        }
        E(graph)$accuracy
    }
    graph
}


edge_summary <- function(graph, accuracy_data, stability_data, statistic = "edge") {
    # adds accuracy and strength data to a table

    result <- tibble()
    result <- accuracy_data[[1]] %>% dplyr::select("Spl.Acc." = value)
    result <- bind_cols(result, accuracy_data[[2]] %>% ungroup())
    linknames <- attributes(E(graph))[1]
    result <- result %>% mutate(link=paste0(node1,'|',node2)) %>% filter(link %in% unlist(linknames))  %>% dplyr::select(-link, -`Spl.Acc.`)
    colnames(result) <- c("Vx.1", "Vx.2", "Avg. Acc.", "L.B.Acc.", "U.B.Acc.")
    stability_cutpoints <- summary(stability_data$bootTable$nPerson)[c(1, 2, 3, 5, 6)] %>%
        unname()
    stability_cutpoints <- tibble(nPerson = stability_cutpoints) %>%
        mutate(nPerson = as.integer(nPerson)) %>%
        mutate(Perc = 1 - nPerson / max(nPerson))
    stability_aggregated <- tibble(stability_data$bootTable) %>%
        left_join(stability_cutpoints) %>%
        group_by(node1, node2, Perc, type) %>%
        summarise(stability = mean(value, na.rm = TRUE)) %>%
        filter(!is.na(node2)) %>%
        filter(!is.na(Perc)) %>%
        filter(!node2 == "") %>%
        ungroup() %>%
        mutate(`Miss.%` = paste0("Miss.%:", round(Perc, 1))) %>%
        dplyr::select(-Perc) %>%
        group_by(node1, node2, type) %>%
        spread(`Miss.%`, stability)

    colnames(stability_aggregated)[1:3] <- c("Vx.1", "Vx.2", "Statistic")

    stability_aggregated <- stability_aggregated %>% filter(Statistic == statistic)
    result <- result %>% left_join(stability_aggregated)
    result <- result %>%
        dplyr::select(-Statistic)
    result %>% dplyr::select("Vx.1", "Vx.2", everything())
}

