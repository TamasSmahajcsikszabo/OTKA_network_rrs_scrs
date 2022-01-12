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
get_bridge_estimate <- function(network, dec = 2, seed=1234) {
    W <- E(network)$strength
    set.seed(seed)
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

    #    }
    if (length(unique(group_estimation$membership)) > 1) {
        bridge_estimate <- networktools::bridge(network, communities = group_estimation$membership)
    } else {
        bridge_estimate <- networktools::bridge(network)
    }

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
    output_tibble
}

network_summary <- function(graph, dec = 2, name='Graph', single_scale=FALSE, seed=1234) {
    summary <- matrix(nrow = vcount(graph))
    summary <- data.frame(summary)
    summary[, 1] <- tibble("Item" = V(graph)$name)
    colnames(summary) <- "Item"
    summary["Subscale"] <- tibble("subscale" = V(graph)$subscale)
    summary["Label"] <- tibble("label" = V(graph)$label)
    summary["Degree"] <- tibble("degree" = degree(graph))
    summary["Strength"] <- tibble("strength" = strength(graph))
    summary["Betweenness"] <- tibble("betweenness" = round(betweenness(graph), dec))
    summary["Closeness"] <- tibble("closeness" = round(closeness(graph), dec))
    influence_df <- networktools::expectedInf(graph, step = c("both"), directed = FALSE)
    summary["1-step Exp. Inf."] <- tibble("step1 Exp. Inf." = (influence_df$step1))
    summary["2-step Exp. Inf."] <- tibble("step2 Exp. Inf." = (influence_df$step2))
    summary <- bind_cols(summary, get_bridge_estimate(graph, seed=seed))
    summary <- summary %>% arrange(Community)
    summary["Title"]  <- c(name, rep("",vcount(graph)-1))
    summary %>% select("Title","Item", "Subscale", "Label", "Community", everything())
}

degree_distribution_summary <- function(graph){
    tibble("Degree" = 0:max(degree(graph), na.rm=TRUE),
           "Degree dist." = degree_distribution(graph, mode='all'))
}

