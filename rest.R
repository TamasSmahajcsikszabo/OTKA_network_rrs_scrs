
### 3.1 Network plots
When estimating networks, we modeled RRS, SCRS and their combined network as Gaussian graph-networks which are networks of partial correlation coefficients [@Epskamp2018]. In many parts of this analysis we also utilized functions from the *qgaph* R package [@Epskamp2012]. Please refer Fig. 1-3.



```{r echo=FALSE, message=FALSE, warning=FALSE, fig.cap = "RRS network with gLASSO"}

RRS_labels <- data_labels[grep("RRS", data_labels$var),]$labels
RRS_labels  <- RRS_labels[1:(length(RRS_labels)-3)]
RSScov2 <- qgraph::cor_auto(RRS_subset)
rss_network <- qgraph(RSScov2, graph = "glasso", tuning = 0.5, layout = "spring", sampleSize = 888, theme = "TeamFortress", details = TRUE, threshold = FALSE, labels = RRS_labels)

#rss_network  <- estimateNetwork(RRS_subset, default="EBICglasso")
```



```{r echo=FALSE, message=FALSE, warning=FALSE, fig.cap = "SCRS network with gLASSO"}

SCRS_labels <- data_labels[grep("SCRS", data_labels$var),]$labels
SCRScov <- qgraph::cor_auto(SCRS_subset_reduced)
scrs_network <- qgraph(SCRScov, graph = "glasso", tuning = 0.5, layout = "spring", sampleSize = 888, theme = "TeamFortress", details = TRUE, threshold = FALSE, labels = SCRS_labels)

#scrs_network  <- estimateNetwork(SCRS_subset, default="EBICglasso")
```




```{r echo=FALSE, message=FALSE, warning=FALSE, fig.cap = "Combined RRS and SCRS network"}

# full data scenario
# MHC not included!
full_data <- dataset[, c(RRS_items, SCRS_items)] 
full_goldbricker <- goldbricker(full_data, p = 0.05, method = "hittner2003", threshold = 0.25, corMin = 0.5, progressbar = TRUE)
saveRDS(full_goldbricker, "output/full_goldbricker.RDS")
full_goldbricker <- readRDS("output/full_goldbricker.RDS")


#full_goldbricker$suggested_reductions
full_reduced <- net_reduce(full_data, full_goldbricker, method = c("PCA", "best_goldbricker"))
full_cov <- qgraph::cor_auto(full_reduced)
full_labels  <- c(RRS_labels, SCRS_labels)
full_network <- qgraph(full_cov, graph = "glasso", tuning = 0.5, layout = "spring", sampleSize = 888, theme = "TeamFortress", details = TRUE, threshold = FALSE, labels = full_labels)

```

### 3.2 Community detection 
Community mapping of RSS and SCRS with the *spin glass algorithm* [@Reichardt2006; @Bernstein2019].


```{r echo=FALSE, message=FALSE, warning=FALSE, fig.cap = "RSS communities identified with parameters *gamma = 0.5, spins = 25, starting temperature = 1, ending temperature = 0.01, cooling factor = 0.99*"}

make_community_graph(rss_network, RSScov2, RRS_labels)
```

RRS items (Fig.4.) are grouped into two communities. In one we find the strongly connected items *alone1* and *alone2* with items *write*, *person* and *analyze*. The rest of the items constitute the other main community.




```{r echo=FALSE, message=FALSE, warning=FALSE, fig.cap = "SCRS communities identified with parameters *gamma = 0.5, spins = 25, starting temperature = 1, ending temperature = 0.01, cooling factor = 0.99*."}

make_community_graph(scrs_network, SCRScov, SCRS_labels)
```

The self-critical rumination items exhibit one large community as a whole (Fig.5.).



```{r echo=FALSE, message=FALSE, warning=FALSE, fig.cap ="RRS and SCRS combined communities identified with parameters *gamma = -1.5, spins = 25, starting temperature = 1, ending temperature = 0.01, cooling factor = 0.99*."}

make_community_graph(full_network, full_cov, full_labels)

```

When all RRS and SCRS items are handled as one covariance and one Gaussian graph-network, 4 items of RRS (*handle*, *problem*, *deserve* and *wish*) are grouped into the community of the SCRS items (Fig.6.).



### 3.3 Expected Influence
In order to estimate the Expected Influence metrics, we used the *expectedInf* function of the *networktools* R package. The influence score are standardized.




```{r echo=FALSE, fig.height=5, fig.width=7, message=FALSE, warning=FALSE, fig.cap= "RRS Expected Influence Measures"}
expected_inf_plot(rss_network, "RRS")
```

Regarding RRS influence estimates (Fig.7.): Items such as *alone1*, *alone2* and *person* are high in terms of expected influence, meaning that these items play more central role in the RRS network, by changing these nodes, the rest of the network changes more likely.




```{r echo=FALSE, fig.height=5, fig.width=7, message=FALSE, warning=FALSE, fig.cap = "SCRS Expected Influence Measures"}
expected_inf_plot(scrs_network, "SCRS")
```

SCRS-items (Fig.8.) such as *habits* and *thoughts* are more influencing nodes than *productive*.


```{r echo=FALSE, fig.height=5, fig.width=7, message=FALSE, warning=FALSE, fig.cap = "Combined RRS and SCRS Expected Influence Measures"}
expected_inf_plot(full_network, "RRS & SCRS")
```

In the combined network scenario (Fig.9.), the RRS item *alone1*, the SCRS items *habits* and *thoughts* maintain their influencing position among the items. Due to the combination of tools, the SCRS item *different* has increased its influence, and the items *mistakes and aspects* have increased their impact as a PCA component. This also means that in the combined network they proved to be too redundant and the *goldbricker* algorithm used their PCA-component instead.


### 3.4 Bridge nodes

Using the *networktools* R library bridge expected influence measure were estimated for the three tools we used in this study. We also added aesthetics to reflect which items belonged to the same community according to our previous clustering steps (Fig. 10-12.).




```{r echo=FALSE, fig.height=4, fig.width=11, message=FALSE, warning=FALSE, fig.cap = "RRS bridge estimates"}
# get the estimate of bridge influence scores
make_bridge_plot(rss_network, "RRS")
```



```{r echo=FALSE, fig.height=4, fig.width=11, message=FALSE, warning=FALSE, fig.cap = "SCRS bridge estimates"}
make_bridge_plot(scrs_network, "SCRS")
```


```{r echo=FALSE, fig.height=4, fig.width=11, message=FALSE, warning=FALSE, fig.cap = "Combined RRS and SCRS bridge estimates"}
# make_bridge_plot(full_network, "RSS & SCRS")
```




### 3.5 Network accuracy and stability
We used the *bootnet* R package to have an estimate of edge-weight accuracy with *non-parametric bootstrap* (2500 times sampling) and also produce measures of centrality stability with *case dropping bootstrap* (to simulate cases where only the subset of the existing data would be available).

Network accuracy plots show pairs of nodes (pairs of items) by ascending order of the edge-weights. 95% confidence intervals and the vertical line of 0 edge-weight are shown separately.



```{r echo=FALSE, fig.height=8, fig.width=10, message=FALSE, warning=FALSE, fig.cap = "Network accuracy estimate for RRS"}
# generate the accuracy estimate and save it
#RRS_network <- estimateNetwork(RRS_subset, default="EBICglasso")
#accuracy_rrs  <- bootnet(RRS_network, nBoots = 2500, cores = 6)
#saveRDS(accuracy_rrs, "../output/accuracy_rrs.RDS")

accuracy_rrs  <- readRDS("output/accuracy_rrs.RDS")
network_accuracy_plot(accuracy_rrs, "RRS")
```
From the pair of items *react* and *alone1* confidence intervals of edge weights start to incorporate the value of 0 and thus we cannot be sure if their weight persist given different samples and varying sample characteristics.




```{r echo=FALSE, fig.height=8, fig.width=10, message=FALSE, warning=FALSE, fig.cap = "Network accuracy estimate for SCRS"}
# generate the accuracy estimate and save it
#SCRS_network <- estimateNetwork(SCRS_subset, default="EBICglasso")
#accuracy_scrs  <- bootnet(SCRS_network, nBoots = 2500, cores = 6)
#saveRDS(accuracy_scrs, "../output/accuracy_scrs.RDS")

accuracy_scrs  <- readRDS("output/accuracy_scrs.RDS")
network_accuracy_plot(accuracy_scrs, "SCRS")
```
Similar to RRS, SCRS item pairs are ordered by their edge-weight estimates. Below the pair of items of *time* and *mistake*, accuracy might not be non-zero under different sampling conditions.
The combined RRS and SCRS network estimates are added to the *Appendix* section of this paper.


Network stability estimates for our three networks were carried out by using the *bootnet* function from the R package bearing the same name. We performed a *case dropping bootstrap* to simulate scenarios where fewer and fewer cases of the original sample would have been available and quantified the stability of the centrality metrics with the *Correlation Stability* (CS) coefficient. When estimating bootstrap, we set the number of bootstraps to 2500.



```{r echo=FALSE, fig.height=8, fig.width=12, message=FALSE, warning=FALSE, fig.cap = "Stability estimates of the RRS network centrality indices"}

#RRS_network <- estimateNetwork(RRS_subset, default="EBICglasso")
#RRS_stability <- bootnet(RRS_network, nBoots = 2500, type="case", nCore = 6, statistics = c("strength", "betweenness", "closeness", "edge", "inStrength", "outStrength"), computeCentrality = TRUE)
#saveRDS(RRS_stability, "../output/RRS_stability.RDS")

#SCRS_network <- estimateNetwork(SCRS_subset, default="EBICglasso")
#SCRS_stability <- bootnet(SCRS_network, nBoots = 2500, type="case", nCore = 6, statistics = c("strength", "betweenness", "closeness", "edge", "inStrength", "outStrength"), computeCentrality = TRUE)
#saveRDS(SCRS_stability, "../output/SCRS_stability.RDS")

#full_network <- estimateNetwork(full_data, default="EBICglasso")
#full_stability <- bootnet(full_network, nBoots = 2500, type="case", nCore = 6, statistics = c("strength", "betweenness", "closeness", "edge", "inStrength", "outStrength"), computeCentrality = TRUE)
#saveRDS(full_stability, "../output/full_stability.RDS")




plot_RRS_stab <- readRDS("output/plot_RRS_stability.RDS")
network_stability_plot(plot_RRS_stab)
```
Edge weights and node strength estimates exhibit strong stability, whereas closeness estimates are more sensitive to sample size drops. Betweenness scores seem to be insensitive to sample sizes and are consistently weak in terms of their low correlations with the original estimates for RRS.





```{r echo=FALSE, fig.height=8, fig.width=12, message=FALSE, warning=FALSE, fig.cap="Stability estimates of the SCRS network centrality indices"}
plot_SCRS_stab <- readRDS("output/plot_SCRS_stability.RDS")
network_stability_plot(plot_SCRS_stab)
```
Edge weights, node strength and closeness estimates exhibit strong stability,while betweenness scores seem to be insensitive to sample sizes and are consistently weak in terms of their low correlations with the original estimates for SCRS, similar to RRS.

```{r echo=FALSE, fig.height=8, fig.width=12, message=FALSE, warning=FALSE}
#stability
RRS_stability <- readRDS("output/RRS_stability.RDS")
SCRS_stability <- readRDS("output/SCRS_stability.RDS")
full_stability <- readRDS("output/full_stability.RDS")



#CS estimates

#RRScorStab <-corStability(RRS_stability)
#saveRDS(RRScorStab, "../output/RRScorStab.RDS")
#SCRScorStab <- corStability(SCRS_stability)
#saveRDS(SCRScorStab, "../output/SCRScorStab.RDS")
#fullcorStab <- corStability(full_stability)
#saveRDS(fullcorStab, "../output/fullcorStab.RDS")

RRScorStab <- readRDS("output/RRScorStab.RDS")
SCRScorStab <- readRDS("output/SCRScorStab.RDS")
fullcorStab <- readRDS("output/fullcorStab.RDS")

CS_estimates <- tibble(RRS = RRScorStab,
       SCRS = SCRScorStab,
       combined = fullcorStab)

rownames(CS_estimates) <- names(RRScorStab)

CS_estimates <- CS_estimates %>% 
  mutate(measure = rownames(CS_estimates)) %>% 
  select(Metric = measure, RRS, SCRS, combined)
knitr::kable(CS_estimates, 
             caption = "Correlation Stability coefficient estimates for RRS, SCRS and their combined network", digits = 4)
```

As shown at [@Epskamp2017], the CS-coefficent is ideally between 0.25 and 0.5, the latter is the cat-off point for the estimate. The results show that RRS and SCRS closeness metrics are stable in the sense that around 28% of the case can be dropped to still maintain 0.7 or higher correlation with their original estimates. The closeness measure (`r round(CS_estimates[2, 4][[1]],4)`) for their combined network is already higher than the 0.5 cut-off. The same applies to all edge and strength estimates, with no regard of the network of RRS, SCRS or their combined variant. Betweenness measures are too low to be considered secure to interpret (all below 0.25) meaning that even slight case drops result in losing their stability.


## 4. Discussion
As [@Jones2018] pointed out, when network modeling is done by using force-directed networks[Fruchterman hivatkozas!!!], node placement is not directly interpretable, as these graph are rendered to be aesthetic. 


```{r echo=FALSE, fig.height=20, fig.width=14, message=FALSE, warning=FALSE, fig.cap="Network accuracy estimate for combined RRS and SCRS"}
# generate the accuracy estimate and save it
#combined_network <- estimateNetwork(full_data, default="EBICglasso")
#accuracy_combined  <- bootnet(combined_network, nBoots = 2500, cores = 6)
#saveRDS(accuracy_combined, "../output/accuracy_combined.RDS")

accuracy_combined  <- readRDS("output/accuracy_combined.RDS")
# network_accuracy_plot(accuracy_combined, "combined RRS and SCRS", save_data=TRUE)
```

This combined accuracy estimates for RRS and SCRS include pairs of items from both questionnaires. The last pair of items where the confidence interval does not include 0 edge-weight is the pair of *handle* from RRS and *mistakes* from SCRS. The majority of the more accurate pair estimates are homogeneous in the sense that they are not RRS-SCRS mixed pairs (they are both either from RRS or SCRS). *Handle* belonged to the community of the SCRS nodes in our community analysis.

