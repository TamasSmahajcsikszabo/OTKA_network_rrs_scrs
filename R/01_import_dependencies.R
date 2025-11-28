#' Attempts to load an R package
#' If not found, tries install the package
#' @param package_name string name of the package
load_library <- function(package_name) {
    tryCatch(
        {
            if (package_name %in% installed.packages()) {
                library(package_name, character.only = TRUE)
            } else {
                install.packages(package_name, dependencies = TRUE, repos = "https://cloud.r-project.org")
            }
        },
        error = function(e) {
            print(e)
        },
        warning = function(w) {
            print(w)
        }
    )
}
options(repos = c(CRAN = "https://cloud.r-project.org"))

load_library("devtools")
load_library("tidyverse")
load_library("networktools")
load_library("qgraph")
load_library("glasso")
load_library("ggraph")
load_library("igraph")
load_library("bootnet")
load_library("statnet")
load_library("modelr")
load_library("scales")
load_library("knitr")
load_library("ggpubr")
load_library("broom")
load_library("ergm")
load_library("fdrtool")
load_library("eigenmodel")
load_library("ROCR")
load_library("bookdown")
load_library("foreign")
load_library("tidyverse")
load_library("haven")
load_library("tibble")
load_library("ergm")
load_library("fdrtool")
load_library("eigenmodel")
load_library("ROCR")
load_library("bookdown")
load_library("extrafont")
load_library("WRS")
load_library("kableExtra")
load_library("ltm")
load_library("pander")

if (!"captioner" %in% installed.packages()) {
    devtools::install_github("adletaw/captioner")
} else {
    library(captioner)
}
if (!"rogme" %in% installed.packages()) {
    devtools::install_github("Grousselet/rogme")
} else {
    library(rogme)
}
if (!"GGMncv" %in% installed.packages()) {
    devtools::install_github("donaldRwilliams/GGMncv")
} else {
    library(GGMncv)
}

if (!"WRS" %in% installed.packages()) {
    # first: install dependent packages
    install.packages(c("MASS", "akima", "robustbase"))

    # second: install suggested packages
    install.packages(c("akima", "cobs", "robust", "mgcv", "scatterplot3d", "quantreg", "rrcov", "lars", "pwr", "trimcluster", "mc2d", "psych", "Rfit", "DepthProc", "class", "fda", "rankFD"))

    # third: install an additional package which provides some C functions
    # install.packages("devtools")
    # NOTE: This seems to be stalled and not functional any more
    # devtools::install_github("mrxiaohe/WRScpp")

    # fourth: install WRS
    devtools::install_github("nicebread/WRS", subdir = "pkg")
}

# if (!"networkAnalysisTools" %in% installed.packages()) {
print("NetworkAnalysisTools not found. Attempting install from GitHub source")
devtools::install_github("TamasSmahajcsikszabo/networkAnalysisTools")
library(networkAnalysisTools)
# } else {
    # library(networkAnalysisTools)
# }
# remotes::install_version("ggplot2", "3.5.2")

print("------- ALL DEPENDENCIES INSTALLED -------")
