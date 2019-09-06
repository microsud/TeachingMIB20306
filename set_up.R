

setup_proteome_shiny <- function(){

.packages = c("flexdashboard",
              "RColorBrewer", 
              "DT", 
              "data.table",
              "reshape", 
              "reshape2", 
              "magrittr", 
              "markdown",
              "ggpubr", 
              "shiny",
              "pheatmap",
              "dplyr", 
              "shinythemes", 
              "devtools", 
              "tidyr",
              "png",
              "ggplot2", 
              "knitr", 
              "plyr", 
              "dplyr",
              "formatR",
              "plotly",
              "stringi",
              "filesstrings",
              "shinycssloaders")

.bioc_packages <- c("DEP",
                    "pathview", 
                    "SummarizedExperiment", 
                    "clusterProfiler",
                    "Biobase",
                    "EnrichmentBrowser",
                    "enrichplot",
                    "EnhancedVolcano",
                    "gage")

# Install CRAN packages (if not already installed)
.inst <- .packages %in% installed.packages()
if(length(.packages[!.inst]) > 0) install.packages(.packages[!.inst])

.inst <- .bioc_packages %in% installed.packages()
if(any(!.inst)) {
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(.bioc_packages[!.inst])
}

message("If there was no error then you are ready to do proteome data analysis")
}


setup_proteome_shiny()