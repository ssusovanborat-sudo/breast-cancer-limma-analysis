# packages.R - install required packages for the pipeline
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

cran_pkgs <- c("tidyverse", "pheatmap", "EnhancedVolcano")
bioc_pkgs  <- c("GEOquery", "limma", "clusterProfiler", "org.Hs.eg.db", "hgu133plus2.db")

install_if_missing <- function(pkgs, bioc = FALSE) {
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      if (bioc) BiocManager::install(p, ask = FALSE, update = FALSE)
      else install.packages(p, repos = "https://cloud.r-project.org")
    }
  }
}

install_if_missing(cran_pkgs, bioc = FALSE)
install_if_missing(bioc_pkgs, bioc = TRUE)

message("All packages installed (or already present).")
