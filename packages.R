# Package installation script for:
# Breast Cancer LIMMA Pipeline
# Author: Susovan Borat

packages <- c(
  "GEOquery","limma","annotate","hgu133plus2.db","clusterProfiler",
  "org.Hs.eg.db","pheatmap","ggplot2","STRINGdb","factoextra"
)

install.packages(setdiff(packages, installed.packages()[,"Package"]), 
                 dependencies=TRUE)
