Breast Cancer Differential Expression Analysis (R + limma)

This repository contains a complete, reproducible pipeline for analyzing breast cancer gene expression data using R, limma, and clusterProfiler.
The workflow includes:

GEO download
Preprocessing & normalization
Differential gene expression (limma)
Volcano plot & heatmap
Probe-to-gene symbol mapping (GPL570)
GO / KEGG enrichment
Export of results

Repository Structure
breast-cancer-limma-analysis/
├── scripts/                # R analysis scripts
├── results/                # DEG tables, enrichment results
├── figures/                # Volcano, heatmap
├── docs/                   # Workflow diagrams
└── README.md               # Project documentation

How to Run
Install packages:
install.packages(c("tidyverse", "pheatmap", "EnhancedVolcano", "GEOquery"))
BiocManager::install(c("limma", "org.Hs.eg.db", "hgu133plus2.db", "clusterProfiler"))


Run pipeline:
source("scripts/run_limma_pipeline.R")

Outputs
DEG_results.csv
GO_BP_enrichment.csv
KEGG_enrichment.csv
volcano.png
heatmap_top50.png
session_info.txt

License

This project is licensed under the MIT License.
