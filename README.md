Overview



This repository provides a fully reproducible pipeline for analyzing breast cancer gene expression data using R and limma.

It is designed for microarray data (Affymetrix GPL570) but can be adapted for other platforms.



The workflow includes:



Data download from GEO/TCGA



Preprocessing \& normalization



Group extraction (Tumor vs Normal)



Limma differential expression analysis



Volcano plot and heatmap visualization



Probe-to-gene mapping for functional enrichment



GO Biological Process (BP) and KEGG pathway analysis



Export of results and reproducibility information



Dataset



GEO Accession: GSE15852



Samples: 43 Normal, 43 Tumor



Platform: GPL570 (Affymetrix Human Genome U133 Plus 2.0 Array)



Note: If using a different dataset, adjust probe-to-gene mapping accordingly.



Repository Structure

breast-cancer-limma-analysis/

│

├── data/       # Raw GEO series matrix or CEL files

├── scripts/    # R scripts for analysis

│   └── run\_limma\_pipeline.R

├── results/    # DEG tables, GO/KEGG enrichment, session info

├── figures/    # Volcano plots, heatmaps

├── docs/       # Flowchart and supplementary materials

└── README.md   # This file



Installation



Install R (>=4.2)



Install required packages:



install.packages(c("tidyverse", "pheatmap", "EnhancedVolcano", "clusterProfiler", "GEOquery"))

BiocManager::install(c("limma", "org.Hs.eg.db", "hgu133plus2.db"))



How to Run the Pipeline

\# Load the pipeline script

source("scripts/run\_limma\_pipeline.R")





The script performs:



Data download and extraction



Differential expression analysis with limma



Volcano and heatmap visualization



Probe-to-gene mapping (GPL570)



GO/KEGG enrichment



Export of all results to results/ and figures/



Outputs



results/DEG\_results.csv — Full DEG table with gene symbols



results/GO\_BP\_enrichment.csv — GO Biological Process enrichment



results/KEGG\_enrichment.csv — KEGG pathway enrichment



results/session\_info.txt — R session info for reproducibility



figures/volcano.png — Volcano plot of DEGs



figures/heatmap\_top50.png — Heatmap of top 50 DEGs



Customization



Replace GSE15852 with your dataset accession



Adjust group extraction if phenotype columns differ



If using a different microarray platform, replace the annotation package (hgu133plus2.db) with the correct one



License



MIT License



Contact



s.susovanborat@gmail.com



GitHub: https://github.com/ssusovanborat-sudo/breast-cancer-limma-analysis



