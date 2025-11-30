Methods

Data Retrieval



Dataset GSE15852 (GPL570 platform) was obtained from GEO using the GEOquery package. The dataset includes 43 normal and 43 tumor breast tissue samples.



Preprocessing



Series matrix values were extracted, log-transformed where required, and inspected for quality control.



Group Assignment



Samples were classified based on histopathological exam:ch1.



"normal" → Normal



Others → Tumor



Differential Expression (limma)



A design matrix was generated, followed by contrasts (Tumor − Normal). The limma workflow (lmFit, contrasts.fit, eBayes) was used to compute moderated statistics.



Probe-to-Gene Mapping



Affymetrix probes (GPL570) were mapped to gene symbols using:



library(hgu133plus2.db)

mapIds(hgu133plus2.db, ...)





Probes without valid annotations were removed.



Visualization



Volcano plots generated using EnhancedVolcano



Heatmaps generated using pheatmap with sample annotations



Functional Enrichment



DEGs with |log2FC| > 1 and FDR < 0.05 were analyzed with:



GO Biological Processes



KEGG pathways



via clusterProfiler.



Reproducibility



Session information was exported via sessionInfo().

