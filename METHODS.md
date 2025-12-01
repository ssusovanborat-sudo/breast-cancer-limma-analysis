1\. Data Collection



Publicly available microarray expression data for breast cancer were retrieved from the Gene Expression Omnibus (GEO) under accession GSE15852, generated on the Affymetrix Human Genome U133 Plus 2.0 Array (GPL570) platform. The dataset includes breast tumor tissue samples and matched or unmatched normal breast tissue samples. The raw CEL files and accompanying metadata were downloaded using the GEOquery R package (v2.70.0).



2\. Preprocessing and Quality Control



Raw expression data in CEL format were imported into R (v4.3.x). The following preprocessing steps were performed using the affy and limma packages:



Background correction



Quantile normalization



Log2 transformation



Probe-level summarization using RMA (Robust Multi-array Average)



Quality control included:



Boxplots of raw vs normalized intensities



Density distribution plots



Sample clustering via PCA



Hierarchical clustering dendrograms



Samples showing extreme deviations were inspected but no samples were removed.



3\. Annotation of Probe IDs



Probe identifiers were mapped to official HGNC gene symbols using the hgu133plus2.db annotation package.

Probes mapping to multiple genes or lacking annotation were excluded.

When multiple probes corresponded to the same gene symbol, the probe with the highest mean expression was retained.



4\. Experimental Design and Contrast Matrix



Samples were manually assigned into two groups based on phenotype metadata:



Tumor (n = XX)



Normal (n = XX)



A design matrix was constructed in limma for the comparison:

Tumor vs Normal

The contrast matrix was defined as:

contrast.matrix <- makeContrasts(Tumor - Normal, levels = design)

5\. Differential Expression Analysis



Differential gene expression was computed using the limma pipeline:



Linear modeling via lmFit()



Empirical Bayes moderation via eBayes()



Extraction of significantly dysregulated genes with:

adj.P.Val < 0.05

|log2FC| ≥ 1

The output included log2 fold changes, raw p-values, adjusted p-values (Benjamini–Hochberg), and statistics.



6\. Visualization



Several standard bioinformatics visualizations were produced:



Volcano Plot



Shows upregulated and downregulated genes based on fold change and adjusted p-values.



Heatmap of Top 50 DEGs



Expression of the top 50 DEGs was visualized using hierarchical clustering (Euclidean distance, complete linkage), scaled by row z-scores.



PCA Plot



Used to assess sample stratification between tumor and normal tissue.



7\. Functional Enrichment Analysis



To determine biological significance of DEGs:



GO Biological Process (BP)



KEGG Pathway analysis



were performed using the clusterProfiler package (v4.x). Significant pathways were identified using:

p.adjust < 0.05

Enriched pathways related to cell cycle regulation, proliferation, DNA replication, and breast cancer–related signaling pathways were highlighted.

