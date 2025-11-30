# =========================================
# Breast Cancer Differential Expression Pipeline
# Using limma
# Author: Susovan Borat
# GitHub: https://github.com/ssusovanborat-sudo/breast-cancer-limma-analysis
# =========================================

# 0. Load Packages
library(GEOquery)
library(limma)
library(tidyverse)
library(EnhancedVolcano)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)

# For Affymetrix GPL570 mapping
library(hgu133plus2.db)

# 1. Download GEO dataset (GSE15852 example)
gse <- getGEO("GSE15852", GSEMatrix = TRUE)
expr <- exprs(gse[[1]])
pdata <- pData(gse[[1]])

# 2. Extract group information
group_raw <- pdata$`histopathological exam:ch1`
group <- ifelse(group_raw == "normal", "Normal", "Tumor")
group <- factor(group)
table(group)

# 3. Design matrix
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# 4. Create contrast
contrast.matrix <- makeContrasts(Tumor - Normal, levels = design)

# 5. limma analysis
fit <- lmFit(expr, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# 6. Extract DEGs
deg <- topTable(fit2, adjust="fdr", number=Inf)
write.csv(deg, "results/DEG_results.csv", row.names=TRUE)

# 7. Volcano plot
EnhancedVolcano(deg,
    lab = rownames(deg),
    x = 'logFC',
    y = 'P.Value',
    pCutoff = 0.05,
    FCcutoff = 1,
    title = "Breast Cancer DEGs",
    subtitle = "Tumor vs Normal"
)
ggsave("figures/volcano.png")

# 8. Heatmap of top 50 genes
top50 <- deg %>% arrange(P.Value) %>% head(50)

# Fix for annotation: rownames must match sample names
annotation_col <- data.frame(Group = group)
rownames(annotation_col) <- colnames(expr)

p <- pheatmap(expr[rownames(top50), ], 
         scale = "row",
         annotation_col = annotation_col,
         show_rownames = TRUE,
         show_colnames = FALSE,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "complete")

# Save heatmap as PNG
png("figures/heatmap_top50.png", width=1000, height=800)
pheatmap(expr[rownames(top50), ], 
         scale = "row",
         annotation_col = annotation_col,
         show_rownames = TRUE,
         show_colnames = FALSE,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "complete")
dev.off()

# 9. Map probes to gene symbols (GPL570)
probe_ids <- rownames(deg)
symbols <- mapIds(hgu133plus2.db, keys=probe_ids,
                  column="SYMBOL", keytype="PROBEID", multiVals="first")

deg$GeneSymbol <- symbols

# Select significant genes with symbols
sig_genes <- deg$GeneSymbol[deg$adj.P.Val < 0.05 & abs(deg$logFC) > 1]
sig_genes <- sig_genes[!is.na(sig_genes)]  # remove unmapped probes

# 10. Functional enrichment (GO + KEGG)
entrez <- bitr(sig_genes, fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Hs.eg.db)

# GO enrichment
go <- enrichGO(entrez$ENTREZID, OrgDb=org.Hs.eg.db,
               ont="BP", pAdjustMethod="fdr", readable=TRUE)
write.csv(as.data.frame(go), "results/GO_BP_enrichment.csv")

# KEGG enrichment
kegg <- enrichKEGG(entrez$ENTREZID, organism="hsa")
write.csv(as.data.frame(kegg), "results/KEGG_enrichment.csv")

# 11. Save R session info for reproducibility
writeLines(capture.output(sessionInfo()), "results/session_info.txt")

