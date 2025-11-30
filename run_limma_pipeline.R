# =========================================
# Breast Cancer Differential Expression Pipeline
# Using limma + clusterProfiler
# Author: Susovan Borat
# =========================================

# 0. Load Packages -----------------------------------------

library(GEOquery)
library(limma)
library(tidyverse)
library(EnhancedVolcano)
library(pheatmap)
library(clusterProfiler)
library(org.Hs.eg.db)
library(hgu133plus2.db)
library(enrichplot)
library(RColorBrewer)
library(ggplot2)

# Create output folders
dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

palette_pub <- brewer.pal(8, "Set2")

# 1. Download GEO dataset ----------------------------------

gse <- getGEO("GSE15852", GSEMatrix = TRUE)
expr <- exprs(gse[[1]])
pdata <- pData(gse[[1]])

# 2. Group Information --------------------------------------

group_raw <- pdata$`histopathological exam:ch1`
group <- ifelse(group_raw == "normal", "Normal", "Tumor")
group <- factor(group)

# 3. PCA Plot ----------------------------------------------

expr_norm <- normalizeBetweenArrays(expr)
pca <- prcomp(t(expr_norm), scale. = TRUE)

pca_df <- data.frame(PC1 = pca$x[,1],
                     PC2 = pca$x[,2],
                     group = group)

p_pca <- ggplot(pca_df, aes(PC1, PC2, color = group)) +
    geom_point(size = 4) +
    scale_color_manual(values = palette_pub) +
    theme_classic(base_size = 16) +
    labs(title = "PCA of Samples (Normalized Expression)")

ggsave("figures/PCA_plot.png", p_pca, width = 7, height = 6, dpi = 300)

# 4. Design matrix -----------------------------------------

design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

# 5. Contrast + limma ---------------------------------------

contrast.matrix <- makeContrasts(Tumor - Normal, levels = design)

fit <- lmFit(expr, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# 6. DEG extraction -----------------------------------------

deg <- topTable(fit2, adjust = "fdr", number = Inf)
write.csv(deg, "results/DEG_results.csv")

# 7. Volcano Plot -------------------------------------------

p_volcano <- EnhancedVolcano(
    deg,
    lab = rownames(deg),
    x = 'logFC',
    y = 'P.Value',
    pCutoff = 0.05,
    FCcutoff = 1,
    title = "Breast Cancer DEGs",
    subtitle = "Tumor vs Normal"
)

ggsave("figures/volcano.png", p_volcano, width = 8, height = 6, dpi = 300)

# 8. Top 50 Heatmap -----------------------------------------

annotation_col <- data.frame(Group = group)
rownames(annotation_col) <- colnames(expr)

top50 <- deg %>% arrange(P.Value) %>% head(50)

png("figures/heatmap_top50.png", width = 1000, height = 800)
pheatmap(expr[rownames(top50), ],
         scale = "row",
         annotation_col = annotation_col,
         show_rownames = TRUE,
         show_colnames = FALSE,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "complete")
dev.off()

# 9. Probe → Gene Symbol Mapping ----------------------------

probe_ids <- rownames(deg)
symbols <- mapIds(hgu133plus2.db, keys = probe_ids,
                  column = "SYMBOL", keytype = "PROBEID", multiVals = "first")

deg$GeneSymbol <- symbols

# 10. Significant Gene List ---------------------------------

sig_deg <- deg %>%
    filter(adj.P.Val < 0.05, abs(logFC) > 1, !is.na(GeneSymbol))

write.csv(sig_deg, "results/Significant_DEGs_with_symbols.csv")

# Export for STRING
write.table(sig_deg$GeneSymbol,
            "results/STRING_gene_list.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# 11. GO & KEGG Enrichment ----------------------------------

entrez <- bitr(sig_deg$GeneSymbol, fromType="SYMBOL",
               toType="ENTREZID", OrgDb=org.Hs.eg.db)

go <- enrichGO(entrez$ENTREZID, OrgDb=org.Hs.eg.db,
               ont="BP", pAdjustMethod="fdr", readable=TRUE)

kegg <- enrichKEGG(entrez$ENTREZID, organism="hsa")

write.csv(as.data.frame(go), "results/GO_BP_enrichment.csv")
write.csv(as.data.frame(kegg), "results/KEGG_enrichment.csv")

# 12. GO Barplot --------------------------------------------

p_go <- barplot(go, showCategory = 20) +
    theme_classic(base_size = 16) +
    ggtitle("GO Biological Processes (Top 20)")

ggsave("figures/GO_barplot.png", p_go, width = 10, height = 7, dpi = 300)

# 13. KEGG Dotplot ------------------------------------------

p_kegg <- dotplot(kegg, showCategory = 20) +
    theme_classic(base_size = 16) +
    ggtitle("KEGG Pathway Enrichment (Top 20)")

ggsave("figures/KEGG_dotplot.png", p_kegg, width = 10, height = 7, dpi = 300)

# 14. Top 20 Up/Downregulated Heatmaps -----------------------

up20 <- deg %>% arrange(desc(logFC)) %>% head(20)
down20 <- deg %>% arrange(logFC) %>% head(20)

# Upregulated
png("figures/upregulated_top20_heatmap.png", width = 1000, height = 800)
pheatmap(expr[rownames(up20), ],
         scale = "row",
         annotation_col = annotation_col)
dev.off()

# Downregulated
png("figures/downregulated_top20_heatmap.png", width = 1000, height = 800)
pheatmap(expr[rownames(down20), ],
         scale = "row",
         annotation_col = annotation_col)
dev.off()

# 15. Session Info ------------------------------------------

writeLines(capture.output(sessionInfo()),
           "results/session_info.txt")
