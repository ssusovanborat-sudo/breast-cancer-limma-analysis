#!/usr/bin/env Rscript
# ============================================================================
# run_limma_pipeline.R
# Breast cancer differential expression pipeline (limma) - robust & parameterized
# Author: Your Name
# Usage:
#   Rscript scripts/run_limma_pipeline.R [GEO_ID] [OUTDIR] [FIGDIR]
# Example:
#   Rscript scripts/run_limma_pipeline.R GSE15852 results figures
# ============================================================================

# ---------------------
# Config / args
# ---------------------
args <- commandArgs(trailingOnly = TRUE)
GEO_ID <- ifelse(length(args) >= 1, args[1], Sys.getenv("GEO_ID", "GSE15852"))
OUTDIR <- ifelse(length(args) >= 2, args[2], Sys.getenv("OUTDIR", "results"))
FIGDIR <- ifelse(length(args) >= 3, args[3], Sys.getenv("FIGDIR", "figures"))

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIGDIR, showWarnings = FALSE, recursive = TRUE)

logfile <- file.path(OUTDIR, "run_log.txt")
log_con <- file(logfile, open = "wt")
sink(log_con, type = "output")
sink(log_con, type = "message")

cat("Running pipeline:", date(), "\n")
cat("GEO_ID:", GEO_ID, "OUTDIR:", OUTDIR, "FIGDIR:", FIGDIR, "\n\n")

# ---------------------
# Load packages
# ---------------------
required_cran <- c("tidyverse", "pheatmap", "EnhancedVolcano")
required_bioc <- c("GEOquery", "limma", "clusterProfiler", "org.Hs.eg.db")

# load function to require packages
require_install <- function(pkgs, bioc = FALSE) {
  for (p in pkgs) {
    if (!suppressWarnings(requireNamespace(p, quietly = TRUE))) {
      if (bioc) BiocManager::install(p, ask = FALSE, update = FALSE)
      else install.packages(p, repos = "https://cloud.r-project.org")
    }
    library(p, character.only = TRUE)
  }
}

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
require_install(required_cran, bioc = FALSE)
require_install(required_bioc, bioc = TRUE)

# Additional annotation packages loaded later as needed

# ---------------------
# Download GEO dataset
# ---------------------
cat("Downloading GEO:", GEO_ID, "\n")
gse <- tryCatch({
  getGEO(GEO_ID, GSEMatrix = TRUE)
}, error = function(e) {
  stop("Failed to download GEO dataset: ", e$message)
})

if (length(gse) < 1) stop("No series matrix returned for ", GEO_ID)
gse_obj <- gse[[1]]

expr <- tryCatch(exprs(gse_obj), error = function(e) stop("No expression matrix found: ", e$message))
meta <- pData(gse_obj)

cat("Expression dims (genes x samples):", dim(expr), "\n")
cat("Sample metadata cols:", paste(colnames(meta), collapse = ", "), "\n")

# ---------------------
# Robust phenotype detection (Tumor vs Normal)
# ---------------------
detect_group_column <- function(meta_df) {
  # common candidate columns
  candidates <- c("histopathological exam:ch1",
                  "characteristics_ch1.1",
                  "characteristics_ch1",
                  "source_name_ch1",
                  "title",
                  "description")
  found <- intersect(candidates, colnames(meta_df))
  # also include any column containing "histopath" or "tissue" or "diagnosis"
  extra <- grep("histopath|tissue|diagnos|patholog", colnames(meta_df), value = TRUE, ignore.case = TRUE)
  found <- unique(c(found, extra))
  if (length(found) == 0) return(NULL)
  return(found[1])
}

colname_used <- detect_group_column(meta)
if (is.null(colname_used)) stop("Could not find phenotype column automatically. Columns available: ", paste(colnames(meta), collapse = ", "))

cat("Using phenotype column:", colname_used, "\n")
group_raw <- as.character(meta[[colname_used]])

# Normalize common forms to Tumor/Normal
group <- ifelse(grepl("normal", group_raw, ignore.case = TRUE), "Normal",
                ifelse(grepl("control", group_raw, ignore.case = TRUE), "Normal",
                       ifelse(grepl("tumor|carcinoma|cancer|invasive|malignan", group_raw, ignore.case = TRUE), "Tumor", NA))))
# fallback: if all NA, try to parse source_name or title for 'normal'
if (all(is.na(group))) {
  group <- ifelse(grepl("normal", meta$source_name_ch1, ignore.case = TRUE), "Normal", "Tumor")
}

if (any(is.na(group))) {
  warning("Some samples could not be assigned group; assigning them to 'Tumor' by default.")
  group[is.na(group)] <- "Tumor"
}

group <- factor(group)
cat("Group counts:\n"); print(table(group))

stopifnot(length(group) == ncol(expr))
if (any(table(group) < 2)) stop("One of the groups has fewer than 2 samples - cannot run limma reliably.")

# ---------------------
# Design and limma
# ---------------------
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)
cat("Design matrix:\n"); print(head(design))

contrast.matrix <- makeContrasts(Tumor - Normal, levels = design)

fit <- lmFit(expr, design)
fit2 <- tryCatch({
  contrasts.fit(fit, contrast.matrix)
}, error = function(e) stop("Error in contrasts.fit: ", e$message))
fit2 <- eBayes(fit2)

deg <- topTable(fit2, number = Inf, sort.by = "P")
deg$ProbeID <- rownames(deg)

# Save full DEG per-probe table
deg_outfile <- file.path(OUTDIR, paste0(GEO_ID, "_DEG_per_probe.csv"))
write.csv(deg, deg_outfile, row.names = FALSE)
cat("Wrote per-probe DEG table to:", deg_outfile, "\n")

# ---------------------
# Optional: detect platform and map probes -> gene symbols
# ---------------------
platform <- annotation(gse_obj)
cat("Platform annotation:", platform, "\n")

# map probe->symbol if probes look like probe IDs (e.g., contain _at or other patterns)
is_probe_like <- any(grepl("_at$|\\d+_at$|switc|AFFX", rownames(expr), ignore.case = TRUE))
deg$GeneSymbol <- NA

if (is_probe_like) {
  # Attempt to load a matching Bioconductor annotation package for GPL570
  if (grepl("GPL570", platform, ignore.case = TRUE) || grepl("hgu133plus2", platform, ignore.case = TRUE)) {
    annotation_pkg <- "hgu133plus2.db"
  } else {
    # generic fallback - try to guess or ask user to install specific annotation
    annotation_pkg <- NA
  }

  if (!is.na(annotation_pkg) && !suppressWarnings(requireNamespace(annotation_pkg, quietly = TRUE))) {
    cat("Installing annotation package:", annotation_pkg, "\n")
    BiocManager::install(annotation_pkg, ask = FALSE, update = FALSE)
  }

  if (!is.na(annotation_pkg) && requireNamespace(annotation_pkg, quietly = TRUE)) {
    library(annotation_pkg, character.only = TRUE)
    symbols <- mapIds(get(annotation_pkg), keys = deg$ProbeID, column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
    deg$GeneSymbol <- symbols
    cat("Mapped probes to gene symbols. Mapped:", sum(!is.na(symbols)), "probes\n")
  } else {
    cat("No specific annotation package available for platform:", platform, "\nAttempting org.Hs.eg.db-based mapping (may fail).\n")
    # try mapping if rownames already look like SYMBOLS
    possible_symbols <- rownames(expr)
    mapped <- mapIds(org.Hs.eg.db, keys = possible_symbols, column = "SYMBOL", keytype = "SYMBOL", multiVals = "first")
    deg$GeneSymbol <- ifelse(!is.na(mapped), mapped, NA)
  }
} else {
  # Assume rownames are gene symbols already
  deg$GeneSymbol <- rownames(expr)
}

# Collapse probes mapping to same gene: choose probe with lowest adj.P.Val
deg_collapsed <- deg %>%
  filter(!is.na(GeneSymbol)) %>%
  as_tibble() %>%
  group_by(GeneSymbol) %>%
  slice_min(order_by = adj.P.Val, n = 1) %>%
  ungroup()

# If no mapped genes found, fall back to deg (probes)
if (nrow(deg_collapsed) == 0) {
  warning("No gene symbols mapped. Continuing with probe-level DEGs.")
  deg_for_enrich <- deg
  rownames(deg_for_enrich) <- deg_for_enrich$ProbeID
} else {
  deg_for_enrich <- as.data.frame(deg_collapsed)
  rownames(deg_for_enrich) <- deg_for_enrich$GeneSymbol
}

# Save collapsed DEG
deg_collapsed_file <- file.path(OUTDIR, paste0(GEO_ID, "_DEG_collapsed.csv"))
write.csv(deg_for_enrich, deg_collapsed_file, row.names = TRUE)
cat("Wrote collapsed DEG table to:", deg_collapsed_file, "\n")

# ---------------------
# Volcano plot (use FDR if available)
# ---------------------
volcano_file <- file.path(FIGDIR, paste0(GEO_ID, "_volcano.png"))
png(volcano_file, width = 1400, height = 1200)
y_col <- ifelse("adj.P.Val" %in% colnames(deg_for_enrich), "adj.P.Val",
                ifelse("P.Value" %in% colnames(deg_for_enrich), "P.Value", colnames(deg_for_enrich)[2]))
EnhancedVolcano(deg_for_enrich,
                lab = rownames(deg_for_enrich),
                x = 'logFC',
                y = y_col,
                pCutoff = 0.05,
                FCcutoff = 1,
                title = paste(GEO_ID, "Tumor vs Normal"))
dev.off()
cat("Wrote volcano to:", volcano_file, "\n")

# ---------------------
# Heatmap (top N genes by adj.P.Val or P.Value)
# ---------------------
topN <- 50
if ("adj.P.Val" %in% colnames(deg_for_enrich)) {
  top_genes <- rownames(head(deg_for_enrich[order(deg_for_enrich$adj.P.Val), ], topN))
} else {
  top_genes <- rownames(head(deg_for_enrich[order(deg_for_enrich$P.Value), ], topN))
}

# create annotation_col with sample names as rownames
annotation_col <- data.frame(Group = group)
rownames(annotation_col) <- colnames(expr)

heat_file <- file.path(FIGDIR, paste0(GEO_ID, "_heatmap_top", topN, ".png"))
png(heat_file, width = 1400, height = 1200)
pheatmap(expr[top_genes, , drop = FALSE],
         scale = "row",
         annotation_col = annotation_col,
         show_rownames = TRUE,
         show_colnames = FALSE,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "complete")
dev.off()
cat("Wrote heatmap to:", heat_file, "\n")

# ---------------------
# Functional enrichment (clusterProfiler)
# ---------------------
# Prepare gene list for enrichment: ENTREZ IDs
sig_genes_symbols <- rownames(deg_for_enrich)[(deg_for_enrich$adj.P.Val < 0.05 & abs(deg_for_enrich$logFC) > 1) |
                                               (deg_for_enrich$P.Value < 0.01 & abs(deg_for_enrich$logFC) > 1)]
sig_genes_symbols <- unique(na.omit(sig_genes_symbols))
cat("Number of significant symbols for enrichment:", length(sig_genes_symbols), "\n")

if (length(sig_genes_symbols) > 0) {
  entrez_df <- bitr(sig_genes_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  if (nrow(entrez_df) > 0) {
    go <- enrichGO(entrez_df$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "fdr", readable = TRUE)
    kegg <- enrichKEGG(entrez_df$ENTREZID, organism = "hsa")
    write.csv(as.data.frame(go), file.path(OUTDIR, paste0(GEO_ID, "_GO_BP_enrichment.csv")), row.names = FALSE)
    write.csv(as.data.frame(kegg), file.path(OUTDIR, paste0(GEO_ID, "_KEGG_enrichment.csv")), row.names = FALSE)
    cat("Wrote enrichment results to:", OUTDIR, "\n")
  } else {
    warning("No ENTREZ IDs mapped from significant symbols; skipping enrichment.")
  }
} else {
  cat("No significant genes passed thresholds - skipping enrichment.\n")
}

# ---------------------
# Session info & close
# ---------------------
writeLines(capture.output(sessionInfo()), file.path(OUTDIR, paste0(GEO_ID, "_session_info.txt")))
cat("Pipeline finished at:", date(), "\n")
sink(type = "message"); sink(type = "output")
close(log_con)
