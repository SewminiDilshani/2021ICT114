# 06_enrichment.R
# Gene Ontology (Biological Process) enrichment of dengue vs control DE genes

library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

res  <- read.csv("results/DE_dengue_vs_control.csv", stringsAsFactors = FALSE)
d    <- readRDS("data/clean_data.rds")
feat <- d$feat

# ---- 1. Add Entrez IDs to the results ----
# If a probe maps to several genes, keep the first one
entrez <- as.character(feat[res$Probe, "ENTREZ_GENE_ID"])
entrez <- sub(" ///.*", "", entrez)
res$Entrez <- entrez

# Remove probes without an Entrez ID
res <- res[!is.na(res$Entrez) & res$Entrez != "", ]

# ---- 2. One probe per gene ----
# res is already sorted by p-value, so keep the first (best) probe of each gene
res <- res[!duplicated(res$Entrez), ]
cat("Genes with an Entrez ID (after removing duplicates):", nrow(res), "\n")

# ---- 3. Gene lists ----
universe  <- res$Entrez                                      # all tested genes
up_genes   <- res$Entrez[res$Significance == "Up in dengue"]
down_genes <- res$Entrez[res$Significance == "Down in dengue"]
cat("Up genes:", length(up_genes), "| Down genes:", length(down_genes), "\n")

# ---- 4. GO enrichment function ----
run_go <- function(genes) {
  ego <- enrichGO(gene          = genes,
                  universe      = universe,
                  OrgDb         = org.Hs.eg.db,
                  keyType       = "ENTREZID",
                  ont           = "BP",          # Biological Process
                  pAdjustMethod = "BH",
                  pvalueCutoff  = 0.05,
                  qvalueCutoff  = 0.05,
                  readable      = TRUE)
  # Merge very similar GO terms
  if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
    ego <- simplify(ego, cutoff = 0.7)
  }
  ego
}

# ---- 5. Run for up and down genes ----
ego_up   <- run_go(up_genes)
ego_down <- run_go(down_genes)

# ---- 6. Save tables, print top terms and draw plots ----
save_and_plot <- function(ego, label, file_tag, fig_name) {
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    cat("No enriched GO terms for", label, "\n")
    return(invisible(NULL))
  }
  tab <- as.data.frame(ego)
  write.csv(tab, paste0("results/GO_", file_tag, ".csv"), row.names = FALSE)
  cat("\nTop GO terms for", label, ":\n")
  print(head(tab[, c("Description", "Count", "p.adjust")], 10))
  
  p <- dotplot(ego, showCategory = 10) +
    ggtitle(paste("GO Biological Process:", label))
  ggsave(paste0("figures/", fig_name), p, width = 8, height = 6, dpi = 300)
}

save_and_plot(ego_up,   "genes up in dengue",   "up",   "fig8_GO_up.png")
save_and_plot(ego_down, "genes down in dengue", "down", "fig9_GO_down.png")