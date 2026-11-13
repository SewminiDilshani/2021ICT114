# 05_limma.R
# Differential expression analysis with limma

library(limma)
library(ggplot2)

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

d    <- readRDS("data/clean_data.rds")
expr <- d$expr      # filtered probes x 37 samples
meta <- d$meta
feat <- d$feat      # probe annotation

# Helper: convert probe IDs to gene symbols
# (if a probe maps to several genes, keep the first one)
gene_of <- function(ids) {
  g <- as.character(feat[ids, "Gene Symbol"])
  g <- sub(" ///.*", "", g)
  ifelse(is.na(g) | g == "", ids, g)
}

# ================= MAIN: Dengue vs Control =================
# Control is the baseline, so positive logFC = higher in dengue
design1 <- model.matrix(~ condition, data = meta)
fit1    <- lmFit(expr, design1)
fit1    <- eBayes(fit1)

res1 <- topTable(fit1, coef = "conditionDengue", number = Inf,
                 adjust.method = "BH", sort.by = "P")
res1$Probe <- rownames(res1)
res1$Gene  <- gene_of(res1$Probe)

# Significant = adjusted p < 0.05 and |log2 fold change| >= 1
res1$Significance <- "Not significant"
res1$Significance[res1$adj.P.Val < 0.05 & res1$logFC >=  1] <- "Up in dengue"
res1$Significance[res1$adj.P.Val < 0.05 & res1$logFC <= -1] <- "Down in dengue"
print(table(res1$Significance))

write.csv(res1, "results/DE_dengue_vs_control.csv", row.names = FALSE)

# ---- Volcano plot ----
top_lab <- head(res1[res1$Significance != "Not significant", ], 15)

p1 <- ggplot(res1, aes(logFC, -log10(P.Value), colour = Significance)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_colour_manual(values = c("Up in dengue" = "firebrick",
                                 "Down in dengue" = "steelblue",
                                 "Not significant" = "grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_text(data = top_lab, aes(label = Gene), size = 3,
            check_overlap = TRUE, vjust = -0.5, colour = "black") +
  labs(title = "Dengue vs Control (whole blood)",
       x = "log2 fold change", y = "-log10 p-value") +
  theme_minimal()
ggsave("figures/fig7_volcano_dengue_vs_control.png", p1,
       width = 7, height = 5.5, dpi = 300)

# ---- Top 10 genes ----
print(head(res1[, c("Gene", "logFC", "AveExpr", "adj.P.Val")], 10))

# ================= EXTRA: DHF vs DF =================
sub         <- meta$group %in% c("DF", "DHF")
meta2       <- meta[sub, ]
meta2$group <- droplevels(meta2$group)
design2     <- model.matrix(~ group, data = meta2)   # DF is the baseline
fit2        <- eBayes(lmFit(expr[, sub], design2))

res2 <- topTable(fit2, coef = "groupDHF", number = Inf,
                 adjust.method = "BH", sort.by = "P")
res2$Probe <- rownames(res2)
res2$Gene  <- gene_of(res2$Probe)

cat("DHF vs DF: probes with adj.P < 0.05 =",
    sum(res2$adj.P.Val < 0.05), "\n")
cat("DHF vs DF: probes with raw P < 0.01 =",
    sum(res2$P.Value < 0.01), "\n")
write.csv(res2, "results/DE_DHF_vs_DF.csv", row.names = FALSE)