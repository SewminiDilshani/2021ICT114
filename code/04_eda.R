# 04_eda.R
# Exploratory data analysis: PCA and heatmaps

library(ggplot2)
library(pheatmap)

dir.create("figures", showWarnings = FALSE)

d    <- readRDS("data/clean_data.rds")
expr <- d$expr      # filtered probes x 37 samples
meta <- d$meta
feat <- d$feat

group_cols <- c(Control = "grey50", DF = "darkorange",
                DHF = "firebrick", Convalescent = "steelblue")

# ---- 1. PCA (Control + DF + DHF) ----
# Use the 2000 most variable probes
top_var <- names(sort(apply(expr, 1, var), decreasing = TRUE))[1:2000]
pca     <- prcomp(t(expr[top_var, ]), scale. = TRUE)
var_exp <- round(100 * summary(pca)$importance[2, 1:2], 1)
cat("PC1 variance:", var_exp[1], "% | PC2 variance:", var_exp[2], "%\n")

pca_df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                     group = meta$group)

p1 <- ggplot(pca_df, aes(PC1, PC2, colour = group)) +
  geom_point(size = 3) +
  scale_colour_manual(values = group_cols) +
  labs(title = "PCA of whole blood expression (top 2000 variable probes)",
       x = paste0("PC1 (", var_exp[1], "%)"),
       y = paste0("PC2 (", var_exp[2], "%)")) +
  theme_minimal()
ggsave("figures/fig3_pca_main.png", p1, width = 7, height = 5, dpi = 300)

# ---- 2. PCA including convalescent samples (same probes) ----
expr_all <- d$expr_all[rownames(expr), ]
meta_all <- d$meta_all
pca2     <- prcomp(t(expr_all[top_var, ]), scale. = TRUE)
var2     <- round(100 * summary(pca2)$importance[2, 1:2], 1)

pca2_df <- data.frame(PC1 = pca2$x[, 1], PC2 = pca2$x[, 2],
                      group = meta_all$group)

p2 <- ggplot(pca2_df, aes(PC1, PC2, colour = group)) +
  geom_point(size = 3) +
  scale_colour_manual(values = group_cols) +
  labs(title = "PCA including convalescent samples",
       x = paste0("PC1 (", var2[1], "%)"),
       y = paste0("PC2 (", var2[2], "%)")) +
  theme_minimal()
ggsave("figures/fig4_pca_all.png", p2, width = 7, height = 5, dpi = 300)

# ---- 3. Sample-to-sample correlation heatmap ----
annot <- data.frame(Group = meta$group, row.names = rownames(meta))
annot_cols <- list(Group = group_cols[c("Control", "DF", "DHF")])

cor_mat <- cor(expr)
png("figures/fig5_sample_correlation.png", width = 1000, height = 900, res = 130)
pheatmap(cor_mat, annotation_col = annot, annotation_row = annot,
         annotation_colors = annot_cols,
         show_rownames = FALSE, show_colnames = FALSE,
         main = "Sample-to-sample correlation")
dev.off()

# ---- 4. Heatmap of the top 50 most variable probes ----
top50  <- top_var[1:50]
labels <- feat[top50, "Gene Symbol"]
labels <- ifelse(is.na(labels) | labels == "", top50, labels)
labels <- make.unique(as.character(labels))

png("figures/fig6_heatmap_top50.png", width = 1100, height = 1100, res = 130)
pheatmap(expr[top50, ], scale = "row",
         labels_row = labels,
         annotation_col = annot, annotation_colors = annot_cols,
         show_colnames = FALSE, fontsize_row = 7,
         main = "Top 50 most variable probes")
dev.off()

cat("EDA figures saved in figures/ folder\n")