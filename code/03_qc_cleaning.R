# 03_qc_cleaning.R
# QC and cleaning of GSE51808 (dengue whole blood expression)

library(GEOquery)
library(Biobase)

dir.create("figures", showWarnings = FALSE)

gse  <- readRDS("data/GSE51808.rds")
expr <- exprs(gse)    # probes x samples
meta <- pData(gse)    # sample metadata
feat <- fData(gse)    # probe annotation

# ---- 1. Clean the group labels ----
meta$group <- factor(meta$`status:ch1`,
                     levels = c("control", "DF", "DHF", "convalescent"),
                     labels = c("Control", "DF", "DHF", "Convalescent"))
print(table(meta$group, useNA = "ifany"))

# ---- 2. Check that sample order matches between matrix and metadata ----
stopifnot(all(colnames(expr) == rownames(meta)))

# ---- 3. Missing values ----
cat("Missing values in expression matrix:", sum(is.na(expr)), "\n")

# ---- 4. Probe annotation columns (to find the gene symbol column) ----
print(colnames(feat))

# ---- 5. QC plots (all 56 samples) ----
cols <- c("grey60", "firebrick", "darkorange", "steelblue")[as.integer(meta$group)]

png("figures/fig1_boxplot_all_samples.png", width = 1400, height = 700)
boxplot(expr, col = cols, las = 2, cex.axis = 0.5, outline = FALSE,
        main = "Expression distribution per sample", ylab = "log2 expression")
legend("topright", legend = levels(meta$group),
       fill = c("grey60", "firebrick", "darkorange", "steelblue"), cex = 0.8)
dev.off()

png("figures/fig2_density_all_samples.png", width = 900, height = 600)
plot(density(expr[, 1]), col = cols[1], main = "Density of expression values",
     xlab = "log2 expression", ylim = c(0, 0.4))
for (i in 2:ncol(expr)) lines(density(expr[, i]), col = cols[i])
legend("topright", legend = levels(meta$group),
       fill = c("grey60", "firebrick", "darkorange", "steelblue"), cex = 0.8)
dev.off()

# ---- 6. Samples for the main analysis: Control + DF + DHF ----
# Convalescent samples come from the same patients as the acute samples,
# so they are excluded from the main comparison.
keep      <- meta$group %in% c("Control", "DF", "DHF")
expr_main <- expr[, keep]
meta_main <- meta[keep, ]
meta_main$group     <- droplevels(meta_main$group)
meta_main$condition <- factor(ifelse(meta_main$group == "Control",
                                     "Control", "Dengue"),
                              levels = c("Control", "Dengue"))
cat("Samples kept:", ncol(expr_main), "\n")
print(table(meta_main$condition))

# ---- 7. Filter low-expression probes ----
# Smallest group has 9 samples, so a probe is kept only if its
# log2 expression is above 3.5 in at least 9 samples.
keep_probes <- rowSums(expr_main > 3.5) >= 9
expr_filt   <- expr_main[keep_probes, ]
feat_filt   <- feat[keep_probes, ]
cat("Probes before filter:", nrow(expr_main),
    "| after filter:", nrow(expr_filt), "\n")

# ---- 8. Save the cleaned data ----
saveRDS(list(expr = expr_filt, meta = meta_main, feat = feat_filt,
             expr_all = expr, meta_all = meta),
        "data/clean_data.rds")