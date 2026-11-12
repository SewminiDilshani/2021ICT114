# 01_download_data.R
# Download dengue patient blood expression dataset GSE51808 from GEO

library(GEOquery)
set.seed(123)

gse <- getGEO("GSE51808", GSEMatrix = TRUE)[[1]]

expr <- exprs(gse)   # probes x samples
meta <- pData(gse)   # sample information

cat("Probes:", nrow(expr), "\n")
cat("Samples:", ncol(expr), "\n")
cat("Platform:", annotation(gse), "\n")
print(colnames(meta))
print(head(meta[, 1:5]))

saveRDS(gse, "data/GSE51808.rds")