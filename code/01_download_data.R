# 01_download_data.R
# Download the dengue patient blood expression dataset GSE51808 from GEO

library(GEOquery)
set.seed(123)   # reproducibility

# Create the data folder if it does not exist
dir.create("data", showWarnings = FALSE)

gse <- getGEO("GSE51808", GSEMatrix = TRUE)[[1]]

expr <- exprs(gse)   # expression matrix: probes x samples
meta <- pData(gse)   # sample metadata

cat("Probes:", nrow(expr), "\n")
cat("Samples:", ncol(expr), "\n")
cat("Platform:", annotation(gse), "\n")
print(colnames(meta))
print(head(meta[, 1:5]))

# Save the raw dataset for the next scripts
saveRDS(gse, "data/GSE51808.rds")