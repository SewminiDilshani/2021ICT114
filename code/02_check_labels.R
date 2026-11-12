# 02_check_labels.R
# Check sample labels before QC and analysis

library(GEOquery)

gse  <- readRDS("data/GSE51808.rds")
meta <- pData(gse)

# group labels walin kiyak thiyenawada
table(meta$`infection:ch1`)
table(meta$`status:ch1`)
table(meta$`tissue:ch1`)

# infection saha status ekata ekka
table(meta$`infection:ch1`, meta$`status:ch1`)

# unique patients kiyak da
length(unique(meta$`subject id:ch1`))

# sample 10k balanna
print(meta[1:10, c("title", "infection:ch1", "status:ch1",
                   "subject id:ch1", "tissue:ch1")])

# expression values log2 da balanna (max ~ 15-16 nam log2)
summary(as.vector(exprs(gse)))