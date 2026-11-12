# 02_check_labels.R
# Check sample labels before QC and analysis

library(GEOquery)

gse  <- readRDS("data/GSE51808.rds")
meta <- pData(gse)

# print() use karanawa, source() walin ewa penna
print(table(meta$`infection:ch1`, useNA = "ifany"))
print(table(meta$`status:ch1`, useNA = "ifany"))
print(table(meta$`tissue:ch1`, useNA = "ifany"))

# infection saha status ekata ekka
print(table(meta$`infection:ch1`, meta$`status:ch1`, useNA = "ifany"))

# unique patients kiyak da
print(length(unique(meta$`subject id:ch1`)))

# titles walin groups balanna (control thiyenawada)
print(table(sub(" Patient.*", "", meta$title)))

# expression values log2 da balanna
print(summary(as.vector(exprs(gse))))