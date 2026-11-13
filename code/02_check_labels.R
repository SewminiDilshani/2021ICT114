# 02_check_labels.R
# Check the sample labels before QC and analysis

library(GEOquery)

gse  <- readRDS("data/GSE51808.rds")
meta <- pData(gse)

# print() is used so the tables are shown when the script is run with source()
print(table(meta$`infection:ch1`, useNA = "ifany"))
print(table(meta$`status:ch1`, useNA = "ifany"))
print(table(meta$`tissue:ch1`, useNA = "ifany"))

# Infection vs status
print(table(meta$`infection:ch1`, meta$`status:ch1`, useNA = "ifany"))

# Number of unique patients/subjects
print(length(unique(meta$`subject id:ch1`)))

# Group names taken from the sample titles (to confirm the controls)
print(table(sub(" Patient.*", "", meta$title)))

# Check the value range (a max of about 15-16 means data is already log2)
print(summary(as.vector(exprs(gse))))