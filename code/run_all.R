# run_all.R
# Run the whole project from start to finish.
# Open the 2021ICT114.Rproj file first so the working directory is the project root.

scripts <- c("code/01_download_data.R",
             "code/02_check_labels.R",
             "code/03_qc_cleaning.R",
             "code/04_eda.R",
             "code/05_limma.R",
             "code/06_enrichment.R",
             "code/07_classification.R")

for (s in scripts) {
  cat("\n========== Running", s, "==========\n")
  source(s, echo = FALSE)
}

cat("\nAll scripts finished without errors.\n")