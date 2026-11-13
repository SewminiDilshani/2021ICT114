# 07_classification.R
# SVM classification: Control vs Dengue (whole blood expression)
# Feature selection (limma) is done INSIDE each cross-validation fold,
# so the test samples never influence which genes are chosen.

library(limma)
library(e1071)
library(ggplot2)

dir.create("results", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)

set.seed(123)   # reproducibility

d    <- readRDS("data/clean_data.rds")
expr <- d$expr          # filtered probes x 37 samples
meta <- d$meta
feat <- d$feat
y    <- meta$condition  # factor: Control, Dengue
stopifnot(all(colnames(expr) == rownames(meta)))

# ---- Settings ----
n_feat <- 50   # number of top probes used by the SVM
k      <- 5    # folds
n_rep  <- 10   # number of times the whole cross-validation is repeated

# ---- Helper functions ----

# Select the top n probes with limma using the training samples only
select_features <- function(x, y, n) {
  design <- model.matrix(~ y)
  fit    <- eBayes(lmFit(x, design))
  tt     <- topTable(fit, coef = 2, number = n, sort.by = "P")
  rownames(tt)
}

# Stratified folds: each fold keeps a similar Control/Dengue ratio
make_folds <- function(y, k) {
  fold <- integer(length(y))
  for (lv in levels(y)) {
    idx       <- which(y == lv)
    fold[idx] <- sample(rep(1:k, length.out = length(idx)))
  }
  fold
}

# AUC from scores (rank-based, equals the Mann-Whitney statistic)
auc <- function(score, pos) {
  r  <- rank(score)
  n1 <- sum(pos)
  n0 <- sum(!pos)
  (sum(r[pos]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

# Points for the ROC curve
roc_points <- function(score, pos) {
  th <- sort(unique(c(-Inf, score, Inf)), decreasing = TRUE)
  data.frame(FPR = sapply(th, function(t) mean(score[!pos] >= t)),
             TPR = sapply(th, function(t) mean(score[pos] >= t)))
}

# ---- 1. Repeated stratified cross-validation ----
all_pred <- list()

for (r in 1:n_rep) {
  fold <- make_folds(y, k)
  
  for (f in 1:k) {
    test  <- which(fold == f)
    train <- which(fold != f)
    
    # Feature selection on the training samples only
    feats <- select_features(expr[, train], y[train], n_feat)
    
    xtr <- t(expr[feats, train])
    xte <- t(expr[feats, test, drop = FALSE])
    
    # Scale test data using the training mean and sd
    mu  <- colMeans(xtr)
    sdv <- apply(xtr, 2, sd)
    sdv[sdv == 0] <- 1
    xtr <- scale(xtr, center = mu, scale = sdv)
    xte <- scale(xte, center = mu, scale = sdv)
    
    # Class weights, because there are fewer controls than patients
    cnt <- table(y[train])
    w   <- setNames(as.numeric(max(cnt) / cnt), names(cnt))
    
    model <- svm(x = xtr, y = y[train], kernel = "linear", cost = 1,
                 class.weights = w, probability = TRUE, scale = FALSE)
    
    pred <- predict(model, xte, probability = TRUE)
    prob <- attr(pred, "probabilities")[, "Dengue"]
    
    all_pred[[length(all_pred) + 1]] <- data.frame(
      rep    = r,
      sample = colnames(expr)[test],
      truth  = as.character(y[test]),
      pred   = as.character(pred),
      prob   = as.numeric(prob),
      stringsAsFactors = FALSE)
  }
}

cv       <- do.call(rbind, all_pred)
cv$truth <- factor(cv$truth, levels = c("Control", "Dengue"))
cv$pred  <- factor(cv$pred,  levels = c("Control", "Dengue"))

# ---- 2. Performance metrics (Dengue = positive class) ----
metrics <- do.call(rbind, lapply(split(cv, cv$rep), function(df) {
  tp <- sum(df$truth == "Dengue"  & df$pred == "Dengue")
  tn <- sum(df$truth == "Control" & df$pred == "Control")
  fp <- sum(df$truth == "Control" & df$pred == "Dengue")
  fn <- sum(df$truth == "Dengue"  & df$pred == "Control")
  data.frame(Accuracy    = (tp + tn) / nrow(df),
             Sensitivity = tp / (tp + fn),
             Specificity = tn / (tn + fp),
             AUC         = auc(df$prob, df$truth == "Dengue"))
}))

summary_tab <- data.frame(Metric = colnames(metrics),
                          Mean   = round(colMeans(metrics), 3),
                          SD     = round(apply(metrics, 2, sd), 3),
                          row.names = NULL)
cat("\nCross-validated performance (", n_rep, "x", k, "-fold CV, top",
    n_feat, "probes ):\n")
print(summary_tab)
write.csv(summary_tab, "results/classification_metrics.csv", row.names = FALSE)

# Baseline: always predicting the larger class (Dengue)
cat("\nBaseline accuracy (always predict Dengue):",
    round(max(table(y)) / length(y), 3), "\n")

# ---- 3. Confusion matrix (average over the repeats; totals = 37 samples) ----
cm <- table(Observed = cv$truth, Predicted = cv$pred) / n_rep
cat("\nConfusion matrix (average per repeat):\n")
print(round(cm, 1))

cm_df <- as.data.frame(cm)
p_cm <- ggplot(cm_df, aes(Predicted, Observed, fill = Freq)) +
  geom_tile(colour = "white") +
  geom_text(aes(label = round(Freq, 1)), size = 6) +
  scale_fill_gradient(low = "white", high = "steelblue") +
  labs(title = "SVM confusion matrix (cross-validation)") +
  theme_minimal() + theme(legend.position = "none")
ggsave("figures/fig11_confusion_matrix.png", p_cm,
       width = 5, height = 4.5, dpi = 300)

# ---- 4. ROC curve (probability averaged over the repeats) ----
avg    <- aggregate(prob ~ sample + truth, data = cv, FUN = mean)
pos    <- avg$truth == "Dengue"
auc_all <- auc(avg$prob, pos)
cat("\nAUC (averaged predictions):", round(auc_all, 3), "\n")

roc_df <- roc_points(avg$prob, pos)
p_roc <- ggplot(roc_df, aes(FPR, TPR)) +
  geom_line(colour = "firebrick", linewidth = 1) +
  geom_abline(linetype = "dashed", colour = "grey50") +
  annotate("text", x = 0.7, y = 0.15,
           label = paste0("AUC = ", round(auc_all, 3))) +
  labs(title = "ROC curve: Dengue vs Control (SVM)",
       x = "False positive rate", y = "True positive rate") +
  theme_minimal()
ggsave("figures/fig10_roc_curve.png", p_roc, width = 5.5, height = 5, dpi = 300)

# ---- 5. Top genes (all 37 samples, for interpretation only) ----
# This list is NOT used to measure performance.
feats_all <- select_features(expr, y, n_feat)
sym <- as.character(feat[feats_all, "Gene Symbol"])
sym <- sub(" ///.*", "", sym)
top_genes <- data.frame(Probe = feats_all,
                        Gene  = ifelse(is.na(sym) | sym == "", feats_all, sym))
write.csv(top_genes, "results/svm_top_genes.csv", row.names = FALSE)
cat("\nTop 15 probes used by the classifier:\n")
print(head(top_genes, 15))