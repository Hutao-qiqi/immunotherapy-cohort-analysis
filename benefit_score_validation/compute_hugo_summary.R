rm(list = ls())

read_tsv <- function(path) {
  read.delim(path, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
}

write_tsv <- function(df, path) {
  write.table(df, file = path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

as_numeric <- function(x) {
  suppressWarnings(as.numeric(x))
}

# Read Hugo predictions
hugo_path <- "HUGO_logistic_prediction.grouped_by_IMvigor210_cutoffs.8q24_PDL1_TMB.tsv"
if (!file.exists(hugo_path)) {
  hugo_path <- "HUGO_logistic_prediction.8q24_PDL1_TMB.tsv"
}
stopifnot(file.exists(hugo_path))

hugo <- read_tsv(hugo_path)
hugo$pred <- as_numeric(hugo$pred)
hugo$CD274_raw <- as_numeric(hugo$CD274_raw)
hugo$TMB <- as_numeric(hugo$TMB)

cat("\n=== HUGO SUMMARY (N=", nrow(hugo), ") ===\n\n", sep = "")

# Pred distribution
cat("--- Pred distribution ---\n")
pred_summary <- summary(hugo$pred)
print(pred_summary)
cat("\nPred quantiles (0%, 25%, 33.3%, 50%, 66.7%, 75%, 100%):\n")
pred_quant <- quantile(hugo$pred, probs = c(0, 0.25, 1/3, 0.5, 2/3, 0.75, 1), na.rm = TRUE)
print(pred_quant)

# CD274 distribution
cat("\n--- CD274_raw distribution ---\n")
cd274_summary <- summary(hugo$CD274_raw)
print(cd274_summary)
cat("\nCD274_raw quantiles:\n")
cd274_quant <- quantile(hugo$CD274_raw, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
print(cd274_quant)

# TMB distribution
cat("\n--- TMB distribution ---\n")
tmb_summary <- summary(hugo$TMB)
print(tmb_summary)

# Group counts
cat("\n--- BenefitGroup counts ---\n")
grp_tab <- table(hugo$BenefitGroup, useNA = "ifany")
print(grp_tab)

# Export detailed table
hugo_detail <- hugo[, c("Cohort", "SampleID", "PatientID", "binaryResponse1", 
                        "eightq24_status", "PDL1_ge5", "CD274_raw", "TMB", 
                        "pred", "BenefitGroup")]
hugo_detail <- hugo_detail[order(hugo_detail$pred), ]
write_tsv(hugo_detail, "HUGO_detailed_predictions_sorted_by_pred.8q24_PDL1_TMB.tsv")
cat("\nWrote: HUGO_detailed_predictions_sorted_by_pred.8q24_PDL1_TMB.tsv\n")

# Summary stats table
stats_df <- data.frame(
  Variable = c("pred", "CD274_raw", "TMB"),
  Min = c(min(hugo$pred, na.rm = TRUE), min(hugo$CD274_raw, na.rm = TRUE), min(hugo$TMB, na.rm = TRUE)),
  Q25 = c(quantile(hugo$pred, 0.25, na.rm = TRUE), quantile(hugo$CD274_raw, 0.25, na.rm = TRUE), quantile(hugo$TMB, 0.25, na.rm = TRUE)),
  Median = c(median(hugo$pred, na.rm = TRUE), median(hugo$CD274_raw, na.rm = TRUE), median(hugo$TMB, na.rm = TRUE)),
  Q75 = c(quantile(hugo$pred, 0.75, na.rm = TRUE), quantile(hugo$CD274_raw, 0.75, na.rm = TRUE), quantile(hugo$TMB, 0.75, na.rm = TRUE)),
  Max = c(max(hugo$pred, na.rm = TRUE), max(hugo$CD274_raw, na.rm = TRUE), max(hugo$TMB, na.rm = TRUE)),
  stringsAsFactors = FALSE
)
write_tsv(stats_df, "HUGO_distribution_summary_stats.8q24_PDL1_TMB.tsv")
cat("Wrote: HUGO_distribution_summary_stats.8q24_PDL1_TMB.tsv\n")
