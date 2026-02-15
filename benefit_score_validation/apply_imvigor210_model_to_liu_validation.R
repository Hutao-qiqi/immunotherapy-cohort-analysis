rm(list = ls())

source(file.path("..", "utils", "rscript_utils.R"))
set_working_dir_to_script()

suppressPackageStartupMessages({
  library(pROC)
})

read_tsv <- function(path) {
  read.delim(path, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
}

write_tsv <- function(df, path) {
  write.table(df, file = path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

as_numeric <- function(x) {
  suppressWarnings(as.numeric(x))
}

make_group <- function(pred, c1, c2) {
  group <- rep(NA_character_, length(pred))
  group[!is.na(pred) & pred <= c1] <- "Low"
  group[!is.na(pred) & pred > c1 & pred <= c2] <- "Inter"
  group[!is.na(pred) & pred > c2] <- "High"
  factor(group, levels = c("Low", "Inter", "High"))
}

resp_by_group <- function(df, cohort_name) {
  tab <- aggregate(binaryResponse1 ~ BenefitGroup, data = df, FUN = function(x) {
    c(n = length(x), responders = sum(x == 1, na.rm = TRUE))
  })
  data.frame(
    Cohort = cohort_name,
    BenefitGroup = tab$BenefitGroup,
    N = tab$binaryResponse1[, "n"],
    Responders = tab$binaryResponse1[, "responders"],
    ResponseRate = tab$binaryResponse1[, "responders"] / tab$binaryResponse1[, "n"],
    stringsAsFactors = FALSE
  )
}

# ------------------------------
# Inputs (all are already in this folder)
# ------------------------------
model_path <- "model_IMvigor210_8q24_PDL1_TMB.rds"
cutoff_path <- "BenefitScore_cutoffs_training_IMvigor210.8q24_PDL1_TMB.tsv"

# Liu TMB in this workspace comes from the SKCM (Liu) table as mutation counts.
# Convert to mut/Mb to match the IMvigor210 model's TMB unit.
# If you have a cohort-specific callable territory size, update this value.
tmb_denominator_mb <- 38

# Liu cohort base (has TMB + eightq24_status)
liu_base_path <- "SKCM_logistic_prediction.8q24_TMB.tsv"

# Provides CD274 TPM proxy for PD-L1 (we will binarize as >=5 TPM)
liu_cd274_path <- "Integrated_melanoma_with_CD274.tsv"

stopifnot(
  file.exists(model_path),
  file.exists(cutoff_path),
  file.exists(liu_base_path),
  file.exists(liu_cd274_path)
)

model <- readRDS(model_path)
cutoff_df <- read_tsv(cutoff_path)
cutoff_low_mid <- as_numeric(cutoff_df$Cutoff_LowMid[1])
cutoff_mid_high <- as_numeric(cutoff_df$Cutoff_MidHigh[1])

# ------------------------------
# Assemble Liu validation dataset
# ------------------------------
liu_base <- read_tsv(liu_base_path)
liu_cd274_all <- read_tsv(liu_cd274_path)
liu_cd274 <- liu_cd274_all[liu_cd274_all$Cohort %in% c("Liu"), c("SampleID", "CD274_raw", "CD274_log2p1")]

liu_base$TMB <- as_numeric(liu_base$TMB)
liu_base$binaryResponse1 <- as_numeric(liu_base$binaryResponse1)
liu_base$eightq24_status <- as_numeric(liu_base$eightq24_status)

liu_base$TMB_mutMb <- liu_base$TMB / tmb_denominator_mb

liu_merged <- merge(
  liu_base,
  liu_cd274,
  by = "SampleID",
  all.x = TRUE,
  sort = FALSE
)

liu_merged$CD274_raw <- as_numeric(liu_merged$CD274_raw)

# PD-L1 proxy: CD274 TPM >= 5 -> PDL1_ge5 = 1
liu_merged$PDL1_ge5 <- NA_integer_
liu_merged$PDL1_ge5[!is.na(liu_merged$CD274_raw) & liu_merged$CD274_raw < 5] <- 0L
liu_merged$PDL1_ge5[!is.na(liu_merged$CD274_raw) & liu_merged$CD274_raw >= 5] <- 1L

liu_df <- data.frame(
  SampleID = liu_merged$SampleID,
  binaryResponse1 = liu_merged$binaryResponse1,
  eightq24_status = liu_merged$eightq24_status,
  PDL1_ge5 = liu_merged$PDL1_ge5,
  TMB = liu_merged$TMB_mutMb,
  BR = liu_merged$BR,
  CD274_raw = liu_merged$CD274_raw,
  TMB_count = liu_merged$TMB,
  stringsAsFactors = FALSE
)

liu_cc <- liu_df[complete.cases(liu_df[, c("binaryResponse1", "eightq24_status", "PDL1_ge5", "TMB")]), ]

# ------------------------------
# Apply fixed IMvigor210 model + cutoffs
# ------------------------------
liu_cc$pred <- as.numeric(predict(model, newdata = liu_cc, type = "response"))
liu_cc$BenefitGroup <- make_group(liu_cc$pred, cutoff_low_mid, cutoff_mid_high)

write_tsv(liu_cc, "LIU_logistic_prediction.8q24_PDL1_TMB.tsv")
write_tsv(liu_cc, "LIU_logistic_prediction.grouped_by_IMvigor210_cutoffs.8q24_PDL1_TMB.tsv")

# ------------------------------
# Response rate summary + counts
# ------------------------------
liu_resp <- resp_by_group(liu_cc, "Liu")
write_tsv(liu_resp, "LIU_response_by_group.8q24_PDL1_TMB.tsv")

# Also export a combined training+validation response-rate table (for reviewer cut-off external validation)
train_grouped_path <- "IMvigor210_logistic_prediction.grouped_by_training_cutoffs.8q24_PDL1_TMB.tsv"
if (file.exists(train_grouped_path)) {
  train_grouped <- read_tsv(train_grouped_path)
  train_grouped$binaryResponse1 <- as_numeric(train_grouped$binaryResponse1)
  train_grouped$BenefitGroup <- factor(train_grouped$BenefitGroup, levels = c("Low", "Inter", "High"))
  train_resp <- resp_by_group(train_grouped, "IMvigor210")
  write_tsv(rbind(train_resp, liu_resp), "ResponseRate_by_group_IMvigor210_train_LIU_validation.8q24_PDL1_TMB.tsv")
}

liu_counts <- data.frame(
  Cohort = "Liu",
  N_total_complete_cases = nrow(liu_cc),
  N_Low = sum(liu_cc$BenefitGroup == "Low", na.rm = TRUE),
  N_Inter = sum(liu_cc$BenefitGroup == "Inter", na.rm = TRUE),
  N_High = sum(liu_cc$BenefitGroup == "High", na.rm = TRUE),
  Cutoff_LowMid = cutoff_low_mid,
  Cutoff_MidHigh = cutoff_mid_high,
  stringsAsFactors = FALSE
)
write_tsv(liu_counts, "LIU_validation_counts_by_group.8q24_PDL1_TMB.tsv")

# ------------------------------
# ROC / AUC
# ------------------------------
roc_liu <- roc(liu_cc$binaryResponse1, liu_cc$pred, quiet = TRUE)
auc_out <- data.frame(
  Cohort = "Liu",
  N = nrow(liu_cc),
  AUC = as.numeric(auc(roc_liu)),
  Features = "8q24_PDL1_TMB (trained on IMvigor210)",
  stringsAsFactors = FALSE
)
write_tsv(auc_out, "AUC_summary.8q24_PDL1_TMB.LIU.tsv")

# Also export a combined training+validation AUC table if training AUC summary exists
train_auc_path <- "AUC_summary.8q24_PDL1_TMB.tsv"
if (file.exists(train_auc_path)) {
  train_auc <- read_tsv(train_auc_path)
  train_auc$Features <- "8q24_PDL1_TMB (trained on IMvigor210)"
  write_tsv(rbind(train_auc, auc_out[, c("Cohort", "N", "AUC", "Features")]), "AUC_summary.8q24_PDL1_TMB.IMvigor210_train_plus_LIU.tsv")
}

cat("Done. Liu validation complete-case N = ", nrow(liu_cc), "\n", sep = "")
cat("Applied IMvigor210 cutoffs: ", format(cutoff_low_mid, digits = 10), ", ", format(cutoff_mid_high, digits = 10), "\n", sep = "")
