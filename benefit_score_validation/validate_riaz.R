rm(list = ls())

source(file.path("..", "utils", "rscript_utils.R"))
set_working_dir_to_script()
project_root <- get_project_root(levels_up = 4)

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
# Inputs
# ------------------------------
model_path <- "model_IMvigor210_8q24_PDL1_TMB.rds"
cutoff_path <- "BenefitScore_cutoffs_training_IMvigor210.8q24_PDL1_TMB.tsv"

annotation_path <- file.path(project_root, "MYC_PVT1_annotation.txt")
ra_pheno_path <- file.path(project_root, "ra_phenodata.txt")
ra_expr_path <- file.path(project_root, "ra_symbol.txt")

# Riaz Mutation Load is mutation counts; convert to mut/Mb (same convention as other cohorts)
tmb_denominator_mb <- 38

stopifnot(
  file.exists(model_path),
  file.exists(cutoff_path),
  file.exists(annotation_path),
  file.exists(ra_pheno_path),
  file.exists(ra_expr_path)
)

model <- readRDS(model_path)
cutoff_df <- read_tsv(cutoff_path)
cutoff_low_mid <- as_numeric(cutoff_df$Cutoff_LowMid[1])
cutoff_mid_high <- as_numeric(cutoff_df$Cutoff_MidHigh[1])

# ------------------------------
# 8q24 status (from annotation)
# ------------------------------
anno <- read_tsv(annotation_path)
if (!all(c("Sample", "MYC_PVT1_Status") %in% names(anno))) {
  stop("MYC_PVT1_annotation.txt must contain columns: Sample, MYC_PVT1_Status")
}
anno$status_bin <- NA_integer_
anno$status_bin[anno$MYC_PVT1_Status == "lo_lo"] <- 0L
anno$status_bin[anno$MYC_PVT1_Status == "hi_hi"] <- 1L

anno_ra_samples <- anno$Sample[grepl("^ra_", anno$Sample)]
match_pre_sample <- function(pt) {
  if (is.na(pt) || pt == "") return(NA_character_)
  pat <- paste0("^ra_", pt, "_Pre")
  hits <- anno_ra_samples[grepl(pat, anno_ra_samples)]
  if (length(hits) == 0) return(NA_character_)
  hits[1]
}

# ------------------------------
# CD274 expression extraction (from ra_symbol.txt)
# ------------------------------
cat("Reading expression matrix header + locating CD274...\n")
expr_lines <- readLines(ra_expr_path, warn = FALSE)
if (length(expr_lines) < 2) stop("ra_symbol.txt looks empty or malformed.")

hdr <- strsplit(expr_lines[1], "\t", fixed = TRUE)[[1]]
if (length(hdr) >= 1 && hdr[1] == "") hdr <- hdr[-1]

idx <- grep("^CD274\\t", expr_lines)
if (length(idx) == 0) stop("Cannot find CD274 in ra_symbol.txt")
if (length(idx) > 1) warning("Multiple CD274 rows found; using the first match.")

cd274_fields <- strsplit(expr_lines[idx[1]], "\t", fixed = TRUE)[[1]]
cd274_vals <- as_numeric(cd274_fields[-1])
if (length(cd274_vals) != length(hdr)) {
  stop("CD274 row length does not match header length: ", length(cd274_vals), " vs ", length(hdr))
}
names(cd274_vals) <- hdr

# ------------------------------
# Assemble Riaz validation dataset
# ------------------------------
ra_pheno <- read_tsv(ra_pheno_path)

if (!all(c("Patient", "Response") %in% names(ra_pheno))) {
  stop("ra_phenodata.txt must contain columns: Patient, Response")
}

mut_col <- grep("Mutation Load", names(ra_pheno), value = TRUE)
if (length(mut_col) == 0) stop("Cannot find 'Mutation Load' column in ra_phenodata.txt")
mut_col <- mut_col[1]

ra_patient <- as.character(ra_pheno[["Patient"]])
ra_resp <- as.character(ra_pheno[["Response"]])
ra_mut <- as_numeric(ra_pheno[[mut_col]])

ra_binary <- rep(NA_integer_, length(ra_resp))
ra_binary[ra_resp %in% c("CR", "PR")] <- 1L
ra_binary[!is.na(ra_resp) & !(ra_resp %in% c("CR", "PR"))] <- 0L

ra_sample <- vapply(ra_patient, match_pre_sample, FUN.VALUE = character(1))
ra_expr_sample <- sub("^ra_", "", ra_sample)

cd274_raw <- unname(cd274_vals[ra_expr_sample])

pd_l1_ge5 <- rep(NA_integer_, length(cd274_raw))
pd_l1_ge5[!is.na(cd274_raw) & cd274_raw < 5] <- 0L
pd_l1_ge5[!is.na(cd274_raw) & cd274_raw >= 5] <- 1L

eightq24_status <- anno$status_bin[match(ra_sample, anno$Sample)]

riaz_df <- data.frame(
  SampleID = ra_sample,
  PatientID = ra_patient,
  Response_raw = ra_resp,
  binaryResponse1 = ra_binary,
  eightq24_status = as_numeric(eightq24_status),
  CD274_raw = as_numeric(cd274_raw),
  PDL1_ge5 = as_numeric(pd_l1_ge5),
  TMB_count = as_numeric(ra_mut),
  stringsAsFactors = FALSE
)

riaz_df$TMB <- riaz_df$TMB_count / tmb_denominator_mb

riaz_cc <- riaz_df[complete.cases(riaz_df[, c("binaryResponse1", "eightq24_status", "PDL1_ge5", "TMB")]), ]

cat("Riaz: total rows in phenodata = ", nrow(riaz_df), "\n", sep = "")
cat("Riaz: complete-case N used in model = ", nrow(riaz_cc), "\n", sep = "")

# ------------------------------
# Apply fixed IMvigor210 model + cutoffs
# ------------------------------
riaz_cc$pred <- as.numeric(predict(model, newdata = riaz_cc, type = "response"))
riaz_cc$BenefitGroup <- make_group(riaz_cc$pred, cutoff_low_mid, cutoff_mid_high)

write_tsv(riaz_cc, "RIAZ_logistic_prediction.8q24_PDL1_TMB.tsv")
write_tsv(riaz_cc, "RIAZ_logistic_prediction.grouped_by_IMvigor210_cutoffs.8q24_PDL1_TMB.tsv")

riaz_for_tests <- riaz_cc
riaz_for_tests$Cohort <- "Riaz"
key_cols <- c("Cohort", "SampleID", "binaryResponse1", "BenefitGroup")
write_tsv(riaz_for_tests[, unique(c(key_cols, names(riaz_for_tests))), drop = FALSE],
          "RIAZ_for_tests.grouped_by_IMvigor210_cutoffs.8q24_PDL1_TMB.tsv")

# ------------------------------
# Response rate summary + counts
# ------------------------------
riaz_resp <- resp_by_group(riaz_cc, "Riaz")
write_tsv(riaz_resp, "RIAZ_response_by_group.8q24_PDL1_TMB.tsv")

riaz_counts <- data.frame(
  Cohort = "Riaz",
  N_total_complete_cases = nrow(riaz_cc),
  N_Low = sum(riaz_cc$BenefitGroup == "Low", na.rm = TRUE),
  N_Inter = sum(riaz_cc$BenefitGroup == "Inter", na.rm = TRUE),
  N_High = sum(riaz_cc$BenefitGroup == "High", na.rm = TRUE),
  Cutoff_LowMid = cutoff_low_mid,
  Cutoff_MidHigh = cutoff_mid_high,
  stringsAsFactors = FALSE
)
write_tsv(riaz_counts, "RIAZ_validation_counts_by_group.8q24_PDL1_TMB.tsv")

# ------------------------------
# ROC / AUC
# ------------------------------
roc_riaz <- roc(riaz_cc$binaryResponse1, riaz_cc$pred, quiet = TRUE)
auc_out <- data.frame(
  Cohort = "Riaz",
  N = nrow(riaz_cc),
  AUC = as.numeric(auc(roc_riaz)),
  Features = "8q24_PDL1_TMB (trained on IMvigor210)",
  stringsAsFactors = FALSE
)
write_tsv(auc_out, "AUC_summary.8q24_PDL1_TMB.RIAZ.tsv")

cat("Done. Riaz validation complete-case N = ", nrow(riaz_cc), "\n", sep = "")
cat("Applied IMvigor210 cutoffs: ", format(cutoff_low_mid, digits = 10), ", ", format(cutoff_mid_high, digits = 10), "\n", sep = "")
