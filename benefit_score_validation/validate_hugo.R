#!/usr/bin/env Rscript

source(file.path("..", "utils", "rscript_utils.R"))
set_working_dir_to_script()
project_root <- get_project_root(levels_up = 4)

suppressPackageStartupMessages({
  library(dplyr)
  library(pROC)
})

# Inputs (Hugo cohort curated locally)
hugo_expr_file <- file.path(project_root, "hugo_symbol_fixed.txt")
hugo_pheno_file <- file.path(project_root, "hugo_phenoData.txt")

# Hugo `TMBIM1` in the expression-like matrix appears to be a mutation-count-like
# quantity (range can be ~20-190), not mut/Mb. Convert to mut/Mb to match the
# scale used in the IMvigor210-trained model.
EXOME_SIZE_MB <- 38

# IMvigor210-trained artifacts (frozen)
model_file <- "model_IMvigor210_8q24_PDL1_TMB.rds"
cutoff_file <- "BenefitScore_cutoffs_training_IMvigor210.8q24_PDL1_TMB.tsv"

stopifnot(file.exists(hugo_expr_file))
stopifnot(file.exists(hugo_pheno_file))
stopifnot(file.exists(model_file))
stopifnot(file.exists(cutoff_file))

model <- readRDS(model_file)
cutoffs <- read.delim(cutoff_file, sep = "\t", stringsAsFactors = FALSE)
cutoff_low_mid <- as.numeric(cutoffs$Cutoff_LowMid[1])
cutoff_mid_high <- as.numeric(cutoffs$Cutoff_MidHigh[1])

pheno <- read.delim(hugo_pheno_file, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

parse_binary_flag <- function(x) {
  # Handle values like 1, 0, "1*", "1^", blanks
  x_chr <- as.character(x)
  x_chr <- trimws(x_chr)
  x_chr <- gsub("[^0-9]", "", x_chr)
  x_chr[x_chr == ""] <- NA_character_
  as.integer(x_chr)
}

required_pheno <- c("Patient ID", "irRECIST", "Biopsy Time", "WES", "RNAseq")
missing_pheno <- setdiff(required_pheno, colnames(pheno))
if (length(missing_pheno) > 0) {
  stop("Missing required columns in hugo_phenoData: ", paste(missing_pheno, collapse = ", "))
}

pheno <- pheno %>%
  transmute(
    PatientID = as.character(`Patient ID`),
    irRECIST = as.character(irRECIST),
    BiopsyTime = as.character(`Biopsy Time`),
    WES = parse_binary_flag(WES),
    RNAseq = parse_binary_flag(RNAseq)
  ) %>%
  filter(!is.na(PatientID) & PatientID != "") %>%
  filter(tolower(BiopsyTime) %in% c("pre-treatment", "pretreatment", "pre_treatment", "pre treatment")) %>%
  filter(WES == 1 & RNAseq == 1) %>%
  distinct(PatientID, .keep_all = TRUE) %>%
  mutate(
    binaryResponse1 = case_when(
      grepl("complete response", irRECIST, ignore.case = TRUE) ~ 1L,
      grepl("partial response", irRECIST, ignore.case = TRUE) ~ 1L,
      TRUE ~ 0L
    )
  )

# Expression matrix: first column is gene symbol (rownames), columns are patients like Pt1, Pt2...
expr <- read.table(hugo_expr_file, header = TRUE, sep = "", row.names = 1,
                   stringsAsFactors = FALSE, check.names = FALSE)

needed_rows <- c("MYC", "PVT1", "CD274", "TMBIM1")
missing_rows <- setdiff(needed_rows, rownames(expr))
if (length(missing_rows) > 0) {
  stop("Missing required rows in hugo_symbol_fixed: ", paste(missing_rows, collapse = ", "))
}

common_patients <- intersect(pheno$PatientID, colnames(expr))
if (length(common_patients) == 0) {
  stop("No overlapping patient IDs between phenoData and expression matrix.")
}

pheno <- pheno %>% filter(PatientID %in% common_patients)

feat <- tibble(
  PatientID = common_patients,
  MYC = as.numeric(expr["MYC", common_patients, drop = TRUE]),
  PVT1 = as.numeric(expr["PVT1", common_patients, drop = TRUE]),
  CD274_raw = as.numeric(expr["CD274", common_patients, drop = TRUE]),
  TMB = as.numeric(expr["TMBIM1", common_patients, drop = TRUE]) / EXOME_SIZE_MB
)

med_myc <- stats::median(feat$MYC, na.rm = TRUE)
med_pvt1 <- stats::median(feat$PVT1, na.rm = TRUE)

feat <- feat %>%
  mutate(
    MYC_PVT1_Status = case_when(
      MYC > med_myc & PVT1 > med_pvt1 ~ "hi_hi",
      MYC < med_myc & PVT1 < med_pvt1 ~ "lo_lo",
      TRUE ~ NA_character_
    ),
    eightq24_status = case_when(
      MYC_PVT1_Status == "hi_hi" ~ 1L,
      MYC_PVT1_Status == "lo_lo" ~ 0L,
      TRUE ~ NA_integer_
    ),
    PDL1_ge5 = ifelse(!is.na(CD274_raw) & CD274_raw >= 5, 1L, 0L)
  )

df <- pheno %>%
  inner_join(feat, by = "PatientID") %>%
  mutate(
    Cohort = "Hugo",
    SampleID = paste0("hugo_", PatientID)
  ) %>%
  select(Cohort, SampleID, PatientID, binaryResponse1, MYC_PVT1_Status, eightq24_status,
         PDL1_ge5, CD274_raw, TMB, irRECIST, BiopsyTime)

# Keep complete cases for model inputs; drop intermediate MYC/PVT1 combinations
df_model <- df %>%
  filter(!is.na(eightq24_status)) %>%
  filter(!is.na(PDL1_ge5) & !is.na(TMB) & !is.na(binaryResponse1))

if (nrow(df_model) == 0) {
  stop("No complete cases for Hugo after filtering; cannot apply model.")
}

df_model$pred <- as.numeric(stats::predict(model, newdata = df_model, type = "response"))

df_model$BenefitGroup <- dplyr::case_when(
  df_model$pred < cutoff_low_mid ~ "Low",
  df_model$pred < cutoff_mid_high ~ "Inter",
  TRUE ~ "High"
)

out_pred <- "HUGO_logistic_prediction.8q24_PDL1_TMB.tsv"
out_grouped <- "HUGO_logistic_prediction.grouped_by_IMvigor210_cutoffs.8q24_PDL1_TMB.tsv"
out_resp <- "HUGO_response_by_group.8q24_PDL1_TMB.tsv"
out_counts <- "HUGO_validation_counts_by_group.8q24_PDL1_TMB.tsv"
out_auc <- "AUC_summary.8q24_PDL1_TMB.HUGO.tsv"

write.table(df_model, out_pred, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(df_model, out_grouped, sep = "\t", quote = FALSE, row.names = FALSE)

resp_by_group <- df_model %>%
  group_by(Cohort, BenefitGroup) %>%
  summarise(
    N = dplyr::n(),
    Responders = sum(binaryResponse1 == 1, na.rm = TRUE),
    ResponseRate = Responders / N,
    .groups = "drop"
  )

counts_by_group <- df_model %>%
  count(Cohort, BenefitGroup, name = "N")

write.table(resp_by_group, out_resp, sep = "\t", quote = FALSE, row.names = FALSE)
write.table(counts_by_group, out_counts, sep = "\t", quote = FALSE, row.names = FALSE)

roc_obj <- pROC::roc(response = df_model$binaryResponse1, predictor = df_model$pred, quiet = TRUE)
auc_val <- as.numeric(pROC::auc(roc_obj))
auc_df <- data.frame(Cohort = "Hugo", N = nrow(df_model), AUC = auc_val, stringsAsFactors = FALSE)
write.table(auc_df, out_auc, sep = "\t", quote = FALSE, row.names = FALSE)

message("Hugo complete-case N = ", nrow(df_model))
message("Hugo AUC = ", signif(auc_val, 6))
message("Wrote: ", out_grouped)
