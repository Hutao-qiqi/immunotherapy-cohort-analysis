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

# Inputs: ONLY RIAZ + LIU + MSK
paths <- c(
  "RIAZ_logistic_prediction.grouped_by_IMvigor210_cutoffs.8q24_PDL1_TMB.tsv",
  "LIU_logistic_prediction.grouped_by_IMvigor210_cutoffs.8q24_PDL1_TMB.tsv",
  "MSK_NSCLC_logistic_prediction.grouped_by_IMvigor210_cutoffs.8q24_PDL1_TMB.tsv"
)

for (p in paths) stopifnot(file.exists(p))

infer_cohort <- function(path) {
  base <- basename(path)
  if (grepl("^RIAZ_", base)) return("Riaz")
  if (grepl("^MSK_NSCLC_", base)) return("MSK_NSCLC")
  if (grepl("^LIU_", base)) return("Liu")
  sub("_logistic_prediction\\.grouped_by_IMvigor210_cutoffs\\.8q24_PDL1_TMB\\.tsv$", "", base)
}

dfs <- list()
for (p in paths) {
  cohort <- infer_cohort(p)
  df <- read_tsv(p)

  need <- c("SampleID", "binaryResponse1", "eightq24_status", "PDL1_ge5", "TMB", "pred", "BenefitGroup")
  miss <- setdiff(need, names(df))
  if (length(miss) > 0) stop("Missing columns in ", p, ": ", paste(miss, collapse = ", "))

  df$binaryResponse1 <- as_numeric(df$binaryResponse1)
  df$eightq24_status <- as_numeric(df$eightq24_status)
  df$PDL1_ge5 <- as_numeric(df$PDL1_ge5)
  df$TMB <- as_numeric(df$TMB)
  df$pred <- as_numeric(df$pred)
  df$BenefitGroup <- factor(df$BenefitGroup, levels = c("Low", "Inter", "High"))
  df$Cohort <- cohort

  key_cols <- c("Cohort", "SampleID", "binaryResponse1", "eightq24_status", "PDL1_ge5", "TMB", "pred", "BenefitGroup")
  extra_cols <- setdiff(names(df), key_cols)
  df <- df[, c(key_cols, extra_cols), drop = FALSE]

  dfs[[cohort]] <- df
}

all_cols <- Reduce(union, lapply(dfs, names))
key_order <- c("Cohort", "SampleID", "binaryResponse1", "eightq24_status", "PDL1_ge5", "TMB", "pred", "BenefitGroup")
tail_cols <- setdiff(all_cols, key_order)
all_cols <- c(key_order, tail_cols)

for (nm in names(dfs)) {
  missing_cols <- setdiff(all_cols, names(dfs[[nm]]))
  if (length(missing_cols) > 0) {
    for (mc in missing_cols) dfs[[nm]][[mc]] <- NA
  }
  dfs[[nm]] <- dfs[[nm]][, all_cols, drop = FALSE]
}

merged <- do.call(rbind, dfs)

out_prefix <- "ValidationMerged_RIAZ_LIU_MSK"
write_tsv(merged, paste0(out_prefix, "_logistic_prediction.grouped_by_IMvigor210_cutoffs.8q24_PDL1_TMB.tsv"))

cohorts <- unique(merged$Cohort)

resp_list <- list()
for (c in cohorts) resp_list[[c]] <- resp_by_group(merged[merged$Cohort == c, ], c)
resp_by_cohort <- do.call(rbind, resp_list)
write_tsv(resp_by_cohort, paste0(out_prefix, "_response_by_group_by_cohort.8q24_PDL1_TMB.tsv"))

resp_all <- resp_by_group(merged, out_prefix)
write_tsv(resp_all, paste0(out_prefix, "_response_by_group.8q24_PDL1_TMB.tsv"))

cc <- merged[complete.cases(merged[, c("binaryResponse1", "pred")]), ]
roc_obj <- roc(cc$binaryResponse1, cc$pred, quiet = TRUE)
auc_out <- data.frame(
  Cohort = out_prefix,
  N = nrow(cc),
  AUC = as.numeric(auc(roc_obj)),
  stringsAsFactors = FALSE
)
write_tsv(auc_out, paste0("AUC_summary.8q24_PDL1_TMB.", out_prefix, ".tsv"))

cat("Merged cohorts: ", paste(cohorts, collapse = ", "), "\n", sep = "")
cat(out_prefix, " N = ", nrow(cc), ", pooled AUC = ", format(auc_out$AUC[1], digits = 10), "\n", sep = "")
