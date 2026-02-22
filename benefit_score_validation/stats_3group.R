#!/usr/bin/env Rscript

source(file.path("..", "utils", "rscript_utils.R"))
set_working_dir_to_script()

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)

in_file <- if (length(args) >= 1 && nzchar(args[1])) {
  args[1]
} else {
  "ValidationMerged_logistic_prediction.grouped_by_IMvigor210_cutoffs.8q24_PDL1_TMB.tsv"
}

out_file <- if (length(args) >= 2 && nzchar(args[2])) {
  args[2]
} else {
  "Low_Inter_High_3group_tests.Response1.ValidationMerged_and_by_cohort.8q24_PDL1_TMB.tsv"
}

merged_label <- if (length(args) >= 3 && nzchar(args[3])) {
  args[3]
} else {
  "ValidationMerged"
}

stopifnot(file.exists(in_file))

df <- read.delim(in_file, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

required_cols <- c("Cohort", "BenefitGroup", "binaryResponse1")
missing_cols <- setdiff(required_cols, colnames(df))
if (length(missing_cols) > 0) {
  stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
}

df <- df %>%
  mutate(
    Cohort = as.character(Cohort),
    BenefitGroup = factor(BenefitGroup, levels = c("Low", "Inter", "High"), ordered = TRUE),
    binaryResponse1 = as.integer(binaryResponse1)
  ) %>%
  filter(!is.na(Cohort), !is.na(BenefitGroup), !is.na(binaryResponse1))

make_3x2_table <- function(dat) {
  tab <- dat %>%
    mutate(Response = ifelse(binaryResponse1 == 1, "Benefit", "No")) %>%
    count(BenefitGroup, Response, name = "n") %>%
    tidyr::complete(BenefitGroup, Response, fill = list(n = 0)) %>%
    tidyr::pivot_wider(names_from = Response, values_from = n, values_fill = 0) %>%
    arrange(BenefitGroup)

  m <- as.matrix(tab[, c("No", "Benefit")])
  rownames(m) <- as.character(tab$BenefitGroup)
  m
}

compute_tests <- function(m) {
  # m: rows Low/Inter/High, cols No/Benefit
  # Fisher exact test for RxC
  fisher_p <- tryCatch(stats::fisher.test(m)$p.value, error = function(e) NA_real_)
  chi_p <- tryCatch(stats::chisq.test(m, correct = FALSE)$p.value, error = function(e) NA_real_)

  # Trend test (Cochran-Armitage) via prop.trend.test
  successes <- m[, "Benefit"]
  totals <- rowSums(m)
  trend_p <- tryCatch(stats::prop.trend.test(successes, totals)$p.value, error = function(e) NA_real_)

  list(
    fisher_p = fisher_p,
    chisq_p = chi_p,
    trend_p = trend_p
  )
}

summarize_rates <- function(m) {
  totals <- rowSums(m)
  responders <- m[, "Benefit"]
  rates <- ifelse(totals > 0, responders / totals, NA_real_)
  data.frame(
    BenefitGroup = factor(rownames(m), levels = c("Low", "Inter", "High"), ordered = TRUE),
    N = as.integer(totals),
    Responders = as.integer(responders),
    ResponseRate = as.numeric(rates),
    stringsAsFactors = FALSE
  )
}

run_one <- function(dat, cohort_label) {
  m <- make_3x2_table(dat)
  tests <- compute_tests(m)
  rates <- summarize_rates(m)

  # flatten output (one row per cohort)
  out <- data.frame(
    Cohort = cohort_label,
    N_Low = rates$N[rates$BenefitGroup == "Low"],
    Benefit_Low = rates$Responders[rates$BenefitGroup == "Low"],
    ResponseRate_Low = rates$ResponseRate[rates$BenefitGroup == "Low"],
    N_Inter = rates$N[rates$BenefitGroup == "Inter"],
    Benefit_Inter = rates$Responders[rates$BenefitGroup == "Inter"],
    ResponseRate_Inter = rates$ResponseRate[rates$BenefitGroup == "Inter"],
    N_High = rates$N[rates$BenefitGroup == "High"],
    Benefit_High = rates$Responders[rates$BenefitGroup == "High"],
    ResponseRate_High = rates$ResponseRate[rates$BenefitGroup == "High"],
    P_Fisher_3x2 = tests$fisher_p,
    P_ChiSquare_3x2 = tests$chisq_p,
    P_Trend_Low_to_High = tests$trend_p,
    stringsAsFactors = FALSE
  )
  out
}

cohorts <- sort(unique(df$Cohort))

res_list <- list()
res_list[[merged_label]] <- run_one(df, merged_label)
for (co in cohorts) {
  dat <- df %>% filter(Cohort == co)
  res_list[[co]] <- run_one(dat, co)
}

res <- bind_rows(res_list)

write.table(res, out_file, sep = "\t", quote = FALSE, row.names = FALSE)
message("Wrote: ", out_file)
