rm(list = ls())

source(file.path("..", "utils", "rscript_utils.R"))
set_working_dir_to_script()

read_tsv <- function(path) {
  read.delim(path, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
}

write_tsv <- function(df, path) {
  write.table(df, file = path, sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
}

as_numeric <- function(x) {
  suppressWarnings(as.numeric(x))
}

args <- commandArgs(trailingOnly = TRUE)

in_path <- if (length(args) >= 1 && nzchar(args[1])) {
  args[1]
} else {
  "ValidationMerged_logistic_prediction.grouped_by_IMvigor210_cutoffs.8q24_PDL1_TMB.tsv"
}

out_path <- if (length(args) >= 2 && nzchar(args[2])) {
  args[2]
} else {
  "Low_vs_High_fisher.Response1.ValidationMerged_and_by_cohort.8q24_PDL1_TMB.tsv"
}

merged_label <- if (length(args) >= 3 && nzchar(args[3])) {
  args[3]
} else {
  "ValidationMerged"
}
stopifnot(file.exists(in_path))

df <- read_tsv(in_path)
df$binaryResponse1 <- as_numeric(df$binaryResponse1)

run_fisher <- function(d, label) {
  d <- d[d$BenefitGroup %in% c("Low", "High"), ]
  d <- d[complete.cases(d[, c("BenefitGroup", "binaryResponse1")]), ]

  # rows: group; cols: response
  tab <- table(
    BenefitGroup = factor(d$BenefitGroup, levels = c("Low", "High")),
    Benefit = factor(d$binaryResponse1, levels = c(0, 1))
  )

  # Ensure 2x2
  if (!all(dim(tab) == c(2, 2))) {
    stop("Unexpected table shape for ", label, ": ", paste(dim(tab), collapse = "x"))
  }

  ft <- fisher.test(tab)

  low_n <- sum(tab["Low", ])
  high_n <- sum(tab["High", ])
  low_benefit <- tab["Low", "1"]
  high_benefit <- tab["High", "1"]

  # If either group is empty, the comparison is not estimable.
  if (low_n == 0 || high_n == 0) {
    ft <- NULL
  }

  data.frame(
    Cohort = label,
    N_Low = low_n,
    Benefit_Low = as.integer(low_benefit),
    ResponseRate_Low = if (low_n > 0) as.numeric(low_benefit) / as.numeric(low_n) else NA_real_,
    N_High = high_n,
    Benefit_High = as.integer(high_benefit),
    ResponseRate_High = if (high_n > 0) as.numeric(high_benefit) / as.numeric(high_n) else NA_real_,
    OR_High_vs_Low = if (!is.null(ft)) unname(ft$estimate) else NA_real_,
    OR_CI95_low = if (!is.null(ft)) unname(ft$conf.int[1]) else NA_real_,
    OR_CI95_high = if (!is.null(ft)) unname(ft$conf.int[2]) else NA_real_,
    P_value = if (!is.null(ft)) unname(ft$p.value) else NA_real_,
    stringsAsFactors = FALSE
  )
}

out <- list()
out[[merged_label]] <- run_fisher(df, merged_label)

if ("Cohort" %in% names(df)) {
  for (c in sort(unique(df$Cohort))) {
    out[[c]] <- run_fisher(df[df$Cohort == c, ], c)
  }
}

res <- do.call(rbind, out)
write_tsv(res, out_path)

cat("Wrote: ", out_path, "\n", sep = "")
