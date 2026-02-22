rm(list = ls())

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

fmt_pct <- function(x, digits = 1) {
  ifelse(is.na(x), "NA", paste0(format(round(100 * x, digits = digits), nsmall = digits), "%"))
}

fmt_num <- function(x, digits = 2) {
  ifelse(is.na(x), "NA", format(round(x, digits = digits), nsmall = digits))
}

standardize_df <- function(df, cohort_name) {
  needed <- c("binaryResponse1", "BenefitGroup", "eightq24_status", "PDL1_ge5", "TMB")
  missing <- setdiff(needed, colnames(df))
  if (length(missing) > 0) {
    stop(paste0("Missing columns in ", cohort_name, ": ", paste(missing, collapse = ", ")))
  }
  out <- df[, needed, drop = FALSE]
  out$binaryResponse1 <- suppressWarnings(as.integer(out$binaryResponse1))
  out$eightq24_status <- suppressWarnings(as.integer(out$eightq24_status))
  out$PDL1_ge5 <- suppressWarnings(as.integer(out$PDL1_ge5))
  out$TMB <- suppressWarnings(as.numeric(out$TMB))
  out$BenefitGroup <- as.character(out$BenefitGroup)
  out$Cohort <- cohort_name
  out
}

make_fraction_facets <- function(df_all, title_txt, subtitle_txt) {
  df_all$BenefitGroup <- factor(df_all$BenefitGroup, levels = c("Low", "Inter", "High"))
  df_all$ResponseBin <- ifelse(df_all$binaryResponse1 == 1, "Benefit", "No benefit")
  df_all$ResponseBin <- factor(df_all$ResponseBin, levels = c("No benefit", "Benefit"))

  ct_all <- as.data.frame(table(df_all$Cohort, df_all$BenefitGroup, df_all$ResponseBin))
  colnames(ct_all) <- c("Cohort", "BenefitGroup", "ResponseBin", "Count")

  n_by_group <- aggregate(Count ~ Cohort + BenefitGroup, data = ct_all, sum)
  colnames(n_by_group) <- c("Cohort", "BenefitGroup", "N")
  agg <- merge(ct_all, n_by_group, by = c("Cohort", "BenefitGroup"), all.x = TRUE)
  agg$Fraction <- ifelse(agg$N == 0, 0, agg$Count / agg$N)

  benefit_df <- subset(agg, ResponseBin == "Benefit")
  no_df <- subset(agg, ResponseBin == "No benefit")[, c("Cohort", "BenefitGroup", "Fraction")]
  colnames(no_df)[3] <- "NoFraction"
  benefit_df <- merge(benefit_df, no_df, by = c("Cohort", "BenefitGroup"), all.x = TRUE)
  benefit_df$label <- ifelse(benefit_df$N == 0, "", paste0(round(100 * benefit_df$Fraction), "%"))
  benefit_df$y <- benefit_df$NoFraction + benefit_df$Fraction / 2

  p <- ggplot(agg, aes(x = BenefitGroup, y = Fraction, fill = ResponseBin)) +
    geom_col(width = 0.85, color = "white") +
    scale_fill_manual(values = c("No benefit" = "#7a2b7f", "Benefit" = "#358043")) +
    scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
    labs(
      title = title_txt,
      subtitle = subtitle_txt,
      x = NULL,
      y = "Fraction",
      fill = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "right",
      strip.text = element_text(face = "bold")
    ) +
    facet_wrap(~Cohort, nrow = 1)

  p <- p + geom_text(
    data = benefit_df,
    aes(x = BenefitGroup, y = y, label = label),
    inherit.aes = FALSE,
    size = 3.6
  )

  list(plot = p, agg = agg)
}

make_cohort_summary_table_plot <- function(df_all) {
  df_all <- df_all[!is.na(df_all$binaryResponse1), ]

  cohorts <- sort(unique(df_all$Cohort))
  rows <- list()
  for (coh in cohorts) {
    d <- df_all[df_all$Cohort == coh, ]
    n_total <- nrow(d)
    orr <- ifelse(n_total == 0, NA_real_, mean(d$binaryResponse1 == 1, na.rm = TRUE))
    amp <- ifelse(n_total == 0, NA_real_, mean(d$eightq24_status == 1, na.rm = TRUE))
    pdl1 <- ifelse(n_total == 0, NA_real_, mean(d$PDL1_ge5 == 1, na.rm = TRUE))
    tmb_med <- ifelse(n_total == 0, NA_real_, median(d$TMB, na.rm = TRUE))
    rows[[length(rows) + 1]] <- data.frame(
      Cohort = coh,
      N = n_total,
      ORR = orr,
      `8q24+` = amp,
      `PDL1_ge5` = pdl1,
      TMB_median = tmb_med,
      check.names = FALSE,
      stringsAsFactors = FALSE
    )
  }
  sum_df <- do.call(rbind, rows)

  long <- data.frame(
    Cohort = rep(sum_df$Cohort, times = 5),
    Metric = rep(c("N", "ORR", "8q24+", "PDL1_ge5", "TMB_median"), each = nrow(sum_df)),
    Value = NA_character_,
    stringsAsFactors = FALSE
  )

  for (i in seq_len(nrow(sum_df))) {
    coh <- sum_df$Cohort[i]
    long$Value[long$Cohort == coh & long$Metric == "N"] <- as.character(sum_df$N[i])
    long$Value[long$Cohort == coh & long$Metric == "ORR"] <- fmt_pct(sum_df$ORR[i], digits = 1)
    long$Value[long$Cohort == coh & long$Metric == "8q24+"] <- fmt_pct(sum_df$`8q24+`[i], digits = 1)
    long$Value[long$Cohort == coh & long$Metric == "PDL1_ge5"] <- fmt_pct(sum_df$PDL1_ge5[i], digits = 1)
    long$Value[long$Cohort == coh & long$Metric == "TMB_median"] <- fmt_num(sum_df$TMB_median[i], digits = 2)
  }

  long$Metric <- factor(long$Metric, levels = c("N", "ORR", "8q24+", "PDL1_ge5", "TMB_median"))
  long$Cohort <- factor(long$Cohort, levels = unique(sum_df$Cohort))

  p <- ggplot(long, aes(x = Cohort, y = Metric)) +
    geom_tile(fill = "grey95", color = "white") +
    geom_text(aes(label = Value), size = 3.6) +
    theme_void(base_size = 11) +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, face = "bold"),
      axis.text.y = element_text(face = "bold"),
      plot.title = element_text(face = "bold")
    ) +
    labs(title = "Cohort characteristics (complete cases)")

  p
}

# Inputs (kept files)
in_msk <- "MSK_NSCLC_logistic_prediction.grouped_by_IMvigor210_cutoffs.8q24_PDL1_TMB.tsv"
in_liu <- "LIU_logistic_prediction.grouped_by_IMvigor210_cutoffs.8q24_PDL1_TMB.tsv"
in_riaz <- "RIAZ_logistic_prediction.grouped_by_IMvigor210_cutoffs.8q24_PDL1_TMB.tsv"
in_hugo <- "HUGO_logistic_prediction.grouped_by_IMvigor210_cutoffs.8q24_PDL1_TMB.tsv"
in_pooled <- "ValidationMerged_RIAZ_LIU_MSK_logistic_prediction.grouped_by_IMvigor210_cutoffs.8q24_PDL1_TMB.tsv"

stopifnot(
  file.exists(in_msk), file.exists(in_liu), file.exists(in_riaz), file.exists(in_hugo), file.exists(in_pooled)
)

msk <- standardize_df(read.delim(in_msk, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE), "MSK")
liu <- standardize_df(read.delim(in_liu, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE), "LIU")
riaz <- standardize_df(read.delim(in_riaz, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE), "RIAZ")
hugo <- standardize_df(read.delim(in_hugo, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE), "HUGO")
pooled <- standardize_df(read.delim(in_pooled, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE), "Pooled(RIAZ+LIU+MSK)")

df_all <- rbind(msk, liu, riaz, pooled, hugo)
df_all <- df_all[!is.na(df_all$BenefitGroup) & !is.na(df_all$binaryResponse1), ]

title_txt <- "Benefit Score external validation (frozen IMvigor210 cut-offs)"
subtitle_txt <- "Stacked bars show response fraction within each benefit group"

facet_obj <- make_fraction_facets(df_all, title_txt = title_txt, subtitle_txt = subtitle_txt)
table_p <- make_cohort_summary_table_plot(df_all)


# Export Panel A (stacked bars per cohort)
out_pdf_A <- file.path("..", "..", "Supplementary_FIG6A.pdf")
out_png_A <- file.path("..", "..", "Supplementary_FIG6A.png")

pdf(out_pdf_A, width = 12.5, height = 3.6)
grid.newpage()
print(facet_obj$plot)
dev.off()

png(out_png_A, width = 12.5, height = 3.6, units = "in", res = 300)
grid.newpage()
print(facet_obj$plot)
dev.off()

# Export Panel B (cohort summary table)
out_pdf_B <- file.path("..", "..", "Supplementary_FIG6B.pdf")
out_png_B <- file.path("..", "..", "Supplementary_FIG6B.png")

pdf(out_pdf_B, width = 12.5, height = 2.8)
grid.newpage()
print(table_p)
dev.off()

png(out_png_B, width = 12.5, height = 2.8, units = "in", res = 300)
grid.newpage()
print(table_p)
dev.off()

cat("Wrote:\n")
cat(out_pdf_A, "\n")
cat(out_png_A, "\n")
cat(out_pdf_B, "\n")
cat(out_png_B, "\n")
