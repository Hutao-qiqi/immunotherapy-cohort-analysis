rm(list = ls())

suppressPackageStartupMessages({
  library(ggplot2)
})

in_path <- "HUGO_logistic_prediction.grouped_by_IMvigor210_cutoffs.8q24_PDL1_TMB.tsv"
stopifnot(file.exists(in_path))

df <- read.delim(in_path, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)

df$BenefitGroup <- factor(df$BenefitGroup, levels = c("Low", "Inter", "High"))
df$ResponseBin <- ifelse(df$binaryResponse1 == 1, "Benefit", "No benefit")
df$ResponseBin <- factor(df$ResponseBin, levels = c("No benefit", "Benefit"))

# Continuous association test in HUGO: glm(binaryResponse1 ~ pred)
glm_df <- df[!is.na(df$binaryResponse1) & !is.na(df$pred), ]
glm_fit <- tryCatch(
  glm(binaryResponse1 ~ pred, data = glm_df, family = binomial(link = "logit")),
  error = function(e) NULL
)

p_pred <- NA_real_
beta_pred <- NA_real_
se_pred <- NA_real_

if (!is.null(glm_fit)) {
  s <- summary(glm_fit)$coefficients
  if ("pred" %in% rownames(s)) {
    beta_pred <- unname(s["pred", "Estimate"])
    se_pred <- unname(s["pred", "Std. Error"])
    p_pred <- unname(s["pred", "Pr(>|z|)"])
  }
}

or_per_0_1 <- ifelse(is.na(beta_pred), NA_real_, exp(beta_pred * 0.1))

subtitle_txt <- paste0(
  "HUGO logistic: binaryResponse1 ~ pred; p(pred)=",
  ifelse(is.na(p_pred), "NA", format(p_pred, digits = 3, scientific = TRUE)),
  ifelse(is.na(or_per_0_1), "", paste0(", OR per 0.1=", format(or_per_0_1, digits = 3)))
)

glm_out <- data.frame(
  Cohort = "HUGO",
  Model = "glm(binaryResponse1 ~ pred, binomial logit)",
  N = nrow(glm_df),
  Beta_pred = beta_pred,
  SE_pred = se_pred,
  P_pred = p_pred,
  OR_per_0.1 = or_per_0_1,
  stringsAsFactors = FALSE
)

write.table(
  glm_out,
  file = "HUGO_glm_binaryResponse1_on_pred.8q24_PDL1_TMB.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE
)

# Counts for plot (all 3 groups)
ct_all <- table(df$BenefitGroup, df$ResponseBin)

# Fractions for stacked bars
agg <- as.data.frame(ct_all)
colnames(agg) <- c("BenefitGroup", "ResponseBin", "Count")

n_by_group <- aggregate(Count ~ BenefitGroup, data = agg, sum)
colnames(n_by_group) <- c("BenefitGroup", "N")
agg <- merge(agg, n_by_group, by = "BenefitGroup", all.x = TRUE)
agg$Fraction <- ifelse(agg$N == 0, 0, agg$Count / agg$N)

# Annotate Benefit% inside the Benefit segment
benefit_df <- subset(agg, ResponseBin == "Benefit")
no_df <- subset(agg, ResponseBin == "No benefit")[, c("BenefitGroup", "Fraction")]
colnames(no_df)[2] <- "NoFraction"
benefit_df <- merge(benefit_df, no_df, by = "BenefitGroup", all.x = TRUE)
benefit_df$label <- ifelse(benefit_df$N == 0, "", paste0(round(100 * benefit_df$Fraction), "%"))
benefit_df$y <- benefit_df$NoFraction + benefit_df$Fraction / 2

p <- ggplot(agg, aes(x = BenefitGroup, y = Fraction, fill = ResponseBin)) +
  geom_col(width = 0.85, color = "white") +
  scale_fill_manual(values = c("No benefit" = "#7a2b7f", "Benefit" = "#358043")) +
  scale_y_continuous(limits = c(0, 1), expand = c(0, 0)) +
  labs(
    title = "Benefit Score (validation: HUGO)",
    subtitle = subtitle_txt,
    x = NULL,
    y = "Fraction",
    fill = NULL
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

p <- p + geom_text(
  data = benefit_df,
  aes(x = BenefitGroup, y = y, label = label),
  inherit.aes = FALSE,
  size = 4
)

ggsave(
  "HUGO_BenefitScore_validation_fraction.8q24_PDL1_TMB.png",
  plot = p,
  width = 4.2,
  height = 3.4,
  dpi = 300
)

ggsave(
  "HUGO_BenefitScore_validation_fraction.8q24_PDL1_TMB.pdf",
  plot = p,
  width = 4.2,
  height = 3.4
)

cat("HUGO counts by group (Benefit/No benefit), all groups:\n")
print(ct_all)

cat("\nHUGO logistic regression: binaryResponse1 ~ pred\n")
if (is.null(glm_fit)) {
  cat("glm failed\n")
} else {
  print(summary(glm_fit)$coefficients)
}
