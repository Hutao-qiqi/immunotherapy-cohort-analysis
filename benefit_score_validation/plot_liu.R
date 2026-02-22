rm(list = ls())

suppressPackageStartupMessages({
  library(pROC)
  library(ggplot2)
})

read_tsv <- function(path) {
  read.delim(path, header = TRUE, sep = "\t", check.names = FALSE, stringsAsFactors = FALSE)
}

liu_pred_path <- "LIU_logistic_prediction.grouped_by_IMvigor210_cutoffs.8q24_PDL1_TMB.tsv"
liu_resp_path <- "LIU_response_by_group.8q24_PDL1_TMB.tsv"

stopifnot(file.exists(liu_pred_path))

df <- read_tsv(liu_pred_path)
df$binaryResponse1 <- as.numeric(df$binaryResponse1)
df$pred <- as.numeric(df$pred)

# ROC
roc_obj <- roc(df$binaryResponse1, df$pred, quiet = TRUE)
auc_val <- as.numeric(auc(roc_obj))

png("LIU_ROC_8q24_PDL1_TMB.png", width = 1200, height = 1000, res = 150)
plot(roc_obj, main = sprintf("Liu validation ROC (AUC=%.3f)", auc_val), col = "blue", lwd = 3)
dev.off()

# Benefit group response barplot
if (file.exists(liu_resp_path)) {
  resp_df <- read_tsv(liu_resp_path)
  resp_df$BenefitGroup <- factor(resp_df$BenefitGroup, levels = c("Low", "Inter", "High"))
  p <- ggplot(resp_df, aes(x = BenefitGroup, y = ResponseRate)) +
    geom_col(fill = c("#d73027", "#fee08b", "#1a9850")) +
    geom_text(aes(label = sprintf("%d/%d", Responders, N)), vjust = -0.5) +
    coord_cartesian(ylim = c(0,1)) +
    labs(title = "Liu: Response rate by BenefitGroup (IMvigor210 cutoffs)", y = "Response rate (CR/PR)", x = NULL) +
    theme_bw(base_size = 12)
  ggsave("LIU_ResponseRate_by_BenefitGroup_8q24_PDL1_TMB.png", plot = p, width = 4, height = 3.2, dpi = 300)
}

cat("Plots written: LIU_ROC_8q24_PDL1_TMB.png, LIU_ResponseRate_by_BenefitGroup_8q24_PDL1_TMB.png\n")
