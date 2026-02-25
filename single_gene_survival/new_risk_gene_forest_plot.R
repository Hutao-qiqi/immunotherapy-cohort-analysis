#!/usr/bin/env Rscript
# ==============================================================================
# new_risk_gene_forest_plot.R
# 三段式森林图（顶刊风格）：基因名 | 森林图 | P值 | HRR柱状图
# 数据来源: hi_hi / lo_lo 组内 Cox 回归 (median-rank split)
# ==============================================================================

# ---- 用户配置 ----
DATA_DIR   <- "E:/data/changyuan/免疫队列/单基因生存分析"
HIHI_FILE  <- file.path(DATA_DIR, "hi_hi_cox_results_median_rank.txt")
LOLO_FILE  <- file.path(DATA_DIR, "lo_lo_cox_results_median_rank.txt")
VI_FILE    <- file.path(DATA_DIR, "cox_interaction_continuous_genomewide.tsv")

# 模式说明:
# - "continuous": 使用连续量交互项 Cox 的 HR/CI/p/VI（与 VI=exp(beta_int) 完全一致）
# - "median_rank": 使用组内 median-rank 二分的组内 Cox（更接近你原始 KM/组内 Cox 口径）
MODE <- "continuous"

TARGET_GENES <- c("AGT", "APOA2", "APOC1", "FGFR4",
                  "LGALS4", "TRIM6", "OASL", "ARG1")

parse_ci <- function(x) {
  if (is.na(x)) return(c(NA_real_, NA_real_))
  s <- as.character(x)
  parts <- strsplit(s, "-", fixed = TRUE)[[1]]
  if (length(parts) != 2) return(c(NA_real_, NA_real_))
  suppressWarnings(c(as.numeric(parts[1]), as.numeric(parts[2])))
}

if (MODE == "continuous") {
  # ---- 连续量交互项 Cox 结果 ----
  if (!file.exists(VI_FILE)) {
    stop(paste0("Missing VI_FILE: ", VI_FILE))
  }

  cont <- read.delim(VI_FILE, stringsAsFactors = FALSE)
  cont <- cont[cont$Gene %in% TARGET_GENES, ]

  # 解析 CI 字符串到 numeric lower/upper
  lo_ci <- t(sapply(cont$CI95_lo_lo, parse_ci))
  hi_ci <- t(sapply(cont$CI95_hi_hi, parse_ci))

  merged <- data.frame(
    Gene = cont$Gene,
    HR_hihi = cont$HR_gene_hi_hi,
    p_value_hihi = cont$p_hi_hi,
    CI_lower_hihi = hi_ci[, 1],
    CI_upper_hihi = hi_ci[, 2],
    HR_lolo = cont$HR_gene_lo_lo,
    p_value_lolo = cont$p_lo_lo,
    CI_lower_lolo = lo_ci[, 1],
    CI_upper_lolo = lo_ci[, 2],
    VI = cont$VI,
    p_interaction = cont$p_interaction,
    stringsAsFactors = FALSE
  )

} else if (MODE == "median_rank") {
  # ---- 组内 Cox (median-rank split) ----
  hihi <- read.delim(HIHI_FILE, stringsAsFactors = FALSE)
  lolo <- read.delim(LOLO_FILE, stringsAsFactors = FALSE)

  hihi <- hihi[hihi$Gene %in% TARGET_GENES, ]
  lolo <- lolo[lolo$Gene %in% TARGET_GENES, ]

  merged <- merge(hihi[, c("Gene", "HR", "p_value", "CI_lower", "CI_upper")],
                  lolo[, c("Gene", "HR", "p_value", "CI_lower", "CI_upper")],
                  by = "Gene", suffixes = c("_hihi", "_lolo"))

  # 使用 binary_interaction_candidates_within_group.tsv 的 VI_from_interaction（更接近 median-rank 的定义）
  vi_path <- file.path(DATA_DIR, "binary_interaction_candidates_within_group.tsv")
  if (!file.exists(vi_path)) stop(paste0("Missing VI file: ", vi_path))
  vi_data <- read.delim(vi_path, stringsAsFactors = FALSE)
  colnames(vi_data)[colnames(vi_data) == "Gene"] <- "Gene"
  vi_data$VI <- vi_data$VI_from_interaction

  merged <- merge(merged, vi_data[, c("Gene", "VI", "p_interaction")], by = "Gene", all.x = TRUE)

} else {
  stop("MODE must be 'continuous' or 'median_rank'")
}

# 按 VI 降序排列（顶部 = 最大 VI）
merged <- merged[order(-merged$VI), ]
rownames(merged) <- NULL

n <- nrow(merged)
genes <- merged$Gene

cat("Genes (sorted by VI descending):\n")
for (i in seq_len(n)) {
  cat(sprintf("  %d. %s  VI=%.3f  p_int=%.4f  hi_hi: HR=%.3f p=%.3f  lo_lo: HR=%.3f p=%.3f\n",
              i, merged$Gene[i], merged$VI[i], merged$p_interaction[i],
              merged$HR_hihi[i], merged$p_value_hihi[i],
              merged$HR_lolo[i], merged$p_value_lolo[i]))
}

# ---- 颜色定义 ----
col_hihi <- "#D95F02"   # 橙色 (hi_hi)
col_lolo <- "#377EB8"   # 蓝色 (lo_lo)

# 基因标签色块: 按 HRR 排名从暖到冷渐变
gene_rect_cols <- colorRampPalette(c("#C0392B", "#E67E22", "#27AE60", "#2980B9"))(n)

# HRR 柱状图: OrRd 渐变
hrr_bar_cols <- colorRampPalette(c("#B03030", "#E07040", "#F0B870"))(n)

# ---- 显著性星号 ----
sig_star <- function(p) {
  ifelse(p < 0.001, "***",
         ifelse(p < 0.01, "**",
                ifelse(p < 0.05, "*", "")))
}

# ---- Y 坐标计算 ----
# 每个基因占 3 个单位: y_hihi (上), y_mid (基因名), y_lolo (下)
spacing  <- 3
y_hihi   <- (n:1) * spacing + 0.45
y_lolo   <- (n:1) * spacing - 0.45
y_mid    <- (n:1) * spacing
ylim_all <- c(0.5, (n + 0.5) * spacing)

# ---- x 范围 ----
all_ci <- c(merged$CI_lower_hihi, merged$CI_upper_hihi,
            merged$CI_lower_lolo, merged$CI_upper_lolo)
xmin <- max(0.1, floor(min(all_ci) * 10) / 10 - 0.05)
xmax <- ceiling(max(all_ci) * 10) / 10 + 0.1
xlim_forest <- c(xmin, xmax)

# ---- 输出文件 ----
timestamp <- as.integer(Sys.time())
out_pdf <- file.path(DATA_DIR, paste0("new_risk_gene_forest_plot_", timestamp, ".pdf"))
out_png <- file.path(DATA_DIR, paste0("new_risk_gene_forest_plot_", timestamp, ".png"))

# ===========================================================================
#  绘图函数
# ===========================================================================
draw_forest <- function() {

  # 4 列布局: 基因名(1.6) | 森林图(5) | P值(1.3) | HRR(2.2)
  layout(matrix(1:4, nrow = 1), widths = c(1.6, 5, 1.3, 2.2))

  # ---------- 交替行阴影 ----------
  shade_col <- "#F5F5F5"

  # ========== Panel 1: 基因名 + 色块 ==========
  par(mar = c(3, 0.5, 2.5, 0), family = "serif")
  plot.new()
  plot.window(xlim = c(0, 1), ylim = ylim_all)

  # 标题
  text(0.5, ylim_all[2] + 0.6, "Gene", font = 2, cex = 1.3, xpd = TRUE, adj = 0.5)

  for (i in seq_len(n)) {
    # 交替阴影
    if (i %% 2 == 1) {
      rect(-0.2, y_mid[i] - spacing / 2, 1.2, y_mid[i] + spacing / 2,
           col = shade_col, border = NA)
    }
    # 色块
    rect(0.02, y_mid[i] - 0.85, 0.12, y_mid[i] + 0.85,
         col = gene_rect_cols[i], border = NA)
    # 基因名（斜体）
    text(0.18, y_mid[i], genes[i], adj = c(0, 0.5), cex = 1.15, font = 3)
  }

  # 分隔线
  for (i in seq_len(n - 1)) {
    sep_y <- (y_mid[i] + y_mid[i + 1]) / 2
    segments(-0.2, sep_y, 1.2, sep_y, col = "#C0C0C0", lwd = 0.6)
  }

  # ========== Panel 2: 森林图 ==========
  par(mar = c(3, 0, 2.5, 0))
  plot.new()
  plot.window(xlim = xlim_forest, ylim = ylim_all)

  # 交替阴影（与 Panel 1 对齐）
  for (i in seq_len(n)) {
    if (i %% 2 == 1) {
      rect(xlim_forest[1] - 0.5, y_mid[i] - spacing / 2,
           xlim_forest[2] + 0.5, y_mid[i] + spacing / 2,
           col = shade_col, border = NA)
    }
  }

  # 分隔线
  for (i in seq_len(n - 1)) {
    sep_y <- (y_mid[i] + y_mid[i + 1]) / 2
    segments(xlim_forest[1], sep_y, xlim_forest[2], sep_y, col = "#C0C0C0", lwd = 0.6)
  }

  # 参考线 HR = 1
  abline(v = 1, lty = 2, col = "grey40", lwd = 1.5)

  # 绘制 CI 线段和菱形点
  for (i in seq_len(n)) {
    # hi_hi (橙)
    segments(merged$CI_lower_hihi[i], y_hihi[i],
             merged$CI_upper_hihi[i], y_hihi[i],
             col = col_hihi, lwd = 2.5)
    points(merged$HR_hihi[i], y_hihi[i], pch = 18, cex = 2.0, col = col_hihi)
    # CI 端帽
    ci_cap <- 0.2
    segments(merged$CI_lower_hihi[i], y_hihi[i] - ci_cap,
             merged$CI_lower_hihi[i], y_hihi[i] + ci_cap, col = col_hihi, lwd = 1.5)
    segments(merged$CI_upper_hihi[i], y_hihi[i] - ci_cap,
             merged$CI_upper_hihi[i], y_hihi[i] + ci_cap, col = col_hihi, lwd = 1.5)
    # 显著性星号
    s <- sig_star(merged$p_value_hihi[i])
    if (nchar(s) > 0) {
      text(merged$CI_upper_hihi[i] + 0.02, y_hihi[i] + 0.15, s,
           col = col_hihi, cex = 1.1, adj = c(0, 0.5), font = 2)
    }

    # lo_lo (蓝)
    segments(merged$CI_lower_lolo[i], y_lolo[i],
             merged$CI_upper_lolo[i], y_lolo[i],
             col = col_lolo, lwd = 2.5)
    points(merged$HR_lolo[i], y_lolo[i], pch = 18, cex = 2.0, col = col_lolo)
    # CI 端帽
    segments(merged$CI_lower_lolo[i], y_lolo[i] - ci_cap,
             merged$CI_lower_lolo[i], y_lolo[i] + ci_cap, col = col_lolo, lwd = 1.5)
    segments(merged$CI_upper_lolo[i], y_lolo[i] - ci_cap,
             merged$CI_upper_lolo[i], y_lolo[i] + ci_cap, col = col_lolo, lwd = 1.5)
    # 显著性星号
    s <- sig_star(merged$p_value_lolo[i])
    if (nchar(s) > 0) {
      text(merged$CI_upper_lolo[i] + 0.02, y_lolo[i] + 0.15, s,
           col = col_lolo, cex = 1.1, adj = c(0, 0.5), font = 2)
    }
  }

  # x 轴
  axis(1, cex.axis = 1.0, col = "grey30", col.axis = "grey20", lwd = 1.5)
  mtext("Hazard ratio", side = 1, line = 2, cex = 1.05, font = 2)

  # 标题
  text(mean(xlim_forest), ylim_all[2] + 0.6, "Forest Plot",
       font = 2, cex = 1.3, xpd = TRUE, adj = 0.5)

  # 图例
  legend("topright",
         legend = c(expression(italic("MYC/PVT1")^"High" ~ "expression group"),
                    expression(italic("MYC/PVT1")^"Low" ~ "expression group")),
         col = c(col_hihi, col_lolo),
         pch = 18, lty = 1, lwd = 2, pt.cex = 1.5,
         cex = 1.56, bty = "n", seg.len = 2.5,
         inset = c(0.01, -0.02), xpd = TRUE)

  # ========== Panel 3: P 值 ==========
  par(mar = c(3, 0, 2.5, 0))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = ylim_all)

  # 交替阴影
  for (i in seq_len(n)) {
    if (i %% 2 == 1) {
      rect(-0.5, y_mid[i] - spacing / 2, 1.5, y_mid[i] + spacing / 2,
           col = shade_col, border = NA)
    }
  }

  # 分隔线
  for (i in seq_len(n - 1)) {
    sep_y <- (y_mid[i] + y_mid[i + 1]) / 2
    segments(-0.5, sep_y, 1.5, sep_y, col = "#C0C0C0", lwd = 0.6)
  }

  # 标题
  text(0.5, ylim_all[2] + 0.6, "P value", font = 2, cex = 1.3, xpd = TRUE, adj = 0.5)

  # P 值文本（各组颜色）
  fmt_p <- function(p) {
    if (p < 0.001) return("< 0.001")
    return(sprintf("%.3f", p))
  }
  for (i in seq_len(n)) {
    text(0.5, y_hihi[i], fmt_p(merged$p_value_hihi[i]),
         col = col_hihi, cex = 1.90, adj = c(0.5, 0.5),
         font = ifelse(merged$p_value_hihi[i] < 0.05, 2, 1))
    text(0.5, y_lolo[i], fmt_p(merged$p_value_lolo[i]),
         col = col_lolo, cex = 1.90, adj = c(0.5, 0.5),
         font = ifelse(merged$p_value_lolo[i] < 0.05, 2, 1))
  }

  # ========== Panel 4: HRR 柱状图 ==========
  par(mar = c(3, 0, 2.5, 1.5))
  hrr_max <- max(merged$VI) * 1.25
  plot.new()
  plot.window(xlim = c(0, hrr_max), ylim = ylim_all)

  # 交替阴影
  for (i in seq_len(n)) {
    if (i %% 2 == 1) {
      rect(-0.5, y_mid[i] - spacing / 2, hrr_max + 0.5, y_mid[i] + spacing / 2,
           col = shade_col, border = NA)
    }
  }

  # 分隔线
  for (i in seq_len(n - 1)) {
    sep_y <- (y_mid[i] + y_mid[i + 1]) / 2
    segments(0, sep_y, hrr_max, sep_y, col = "#C0C0C0", lwd = 0.6)
  }

  # 标题
  text(hrr_max / 2, ylim_all[2] + 0.6, "VI (exp(\u03b2_int))", font = 2, cex = 1.3, xpd = TRUE, adj = 0.5)

  # 水平柱状图
  bar_h <- 1.2
  for (i in seq_len(n)) {
    rect(0, y_mid[i] - bar_h / 2, merged$VI[i], y_mid[i] + bar_h / 2,
         col = hrr_bar_cols[i], border = NA)
    text(merged$VI[i] + 0.01, y_mid[i], sprintf("%.2f", merged$VI[i]),
         adj = c(0, 0.5), cex = 1.8, col = "grey20", xpd = TRUE)
  }

  # x 轴
  axis(1, cex.axis = 1.0, col = "grey30", col.axis = "grey20", lwd = 1.5)

  # ---------- 全图顶线和底线 ----------
  # (通过在各 panel 添加顶/底边框实现)

}

# ---- 输出 PDF ----
pdf(out_pdf, width = 13, height = max(5.5, n * 0.75))
par(oma = c(0, 0, 0.5, 0))  # 外边距留出顶线空间
draw_forest()
dev.off()
cat("PDF saved:", out_pdf, "\n")

# ---- 输出 PNG ----
png(out_png, width = 3900, height = max(1650, n * 225), res = 300)
par(oma = c(0, 0, 0.5, 0))
draw_forest()
dev.off()
cat("PNG saved:", out_png, "\n")
cat("Done!\n")
