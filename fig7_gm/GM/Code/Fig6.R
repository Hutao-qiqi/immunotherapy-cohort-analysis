# 清空环境
rm(list = ls())
gc()

# Set working directory to GM root (script_dir/..), avoid absolute paths.
cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- cmd[grepl("^--file=", cmd)]
if (length(file_arg) > 0) {
  script_path <- sub("^--file=", "", file_arg[1])
  script_dir <- normalizePath(dirname(script_path), winslash = "/", mustWork = FALSE)
  gm_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = FALSE)
  setwd(gm_root)
}

####== 1. 加载库 ==####
library(dplyr); library(ggsci); library(ComplexHeatmap); library(circlize); library(stringr)

####== 2. 加载数据并预处理 ==####
load('data/CRISPR/crispr.Rdata')

####== 3. 数据处理 ==####
df <- crispr
df$mean <- apply(df, 1, function(z) { mean(z, na.rm = TRUE) })
df_filtered <- df[!is.na(df$mean), ]
df_sorted <- df_filtered[order(df_filtered$mean), ]
df_subset <- df_sorted[c(1:8, (nrow(df_sorted) - 7):nrow(df_sorted)), ]
plot_matrix <- as.matrix(df_subset[, -ncol(df_subset)])

####== 4. 创建顶部注释 ==####
cohort_names <- colnames(plot_matrix)
anno_df <- data.frame(row.names = cohort_names)
anno_df$cancer <- case_when(
  str_detect(cohort_names, "Renal") ~ "Renal carcinoma",
  str_detect(cohort_names, "Melanoma") ~ "Melanoma",
  str_detect(cohort_names, "Colon") ~ "Colon adenocarcinoma",
  str_detect(cohort_names, "Breast") ~ "Breast carcinoma",
  TRUE ~ "Unknown"
)
anno_df$cohort <- cohort_names

# --- 定义颜色 ---
cancer_colors <- c("Renal carcinoma" = "black", "Melanoma" = "#666666", "Colon adenocarcinoma" = "#B59D77", "Breast carcinoma" = "#E69F00")
unique_cohorts <- unique(anno_df$cohort)
cohort_colors <- setNames(rand_color(length(unique_cohorts)), unique_cohorts)

# --- 创建注释对象 (已修正图例边框问题) ---
ha_top <- HeatmapAnnotation(
  cancer = anno_df$cancer,
  cohort = anno_df$cohort,
  col = list(cancer = cancer_colors, cohort = cohort_colors),
  show_legend = c(cancer = TRUE, cohort = TRUE), 
  annotation_legend_param = list(
    cancer = list(
      title = "Cancer",
      title_gp = gpar(fontface = "bold"),
      border = "black" # <-- 使用 border 为离散图例添加边框
    ), 
    cohort = list(
      title = "Cohort",
      ncol = 2,
      labels_gp = gpar(fontsize = 8),
      title_gp = gpar(fontface = "bold"),
      border = "black" # <-- 使用 border 为离散图例添加边框
    )
  ),
  annotation_name_gp = gpar(fontsize = 10),
  gap = unit(1, 'mm'),
  gp = gpar(col = "black", lwd = 1) 
)

####== 5. 绘制热图 ==####
main_col_fun <- colorRamp2(c(-2, 0, 2), c("#440154", "#F5F5F5", "#FDE725"))
row_split_vec <- c(rep("Immune resistant", 8), rep("Immune sensitive", 8))

ht <- Heatmap(plot_matrix, 
              name = "z scores", col = main_col_fun, top_annotation = ha_top,
              cluster_rows = FALSE, cluster_columns = FALSE, row_names_side = 'left',
              row_split = row_split_vec, row_title_rot = 90, row_title_gp = gpar(fontsize = 12),
              row_gap = unit(3, "mm"), show_column_names = FALSE, na_col = 'white',
              rect_gp = gpar(col = "black", lwd = 1), 
              width = ncol(plot_matrix) * unit(5, "mm"), 
              height = nrow(plot_matrix) * unit(5, "mm"),
              heatmap_legend_param = list(
                title = "z scores",
                title_gp = gpar(fontface = "bold"),
                grid_border = "black" # <-- 为连续图例添加边框
              )
)

####== 6. 保存为PDF文件 ==####
# --- 步骤 1: 开启PDF设备 ---
# 定义文件名和尺寸（单位：英寸）
# 您可能需要反复调整 width 和 height 来获得最佳布局
pdf("heatmap_output.pdf", width = 12, height = 8)

# --- 步骤 2: 绘制图形 ---
# 这会将图形内容输出到PDF文件而不是R的绘图窗口
draw(ht, 
     heatmap_legend_side = "right", 
     annotation_legend_side = "right",
     legend_gap = unit(10, "mm")
)

# --- 步骤 3: 关闭PDF设备，完成保存 ---
# 这是非常重要的一步，否则文件会损坏或为空
dev.off()

####=== 将图例移至底部的最终代码 ===####

# --- 1. 确保已加载所有必需的库 ---
library(ComplexHeatmap)
library(circlize)

# --- 2. 数据准备 ---
# (这部分代码保持不变)
selected_genes <- c(
  "ADAM10", "MMP1", "ADH1C", "CD44", "MAP4K5", "CD40LG", "SCN7A", 
  "CYP4F8", "TUBB4B", "CYP2R1", "HCRTR2", "CBR1", "FGF16", "CDH3", 
  "CFTR", "AKR1A1", "CA2", "IL3RA", "PGF", "SQLE", "MMP2", "ERBB2"
)
# 假设 crispr 数据框已存在
df_subset <- crispr[selected_genes, ]
plot_matrix <- t(df_subset)
mean_scores <- apply(df_subset, 1, mean, na.rm = TRUE)

# --- 3. 创建顶部注释和主热图对象 ---
# (这部分代码保持不变)
safe_font <- "sans" 
col_fun <- colorRamp2(c(-2, 0, 2), c("#440154", "grey", "#FDE725"))
ha_top <- HeatmapAnnotation(
  `mean Z Scores` = mean_scores,
  col = list(`mean Z Scores` = col_fun),
  gp = gpar(col = "black"),
  show_legend = FALSE,
  annotation_name_side = "right",
  annotation_name_gp = gpar(fontfamily = safe_font, fontface = "bold", fontsize = 10),
  gap = unit(4, 'mm')
)
ht <- Heatmap(plot_matrix,
              name = "Z scores", 
              col = col_fun,
              top_annotation = ha_top,
              cluster_rows = FALSE, 
              cluster_columns = FALSE,
              row_names_side = 'right',
              column_names_side = 'top',
              column_names_rot = 45,
              row_names_gp = gpar(fontfamily = safe_font, fontface = "bold", fontsize = 10),
              column_names_gp = gpar(fontfamily = safe_font, fontface = "bold", fontsize = 10),
              rect_gp = gpar(col = "black", lwd = 1),
              na_col = 'white',
              # 我们在这里配置图例的样式，例如方向
              heatmap_legend_param = list(
                title = "Z scores",
                title_gp = gpar(fontfamily = safe_font, fontface = "bold"),
                labels_gp = gpar(fontfamily = safe_font, fontface = "bold"),
                legend_direction = "horizontal", # 设置为水平方向，适合放在底部
                legend_width = unit(6, "cm"),   # 可以调整图例条的长度
                title_position = "topcenter"
              )
)

# --- 4. 绘制图形，并明确指定图例位置 ---

# 在绘图前，请确保您的RStudio Plots窗格足够高，以便能看到下方的图例
draw(
  ht, 
  # --- 这是最关键的指令 ---
  # 明确告诉 R 将热图的主图例 ("Z scores") 放置在底部
  heatmap_legend_side = "bottom"
)
