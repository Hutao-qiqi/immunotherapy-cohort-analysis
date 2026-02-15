# Title: 3D PCA Plot for Gene Expression Data
# Description: This script performs PCA on gene expression data and
#              visualizes the first three principal components in a 3D scatter plot.
#              The points are colored by batch, following the style of the provided tutorial.

# --- 1. Setup ---

# Set working directory to script location, so relative paths resolve.
cmd <- commandArgs(trailingOnly = FALSE)
file_arg <- cmd[grepl("^--file=", cmd)]
if (length(file_arg) > 0) {
  script_path <- sub("^--file=", "", file_arg[1])
  script_dir <- normalizePath(dirname(script_path), winslash = "/", mustWork = FALSE)
  setwd(script_dir)
}

# Load dependency (no auto-install in repo scripts).
if (!requireNamespace("plot3D", quietly = TRUE)) {
  stop("Missing R package: plot3D. Install it first with install.packages('plot3D')")
}
suppressPackageStartupMessages(library(plot3D))

cat("Successfully loaded plot3D package.\n")

# --- 2. Load Data ---

# Define file paths
expression_file_after <- "combined_expression_combat_corrected.txt"
expression_file_before <- "combined_expression_before_combat.txt"
batch_file <- "batch_info.txt"

# Check if files exist
if (!file.exists(expression_file_after)) {
  stop("Error: Corrected expression data file not found: ", expression_file_after)
}
if (!file.exists(expression_file_before)) {
  stop("Error: Original expression data file not found: ", expression_file_before)
}
if (!file.exists(batch_file)) {
  stop("Error: Batch info file not found: ", batch_file)
}

# Load the batch-corrected expression data
cat("Loading batch-corrected expression data from '", expression_file_after, "'...\n", sep="")
expr_data_after <- read.delim(expression_file_after, row.names = 1, check.names = FALSE)

# Load the original (before correction) expression data
cat("Loading original expression data from '", expression_file_before, "'...\n", sep="")
expr_data_before <- read.delim(expression_file_before, row.names = 1, check.names = FALSE)

# Load the batch information for each sample
cat("Loading batch information from '", batch_file, "'...\n", sep="")
batch_info <- read.delim(batch_file, row.names = 1, check.names = FALSE)

# --- 3. Data Preparation ---

# Ensure that the samples in the expression data and batch info are in the same order.
cat("Aligning samples between expression data and batch info...\n")
common_samples <- intersect(colnames(expr_data_after), rownames(batch_info))
if (length(common_samples) < 3) {
    stop("Error: Not enough common samples found.")
}
expr_data_after <- expr_data_after[, common_samples, drop = FALSE]
expr_data_before <- expr_data_before[, common_samples, drop = FALSE]
batch_info_aligned <- batch_info[common_samples, , drop = FALSE]
batch_labels <- factor(batch_info_aligned$Batch)

# Transpose the expression matrices for PCA
data_for_pca_after <- t(expr_data_after)
data_for_pca_before <- t(expr_data_before)
cat(sprintf("Prepared data for PCA with %d samples and %d genes.\n", nrow(data_for_pca_after), ncol(data_for_pca_after)))

# --- 4. Principal Component Analysis (PCA) ---

# Perform PCA on both datasets
cat("Running PCA on 'before' and 'after' data...\n")
pca_res_before <- prcomp(data_for_pca_before, center = TRUE, scale. = TRUE)
pca_res_after <- prcomp(data_for_pca_after, center = TRUE, scale. = TRUE)

# Extract variance explained for both
variance_explained_before <- summary(pca_res_before)$importance[2, 1:3] * 100
variance_explained_after <- summary(pca_res_after)$importance[2, 1:3] * 100

# Create data frames with the PCA scores
pca_scores_before <- as.data.frame(pca_res_before$x)
pca_scores_after <- as.data.frame(pca_res_after$x)

# --- 5. Plotting Preparation ---

# Get unique batch names and assign colors.
unique_batches <- levels(batch_labels)
num_batches <- length(unique_batches)

# Define a color palette similar to the Python script (Nature journal style).
nature_colors <- c('#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f')
if (num_batches > length(nature_colors)) {
    # If more batches than colors, use a different palette
    batch_colors <- rainbow(num_batches)
} else {
    batch_colors <- nature_colors[1:num_batches]
}
names(batch_colors) <- unique_batches

# Define prettier names for the legend, as in the Python script.
batch_name_mapping <- c(
  'IMvigor210' = 'Mariathasan et al.',
  'liu' = 'Liu et al.',
  'ra' = 'Riaz et al.',
  'hugo' = 'Hugo et al.'
)
legend_names <- batch_name_mapping[unique_batches]
# Fallback for any names not in the mapping
legend_names[is.na(legend_names)] <- unique_batches[is.na(legend_names)]


# --- 6. Generate 3D PCA Plots ---

output_file <- "3d_pca_comparison_plot.pdf"
cat("Generating 3D PCA comparison plot and saving to '", output_file, "'...\n", sep="")

# Open a standard PDF device
pdf(output_file, width = 16, height = 8)

# Set up the layout to have two plots side-by-side
par(mfrow = c(1, 2), mar = c(5.1, 4.1, 4.1, 2.1))

# --- Plot 1: Before ComBat ---
scatter3D(
  x = pca_scores_before$PC1, y = pca_scores_before$PC2, z = pca_scores_before$PC3,
  col = "black", pch = 19, cex = 1.3,
  ticktype = "detailed", bty = "b2", theta = 40, phi = 20, d = 3,
  main = "Before ComBat Correction",
  xlab = sprintf("PC1 (%.1f%% variance)", variance_explained_before[1]),
  ylab = sprintf("PC2 (%.1f%% variance)", variance_explained_before[2]),
  zlab = sprintf("PC3 (%.1f%% variance)", variance_explained_before[3]),
  font.main = 2, # 2 = bold
  font.lab = 2,  # 2 = bold
  colkey = FALSE
)
scatter3D(
  x = pca_scores_before$PC1, y = pca_scores_before$PC2, z = pca_scores_before$PC3,
  colvar = as.numeric(batch_labels), col = batch_colors,
  pch = 19, cex = 1.2, add = TRUE, colkey = FALSE
)

# --- Plot 2: After ComBat ---
scatter3D(
  x = pca_scores_after$PC1, y = pca_scores_after$PC2, z = pca_scores_after$PC3,
  col = "black", pch = 19, cex = 1.3,
  ticktype = "detailed", bty = "b2", theta = 40, phi = 20, d = 3,
  main = "After ComBat Correction",
  xlab = sprintf("PC1 (%.1f%% variance)", variance_explained_after[1]),
  ylab = sprintf("PC2 (%.1f%% variance)", variance_explained_after[2]),
  zlab = sprintf("PC3 (%.1f%% variance)", variance_explained_after[3]),
  font.main = 2, # 2 = bold
  font.lab = 2,  # 2 = bold
  colkey = FALSE
)
scatter3D(
  x = pca_scores_after$PC1, y = pca_scores_after$PC2, z = pca_scores_after$PC3,
  colvar = as.numeric(batch_labels), col = batch_colors,
  pch = 19, cex = 1.2, add = TRUE, colkey = FALSE
)

# Add a shared legend in the outer margin
par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
legend(
  "right",
  legend = legend_names,
  col = batch_colors,
  pch = 19,
  cex = 1.2,
  title = "Cohort",
  pt.cex = 2,
  bty = "n",
  text.font = 2
)


# Close the PDF device.
dev.off()

cat("Done. Plot saved to '", output_file, "'.\n", sep="")
