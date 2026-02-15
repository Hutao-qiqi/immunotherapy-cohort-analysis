# --- Load required libraries ---
# If a library is not installed, the script will prompt for installation.
if (!requireNamespace("tidyverse", quietly = TRUE)) install.packages("tidyverse")
if (!requireNamespace("survival", quietly = TRUE)) install.packages("survival")
if (!requireNamespace("survminer", quietly = TRUE)) install.packages("survminer")
if (!requireNamespace("tictoc", quietly = TRUE)) install.packages("tictoc")
if (!requireNamespace("progress", quietly = TRUE)) install.packages("progress")
if (!requireNamespace("broom", quietly = TRUE)) install.packages("broom")

library(tidyverse)
library(survival)
library(survminer)
library(tictoc)
library(progress)
library(broom)


#' Calculate p-value for a specific gene within a specific group
#'
#' This function calculates the log-rank p-value for high vs. low gene expression
#' within a defined subgroup (e.g., 'hi_hi' or 'lo_lo').
#'
#' @param data The main dataframe containing survival and expression data.
#' @param gene_name The name of the gene column to analyze.
#' @param group_name The name of the subgroup in 'MYC_PVT1_Status'.
#' @return A list containing the p-value and sample sizes for high/low groups.
calculate_p_value_for_group <- function(data, gene_name, group_name) {
  # Filter data for the specific group
  group_data <- data %>%
    filter(MYC_PVT1_Status == group_name)
  
  # Ensure group has data and the gene column exists and is not all NA
  if (nrow(group_data) == 0 || !gene_name %in% names(group_data) || all(is.na(group_data[[gene_name]]))) {
    return(list(p_value = NA_real_, n_high = 0, n_low = 0))
  }
  
  # Ensure valid survival data
  group_data <- group_data %>%
    drop_na(OS_months, OS_event)
  
  if (nrow(group_data) == 0) {
    return(list(p_value = NA_real_, n_high = 0, n_low = 0))
  }
  
  # Calculate median expression for this specific group
  median_expression <- median(group_data[[gene_name]], na.rm = TRUE)
  
  # If median is NA (e.g., all values were NA), can't proceed
  if (is.na(median_expression)) {
    return(list(p_value = NA_real_, n_high = 0, n_low = 0))
  }
  
  # Create expression group based on the median
  group_data <- group_data %>%
    mutate(Expression_Group = if_else(.data[[gene_name]] >= median_expression, "high", "low"))
  
  # Check if both high and low groups exist and have events
  # `survdiff` requires at least two groups.
  if (n_distinct(group_data$Expression_Group) < 2) {
    return(list(p_value = NA_real_,
                n_high = sum(group_data$Expression_Group == "high"),
                n_low = sum(group_data$Expression_Group == "low")))
  }
  
  # Perform log-rank test
  # The `tryCatch` block handles cases where the test fails for mathematical reasons
  # (e.g., no events in one group, which survdiff can sometimes handle, but can be problematic)
  logrank_result <- tryCatch({
    survdiff(Surv(OS_months, OS_event) ~ Expression_Group, data = group_data)
  }, error = function(e) {
    NULL
  })
  
  if (is.null(logrank_result)) {
    return(list(p_value = NA_real_,
                n_high = sum(group_data$Expression_Group == "high"),
                n_low = sum(group_data$Expression_Group == "low")))
  }
  
  # Extract p-value using the broom package for a clean result
  p_value <- tidy(logrank_result)$p.value
  
  return(list(
    p_value = p_value,
    n_high = sum(group_data$Expression_Group == "high"),
    n_low = sum(group_data$Expression_Group == "low")
  ))
}


#' Plot combined survival curves for a gene
#'
#' Plots four survival curves on a single plot:
#' hi_hi (high/low expr), lo_lo (high/low expr)
#'
#' @param data The main dataframe.
#' @param gene_name The gene to plot.
#' @param p_hi_hi The pre-calculated p-value for the 'hi_hi' group.
#' @param p_lo_lo The pre-calculated p-value for the 'lo_lo' group.
#' @return A ggplot object.
plot_combined_survival_curves <- function(data, gene_name, p_hi_hi, p_lo_lo) {
  # Prepare data for plotting
  plot_data <- data %>%
    filter(MYC_PVT1_Status %in% c("hi_hi", "lo_lo")) %>%
    # CRITICAL: Calculate median and expression group *within* each MYC_PVT1_Status
    group_by(MYC_PVT1_Status) %>%
    mutate(
      median_expr = median(.data[[gene_name]], na.rm = TRUE),
      Expression_Group = if_else(.data[[gene_name]] >= median_expr, "high", "low")
    ) %>%
    ungroup() %>%
    # Create a combined group for plotting aesthetics
    mutate(
      Plot_Group = str_c(MYC_PVT1_Status, " & ", gene_name, "_", Expression_Group)
    )
  
  # Fit a single survival model with the combined group
  fit <- survfit(Surv(OS_months, OS_event) ~ Plot_Group, data = plot_data)
  
  # Create the plot using survminer
  p <- ggsurvplot(
    fit,
    data = plot_data,
    conf.int = FALSE,
    pval = FALSE, # We will add custom p-value annotations
    palette = c("red", "blue", "darkorange", "deepskyblue"), # Corresponds to hi_hi_high, hi_hi_low, lo_lo_high, lo_lo_low
    linetype = c("solid", "solid", "dashed", "dashed"),
    legend = "top.right",
    legend.title = "Groups",
    ggtheme = theme_survminer() # Or theme_bw()
  )
  
  # Add custom annotations and labels
  p_text <- sprintf("p-value (hi_hi) = %.4f\np-value (lo_lo) = %.4f", p_hi_hi, p_lo_lo)
  
  p$plot <- p$plot +
    labs(
      title = paste("Survival Analysis for", gene_name),
      x = "Overall survival time (Months)",
      y = "Probability of Survival"
    ) +
    annotate(
      "text", x = max(plot_data$OS_months, na.rm=T) * 0.1, y = 0.15,
      label = p_text, hjust = 0,
      size = 4,
      fontface = "bold"
    ) +
    coord_cartesian(ylim = c(0, 1.05))
  
  return(p)
}


# --- Main execution block ---
main <- function() {
  # --- System information output ---
  cat("--- System Information ---\n")
  cat(sprintf("Detected CPU logical cores (threads): %d\n", parallel::detectCores()))
  cat("Mode: Single-threaded sequential execution\n")
  cat("--------------------------\n\n")
  
  # --- File path configuration ---
  expression_file <- 'combined_expression_combat_corrected.txt'
  survival_file <- 'updated_survival_data.txt'
  annotation_file <- 'MYC_PVT1_annotation.txt'
  output_dir <- 'survival_plots_filtered'
  results_file <- 'filtered_genes_results.tsv' # Use .tsv for tab-separated
  
  # --- Timer initialization ---
  tictoc::tic("Total Runtime")
  tictoc::tic("Step 1")
  
  # --- Step 1: Load and preprocess data ---
  cat("Step 1/5: Starting data loading and preprocessing...\n")
  dir.create(output_dir, showWarnings = FALSE)
  
  tryCatch({
    cat("  - Loading expression file...\n")
    # read_tsv is fast and robust for tab or space separated files
    df_expr_raw <- read_tsv(expression_file, col_types = cols(.default = "d", Hugo_Symbol = "c")) %>%
      # Ensure unique gene names before setting as rownames for transpose
      distinct(Hugo_Symbol, .keep_all = TRUE) %>% 
      column_to_rownames(var = "Hugo_Symbol")
    
    cat(sprintf("    > Successfully loaded %d genes, %d samples.\n", nrow(df_expr_raw), ncol(df_expr_raw)))
    
    # Transpose data so samples are rows and genes are columns
    df_expr <- as.data.frame(t(df_expr_raw)) %>%
      rownames_to_column(var = "Sample_ID")
    
    cat("  - Loading survival file...\n")
    df_surv <- read_tsv(survival_file)
    
    cat("  - Loading annotation file...\n")
    df_annot <- read_tsv(annotation_file) %>%
      rename(Sample_ID = Sample)
    
    # --- Merge data ---
    cat("  - Merging data tables...\n")
    full_data <- df_surv %>%
      inner_join(df_annot, by = "Sample_ID") %>%
      inner_join(df_expr, by = "Sample_ID")
    
    cat(sprintf("    > Data merge complete. %d samples for analysis.\n", nrow(full_data)))
    
    # --- Clean merged data ---
    base_data <- full_data %>%
      mutate(
        OS_months = as.numeric(OS_months),
        OS_event = as.numeric(OS_event)
      ) %>%
      drop_na(OS_months, OS_event, MYC_PVT1_Status)
    
  }, error = function(e) {
    stop(paste("Error during data loading:", e$message))
  })
  
  cat(sprintf("Step 1/5: Data loading and preprocessing complete. "), "Time elapsed: ", unlist(tictoc::toc(quiet=T)), " seconds\n\n")
  
  # Steps 2 & 3 combined: Gene screening
  tictoc::tic("Step 2/3")
  all_genes <- names(df_expr_raw)
  cat(sprintf("Step 2-3/5: Starting single-threaded scan of %d genes...\n", length(all_genes)))
  
  # Initialize progress bar
  pb <- progress_bar$new(
    format = "Gene screening [:bar] :percent in :elapsed | ETA: :eta",
    total = length(all_genes), clear = FALSE, width = 80
  )
  
  # Loop through all genes to find those meeting the criteria
  filtered_gene_results <- list()
  for (gene in all_genes) {
    pb$tick() # Update progress bar
    
    # Skip if gene column doesn't exist for some reason
    if (!gene %in% names(base_data)) next
    
    res_hi_hi <- calculate_p_value_for_group(base_data, gene, "hi_hi")
    res_lo_lo <- calculate_p_value_for_group(base_data, gene, "lo_lo")
    
    # Check if the conditions are met (p-values must not be NA)
    if (!is.na(res_hi_hi$p_value) && !is.na(res_lo_lo$p_value)) {
      if (res_hi_hi$p_value < 0.05 && res_lo_lo$p_value > 0.05) {
        # If conditions met, store the result
        filtered_gene_results[[gene]] <- tibble(
          Gene = gene,
          p_value_hi_hi = res_hi_hi$p_value,
          p_value_lo_lo = res_lo_lo$p_value,
          hi_hi_high_n = res_hi_hi$n_high,
          hi_hi_low_n = res_hi_hi$n_low,
          lo_lo_high_n = res_lo_lo$n_high,
          lo_lo_low_n = res_lo_lo$n_low
        )
      }
    }
  }
  
  if (length(filtered_gene_results) == 0) {
    cat("\nStep 2-3/5: Scan complete. No genes met the specified criteria. Exiting.\n")
    return()
  }
  
  cat(sprintf("\nStep 2-3/5: Scan complete. "), "Time elapsed: ", unlist(tictoc::toc(quiet=T)), " seconds\n\n")
  
  # --- Step 4: Save results to a file ---
  tictoc::tic("Step 4")
  cat("Step 4/5: Saving filtered results...\n")
  
  # Combine list of tibbles into a single dataframe
  results_df <- bind_rows(filtered_gene_results) %>%
    arrange(p_value_hi_hi) # Sort by p-value
  
  write_tsv(results_df, results_file)
  
  cat(sprintf("  - Found %d matching genes. Results saved to '%s'.\n", nrow(results_df), results_file))
  cat(sprintf("Step 4/5: Results saved. "), "Time elapsed: ", unlist(tictoc::toc(quiet=T)), " seconds\n\n")
  
  # --- Step 5: Generate plots for filtered genes ---
  tictoc::tic("Step 5")
  cat(sprintf("Step 5/5: Starting to plot for %d genes...\n", nrow(results_df)))
  
  # Initialize progress bar for plotting
  pb_plot <- progress_bar$new(
    format = "Plotting [:bar] :percent in :elapsed | ETA: :eta",
    total = nrow(results_df), clear = FALSE, width = 80
  )
  
  # Use pwalk for a clean loop over the results dataframe rows
  pwalk(results_df, function(Gene, p_value_hi_hi, p_value_lo_lo, ...) {
    pb_plot$tick() # Update progress bar
    
    # Generate the plot
    survival_plot <- plot_combined_survival_curves(base_data, Gene, p_value_hi_hi, p_value_lo_lo)
    
    # Save the plot
    output_path <- file.path(output_dir, paste0("survival_curve_", Gene, ".png"))
    ggsave(output_path, plot = survival_plot, width = 10, height = 7, dpi = 300)
  })
  
  cat(sprintf("\nStep 5/5: Plotting complete. "), "Time elapsed: ", unlist(tictoc::toc(quiet=T)), " seconds\n")
  cat(sprintf("\nAll images saved to '%s' directory.\n", output_dir))
  
  cat("--- ")
  tictoc::toc()
  cat(" ---\n")
}

# --- Run the main function ---
main()