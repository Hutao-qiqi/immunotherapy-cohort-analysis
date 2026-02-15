import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys

def analyze_distribution(data, gene_name, group_name):
    """Analyzes and prints the distribution statistics for a given gene and group."""
    print(f"\n--- Analyzing '{gene_name}' expression in group: '{group_name}' ---")

    group_data = data[data['MYC_PVT1_Status'] == group_name][gene_name].dropna()

    if group_data.empty:
        print("No data available for this group.")
        return None

    # Calculate statistics
    median_val = group_data.median()
    mean_val = group_data.mean()
    min_val = group_data.min()
    max_val = group_data.max()
    total_count = len(group_data)

    # Key analysis: count values relative to the median
    count_less = (group_data < median_val).sum()
    count_equal = (group_data == median_val).sum()
    count_greater = (group_data > median_val).sum()
    
    # The 'high' group in the original script was >= median
    high_group_count = count_equal + count_greater
    low_group_count = count_less

    print(f"Total non-null samples: {total_count}")
    print(f"Statistics:")
    print(f"  - Min:    {min_val:.4f}")
    print(f"  - Mean:   {mean_val:.4f}")
    print(f"  - Median: {median_val:.4f}")
    print(f"  - Max:    {max_val:.4f}")
    print("\nVerification of group split (using median):")
    print(f"  - Value of Median: {median_val:.4f}")
    print(f"  - Samples with value < median: {count_less} (This would be the 'low' group)")
    print(f"  - Samples with value = median: {count_equal}")
    print(f"  - Samples with value > median: {count_greater}")
    print(f"  - Resulting 'high' group size (value >= median): {high_group_count}")
    print(f"  - Resulting 'low' group size (value < median):  {low_group_count}")
    
    if count_equal > 1:
        print("\n*** Hypothesis Confirmed ***")
        print(f"There are {count_equal} samples with expression exactly at the median value.")
        print("These are all assigned to the 'high' group by the '>=' logic, causing the imbalance.")

    return group_data

def main():
    # --- Configuration ---
    expression_file = 'combined_expression_combat_corrected.txt'
    annotation_file = 'MYC_PVT1_annotation.txt'
    
    # Use 'IL5' by default, or take a gene name from the command line
    if len(sys.argv) > 1:
        gene_to_check = sys.argv[1]
    else:
        gene_to_check = 'IL5'
    
    print(f"Attempting to analyze gene: '{gene_to_check}'")

    # --- Load Data ---
    print("Loading data...")
    try:
        df_expr_raw = pd.read_csv(expression_file, sep='\s+', index_col=0, engine='python')
        df_expr = df_expr_raw.T
        df_annot = pd.read_csv(annotation_file, sep='\s+', engine='python')
    except FileNotFoundError as e:
        print(f"Error loading file: {e}")
        return

    # Check if gene exists
    if gene_to_check not in df_expr.columns:
        print(f"\nError: Gene '{gene_to_check}' not found in the expression file '{expression_file}'.")
        print("Please check the gene name. Available genes include:", df_expr.columns.tolist()[:10], "...")
        return

    # --- Merge and Prepare Data ---
    df_annot.rename(columns={'Sample': 'Sample_ID'}, inplace=True)
    merged_data = pd.merge(df_annot, df_expr[[gene_to_check]], left_on='Sample_ID', right_index=True, how='inner')
    merged_data[gene_to_check] = pd.to_numeric(merged_data[gene_to_check], errors='coerce')

    # --- Analyze and Plot ---
    hi_hi_data = analyze_distribution(merged_data, gene_to_check, 'hi_hi')
    lo_lo_data = analyze_distribution(merged_data, gene_to_check, 'lo_lo')

    # Plotting
    print("\nGenerating distribution plots...")
    fig, axes = plt.subplots(1, 2, figsize=(18, 7), sharey=True)
    fig.suptitle(f"Distribution of '{gene_to_check}' Expression", fontsize=18)
    
    sns.set_style("whitegrid")

    # Plot for hi_hi group
    if hi_hi_data is not None and not hi_hi_data.empty:
        sns.histplot(hi_hi_data, ax=axes[0], kde=True, bins=30, color='skyblue')
        median_val = hi_hi_data.median()
        axes[0].axvline(median_val, color='r', linestyle='--', lw=2, label=f'Median: {median_val:.4f}')
        axes[0].set_title('hi_hi Group', fontsize=14)
        axes[0].legend()
    else:
        axes[0].text(0.5, 0.5, "No data for hi_hi group", ha='center', va='center')


    # Plot for lo_lo group
    if lo_lo_data is not None and not lo_lo_data.empty:
        sns.histplot(lo_lo_data, ax=axes[1], kde=True, bins=30, color='salmon')
        median_val = lo_lo_data.median()
        axes[1].axvline(median_val, color='r', linestyle='--', lw=2, label=f'Median: {median_val:.4f}')
        axes[1].set_title('lo_lo Group', fontsize=14)
        axes[1].legend()
    else:
        axes[1].text(0.5, 0.5, "No data for lo_lo group", ha='center', va='center')


    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plot_filename = f'distribution_of_{gene_to_check}.png'
    plt.savefig(plot_filename)
    print(f"Plot saved to '{plot_filename}'")
    
    # Try to show the plot, but don't fail if there's no GUI
    try:
        plt.show()
    except Exception as e:
        print(f"\nCould not display plot interactively (e.g., on a server without a GUI). Image is saved to file.")
        print(f"Error details: {e}")


if __name__ == '__main__':
    main() 