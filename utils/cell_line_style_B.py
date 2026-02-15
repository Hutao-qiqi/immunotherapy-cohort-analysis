import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import warnings
from pathlib import Path

# --- 1. Global plotting parameters ---
warnings.filterwarnings('ignore')
plt.rcParams.update({
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'font.family': 'Arial',
    'font.size': 10,
    'axes.linewidth': 0,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.2,
    'pdf.fonttype': 42,
})

# --- 2. Data loading ---
def _find_default_input_file(filename: str = "996celllines.txt") -> Path:
    """Find the default input file by searching parent directories."""
    script_path = Path(__file__).resolve()
    for parent in [script_path.parent] + list(script_path.parents):
        candidate = parent / filename
        if candidate.exists():
            return candidate
    return Path(filename)


def load_data(input_file: Path):
    """Load and clean input data."""
    try:
        df = pd.read_csv(input_file, sep='\t')
        print(f"Loaded input data: {len(df)} rows")
    except FileNotFoundError:
        print(f"ERROR: input file not found: {input_file}")
        print("Place the file next to the script, or keep it in the project root.")
        return None
    
    # Standardize column names
    df.columns = df.columns.str.strip().str.lower()
    
    # Check required columns
    if 'site_primary' not in df.columns or 'tcga_code' not in df.columns:
        print("ERROR: missing required columns 'Site_Primary' or 'tcga_code'.")
        return None

    df['site_primary'] = df['site_primary'].fillna('Unknown')
    df['tcga_code'] = df['tcga_code'].fillna('Unknown')
    
    df_clean = df[(df['site_primary'] != 'Unknown') & (df['tcga_code'] != 'Unknown')].copy()
    print(f"Rows used for plotting after cleaning: {len(df_clean)}")
    
    return df_clean

# --- 3. Style definitions ---
def get_style_definitions():
    """Define color and display name mappings."""
    organ_colors = {
        'lung': '#FDB462', 'haematopoietic_and_lymphoid_tissue': '#FB8072', 'skin': '#E5E5E5',
        'colon': '#8DD3C7', 'breast': '#FFFFB3', 'ovary': '#7FC97F', 'central_nervous_system': '#BEBADA',
        'pancreas': '#BEAED4', 'stomach': '#FCCDE5', 'kidney': '#DEB887', 'upper_aerodigestive_tract': '#E0E0E0',
        'endometrium': '#F0F8FF', 'liver': '#80B1D3', 'oesophagus': '#FBCDE5', 'urinary_tract': '#FFE4E1',
        'bone': '#F5F5DC', 'soft_tissue': '#B3DE69', 'autonomic_ganglia': '#D9D9D9', 'thyroid': '#E6E6FA',
        'pleura': '#DCDCDC', 'large_intestine': '#8DD3C7', 'prostate': '#FFFF99', 'salivary_gland': '#F0F0F0',
        'biliary_tract': '#F8F8FF', 'small_intestine': '#8DD3C7', 'Unknown': '#F5F5F5'
    }
    tcga_colors = {
        'ALL': '#A50F15', 'BLCA': '#FFE4E1', 'BRCA': '#FFFFD4', 'COAD/READ': '#74C476', 'DLBC': '#EF3B2C',
        'ESCA': '#D4B9DA', 'GBM': '#C7E9B4', 'HNSC': '#DDDDDD', 'KIRC': '#CDAA7D', 'LAML': '#FB6A4A',
        'LGG': '#7FCDBB', 'LIHC': '#6BAED6', 'LUAD': '#FD8D3C', 'LUSC': '#FC4E2A', 'MESO': '#CCCCCC',
        'MM': '#CB181D', 'OV': '#2171B5', 'PAAD': '#084594', 'SCLC': '#E31A1C', 'SKCM': '#E5E5E5',
        'STAD': '#F1EEF6', 'UCEC': '#A6BDDB', 'LCML': '#696969', 'PRAD': '#FFEDA0', 'SARC': '#A1DAB4',
        'THCA': '#D8BFD8', 'PCPG': '#D9D9D9', 'ACC': '#F0F0F0', 'COAD': '#74C476', 'READ': '#41AB5D',
        'KIRP': '#D2691E', 'KICH': '#A0522D', 'CHOL': '#4292C6', 'UCS': '#E0E0E0', 'TGCT': '#D8D8D8', 
        'THYM': '#D0D0D0', 'UVM': '#90EE90', 'CHRD': '#B3DE69', 'CLL': '#FFFFFF', 'MB': '#FFFFFF',
        'NB': '#FFFFFF', 'Unknown': '#F5F5F5'
    }
    organ_name_map = {
        'haematopoietic_and_lymphoid_tissue': 'Hematopoietic And Lymphoid Tissue',
        'central_nervous_system': 'Brain', 'large_intestine': 'Colon',
        'upper_aerodigestive_tract': 'Upper Aerodigestive Tract', 'urinary_tract': 'Bladder',
    }
    return organ_colors, tcga_colors, organ_name_map

# --- 4. Main plotting function (legend columns swapped) ---
def create_final_chart(df):
    """Create the chart with the TYPE legend on the left."""
    organ_colors, tcga_colors, organ_name_map = get_style_definitions()

    # --- Data preparation ---
    df['site_primary'] = df['site_primary'].str.lower().str.replace(' ', '_')
    outer_data = df.groupby('site_primary').size().sort_values(ascending=False)
    
    inner_data, inner_labels, inner_colors = [], [], []
    for site in outer_data.index:
        subset = df[df['site_primary'] == site]
        tcga_counts = subset.groupby('tcga_code').size().sort_values(ascending=False)
        for tcga_code, count in tcga_counts.items():
            inner_data.append(count)
            inner_labels.append(tcga_code)
            inner_colors.append(tcga_colors.get(tcga_code, '#F5F5F5'))

    outer_colors_list = [organ_colors.get(organ, '#F5F5F5') for organ in outer_data.index]

    # --- Plot area setup ---
    fig = plt.figure(figsize=(16, 10))
    gs = fig.add_gridspec(1, 2, width_ratios=[2.5, 1.5], wspace=0.1)
    ax_main = fig.add_subplot(gs[0])
    ax_legend = fig.add_subplot(gs[1])
    
    # --- Draw chart ---
    wedgeprops_outer = dict(width=0.3, edgecolor='white', linewidth=1)
    wedgeprops_inner = dict(width=0.25, edgecolor='white', linewidth=0.6)
    
    ax_main.pie(outer_data, radius=1.0, colors=outer_colors_list, wedgeprops=wedgeprops_outer, startangle=90, counterclock=False)
    ax_main.pie(inner_data, radius=0.7, colors=inner_colors, wedgeprops=wedgeprops_inner, startangle=90, counterclock=False)
    
    ax_main.set_aspect('equal')
    ax_main.text(0, 0, 'Cancer Cell\nLine', ha='center', va='center', fontsize=28, fontweight='bold')
    
    # --- Draw legends (TYPE first) ---
    ax_legend.axis('off')

    # Column 1: TYPE (left)
    ax_legend.text(0.0, 1.0, 'TYPE', fontsize=16, fontweight='bold', transform=ax_legend.transAxes)
    y_pos_type = 0.95
    tcga_legend_items = sorted([code for code in df['tcga_code'].unique() if code != 'Unknown'])
    for tcga in tcga_legend_items:
        color = tcga_colors.get(tcga, '#F5F5F5')
        rect = mpatches.Rectangle((0.0, y_pos_type - 0.015), 0.1, 0.03, facecolor=color, edgecolor='black', linewidth=0.5, transform=ax_legend.transAxes)
        ax_legend.add_patch(rect)
        ax_legend.text(0.15, y_pos_type, tcga, fontsize=12, va='center', transform=ax_legend.transAxes)
        y_pos_type -= 0.045
        if y_pos_type < -0.05:
            break  # stop if there are too many items

    # Column 2: ORGAN (right)
    ax_legend.text(0.5, 1.0, 'ORGAN', fontsize=16, fontweight='bold', transform=ax_legend.transAxes)
    y_pos_organ = 0.95
    organ_legend_items = outer_data.index.tolist()
    for organ in organ_legend_items:
        color = organ_colors.get(organ, '#F5F5F5')
        rect = mpatches.Rectangle((0.5, y_pos_organ - 0.015), 0.1, 0.03, facecolor=color, edgecolor='black', linewidth=0.5, transform=ax_legend.transAxes)
        ax_legend.add_patch(rect)
        display_name = organ_name_map.get(organ, organ.replace('_', ' ').title())
        ax_legend.text(0.65, y_pos_organ, display_name, fontsize=12, va='center', transform=ax_legend.transAxes)
        y_pos_organ -= 0.045
        if y_pos_organ < -0.05:
            break  # stop if there are too many items

    plt.tight_layout(pad=1.0, w_pad=1.0)
    return fig

# --- 5. Entry point ---
def main():
    """Run end-to-end generation."""
    print("--- Generating cell line distribution plot (legend columns swapped) ---")

    input_file = _find_default_input_file("996celllines.txt")
    
    df = load_data(input_file)
    if df is None:
        print("--- Aborted due to input load failure ---")
        return
    
    print("Rendering figure...")
    fig = create_final_chart(df)
    
    # Save output
    output_dir = Path(__file__).resolve().parent / "outputs"
    output_dir.mkdir(parents=True, exist_ok=True)
    output_filename = output_dir / "Cell_Line_Distribution_Legend_Swapped.png"
    try:
        fig.savefig(output_filename, dpi=300)
        print(f"Saved: {output_filename}")
    except Exception as e:
        print(f"ERROR: failed to save figure: {e}")
    
    # Display
    plt.show()
    print("--- Done ---")

if __name__ == "__main__":
    main()