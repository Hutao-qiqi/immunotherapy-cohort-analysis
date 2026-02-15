import plotly.graph_objects as go

# --- 1. LABELS (Final hierarchical structure) ---
labels = [
    # --- Stage 0: Cancer Types (Indices 0-2) ---
    "<b>Melanoma</b>",
    "<b>Urothelial Carcinoma</b>",
    "<b>NSCLC</b>",

    # --- Stage 1: Source Cohorts (Indices 3-7) ---
    "Riaz et al. (n=109)",
    "Hugo et al. (n=28)",
    "Liu et al. (n=121)",
    "IMvigor210 (n=348)",
    "Rizvi et al (n=240)",

    # --- Stage 2: Aggregation & Filtering (Indices 8-9) ---
    "<b>Combined Analysis Cohort (n=650)</b>",
    "Excluded (n=196)",

    # --- Stage 3: Final Stratification (Indices 10-11) ---
    "<b>8q24-Amplified (n=157)</b>",
    "8q24-Neutral/Deleted (n=493)"
]

# --- 2. PALETTE (Updated for cancer types) ---
palette = {
    'Hero_Gold':    '#E4CD87', 'Control_Terra':'#C69287',
    'Combined_Terra': '#C69287', 'Excluded_Grey':'#D3D3D3',
    'Melanoma_Family': '#FAE5B8', 'Urothelial_Family': '#E79A90',
    'NSCLC_Family': '#a9c3e2', 'Text_Strong':  '#212121'
}

node_colors = [
    palette['Melanoma_Family'], palette['Urothelial_Family'], palette['NSCLC_Family'],
    palette['Melanoma_Family'], palette['Melanoma_Family'], palette['Melanoma_Family'],
    palette['Urothelial_Family'], palette['NSCLC_Family'],
    palette['Combined_Terra'], palette['Excluded_Grey'],
    palette['Hero_Gold'], palette['Control_Terra']
]

# --- 3. DATA FLOW (Final hierarchical logic) ---
source_nodes = [
    0, 0, 0, 1, 2, # Stage 0 -> 1
    3, 3, 4, 4, 5, 5, 6, 6, 7, 7, # Stage 1 -> 2
    8, 8 # Stage 2 -> 3
]
target_nodes = [
    3, 4, 5, 6, 7, # Stage 0 -> 1
    8, 9, 8, 9, 8, 9, 8, 9, 8, 9, # Stage 1 -> 2
    10, 11 # Stage 2 -> 3
]
values = [
    109, 28, 121, 348, 240, # Stage 0 -> 1
    24, 85, 26, 2, 109, 12, 298, 50, 193, 47, # Stage 1 -> 2
    157, 493 # Stage 2 -> 3
]

# --- 4. LINK COLORS ---
link_colors_rgba = [
    # Stage 0 -> 1
    'rgba(250, 229, 184, 0.7)', 'rgba(250, 229, 184, 0.7)', 'rgba(250, 229, 184, 0.7)',
    'rgba(231, 154, 144, 0.7)', 'rgba(169, 195, 226, 0.7)',
    # Stage 1 -> 2
    'rgba(250, 229, 184, 0.7)', 'rgba(211, 211, 211, 0.7)', # Riaz
    'rgba(250, 229, 184, 0.7)', 'rgba(211, 211, 211, 0.7)', # Hugo
    'rgba(250, 229, 184, 0.7)', 'rgba(211, 211, 211, 0.7)', # Liu
    'rgba(231, 154, 144, 0.7)', 'rgba(211, 211, 211, 0.7)', # IMvigor
    'rgba(169, 195, 226, 0.7)', 'rgba(211, 211, 211, 0.7)', # MSKCC
    # Stage 2 -> 3
    'rgba(228, 205, 135, 0.9)', 'rgba(198, 146, 135, 0.7)'
]

# --- 5. FIGURE CONSTRUCTION WITH REFINED LAYOUT ---
fig = go.Figure(data=[go.Sankey(
    arrangement='snap',
    node=dict(
        pad=15,
        thickness=20,
        line=dict(color=palette['Text_Strong'], width=0.5),
        label=[f"<b>{l}</b>" for l in labels],
        color=node_colors,
        # --- KEY MODIFICATION: Manual X and Y for compact grouping ---
        x=[0.01, 0.01, 0.01,  # Cancer types
           0.25, 0.25, 0.25, 0.25, 0.25, # Cohorts
           0.6, 0.6,      # Aggregation
           1, 1],          # Final
        y=[0.85, 0.4, 0.1,    # Cancer types (spread out)
           0.95, 0.85, 0.75, # Melanoma cohorts (grouped)
           0.4,              # Urothelial cohort
           0.1,              # NSCLC cohort
           0.4, 0.9,         # Aggregation (spread)
           0.3, 0.6]         # Final (spread)
    ),
    link=dict(
        source=source_nodes,
        target=target_nodes,
        value=values,
        color=link_colors_rgba
    )
)])

# --- 6. VISUAL OPTIMIZATION ---
fig.update_layout(
    title=dict(
        text="<b>Figure 1. Patient Flow and Biomarker Stratification of a Pan-Cancer ICB Cohort</b>",
        font=dict(family="Arial, sans-serif", size=22, color=palette['Text_Strong']),
        x=0.5
    ),
    font=dict(
        family="Arial, sans-serif",
        size=14,
        color=palette['Text_Strong']
    ),
    plot_bgcolor='white',
    paper_bgcolor='#FDFBF6',
    margin=dict(l=40, r=40, t=80, b=40)
)

fig.show()