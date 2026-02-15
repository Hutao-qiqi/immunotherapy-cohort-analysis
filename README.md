# Immunotherapy cohort analysis

This repository contains analysis scripts used for a multi-cohort immunotherapy study, including:

- Benefit score validation (IMvigor210-trained logistic model applied to external cohorts)
- Pan-cancer immune analyses and visualizations (bubble heatmap, trend line plot, percent barplot)
- Additional supporting scripts that are part of the full analysis pipeline (organized by function)

## Repository structure

- `benefit_score_validation/`: model application, validation stats, and figures
- `bubble_heatmap/`, `trend_line_plot/`, `percent_barplot/`: visualization modules
- `wilcoxon_subtype/`, `cox_uni_multi/`: statistical analysis modules
- `utils/`: shared utilities
- `data_processing/`: general data processing helpers (merging clinical info, intersections, filtering, etc.)
- `single_gene_survival/`: single-gene survival analyses and related plots
- `batch_correction/`: batch correction utilities
- `survival_curves/`: survival curve plotting utilities
- `fig7_gm/`: scripts used for Figure 7 (GM module)

## Notes

- Scripts are organized to avoid hard-coded absolute paths where possible.
- Some scripts expect specific input files that are not included in this repo.
