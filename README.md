# Immunotherapy cohort analysis

This repository contains analysis scripts used for a multi-cohort immunotherapy study, including:

- Benefit score validation (IMvigor210-trained logistic model applied to external cohorts)
- Pan-cancer immune analyses and visualizations (bubble heatmap, trend line plot, percent barplot)
- Additional supporting scripts used during the project (see `extra/`)

## Repository structure

- `benefit_score_validation/`: model application, validation stats, and figures
- `bubble_heatmap/`, `trend_line_plot/`, `percent_barplot/`: visualization modules
- `wilcoxon_subtype/`, `cox_uni_multi/`: statistical analysis modules
- `utils/`: shared utilities
- `extra/`: other scripts from the project workspace (single-gene survival, batch correction, survival curves, fig7/GM, etc.)

## Notes

- Scripts are organized to avoid hard-coded absolute paths where possible.
- Some scripts expect specific input files that are not included in this repo.
