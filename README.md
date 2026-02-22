# 8q24 amplification and immunotherapy resistance — analysis code

Analysis scripts for the manuscript:  
**"Pan-cancer 8q24 amplification predicts primary immunotherapy resistance and therapeutic vulnerabilities"**  
*(iScience, ISCIENCE-D-25-20149)*

---

## Repository structure

```
benefit_score_validation/   # Benefit score model: training (IMvigor210) → validation (RIAZ/LIU/MSK)
  validate_hugo.R             apply frozen cut-offs to HUGO cohort
  validate_liu.R              apply frozen cut-offs to LIU cohort
  validate_riaz.R             apply frozen cut-offs to RIAZ cohort
  validate_riaz_excl_ne.R    RIAZ with NE samples excluded
  merge_cohorts.R             pool RIAZ+LIU+MSK for trend test
  hugo_summary.R              HUGO summary stats
  stats_3group.R              Cochran-Armitage trend test (Low/Inter/High)
  stats_fisher.R              Low vs High Fisher tests
  plot_hugo.R / plot_liu.R / plot_msk.R / plot_merged.R   fraction barplots
  fig_supp6.R                 Supplementary Figure 6

single_gene_survival/       # Vulnerability Index & interaction-term Cox analysis
  cox_interaction_term.py     Cox model with Gene×8q24_status interaction term
  single_gene_cox_analysis.py genome-wide single-gene Cox scan
  filter_hiHiRisk_loLoNoRisk.py  filter by hi_hi risk / lo_lo no-risk criteria
  forest_plot_3panel.py       3-panel forest plot (main figure)
  final_combined_forest_plot.py  combined forest plot
  forest_plot_like.py / forest_plot_like.R
  forest.R                    forest plot utilities
  cox_main.py                 Cox analysis entry point
  enrichment_analysis.py / enrichment_main.py  pathway enrichment
  draw_waterfall_plot.py
  query_interaction_for_hrr_genes.py
  interaction_for_hrr_genes_hilo_within_status.py
  verify_epo_analysis.py
  cell_lines/                 cell-line-level analyses

bubble_heatmap/             # Bubble heatmap (Fig 3, immune infiltration)
  bubble_heatmap.R
  statistical_tests.R

percent_barplot/            # Percentage barplots
  percent_barplot.R

trend_line_plot/            # TIL map trend line plots
  trend_line_plot.R

survival_curves/            # Kaplan-Meier survival curves
  plot_survival.py
  batch_survival_curves_filtered.R
  plot_survival_parallel.py
  check_distribution.py

batch_correction/           # ComBat batch correction for multi-cohort integration
  batch_correction_combat.py
  plot_3d_pca.R

cox_uni_multi/              # Univariate & multivariate Cox PH (Supplementary)
  cox_uni_multi_ph.R

wilcoxon_subtype/           # Wilcoxon rank-sum tests across subtypes
  wilcox_ranksum_subtype.R

data_processing/            # Data loading, merging, filtering utilities
  merge_clinical.py
  preprocess.py
  filter_data.py
  intersect.py
  read_data.R
  diagnosis.py
  scan_for_tpm.py

utils/                      # Shared utilities
  cell_line_style_B.py
  rscript_utils.R

fig7_gm/                    # Figure 7 GM module (gene-module analysis)
  GM/Code/                  R scripts for GM figure panels
  GM/data/sig/              signature gene sets
```

---

## Key analysis notes

- **Benefit score cut-offs** are defined once in the IMvigor210 training set (tertiles of predicted probability: Low/Mid = 0.1266, Mid/High = 0.2660) and applied unchanged to all external cohorts.
- **Vulnerability Index (VI)** = exp(β_interaction) from a Cox model with an interaction term (Gene × 8q24_status), not a simple HR ratio.
- **CRISPR integration** uses a two-layer strategy: (1) DepMap 8q24-specific dependencies; (2) external immune co-culture screens for immune relevance.

## Notes

- Scripts avoid hard-coded absolute paths where possible; input file paths may need to be adjusted.
- Some input files (raw expression matrices, clinical tables) are not included due to data sharing restrictions.
