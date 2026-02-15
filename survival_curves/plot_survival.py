import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
import os
from tqdm import tqdm
import time
from multiprocessing import Pool, cpu_count
from functools import partial
import math
import matplotlib

# 全局字体设置：Arial 加粗
matplotlib.rcParams['font.family'] = 'Arial'
matplotlib.rcParams['font.weight'] = 'bold'


def calculate_p_value_for_group(data, gene_name, group_name, strategy='median_rank'):
    """
    为指定分组（hi_hi 或 lo_lo）按中位数排名将基因表达分为高低两组，
    返回：log-rank p值、高/低两组的中位生存期、两组样本量
    """
    group_data = data[data['MYC_PVT1_Status'] == group_name].copy()

    if group_data.empty or group_data[gene_name].isnull().all() or len(group_data) < 4:
        return np.nan, np.nan, np.nan, 0, 0

    group_data.dropna(subset=['OS_months', 'OS_event', gene_name], inplace=True)
    if len(group_data) < 4:
        return np.nan, np.nan, np.nan, 0, 0

    # 仅使用中位数排名策略
    sorted_data = group_data.sort_values(by=gene_name)
    n_samples = len(sorted_data)
    n_half = n_samples // 2
    low_group = sorted_data.iloc[:n_half]
    high_group = sorted_data.iloc[-n_half:]

    if (high_group.empty or low_group.empty or
        high_group['OS_event'].sum() == 0 or low_group['OS_event'].sum() == 0):
        return np.nan, np.nan, np.nan, len(high_group), len(low_group)

    results = logrank_test(
        high_group['OS_months'], low_group['OS_months'],
        event_observed_A=high_group['OS_event'], event_observed_B=low_group['OS_event']
    )
    p_value = results.p_value

    kmf_high = KaplanMeierFitter().fit(high_group['OS_months'], high_group['OS_event'])
    median_survival_high = kmf_high.median_survival_time_

    kmf_low = KaplanMeierFitter().fit(low_group['OS_months'], low_group['OS_event'])
    median_survival_low = kmf_low.median_survival_time_

    return p_value, median_survival_high, median_survival_low, len(high_group), len(low_group)


def process_gene(gene, base_data, strategy):
    """对单个基因执行log-rank检验并返回汇总结果"""
    p_hi_hi, med_surv_h_hi, med_surv_l_hi, n_h_hi, n_l_hi = calculate_p_value_for_group(base_data, gene, 'hi_hi', strategy=strategy)
    p_lo_lo, med_surv_h_lo, med_surv_l_lo, n_h_lo, n_l_lo = calculate_p_value_for_group(base_data, gene, 'lo_lo', strategy=strategy)

    return {
        'Gene': gene,
        'p_value_hi_hi': p_hi_hi,
        'median_survival_hi_hi_high': med_surv_h_hi,
        'median_survival_hi_hi_low': med_surv_l_hi,
        'p_value_lo_lo': p_lo_lo,
        'median_survival_lo_lo_high': med_surv_h_lo,
        'median_survival_lo_lo_low': med_surv_l_lo,
        'hi_hi_high_n': n_h_hi,
        'hi_hi_low_n': n_l_hi,
        'lo_lo_high_n': n_h_lo,
        'lo_lo_low_n': n_l_lo
    }


def plot_combined_survival_curves(data, gene_name, p_hi_hi, p_lo_lo, ax, strategy):
    """在同一坐标轴上绘制 hi_hi 与 lo_lo 两组的高低表达四条生存曲线"""
    original_fontsize = 12
    scaled_fontsize = original_fontsize * 1.5  # 放大1.5倍

    kmf = KaplanMeierFitter()

    # 调整配色：hi_hi更深、lo_lo更浅
    colors = {'hi_hi': ('#8B0000', '#00008B'), 'lo_lo': ('#FFB74D', '#90CAF9')}
    styles = {'hi_hi': '-', 'lo_lo': '--'}

    for group_status in ['hi_hi', 'lo_lo']:
        group_data = data[data['MYC_PVT1_Status'] == group_status].copy()
        if len(group_data) < 4:
            continue

        sorted_data = group_data.sort_values(by=gene_name)
        n_samples = len(sorted_data)
        n_half = n_samples // 2
        low_group = sorted_data.iloc[:n_half]
        high_group = sorted_data.iloc[-n_half:]

        if not high_group.empty:
            label = f'{group_status} & {gene_name}_high (n={len(high_group)})'
            kmf.fit(high_group['OS_months'], high_group['OS_event'], label=label)
            kmf.plot_survival_function(ax=ax, color=colors[group_status][0], linestyle=styles[group_status], ci_show=False, linewidth=3)

        if not low_group.empty:
            label = f'{group_status} & {gene_name}_low (n={len(low_group)})'
            kmf.fit(low_group['OS_months'], low_group['OS_event'], label=label)
            kmf.plot_survival_function(ax=ax, color=colors[group_status][1], linestyle=styles[group_status], ci_show=False, linewidth=3)

    p_hi_hi_text = f'{p_hi_hi:.4f}' if pd.notna(p_hi_hi) else 'N/A'
    p_lo_lo_text = f'{p_lo_lo:.4f}' if pd.notna(p_lo_lo) else 'N/A'
    p_text = f'p-value (hi_hi) = {p_hi_hi_text}\np-value (lo_lo) = {p_lo_lo_text}'
    ax.text(0.1, 0.1, p_text, transform=ax.transAxes, fontsize=scaled_fontsize,
            bbox=dict(boxstyle='round,pad=0.5', fc='wheat', alpha=0.5))

    ax.set_title(f'Survival Analysis for {gene_name} (Strategy: {strategy})', fontsize=scaled_fontsize, fontweight='bold')
    ax.set_xlabel('Overall survival time (Months)', fontsize=scaled_fontsize, fontweight='bold')
    ax.set_ylabel('Probability of Survival', fontsize=scaled_fontsize, fontweight='bold')
    ax.legend(loc='upper right', fontsize=scaled_fontsize, prop={'weight': 'bold'})
    ax.set_ylim(0, 1.05)
    ax.tick_params(axis='both', which='major', labelsize=scaled_fontsize, width=1.6, length=6)

    # 加粗坐标轴边框
    for spine in ax.spines.values():
        spine.set_linewidth(1.8)


def main():
    # --- 策略配置 ---
    print("--- System and Strategy Configuration ---")
    strategy = 'median_rank'  # 只使用中位数排名策略
    print(f"Strategy to run: {strategy}")
    print("----------------------------------------\n")

    # --- 文件路径配置 ---
    expression_file = 'combined_expression_combat_corrected.txt'
    survival_file = 'updated_survival_data.txt'
    annotation_file = 'MYC_PVT1_annotation.txt'

    # --- 加载一次数据 ---
    print("--- Data Loading (once) ---")
    try:
        df_expr_raw = pd.read_csv(expression_file, sep='\s+', index_col=0, engine='python')
        df_expr = df_expr_raw.T
        df_surv = pd.read_csv(survival_file, sep='\s+', engine='python')
        df_annot = pd.read_csv(annotation_file, sep='\s+', engine='python')
    except FileNotFoundError as e:
        print(f"Error loading file: {e}")
        return

    df_annot.rename(columns={'Sample': 'Sample_ID'}, inplace=True)
    merged_data = pd.merge(df_surv, df_annot, on='Sample_ID', how='inner')
    merged_data.set_index('Sample_ID', inplace=True)
    full_data = pd.merge(merged_data, df_expr, left_index=True, right_index=True, how='inner')
    full_data['OS_months'] = pd.to_numeric(full_data['OS_months'], errors='coerce')
    full_data['OS_event'] = pd.to_numeric(full_data['OS_event'], errors='coerce')
    base_data = full_data.dropna(subset=['OS_months', 'OS_event', 'MYC_PVT1_Status']).copy()

    # --- 指定要处理的8个基因 ---
    target_genes = ['ARG1', 'SOX9', 'CDK6', 'HAPLN1', 'ZMYND8', 'FOLH1', 'F2RL1', 'ZHX1']

    for gene in tqdm(target_genes, desc="Converting gene columns to numeric"):
        if gene in base_data.columns:
            base_data[gene] = pd.to_numeric(base_data[gene], errors='coerce')
    print("----------------------------------------\n")

    # --- 对指定基因执行分析 ---
    print(f"########## Running Analysis for {len(target_genes)} Target Genes ##########")
    total_start_time = time.time()

    output_dir = f'survival_plots_{strategy}'
    os.makedirs(output_dir, exist_ok=True)

    # 处理与绘图
    results = []
    for gene in tqdm(target_genes, desc=f"Processing {len(target_genes)} genes"):
        res = process_gene(gene, base_data=base_data, strategy=strategy)
        results.append(res)

    results_df = pd.DataFrame(results)

    for _, gene_info in tqdm(results_df.iterrows(), total=len(results_df), desc=f"Plotting ({strategy})"):
        fig, ax = plt.subplots(figsize=(12, 8))
        plot_combined_survival_curves(base_data, gene_info['Gene'], gene_info['p_value_hi_hi'], gene_info['p_value_lo_lo'], ax, strategy)
        plt.tight_layout()
        output_path = os.path.join(output_dir, f"survival_curve_{gene_info['Gene']}.pdf")
        plt.savefig(output_path, format='pdf', bbox_inches='tight')
        plt.close(fig)

    print(f"\n########## Finished Analysis. Total time: {time.time() - total_start_time:.2f}s ##########\n")


if __name__ == '__main__':
    main() 