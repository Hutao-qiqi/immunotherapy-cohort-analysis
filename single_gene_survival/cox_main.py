#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Single-gene survival analysis.

Runs Cox regression for each gene within MYC_PVT1 status groups (hi_hi / lo_lo),
and exports per-group results.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import mannwhitneyu
import seaborn as sns
from lifelines import CoxPHFitter
from lifelines.statistics import logrank_test
import warnings
warnings.filterwarnings('ignore')
import argparse
from pathlib import Path
from joblib import Parallel, delayed
from tqdm import tqdm

# 设置matplotlib的字体和样式
plt.rcParams['font.size'] = 12
plt.rcParams['font.family'] = 'Arial'
plt.rcParams['axes.linewidth'] = 1.5
plt.rcParams['axes.spines.top'] = False
plt.rcParams['axes.spines.right'] = False
plt.rcParams['xtick.major.size'] = 6
plt.rcParams['ytick.major.size'] = 6
plt.rcParams['xtick.major.width'] = 1.5
plt.rcParams['ytick.major.width'] = 1.5


def _default_base_dir() -> Path:
    return Path(__file__).resolve().parent


def _default_survival_file(base_dir: Path) -> Path:
    # Prefer survival_curves/ if present in this repo layout
    candidate = (base_dir / ".." / "survival_curves" / "updated_survival_data.txt").resolve()
    return candidate


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Single-gene Cox survival analysis (hi_hi / lo_lo)")
    parser.add_argument(
        "--base-dir",
        type=Path,
        default=_default_base_dir(),
        help="Directory containing combined_expression_combat_corrected.txt and MYC_PVT1_annotation.txt",
    )
    parser.add_argument(
        "--survival-file",
        type=Path,
        default=None,
        help="Path to updated_survival_data.txt (default: ../survival_curves/updated_survival_data.txt)",
    )
    parser.add_argument(
        "--strategy",
        type=str,
        default="median_rank",
        choices=["median", "median_rank", "quartiles"],
        help="Expression grouping strategy",
    )
    parser.add_argument(
        "--n-jobs",
        type=int,
        default=0,
        help="Parallel jobs (0 = use a safe default)",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=_default_base_dir() / "outputs",
        help="Output directory for result tables",
    )
    return parser.parse_args()

def load_and_prepare_data(base_dir: Path, survival_file: Path):
    """
    加载并准备数据，确保样本ID匹配
    """
    print("正在加载数据并检查样本ID匹配...")

    # 1. 加载表达数据
    try:
        expression_data = pd.read_csv(base_dir / 'combined_expression_combat_corrected.txt', sep='\t', index_col=0)
        print(f"✓ 表达数据加载成功: {expression_data.shape}")
    except FileNotFoundError:
        print("❌ 未找到表达数据文件")
        return None, None, None
    
    # 2. 加载分组信息
    try:
        annotation_df = pd.read_csv(base_dir / 'MYC_PVT1_annotation.txt', sep='\t')
        group_info = pd.Series(annotation_df['MYC_PVT1_Status'].values, 
                              index=annotation_df['Sample'].values)
        print(f"✓ 分组信息加载成功: {len(group_info)} 个样本")
        print(f"  分组分布: {group_info.value_counts().to_dict()}")
    except FileNotFoundError:
        print("❌ 未找到分组信息文件")
        return None, None, None
    
    # 3. 加载生存数据
    try:
        survival_data = pd.read_csv(survival_file, sep='\s+', engine='python')
        survival_data = survival_data.set_index('Sample_ID')
        print(f"✓ 生存数据加载成功: {len(survival_data)} 个样本")
        if 'Dataset' in survival_data.columns:
            print(f"  数据集分布: {survival_data['Dataset'].value_counts().to_dict()}")
    except FileNotFoundError:
        print("❌ 未找到生存数据文件")
        return None, None, None
    
    # 4. 直接匹配三个数据集（现在都是带前缀格式）
    print("  检查三个数据集的样本匹配...")
    
    expr_samples = list(expression_data.columns)
    group_samples = list(group_info.index)
    surv_samples = list(survival_data.index)
    
    # 找到三个数据集的共同样本
    expr_samples_set = set(expr_samples)
    group_samples_set = set(group_samples)
    surv_samples_set = set(surv_samples)
    
    common_samples = expr_samples_set & group_samples_set & surv_samples_set
    print(f"✓ 三个数据集共同样本数: {len(common_samples)}")
    
    if len(common_samples) < 50:
        print("❌ 共同样本数太少，无法进行有效分析")
        return None, None, None
    
    # 5. 筛选共同样本的数据
    expression_data_filtered = expression_data[list(common_samples)]
    group_info_filtered = group_info[list(common_samples)]
    survival_data_filtered = survival_data.loc[list(common_samples)]
    
    # 6. 检查生存数据完整性并进一步筛选
    survival_complete = survival_data_filtered.dropna(subset=['OS_months', 'OS_event'])
    print(f"  完整生存数据: {len(survival_complete)} 个样本")
    
    # 7. 确保三个数据集样本ID完全一致（去重）
    final_samples = survival_complete.index.unique().tolist()
    expression_data_final = expression_data_filtered[final_samples]
    group_info_final = group_info_filtered[final_samples]
    survival_complete_final = survival_complete.loc[final_samples]
    
    # 8. 检查各组样本数
    group_counts = group_info_final.value_counts()
    print(f"  hi_hi组: {group_counts.get('hi_hi', 0)} 个样本")
    print(f"  lo_lo组: {group_counts.get('lo_lo', 0)} 个样本")
    
    return expression_data_final, group_info_final, survival_complete_final

def perform_cox_analysis(expression_data, survival_data, gene_name, group_info, strategy='median'):
    """
    对单个基因进行Cox回归分析 (支持不同分组策略)
    strategy: 'median', 'median_rank', 'quartiles'
    """
    try:
        common_samples = list(expression_data.columns)
        
        if len(common_samples) < 10:
            return None, None
        
        if gene_name not in expression_data.index:
            return None, None
            
        gene_expr = expression_data.loc[gene_name, common_samples]
        
        results = {}
        
        for group in ['hi_hi', 'lo_lo']:
            group_samples = [s for s in common_samples if group_info[s] == group]
            
            if len(group_samples) < 10:
                continue
            
            analysis_df = pd.DataFrame({
                'expression': gene_expr[group_samples],
                'OS_time': survival_data.loc[group_samples, 'OS_months'],
                'OS_event': survival_data.loc[group_samples, 'OS_event']
            }).dropna()
            
            if len(analysis_df) < 10:
                continue

            # --- 核心分组逻辑 ---
            if strategy == 'median':
                median_expr = analysis_df['expression'].median()
                if pd.isna(median_expr) or analysis_df['expression'].std() == 0: continue
                analysis_df = analysis_df[analysis_df['expression'] != median_expr]
                analysis_df['Expression_Group'] = (analysis_df['expression'] > median_expr).astype(int)
            
            elif strategy == 'median_rank':
                sorted_df = analysis_df.sort_values('expression').reset_index()
                half_point = len(sorted_df) // 2
                low_group_indices = sorted_df.iloc[:half_point].index
                high_group_indices = sorted_df.iloc[half_point:].index
                sorted_df['Expression_Group'] = 0
                sorted_df.loc[high_group_indices, 'Expression_Group'] = 1
                analysis_df = sorted_df.set_index('index')

            elif strategy == 'quartiles':
                q1 = analysis_df['expression'].quantile(0.25)
                q3 = analysis_df['expression'].quantile(0.75)
                if q1 == q3: continue
                analysis_df = analysis_df[(analysis_df['expression'] < q1) | (analysis_df['expression'] > q3)]
                analysis_df['Expression_Group'] = (analysis_df['expression'] > q3).astype(int)
            
            else:
                raise ValueError(f"未知策略: {strategy}")

            high_group = analysis_df[analysis_df['Expression_Group'] == 1]
            low_group = analysis_df[analysis_df['Expression_Group'] == 0]

            if high_group.empty or low_group.empty or high_group['OS_event'].sum() == 0 or low_group['OS_event'].sum() == 0:
                continue
            
            analysis_df = analysis_df[analysis_df['OS_time'] > 0]
            
            if len(analysis_df) < 10:
                continue
            
            try:
                cph = CoxPHFitter()
                cph.fit(analysis_df[['OS_time', 'OS_event', 'Expression_Group']], 
                       duration_col='OS_time', event_col='OS_event')
                
                hr = cph.hazard_ratios_['Expression_Group']
                p_value = cph.summary.loc['Expression_Group', 'p']
                ci_lower = cph.summary.loc['Expression_Group', 'exp(coef) lower 95%']
                ci_upper = cph.summary.loc['Expression_Group', 'exp(coef) upper 95%']
                
                results[group] = {
                    'HR': hr,
                    'p_value': p_value,
                    'CI_lower': ci_lower,
                    'CI_upper': ci_upper,
                    'n_samples': len(analysis_df),
                    'n_high_group': len(high_group),
                    'n_low_group': len(low_group)
                }
            except Exception:
                continue
        
        return results.get('hi_hi'), results.get('lo_lo')
    
    except Exception:
        return None, None


def main():
    """
    主函数 (仅执行分析)
    """
    args = parse_args()
    base_dir = args.base_dir.resolve()
    survival_file = (args.survival_file.resolve() if args.survival_file else _default_survival_file(base_dir))
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 60)
    print("单基因生存分析 - V4 (仅分析)")
    print("=" * 60)

    strategy = args.strategy
    
    # 1. 加载和准备数据 (只需加载一次)
    print("\n1. 加载和准备数据...")
    expression_data, group_info, survival_data = load_and_prepare_data(base_dir=base_dir, survival_file=survival_file)
    
    if expression_data is None:
        print("❌ 数据加载失败，程序退出")
        return

    print("\n" + "=" * 60)
    print(f"开始执行策略: {strategy}")
    print("=" * 60)

    # 2. 进行单基因Cox分析
    default_jobs = 8
    n_jobs = args.n_jobs if args.n_jobs and args.n_jobs > 0 else default_jobs
    print(f"\n2. 使用 {n_jobs} 个并行任务进行单基因Cox分析 (策略: {strategy})...")

    genes_to_analyze = expression_data.index[:]

    results = Parallel(n_jobs=n_jobs)(
        delayed(perform_cox_analysis)(expression_data, survival_data, gene, group_info, strategy)
        for gene in tqdm(genes_to_analyze, desc=f"分析进度 ({strategy})")
    )

    # 处理结果
    results_hi_hi, results_lo_lo = {}, {}
    for gene, (result_hi, result_lo) in zip(genes_to_analyze, results):
        if result_hi:
            results_hi_hi[gene] = result_hi
        if result_lo:
            results_lo_lo[gene] = result_lo

    print(f"\n分析完成 (策略: {strategy})!")
    print(f"✓ hi_hi组有效结果: {len(results_hi_hi)} 个基因")
    print(f"✓ lo_lo组有效结果: {len(results_lo_lo)} 个基因")

    # 3. 保存完整结果 (根据策略命名)
    print("\n3. 保存完整结果...")
    output_files = []

    if results_hi_hi:
        df_hi_hi = pd.DataFrame.from_dict(results_hi_hi, orient='index')
        df_hi_hi.index.name = 'Gene'
        out_hi = output_dir / f'hi_hi_cox_results_{strategy}.txt'
        df_hi_hi.to_csv(out_hi, sep='\t')
        output_files.append(str(out_hi))

    if results_lo_lo:
        df_lo_lo = pd.DataFrame.from_dict(results_lo_lo, orient='index')
        df_lo_lo.index.name = 'Gene'
        out_lo = output_dir / f'lo_lo_cox_results_{strategy}.txt'
        df_lo_lo.to_csv(out_lo, sep='\t')
        output_files.append(str(out_lo))

    print(f"\n策略 '{strategy}' 执行完毕。输出文件:")
    for file in output_files:
        print(f"  • {file}")
    
    print("\n" + "=" * 60)
    print("✅ 所有分析任务完成！")
    print("\n请现在运行 'draw_waterfall_plot.py' 脚本来生成瀑布图。")
    print("=" * 60)

if __name__ == "__main__":
    main() 