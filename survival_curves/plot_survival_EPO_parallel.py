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

def calculate_p_value_for_group(data, gene_name, group_name):
    """为指定的分组计算基因高低表达生存差异的p-value和样本数"""
    group_data = data[data['MYC_PVT1_Status'] == group_name].copy()

    # 确保组内有数据且基因表达值不全为空
    if group_data.empty or group_data[gene_name].isnull().all():
        return np.nan, 0, 0

    # 确保生存数据有效
    group_data.dropna(subset=['OS_months', 'OS_event'], inplace=True)
    if group_data.empty:
        return np.nan, 0, 0

    median_expression = group_data[gene_name].median()
    
    # 如果中位数为 NaN (例如，列中所有值都为 NaN)，则无法分组
    if pd.isna(median_expression):
        return np.nan, 0, 0
        
    group_data['Expression_Group'] = np.where(group_data[gene_name] >= median_expression, 'high', 'low')

    high_group = group_data[group_data['Expression_Group'] == 'high']
    low_group = group_data[group_data['Expression_Group'] == 'low']

    # 确保高表达和低表达组都有事件发生，否则logrank检验会出错
    if (high_group.empty or low_group.empty or
        high_group['OS_event'].sum() == 0 or low_group['OS_event'].sum() == 0):
        return np.nan, len(high_group), len(low_group)

    results = logrank_test(
        high_group['OS_months'], low_group['OS_months'],
        event_observed_A=high_group['OS_event'], event_observed_B=low_group['OS_event']
    )
    return results.p_value, len(high_group), len(low_group)


def process_gene(gene, base_data):
    """
    对单个基因执行完整的p-value计算。
    此函数设计为可被并行处理器调用。
    """
    p_hi_hi, n_hi_high, n_hi_low = calculate_p_value_for_group(base_data, gene, 'hi_hi')
    p_lo_lo, n_lo_high, n_lo_low = calculate_p_value_for_group(base_data, gene, 'lo_lo')
    
    return {
        'Gene': gene,
        'p_value_hi_hi': p_hi_hi,
        'p_value_lo_lo': p_lo_lo,
        'hi_hi_high_n': n_hi_high,
        'hi_hi_low_n': n_hi_low,
        'lo_lo_high_n': n_lo_high,
        'lo_lo_low_n': n_lo_low
    }


def plot_combined_survival_curves(data, gene_name, p_hi_hi, p_lo_lo, ax):
    """在同一个坐标轴上绘制四条生存曲线"""
    kmf = KaplanMeierFitter()

    # 定义颜色和线型
    colors = {'hi_hi': ('red', 'blue'), 'lo_lo': ('darkorange', 'deepskyblue')}
    styles = {'hi_hi': '-', 'lo_lo': '--'}

    for group_status in ['hi_hi', 'lo_lo']:
        group_data = data[data['MYC_PVT1_Status'] == group_status].copy()
        if group_data.empty or group_data[gene_name].isnull().all():
            continue

        median_expression = group_data[gene_name].median()
        group_data['Expression_Group'] = np.where(group_data[gene_name] >= median_expression, 'high', 'low')

        # 绘制高表达曲线
        high_df = group_data[group_data['Expression_Group'] == 'high']
        if not high_df.empty:
            kmf.fit(high_df['OS_months'], high_df['OS_event'], label=f'{group_status} & {gene_name}_high (n={len(high_df)})')
            kmf.plot_survival_function(ax=ax, color=colors[group_status][0], linestyle=styles[group_status], ci_show=False)

        # 绘制低表达曲线
        low_df = group_data[group_data['Expression_Group'] == 'low']
        if not low_df.empty:
            kmf.fit(low_df['OS_months'], low_df['OS_event'], label=f'{group_status} & {gene_name}_low (n={len(low_df)})')
            kmf.plot_survival_function(ax=ax, color=colors[group_status][1], linestyle=styles[group_status], ci_show=False)
    
    # 标注 p-values, 优雅地处理NaN值
    p_hi_hi_text = f'{p_hi_hi:.4f}' if pd.notna(p_hi_hi) else 'N/A'
    p_lo_lo_text = f'{p_lo_lo:.4f}' if pd.notna(p_lo_lo) else 'N/A'
    p_text = f'p-value (hi_hi) = {p_hi_hi_text}\np-value (lo_lo) = {p_lo_lo_text}'
    ax.text(0.1, 0.1, p_text, transform=ax.transAxes, fontsize=12,
            bbox=dict(boxstyle='round,pad=0.5', fc='wheat', alpha=0.5))

    # 设置图表样式
    ax.set_title(f'Survival Analysis for {gene_name}')
    ax.set_xlabel('Overall survival time (Months)')
    ax.set_ylabel('Probability of Survival')
    ax.legend(loc='upper right')
    ax.set_ylim(0, 1.05)


def main():
    # --- 系统信息输出 ---
    print("--- 系统信息 ---")
    try:
        cores_to_use = cpu_count()
        print(f"检测到 CPU 逻辑核心 (线程): {cores_to_use} 个")
        print(f"模式: 多核心并行计算 (专为Linux服务器优化)")
    except NotImplementedError:
        cores_to_use = 1
        print("无法确定CPU核心数，将使用单线程模式。")
    print("------------------\n")


    # --- 文件路径配置 ---
    expression_file = 'combined_expression_combat_corrected.txt'
    survival_file = 'updated_survival_data.txt'
    annotation_file = 'MYC_PVT1_annotation.txt'
    output_dir = 'survival_plots_filtered'
    results_file = 'all_genes_analysis_results.txt'

    # --- 计时器初始化 ---
    total_start_time = time.time()
    step_start_time = time.time()

    # --- 加载和预处理数据 ---
    print("步骤 1/5: 开始加载和预处理数据...")
    os.makedirs(output_dir, exist_ok=True)

    try:
        print("  - 正在加载表达文件...")
        df_expr_raw = pd.read_csv(expression_file, sep='\s+', index_col=0, engine='python')
        print(f"    > 成功加载 {df_expr_raw.shape[0]} 个基因, {df_expr_raw.shape[1]} 个样本。")
    except FileNotFoundError:
        print(f"错误: 表达文件 '{expression_file}' 未找到。")
        return
    df_expr = df_expr_raw.T
    
    try:
        print("  - 正在加载生存文件...")
        df_surv = pd.read_csv(survival_file, sep='\s+', engine='python')
    except FileNotFoundError:
        print(f"错误: 生存文件 '{survival_file}' 未找到。")
        return

    try:
        print("  - 正在加载注释文件...")
        df_annot = pd.read_csv(annotation_file, sep='\s+', engine='python')
    except FileNotFoundError:
        print(f"错误: 注释文件 '{annotation_file}' 未找到。")
        return

    # --- 合并数据 ---
    print("  - 正在合并数据表...")
    df_annot.rename(columns={'Sample': 'Sample_ID'}, inplace=True)
    merged_data = pd.merge(df_surv, df_annot, on='Sample_ID', how='inner')
    merged_data.set_index('Sample_ID', inplace=True)
    full_data = pd.merge(merged_data, df_expr, left_index=True, right_index=True, how='inner')
    print(f"    > 数据合并完成。总共 {len(full_data)} 个样本进入分析。")
    
    full_data['OS_months'] = pd.to_numeric(full_data['OS_months'], errors='coerce')
    full_data['OS_event'] = pd.to_numeric(full_data['OS_event'], errors='coerce')
    base_data = full_data.dropna(subset=['OS_months', 'OS_event', 'MYC_PVT1_Status']).copy()
    print(f"步骤 1/5: 数据加载和预处理完成。耗时: {time.time() - step_start_time:.2f} 秒。\n")
    step_start_time = time.time()

    # --- 筛选符合条件的基因 (并行处理) ---
    all_genes = df_expr.columns.tolist()

    # 优化：一次性转换所有基因列的数据类型
    print(f"步骤 2/5: 预转换 {len(all_genes)} 个基因列的数据类型...")
    for gene in tqdm(all_genes, desc="数据类型转换"):
        if gene in base_data.columns:
            base_data[gene] = pd.to_numeric(base_data[gene], errors='coerce')
    print(f"步骤 2/5: 数据类型转换完成。耗时: {time.time() - step_start_time:.2f} 秒。\n")
    step_start_time = time.time()

    # --- 并行计算 ---
    # 根据核心数和任务总数计算合适的分块大小，以提高效率
    chunk_size = 1
    if cores_to_use > 1:
        chunk_size = math.ceil(len(all_genes) / (cores_to_use * 4))
    
    print(f"步骤 3/5: 开始使用 {cores_to_use} 个核心并行扫描 {len(all_genes)} 个基因 (分块大小: {chunk_size})...")
    
    # 使用 functools.partial 包装 process_gene 函数，以固定 base_data 参数
    process_func = partial(process_gene, base_data=base_data)

    all_gene_data = []
    # 使用 with 语句确保进程池被正确关闭
    with Pool(processes=cores_to_use) as pool:
        # 使用 imap_unordered 以在结果可用时立即处理它们，并与 tqdm 集成以显示进度条
        results_iterator = pool.imap_unordered(process_func, all_genes, chunksize=chunk_size)
        
        for result in tqdm(results_iterator, total=len(all_genes), desc="基因筛选进度"):
            if result:
                all_gene_data.append(result)

    print(f"\n步骤 3/5: 扫描完成。耗时: {time.time() - step_start_time:.2f} 秒。\n")
    step_start_time = time.time()

    # --- 保存结果到文件 ---
    print(f"步骤 4/5: 正在保存所有基因的分析结果...")
    all_results_df = pd.DataFrame(all_gene_data)
    # 按 p_value_hi_hi 排序，方便查看最显著的结果
    all_results_df.sort_values(by='p_value_hi_hi', inplace=True)
    all_results_df.to_csv(results_file, sep='\t', index=False, float_format='%.4f')
    print(f"  - 分析了 {len(all_results_df)} 个基因。结果已保存到 '{results_file}'。")
    print(f"步骤 4/5: 结果保存完成。耗时: {time.time() - step_start_time:.2f} 秒。\n")
    step_start_time = time.time()

    # --- 仅为 'EPO' 基因绘图 ---
    print(f"步骤 5/5: 开始为特定基因 'EPO' 绘图...")
    gene_to_plot = 'EPO'

    # 检查 'EPO' 是否在分析结果中
    if 'Gene' in all_results_df.columns:
        gene_info_row = all_results_df[all_results_df['Gene'] == gene_to_plot]

        if not gene_info_row.empty:
            gene_info = gene_info_row.iloc[0]
            
            # 检查p值是否有效，再进行绘图
            if pd.notna(gene_info['p_value_hi_hi']) or pd.notna(gene_info['p_value_lo_lo']):
                fig, ax = plt.subplots(figsize=(12, 8))
                
                plot_combined_survival_curves(base_data, gene_to_plot, gene_info['p_value_hi_hi'], gene_info['p_value_lo_lo'], ax)

                plt.tight_layout()
                output_path = os.path.join(output_dir, f'survival_curve_{gene_to_plot}.png')
                plt.savefig(output_path)
                plt.close(fig)
                print(f"  - 成功为基因 '{gene_to_plot}' 创建并保存生存曲线图。")
            else:
                print(f"  - 未能为 '{gene_to_plot}' 绘图，因为无法为其计算有效的p-value。")
        else:
            # 这种情况理论上不应该发生，因为我们分析了所有基因
            print(f"  - 奇怪的错误：未能为 '{gene_to_plot}' 绘图，因为它不在最终的分析结果中。")
    else:
        print(f"  - 基因分析结果为空，无法为 '{gene_to_plot}' 绘图。")

    print(f"步骤 5/5: 绘图完成。耗时: {time.time() - step_start_time:.2f} 秒。")
    print(f"\n任务完成。总运行时间: {time.time() - total_start_time:.2f} 秒 ---")


if __name__ == '__main__':
    # 在Windows上, multiprocessing 需要这个保护
    # 在Linux上虽然不是必须的，但这是最佳实践
    main() 