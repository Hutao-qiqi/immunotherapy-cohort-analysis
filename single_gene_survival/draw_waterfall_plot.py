#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
瀑布图绘制脚本
- 从Cox分析结果文件读取数据
- 绘制符合出版要求的瀑布图
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter
import os
import time

def create_waterfall_plot(results_df, group_name, strategy):
    """
    创建符合出版要求的瀑布图 (参考IMvigor210风格)
    """
    print(f"为 {group_name} 组创建瀑布图...")

    # 1. 准备数据
    df_plot = results_df.sort_values('HR').reset_index(drop=True)
    df_plot['rank'] = df_plot.index
    
    # 计算用于气泡大小的-log10(p)
    df_plot['-log10(p)'] = -np.log10(df_plot['p_value'].astype(float).clip(1e-300))

    # 为所有基因设置气泡大小
    min_size, max_size = 20, 2500  # 调整尺寸范围以增强对比度
    p_log_values = df_plot['-log10(p)']
    size_range = np.clip(p_log_values, 0.5, 7)  # 限制在0.5-7范围内, 适配新的p值范围
    # 使用更高次幂(4)来急剧拉大p值大小的梯度
    power_exponent = 4
    norm_factor = np.power(7 / 0.5, power_exponent) - 1
    df_plot['size'] = (min_size + (np.power(size_range/0.5, power_exponent) - 1) * (max_size - min_size) / norm_factor) * 1.3

 
    
    # 3. 识别要标注的基因
    # 选择HR最高和最低的各4个基因（不考虑p值）
    if len(df_plot) >= 8:
        genes_to_label = pd.concat([
            df_plot.nsmallest(4, 'HR'),  # HR最低的4个
            df_plot.nlargest(4, 'HR')    # HR最高的4个
        ]).drop_duplicates()
    else:
        genes_to_label = df_plot  # 如果基因不足8个，则全部标注
    
    labeled_genes_set = set(genes_to_label['Gene'])
    print(f"  将标注 {len(labeled_genes_set)} 个基因")

    # 4. 创建图形
    fig, ax = plt.subplots(figsize=(9.89, 5.02)) # 宽高比更新为 98.874 : 50.173
    plt.subplots_adjust(right=0.75)  # 为右侧图例留出更多空间(25%)

    # 5. 绘图元素
    # a) 绘制所有未被标注的基因点 (灰色空心圈)
    other_genes_df = df_plot[~df_plot['Gene'].isin(labeled_genes_set)]
    ax.scatter(other_genes_df['rank'], other_genes_df['HR'], s=other_genes_df['size'], 
               facecolors='none', edgecolors='grey', linewidth=0.8, alpha=0.5, zorder=1)

    # b) 绘制被标注的基因点 (红蓝色区分，大小根据p值变化)
    if not genes_to_label.empty:
        # 按HR值分组
        red_genes = genes_to_label[genes_to_label['HR'] > 1]
        blue_genes = genes_to_label[genes_to_label['HR'] < 1]
        
        # 绘制红色和蓝色气泡
        if not red_genes.empty:
            red_df = df_plot[df_plot['Gene'].isin(red_genes['Gene'])]
            ax.scatter(red_df['rank'], red_df['HR'], s=red_df['size'] * 2,
                      facecolors='none', edgecolors='red', linewidth=1.5, zorder=4, alpha=0.9)
        
        if not blue_genes.empty:
            blue_df = df_plot[df_plot['Gene'].isin(blue_genes['Gene'])]
            ax.scatter(blue_df['rank'], blue_df['HR'], s=blue_df['size'] * 2,
                      facecolors='none', edgecolors='blue', linewidth=1.5, zorder=4, alpha=0.9)

    # 6. 标注基因名称 (优化布局避免重叠)
    if not genes_to_label.empty:
        # 高HR基因 (右侧)
        high_hr_genes = genes_to_label[genes_to_label['HR'] > 1].sort_values('HR', ascending=False)
        
        # 优化标签位置 - 右侧
        if len(high_hr_genes) > 0:
            log_hr_values = np.log10(high_hr_genes['HR'].values)
            y_text_min, y_text_max = 1.2, 2.4 # 标签Y轴范围

            # 将log(HR)值映射到标签的Y轴范围
            if log_hr_values.max() > log_hr_values.min():
                y_positions = y_text_min + (log_hr_values - log_hr_values.min()) / (log_hr_values.max() - log_hr_values.min()) * (y_text_max - y_text_min)
            else:
                y_positions = np.linspace(y_text_max, y_text_min, len(high_hr_genes))

            # 排序是降序的，所以y_positions也应该是降序
            y_positions = np.sort(y_positions)[::-1]

            # 防止标签重叠的斥力逻辑
            min_gap = 0.08 # 最小垂直间距
            for i in range(len(y_positions) - 1):
                if y_positions[i+1] > y_positions[i] - min_gap:
                    y_positions[i+1] = y_positions[i] - min_gap

            for i, (_, row) in enumerate(high_hr_genes.iterrows()):
                ax.annotate(row['Gene'], 
                           xy=(row['rank'], row['HR']),
                           xytext=(16500, y_positions[i]), 
                           textcoords='data',
                           fontsize=9, color='red', fontweight='bold', # 减小字号
                           ha='right', va='center',
                           arrowprops=dict(arrowstyle='-', color='red', lw=1.0),
                           zorder=10)
        
        # 低HR基因 (左侧)
        low_hr_genes = genes_to_label[genes_to_label['HR'] < 1].sort_values('HR')
        
        # 优化标签位置 - 左侧
        if len(low_hr_genes) > 0:
            log_hr_values = np.log10(low_hr_genes['HR'].values)
            y_text_min, y_text_max = 0.3, 0.5 # 标签Y轴范围
            
            if log_hr_values.max() > log_hr_values.min():
                y_positions = y_text_min + (log_hr_values - log_hr_values.min()) / (log_hr_values.max() - log_hr_values.min()) * (y_text_max - y_text_min)
            else:
                y_positions = np.linspace(y_text_min, y_text_max, len(low_hr_genes))

            # 排序是升序的, y_positions也应该是升序
            y_positions = np.sort(y_positions)

            # 防止标签重叠的斥力逻辑
            min_gap = 0.05 # 最小垂直间距
            for i in range(1, len(y_positions)):
                if y_positions[i] < y_positions[i-1] + min_gap:
                    y_positions[i] = y_positions[i-1] + min_gap

            for i, (_, row) in enumerate(low_hr_genes.iterrows()):
                ax.annotate(row['Gene'], 
                           xy=(row['rank'], row['HR']),
                           xytext=(1000, y_positions[i]), 
                           textcoords='data',
                           fontsize=9, color='blue', fontweight='bold', # 减小字号
                           ha='left', va='center',
                           arrowprops=dict(arrowstyle='-', color='blue', lw=1.0),
                           zorder=10)

    # 7. 设置坐标轴、标题和网格
    ax.axhline(y=1, color='black', linestyle='--', alpha=0.7, linewidth=1.0)
    
    n_total = len(df_plot)
    n_sig = len(df_plot[df_plot['p_value'] < 0.05])
    ax.set_title(f'{group_name.replace("_", " ")} Group: Hazard Ratios (n={n_total}, sig={n_sig})', fontsize=12, fontweight='bold')
    ax.set_xlabel('Gene Rank (sorted by HR)', fontsize=10)
    ax.set_ylabel('Hazard Ratio (HR)', fontsize=10)
    ax.set_yscale('log')
    
    # 设置Y轴范围和刻度，为大的圆圈留出更多空间
    y_min = max(0.1, df_plot['HR'].min() * 0.7)
    y_max = min(3.5, df_plot['HR'].max() * 1.5)
    ax.set_ylim(y_min, y_max)
    
    # 强制Y轴使用小数格式，不使用科学计数法
    from matplotlib.ticker import ScalarFormatter
    formatter = ScalarFormatter(useOffset=False)
    formatter.set_scientific(False)
    ax.yaxis.set_major_formatter(formatter)
    
    # 设置X轴范围
    ax.set_xlim(-1000, 18000)
    
    ax.tick_params(axis='both', which='major', labelsize=8)
    ax.grid(True, which="both", ls="--", linewidth=0.5, alpha=0.3)
    
    # 恢复所有边框并设置样式
    ax.spines['top'].set_visible(True)
    ax.spines['right'].set_visible(True)
    
    # 设置所有边框为黑色细线
    for spine in ['top', 'bottom', 'left', 'right']:
        ax.spines[spine].set_linewidth(1.0)
        ax.spines[spine].set_color('black')
    
    # 8. 创建图例 - 放在图表右侧
    # 基因类型图例
    legend_elements = [
        Line2D([0], [0], marker='o', color='w', label='Genes with extreme HR (HR > 1)',
               markerfacecolor='none', markeredgecolor='red', markersize=10, markeredgewidth=1.5, alpha=0.9),
        Line2D([0], [0], marker='o', color='w', label='Genes with extreme HR (HR < 1)',
               markerfacecolor='none', markeredgecolor='blue', markersize=10, markeredgewidth=1.5, alpha=0.9),
        Line2D([0], [0], color='black', linestyle='--', lw=1.0, label='HR = 1')
    ]
    
    # 创建p值图例 - 更明显的大小差异和更多梯度
    p_values = [0.05, 0.01, 1e-3, 1e-4, 1e-5]
    p_labels = ['p < 0.05', 'p < 0.01', 'p < 0.001', 'p < 10⁻⁴', 'p < 10⁻⁵']
    p_log_values = -np.log10(p_values)
    # 使用与主图完全相同的计算方式，移除导致失真的乘数
    size_range = np.clip(p_log_values, 0.5, 7)
    sizes = (min_size + (np.power(size_range/0.5, power_exponent) - 1) * (max_size - min_size) / norm_factor) * 1.3
    
    # 创建两个图例并放置在图表右侧
    legend1 = ax.legend(handles=legend_elements, loc='upper left', 
                       bbox_to_anchor=(1.02, 0.8), fontsize=8, title='Gene Type')
    ax.add_artist(legend1)
    
    # 创建p值图例 - 更清晰的标题和标签
    size_handles = [plt.scatter([], [], s=s, color='black', alpha=0.7) for s in sizes]
    legend2 = ax.legend(size_handles, p_labels, 
                       loc='upper left', bbox_to_anchor=(1.02, 0.6), 
                       fontsize=8, title='Significance (-log10p)', scatterpoints=1, labelspacing=2.5, borderpad=1.2)
    
    # 9. 调整布局并保存
    output_dir = f"survival_plots_{strategy}"
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    timestamp = int(time.time())
    filename = os.path.join(output_dir, f'{group_name}_waterfall_plot_{timestamp}.pdf')
    fig.savefig(filename, dpi=300, format='pdf') # 暂时移除 bbox_inches='tight'
    print(f"✓ 已保存瀑布图: {filename}")
    plt.close(fig)

def main():
    """
    主函数：加载数据并调用绘图函数
    """
    # 切换到脚本所在目录，方便文件读取
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    
    strategy = 'median_rank'
    
    # 定义要处理的文件和组名
    files_to_process = {
        'hi_hi': f'hi_hi_cox_results_{strategy}.txt',
        'lo_lo': f'lo_lo_cox_results_{strategy}.txt'
    }
    
    print("=" * 60)
    print("开始绘制瀑布图...")
    
    for group, filename in files_to_process.items():
        if os.path.exists(filename):
            print(f"\n正在加载文件: {filename}")
            results_df = pd.read_csv(filename, sep='\t')
            
            # 兼容处理：如果'Gene'列不存在，则假定基因名在第一列
            if 'Gene' not in results_df.columns:
                print("  检测到旧版结果文件格式，正在转换...")
                results_df.rename(columns={results_df.columns[0]: 'Gene'}, inplace=True)

            create_waterfall_plot(results_df, group, strategy)
        else:
            print(f"❌ 警告: 未找到结果文件 {filename}，跳过该组。")
            print("  请先运行 '单基因生存分析.py' 以生成结果文件。")
            
    print("\n" + "=" * 60)
    print("所有绘图任务完成！")
    print("=" * 60)

if __name__ == "__main__":
    main() 