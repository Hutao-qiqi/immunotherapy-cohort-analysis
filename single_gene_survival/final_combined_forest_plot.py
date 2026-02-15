#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.legend_handler import HandlerLine2D
from matplotlib.lines import Line2D
import matplotlib.ticker as ticker # 导入 ticker 模块
import warnings
import time
import argparse
from pathlib import Path

warnings.filterwarnings('ignore')

# --- user config (defaults) ---
DEFAULT_DATA_DIR = Path(__file__).resolve().parent
DATA_DIR = DEFAULT_DATA_DIR
OUTPUT_DIR = DEFAULT_DATA_DIR / "outputs"

HIHI_FILE = 'hi_hi_cox_results_median_rank.txt'
LOLO_FILE = 'lo_lo_cox_results_median_rank.txt'
INTERACTION_FILE = 'cox_interaction_hiHiRisk_loLoNoRisk.tsv'

# 指定需要展示的基因列表
SPECIFIC_GENE_LIST = [
    # Interaction-significant (q<0.05) and immune / immunotherapy-adjacent candidates
    # Selected from cox_interaction_hiHiRisk_loLoNoRisk.tsv (Top8 by q_interaction)
    'CDHR2', 'PROZ', 'SERPINA1', 'VTN', 'CLEC1B', 'CEACAM5', 'HRG', 'FCN3'
]
# --- 配置区结束 ---


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Combined forest plot (hi_hi vs lo_lo + interaction)")
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=DEFAULT_DATA_DIR,
        help="Directory containing input result files",
    )
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=DEFAULT_DATA_DIR / "outputs",
        help="Output directory (default: ./outputs)",
    )
    return parser.parse_args()

# 尝试应用SciencePlots样式
try:
    import scienceplots
    plt.style.use(['science', 'no-latex'])
    print("✓ 成功应用SciencePlots 'science'风格")
except ImportError:
    print("⚠ SciencePlots包未安装，将使用默认Matplotlib样式。")

# 设置全局字体和DPI
plt.rcParams['font.family'] = 'serif'
plt.rcParams['font.serif'] = ['Times New Roman']
plt.rcParams['axes.unicode_minus'] = False
plt.rcParams['figure.dpi'] = 300
plt.rcParams['savefig.dpi'] = 300

# --- 自定义图例处理器 ---
class HandlerErrorBar(HandlerLine2D):
    """
    自定义图例处理器，用于创建“方块+误差棒”样式的图例。
    """
    def create_artists(self, legend, orig_handle,
                       xdescent, ydescent, width, height, fontsize, trans):
        # 绘制误差棒（横线）
        line = Line2D([xdescent, xdescent + width],
                      [ydescent + height / 2, ydescent + height / 2],
                      color=orig_handle.get_color(), lw=orig_handle.get_linewidth())
        
        # 绘制HR点（方块）
        size = height * 0.8
        marker_x = xdescent + 0.5 * width
        marker_y = ydescent + 0.5 * height
        
        square = plt.Rectangle([marker_x - size / 2, marker_y - size / 2], size, size,
                               facecolor=orig_handle.get_color(),
                               transform=trans,
                               edgecolor='white', # 添加白色描边，使其更清晰
                               lw=0.5)

        return [line, square]

def load_and_process_data():
    """数据加载与处理"""
    print("\n" + "=" * 70 + "\n1. 数据加载与处理...\n" + "=" * 70)
    try:
        hihi_data = pd.read_csv(DATA_DIR / HIHI_FILE, sep='\t', index_col=0)
        lolo_data = pd.read_csv(DATA_DIR / LOLO_FILE, sep='\t', index_col=0)
        interaction = pd.read_csv(DATA_DIR / INTERACTION_FILE, sep='\t')
    except FileNotFoundError as e:
        print(f"✗ 错误：找不到数据文件 {e.filename}。请检查文件名和路径。")
        return None, None, None

    if "Gene" not in interaction.columns:
        print("✗ 错误：交互项结果文件缺少 'Gene' 列。")
        return None, None, None

    interaction = interaction.set_index("Gene", drop=False)

    common_genes = [
        g
        for g in SPECIFIC_GENE_LIST
        if (g in hihi_data.index) and (g in lolo_data.index) and (g in interaction.index)
    ]
    print(f"   - 找到 {len(common_genes)} 个共同基因用于分析。")
    if not common_genes: return None, None, None

    hihi_final = hihi_data.loc[common_genes]
    lolo_final = lolo_data.loc[common_genes]
    
    vi_df = interaction.loc[common_genes].copy()
    if "beta_interaction" not in vi_df.columns or "p_interaction" not in vi_df.columns:
        print("✗ 错误：交互项结果文件缺少 beta_interaction / p_interaction 列。")
        return None, None, None

    vi_df["VI"] = np.exp(vi_df["beta_interaction"].astype(float))
    if "q_interaction" in vi_df.columns:
        vi_df = vi_df.sort_values(["q_interaction", "p_interaction"], ascending=[True, True])
    else:
        vi_df = vi_df.sort_values(["p_interaction"], ascending=[True])

    ordered_genes = vi_df["Gene"].tolist()
    
    final_data_parts = []
    for gene in ordered_genes:
        for group_df, group_name in [(lolo_final, 'lolo'), (hihi_final, 'hihi')]:
            row = group_df.loc[[gene]].copy()
            row['Gene'] = gene
            row['Group'] = group_name
            row['p_interaction'] = float(vi_df.loc[gene, 'p_interaction'])
            row['q_interaction'] = float(vi_df.loc[gene, 'q_interaction']) if 'q_interaction' in vi_df.columns else np.nan
            row['VI'] = float(vi_df.loc[gene, 'VI'])
            final_data_parts.append(row)
    
    combined_data = pd.concat(final_data_parts).reset_index(drop=True)
    return combined_data, ordered_genes, vi_df[["Gene", "VI", "p_interaction"] + (["q_interaction"] if "q_interaction" in vi_df.columns else [])]


def create_combined_forest_plot(data, gene_list, vi_df):
    """创建整合了森林图和HRR柱状图的四栏式组合图"""
    print("\n" + "=" * 70 + "\n2. 创建整合式森林图...\n" + "=" * 70)
    
    # 根据基因数量动态计算高度
    fig_height = len(gene_list) * 0.9
    # 设置宽度为高度的3倍，以自动适应3:1的长宽比
    fig_width = fig_height * 3.7
    
    fig, axes = plt.subplots(
        ncols=4, figsize=(fig_width, fig_height),
        # 将wspace设置为0，彻底移除列之间的白色间隙
        gridspec_kw={'width_ratios': [0.7, 2.5, 1.5, 1.2], 'wspace': 0}
    )
    ax1, ax2, ax3, ax4 = axes
    
    total_rows = len(data)
    y_positions = np.arange(total_rows)[::-1]
    
    # 使用标准、柔和的灰色作为阴影颜色
    colors = {'hihi': '#D95F02', 'lolo': '#377eb8', 'shading': '#F0F0F0'}
    
    # --- 通用设置 ---
    # 彻底重写阴影逻辑，确保每个交替基因的两行都被完整覆盖
    print("   - 应用新的 '有阴影-无阴影' 交替样式...")
    for i, gene in enumerate(gene_list):
        if i % 2 == 0:  # 为索引为偶数的基因行 (第1, 3, 5...个) 添加阴影
            # 计算能覆盖两行的准确Y轴范围
            y_min = total_rows - (2 * i) - 2.5 
            y_max = total_rows - (2 * i) - 0.5
            for ax in axes:
                ax.axhspan(y_min, y_max, color=colors['shading'], zorder=0, edgecolor='none')
    
    separator_positions = [total_rows - (j * 2) - 0.5 for j in range(1, len(gene_list))]
    for sep_pos in separator_positions:
        for ax in axes:
            ax.axhline(y=sep_pos, color='#B0B0B0', linestyle='-', linewidth=0.7, alpha=0.8)

    # --- 列1: 基因名 ---
    ax1.axis('off')
    for i, gene in enumerate(gene_list):
        ax1.text(0.1, total_rows - (2 * i) - 1.5, gene, fontsize=14, fontweight='bold', ha='left', va='center')
    
    # --- 列2: 森林图 ---
    ax2.set_xscale('log')
    ax2.axvline(x=1, color='black', linestyle='--', alpha=0.8, linewidth=1.2)
    for i, (_, row) in enumerate(data.iterrows()):
        ax2.plot([row['CI_lower'], row['CI_upper']], [y_positions[i], y_positions[i]], 
                 color=colors[row['Group']], linewidth=2, zorder=3, solid_capstyle='butt')
        ax2.scatter(row['HR'], y_positions[i], color=colors[row['Group']], marker='s', s=50, zorder=5, 
                    edgecolor='white', linewidth=1)
        
        # 添加交互项P值显著性星号（每个基因仅显示一次）
        p_val = row['p_interaction']
        star = ''
        if p_val < 0.001:
            star = '***'
        elif p_val < 0.01:
            star = '**'
        elif p_val < 0.05:
            star = '*'
        
        if star and row['Group'] == 'hihi':
            # 将星号放在HR点的右上方
            ax2.text(row['HR'], y_positions[i] + 0.2, star, 
                     ha='center', va='bottom', color=colors[row['Group']], 
                     fontsize=12, fontweight='bold')

    ax2.set_xlabel('Hazard Ratio (HR)', fontsize=14, fontweight='bold')
    
    # 修改X轴刻度和范围（避免新基因超出范围被截断）
    x_ticks = [0.6, 1, 2, 3]
    ax2.set_xticks(x_ticks)
    ax2.set_xticklabels([str(x) for x in x_ticks], fontsize=12)
    ax2.xaxis.set_minor_locator(ticker.NullLocator())

    x_left = 0.6
    x_right = max(3.0, float(np.nanmax(data['CI_upper'])) * 1.05)
    ax2.set_xlim(x_left, x_right)

    # --- 列3: HR和P值文本 ---
    ax3.axis('off')
    p_formatter = lambda p: "< 0.001" if p < 0.001 else f"{p:.3f}"
    for i, (_, row) in enumerate(data.iterrows()):
        hr_ci_text = f"{row['HR']:.2f} ({row['CI_lower']:.2f}, {row['CI_upper']:.2f})"
        p_text = p_formatter(row['p_interaction']) if row['Group'] == 'hihi' else ""
        ax3.text(0.05, y_positions[i], hr_ci_text, fontsize=12, ha='left', va='center')
        ax3.text(0.95, y_positions[i], p_text, fontsize=12, ha='right', va='center', 
                 fontweight='bold' if (row['Group'] == 'hihi' and row['p_interaction'] < 0.05) else 'normal')

    # --- 列4: Vulnerability Index (from interaction) 柱状图 ---
    vi_plot_df = vi_df.set_index('Gene').loc[gene_list].reset_index()
    vi_values = vi_plot_df['VI']
    vi_y_pos = [total_rows - (2 * i) - 1.5 for i in range(len(gene_list))]
    
    norm = plt.Normalize(vi_values.min() * 0.9, vi_values.max() * 1.1)
    cmap = plt.get_cmap('OrRd')
    ax4.barh(vi_y_pos, vi_values, color=cmap(norm(vi_values)), height=1.6)
    
    for i, vi in enumerate(vi_values):
        ax4.text(vi + 0.03, vi_y_pos[i], f'{vi:.2f}', va='center', ha='left', fontsize=11, fontweight='medium')
    ax4.set_xlabel('Vulnerability Index', fontsize=14, fontweight='bold')
    
    # --- 统一Y轴和边框 ---
    for ax in axes:
        ax.set_ylim(-0.8, total_rows) # 增大顶部空间，为列名和顶线留出间距
        ax.tick_params(axis='y', which='both', length=0)
        for spine_pos in ['top', 'right', 'left']: ax.spines[spine_pos].set_visible(False)
    for ax in [ax1, ax3]: ax.spines['bottom'].set_visible(False)
    for ax in [ax2, ax4]: ax.spines['bottom'].set_linewidth(1.5)
    
    ax2.tick_params(axis='y', labelleft=False)
    ax4.tick_params(axis='y', labelleft=False)
    
    # --- 添加表头和整体线条 ---
    fig.canvas.draw()
    
    # 将列名放置在黑色顶线的上方
    # 使用 transAxes 坐标系 (y=1.02) 将标题放置在子图区域的正上方
    print("   - 将列名移动到顶线上方...")
    ax1.text(0.1, 1.02, 'Gene', transform=ax1.transAxes, fontsize=16, fontweight='bold', ha='left', va='bottom')
    ax3.text(0.05, 1.02, 'HR (95% CI)', transform=ax3.transAxes, fontsize=14, fontweight='bold', ha='left', va='bottom')
    ax3.text(0.95, 1.02, 'P (interaction)', transform=ax3.transAxes, fontsize=14, fontweight='bold', ha='right', va='bottom')

    ax1_pos, ax4_pos = ax1.get_position(), ax4.get_position()
    line_x_start, line_x_end = ax1_pos.x0, ax4_pos.x1
    # 顶线作为标题和内容的分隔
    fig.add_artist(plt.Line2D([line_x_start, line_x_end], [ax1_pos.y1, ax1_pos.y1], color='black', linewidth=2))
    # 底线
    fig.add_artist(plt.Line2D([line_x_start, line_x_end], [ax1_pos.y0, ax1_pos.y0], color='black', linewidth=1.5))
    
    # --- 添加自定义图例 ---
    legend_elements = [
        Line2D([0], [0], color=colors['hihi'], lw=2, label='MYC_PVT1High_expression group'),
        Line2D([0], [0], color=colors['lolo'], lw=2, label='MYC_PVT1Low_expression group')
    ]
    
    fig.legend(
        handles=legend_elements,
        handler_map={Line2D: HandlerErrorBar()}, # 应用自定义处理器
        loc='upper right',
        bbox_to_anchor=(line_x_end - 0.01, ax1_pos.y1 + 0.12),
        frameon=False,
        fontsize=12,
        handlelength=2.0 # 调整图例中图形的长度
    )
    
    return fig

def main():
    """主函数"""
    global DATA_DIR, OUTPUT_DIR
    args = parse_args()
    DATA_DIR = args.data_dir.resolve()
    OUTPUT_DIR = args.out_dir.resolve()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    print(f"✓ Using data dir: {DATA_DIR}")
    print(f"✓ Using output dir: {OUTPUT_DIR}")

    data, gene_list, hrr_df = load_and_process_data()
    if data is None:
        print("✗ 数据处理失败，程序终止。")
        return

    fig = create_combined_forest_plot(data, gene_list, hrr_df)
    
    timestamp = int(time.time())
    output_base = f'Final_Combined_Forest_Plot_{timestamp}'
    plt.savefig(OUTPUT_DIR / f'{output_base}.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(OUTPUT_DIR / f'{output_base}.pdf', bbox_inches='tight', facecolor='white')
    
    print("\n" + "=" * 70 + "\n✓ 最终版整合式森林图生成成功！\n" + "=" * 70)
    print(f"   - 文件已保存为 {output_base}.png 和 {output_base}.pdf")
    
    plt.show()

if __name__ == "__main__":
    main()