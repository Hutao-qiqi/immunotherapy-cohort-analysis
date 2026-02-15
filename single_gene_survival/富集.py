import pandas as pd
import numpy as np
import gseapy as gp
from gseapy.plot import barplot, dotplot
from mygene import MyGeneInfo
import matplotlib.pyplot as plt
import seaborn as sns
import os
from gseapy import Msigdb
import re

# 设置工作目录
os.chdir("E:/data/changyuan/免疫队列/单基因生存分析")

# 检查可用的基因集
def check_available_gene_sets():
    """检查当前可用的基因集"""
    try:
        from gseapy import get_library_name
        available_libs = get_library_name()
        print("可用的基因集库：")
        # 筛选出可能包含免疫相关信息的基因集
        immune_related = [lib for lib in available_libs if any(keyword in lib.lower() 
                         for keyword in ['immun', 'hallmark', 'reactome', 'kegg', 'go', 'biocarta', 'msig'])]
        for i, lib in enumerate(immune_related[:20]):  # 只显示前20个
            print(f"{i+1}. {lib}")
        return immune_related
    except Exception as e:
        print(f"无法获取基因集列表: {e}")
        return []

# 检查可用基因集
available_gene_sets = check_available_gene_sets()

# 1. 读取基因列表 (使用您提供的列表)
gene_symbols = [
    'MIR548I1', 'SLC30A5', 'FOLH1', 'DNAJB8', 'SLC5A12', 'RPS11', 'RPL23',
    'MTRNR2L3', 'TES', 'RPL22L1', 'OR56B4', 'FIS1', 'FAM81B', 'HORMAD2',
    'OLA1', 'TMEM92', 'RPS4X', 'NXPH1', 'MRGPRX4', 'MCMDC2', 'MAMSTR',
    'ZMYND8', 'C10orf126', 'OR6V1', 'GFM2', 'GALNT4', 'REG3G', 'PSG11',
    'RPS28', 'COMMD6', 'SLC6A10P', 'SIX6', 'IMMP2L', 'COL20A1', 'SPACA5',
    'OR4F4', 'ANKRD62P1-PARP4P3', 'GNRHR', 'C12orf56', 'F2RL1', 'EIF3E',
    'SFTPC', 'PTCD1', 'MIR941-2', 'GLIPR1L1', 'COX6A2', 'TPT1', 'MDFI',
    'ARG1', 'IL13RA2', 'NAA11', 'RPL13A', 'RPS3', 'MYADML2', 'POLR1D',
    'C7orf33', 'OR2C3', 'SYNPR', 'RSPH3', 'MAGED4B', 'GSTO2', 'RPL18',
    'RPL13AP5', 'FLG2', 'VWA3A', 'USP17L15', 'FSIP2', 'CA6', 'FCHSD1',
    'TRIM51HP', 'LPIN3', 'CFHR4', 'POU5F1B', 'PLA2G4C-AS1', 'KCTD14',
    'STARD13-AS', 'ATP13A4-AS1', 'ZHX1', 'RPL31', 'SOX9', 'PRICKLE2-AS3',
    'POLM', 'RNVU1-7', 'LCE1C', 'RPL27', 'HAPLN1', 'UBE2H', 'APOC2',
    'LCE5A', 'RPS19', 'FAM138D', 'LY75-CD302', 'LINC00303', 'CDK6', 'IBTK',
    'OR51I2', 'CKMT1A', 'TMCO2', 'OR52E8', 'RPL31P11', 'TMEM161B',
    'PRSS37', 'ASB4', 'SERPINB5', 'CXorf49B', 'DBX1', 'THEM5', 'DPY19L4',
    'RPL13AP6', 'FBP2', 'NKX2-5', 'ANP32AP1', 'RPAP1', 'THUMPD2',
    'OR52D1', 'OR9G1', 'BAAT', 'CTCFL', 'ZHX1-C8orf76', 'KRT10', 'UQCRB',
    'PRSS48', 'KLF3', 'FAM47E-STBD1', 'IGIP', 'SCARNA13', 'MGST2',
    'NIPAL2', 'ASAH2', 'OR52B6', 'NDUFA5', 'SYPL1', 'EDN3', 'PCDHB14',
    'OR2M2', 'FAM217A', 'FAM135B', 'RGS21', 'SFT2D1', 'TCEAL8', 'SSX1',
    'RAG2', 'OR3A3', 'KIAA0087', 'CTAGE8', 'OR2AT4', 'OR1C1', 'LPAR6',
    'SEZ6', 'ZNF561', 'GALNT14', 'ST13P4', 'OR5J2', 'TNNI3K', 'DIS3',
    'KCTD4', 'HHATL', 'UBR5', 'LMBR1', 'DLEU7', 'LRRTM1', 'AFF2', 'CCDC54',
    'ADAM21P1', 'NACA2', 'GUSBP3', 'FGF5', 'RPS5', 'CHRNA7', 'LRPPRC'
]

print(f"\n总基因数: {len(gene_symbols)}")
print("前10个基因符号示例:", gene_symbols[:10])

# 2. 基因符号转换为Entrez ID
mg = MyGeneInfo()
gene_info = mg.querymany(gene_symbols, scopes='symbol', fields='entrezgene', species='human')

# 创建两个列表分别存储转换成功和失败的基因
converted_genes = []
failed_genes = []

for g in gene_info:
    if 'entrezgene' in g:
        converted_genes.append({'symbol': g['query'], 'entrez': str(g['entrezgene'])})
    else:
        failed_genes.append(g['query'])

print(f"\n成功转换的基因数: {len(converted_genes)}")
print(f"转换失败的基因数: {len(failed_genes)}")
print("\n转换失败的基因符号:", failed_genes[:10], "..." if len(failed_genes) > 10 else "")

# 获取转换成功的Entrez ID列表
gene_entrez = [g['entrez'] for g in converted_genes]
symbol_to_entrez = {g['symbol']: g['entrez'] for g in converted_genes}

if len(gene_entrez) == 0:
    print("\n错误：没有成功转换的基因，无法进行富集分析")
    exit()

# 3. GO富集分析
print("\n开始GO富集分析...")
go_bp = gp.enrichr(gene_list=gene_entrez,
                   gene_sets=['GO_Biological_Process_2023'],
                   organism='human',
                   outdir=None,
                   cutoff=0.1)

# 4. KEGG富集分析
print("\n开始KEGG富集分析...")
kegg = gp.enrichr(gene_list=gene_entrez,
                  gene_sets=['KEGG_2021_Human'],
                  organism='human',
                  outdir=None,
                  cutoff=0.1)

# 5. 免疫治疗相关基因集分析
print("\n开始免疫治疗相关基因集分析...")

# 从可用基因集中选择免疫相关的
immune_candidates = []
for gene_set in available_gene_sets:
    if any(keyword in gene_set.lower() for keyword in ['hallmark', 'immun', 'reactome']):
        immune_candidates.append(gene_set)

print(f"找到 {len(immune_candidates)} 个可能的免疫相关基因集")

# 尝试使用找到的免疫相关基因集
immune_enrich = None
successful_geneset = None

for geneset in immune_candidates[:3]:  # 尝试前3个
    try:
        print(f"尝试使用基因集: {geneset}")
        immune_enrich = gp.enrichr(gene_list=gene_entrez,
                                  gene_sets=[geneset],
                                  organism='human',
                            outdir=None,
                            cutoff=0.2)
        successful_geneset = geneset
        print(f"成功使用基因集: {geneset}")
        break
    except Exception as e:
        print(f"基因集 {geneset} 分析失败: {e}")
        continue

# 如果所有免疫基因集都失败，跳过免疫分析
if immune_enrich is None:
    print("所有免疫基因集分析都失败，跳过免疫分析部分")
    # 创建一个空的结果对象
    class EmptyResult:
        def __init__(self):
            self.results = None
    immune_enrich = EmptyResult()
    successful_geneset = "None"

# 6. 结果可视化 (修正缩进)
# GO富集结果
if go_bp.results is not None and not go_bp.results.empty:
    dotplot(go_bp.results,
            title='GO Biological Process Enrichment',
            cmap='viridis',
            size=20,
            top_term=20,
            figsize=(12, 8))
    plt.tight_layout()
    plt.savefig('GO_enrichment.png', dpi=300, bbox_inches='tight')
    plt.close()

# KEGG通路富集
if kegg.results is not None and not kegg.results.empty:
    dotplot(kegg.results,
            title='KEGG Pathway Enrichment',
            cmap='RdYlBu_r',
            size=20,
            top_term=20,
            figsize=(12, 8))
    plt.tight_layout()
    plt.savefig('KEGG_enrichment.png', dpi=300, bbox_inches='tight')
    plt.close()

# 免疫治疗相关基因集富集
if immune_enrich.results is not None and not immune_enrich.results.empty:
    dotplot(immune_enrich.results,
            title='Immunologic Signatures Enrichment',
            cmap='plasma',
            size=20,
            top_term=20,
            figsize=(12, 8))
    plt.tight_layout()
    plt.savefig('Immune_enrichment.png', dpi=300, bbox_inches='tight')
    plt.close()

# 保存结果到CSV文件
print("\n保存富集分析结果...")
if go_bp.results is not None:
    go_bp.results.to_csv("GO_enrichment_results.csv", index=False)
    print(f"GO富集分析结果已保存 (共{len(go_bp.results)}个条目)")

if kegg.results is not None:
    kegg.results.to_csv("KEGG_enrichment_results.csv", index=False)
    print(f"KEGG富集分析结果已保存 (共{len(kegg.results)}个条目)")

if immune_enrich.results is not None:
    immune_enrich.results.to_csv("Immune_enrichment_results.csv", index=False)
    print(f"免疫富集分析结果已保存 (共{len(immune_enrich.results)}个条目)")
    print(f"使用的基因集: {successful_geneset}")
else:
    print("免疫富集分析无结果")

print("\n富集分析完成！")
print(f"总共分析了 {len(gene_entrez)} 个基因")
print(f"其中 {len(converted_genes)} 个基因成功转换为Entrez ID")

# 7. 识别关键免疫治疗基因
print("\n分析免疫治疗相关基因...")

# 从免疫富集结果中提取显著富集的基因
immune_therapy_genes = []
if immune_enrich.results is not None and not immune_enrich.results.empty:
    # 获取显著富集的通路（p < 0.05）
    significant_pathways = immune_enrich.results[immune_enrich.results['Adjusted P-value'] < 0.05]
    
    # 从这些通路中提取基因
    for _, row in significant_pathways.iterrows():
        pathway_name = row['Term']
        genes = row['Genes'].split(';')
    for gene in genes:
            if gene in gene_entrez:  # 确保基因在我们的输入列表中
                # 获取基因符号（如果有）
                symbol = next((g['symbol'] for g in converted_genes if g['entrez'] == gene), gene)
            immune_therapy_genes.append({
                'Gene_Symbol': symbol,
                'Entrez_ID': gene,
                    'Pathway': pathway_name,
                    'Adjusted_P_value': row['Adjusted P-value']
            })

# 转换为DataFrame并去重
if immune_therapy_genes:
    immune_summary = pd.DataFrame(immune_therapy_genes).drop_duplicates()
    # 按基因分组，合并通路信息
    immune_summary = immune_summary.groupby(['Gene_Symbol', 'Entrez_ID']).agg({
        'Pathway': lambda x: '; '.join(set(x)),
        'Adjusted_P_value': 'min'  # 取最显著的p值
    }).reset_index()
    # 按p值排序
    immune_summary = immune_summary.sort_values('Adjusted_P_value')
else:
    immune_summary = pd.DataFrame(columns=['Gene_Symbol', 'Entrez_ID', 'Pathway', 'Adjusted_P_value'])

# 8. 输出结果
print("\n=== 潜在的免疫治疗相关基因 ===")
if not immune_summary.empty:
    print(f"找到 {len(immune_summary)} 个显著富集的免疫相关基因")
    print("\n前10个最显著的基因：")
    print(immune_summary.head(10))
else:
    print("未找到显著富集的免疫相关基因")

# 保存免疫治疗基因结果
immune_summary.to_csv("immunotherapy_genes.csv", index=False)
print("\n免疫治疗基因结果已保存到 immunotherapy_genes.csv")

# 9. 高级可视化 - 基因-通路网络 (修正缩进)
if not immune_summary.empty:
    print("\n创建基因-通路网络图...")
    # 创建网络数据
    network_data = []
    for _, row in immune_summary.iterrows():
        pathways = row['Pathway'].split('; ')
        for pathway in pathways:
            network_data.append({
                'Gene': row['Gene_Symbol'],
                'Pathway': pathway,
                'P_value': row['Adjusted_P_value']
            })

    network_df = pd.DataFrame(network_data)

    # 绘制热图
    plt.figure(figsize=(15, 10))

    # 创建交叉表
    cross_tab = pd.crosstab(network_df['Gene'], network_df['Pathway'])

    # 使用p值作为颜色深度
    p_value_matrix = pd.pivot_table(
        network_df,
        values='P_value',
        index='Gene',
        columns='Pathway',
        aggfunc='min'  # 使用最小p值
    ).fillna(1)  # 填充缺失值为1

    # 绘制热图
    sns.heatmap(
        -np.log10(p_value_matrix),  # 转换p值为-log10(p)
        cmap="YlOrRd",
        cbar_kws={'label': '-log10(Adjusted P-value)'},
        xticklabels=True,
        yticklabels=True
    )

    plt.title("Gene-Pathway Network in Cancer Immunotherapy")
    plt.xticks(rotation=45, ha='right')
    plt.tight_layout()
    plt.savefig('gene_pathway_network.png', dpi=300, bbox_inches='tight')
    plt.close()
    print("网络图已保存为 gene_pathway_network.png")
else:
    print("\n由于没有显著富集的基因，跳过网络图绘制")
