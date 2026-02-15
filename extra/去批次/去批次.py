import pandas as pd
import numpy as np
from combat.pycombat import pycombat
from sklearn.decomposition import PCA
from sklearn.cluster import AgglomerativeClustering
from sklearn.metrics import silhouette_score, adjusted_rand_score, normalized_mutual_info_score
from sklearn.metrics import calinski_harabasz_score, davies_bouldin_score
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.patches import Ellipse
from scipy.stats import chi2
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.spatial.distance import pdist, squareform
import warnings
from sklearn.preprocessing import StandardScaler
warnings.filterwarnings('ignore')

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

def confidence_ellipse(x, y, ax, n_std=1.0, facecolor='none', edgecolor='red', linewidth=2, alpha=0.5):
    """
    绘制置信椭圆
    
    Parameters:
    -----------
    x, y : array-like, shape (n, )
        数据点的x, y坐标
    ax : matplotlib.axes.Axes
        绘图的axes对象
    n_std : float
        标准差的倍数，用于确定椭圆的大小（默认1.0，原来是2.0）
    """
    if x.size != y.size:
        raise ValueError("x and y must be the same size")
    
    # 计算协方差矩阵
    cov = np.cov(x, y)
    
    # 计算椭圆的参数
    eigenvals, eigenvecs = np.linalg.eigh(cov)
    order = eigenvals.argsort()[::-1]
    eigenvals, eigenvecs = eigenvals[order], eigenvecs[:, order]
    
    # 计算椭圆的角度
    angle = np.degrees(np.arctan2(*eigenvecs[:, 0][::-1]))
    
    # 计算椭圆的宽度和高度
    width, height = 2 * n_std * np.sqrt(eigenvals)
    
    # 绘制椭圆
    ellip = Ellipse(xy=(np.mean(x), np.mean(y)), width=width, height=height, angle=angle,
                    facecolor=facecolor, edgecolor=edgecolor, linewidth=linewidth, alpha=alpha)
    
    return ax.add_patch(ellip)

def plot_pca(data, batch_info, title, ax):
    """
    绘制发表级别的PCA散点图
    
    Parameters:
    -----------
    data: 表达矩阵 (基因 x 样本)
    batch_info: 包含批次信息的Series
    title: 图表标题
    ax: matplotlib的Axes对象
    """
    # Sklearn的PCA期望样本作为行，因此需要转置数据
    pca = PCA(n_components=2)
    principal_components = pca.fit_transform(data.T)

    # 创建一个用于绘图的DataFrame
    pca_df = pd.DataFrame(data=principal_components,
                          columns=['Dim1', 'Dim2'],
                          index=data.columns)
    pca_df['Group'] = batch_info.values

    # 获取每个主成分解释的方差百分比
    explained_variance = pca.explained_variance_ratio_ * 100

    # Nature期刊常用配色
    nature_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f']
    
    # 专用于MYC/PVT1分组的配色
    myc_pvt1_colors = {
        'hi_hi': '#d62728',  # 红色 - 高表达
        'hi_lo': '#ff7f0e',  # 橙色
        'lo_hi': '#2ca02c',  # 绿色
        'lo_lo': '#1f77b4'   # 蓝色 - 低表达
    }
    
    unique_groups = pca_df['Group'].unique()
    
    # 判断是否为MYC/PVT1分组
    is_myc_pvt1 = set(unique_groups).issubset({'hi_hi', 'hi_lo', 'lo_hi', 'lo_lo'})
    
    # 选择配色方案
    if is_myc_pvt1:
        # 如果是MYC/PVT1分组，使用专用配色
        colors = [myc_pvt1_colors.get(group, nature_colors[i]) for i, group in enumerate(unique_groups)]
    else:
        colors = [nature_colors[i % len(nature_colors)] for i in range(len(unique_groups))]

    # 批次名称映射
    batch_name_mapping = {
        'IMvigor210': 'Mariathasan et al.',
        'liu': 'Liu et al.',
        'ra': 'Riaz et al.',
        'hugo': 'Hugo et al.'
    }
    
    # 绘制散点图
    for i, group in enumerate(unique_groups):
        group_data = pca_df[pca_df['Group'] == group]
        # 使用映射后的批次名称
        display_name = batch_name_mapping.get(group, group)
        ax.scatter(group_data['Dim1'], group_data['Dim2'],
                  c=colors[i],
                  s=60,
                  alpha=0.8,
                  label=display_name,
                  edgecolors='white',
                  linewidth=0.5)
        
        # 只为cohort分组添加置信椭圆
        if not is_myc_pvt1 and len(group_data) >= 3:
            try:
                confidence_ellipse(
                    group_data['Dim1'], group_data['Dim2'],
                    ax, n_std=2.0,  # 2个标准差
                    edgecolor=colors[i],
                    facecolor=colors[i],
                    alpha=0.15,
                    linewidth=2
                )
            except Exception as e:
                print(f"跳过组 {group} 的椭圆绘制: {str(e)}")

    # 设置字体和样式
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['pdf.fonttype'] = 42
    
    # 设置标题和轴标签
    ax.set_title(title, fontsize=16, fontweight='bold', pad=20)
    ax.set_xlabel(f'PC1 ({explained_variance[0]:.1f}%)', fontsize=14, fontweight='bold')
    ax.set_ylabel(f'PC2 ({explained_variance[1]:.1f}%)', fontsize=14, fontweight='bold')
    
    # 显示完整的矩形边框
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(1.5)
    
    # 设置刻度标签字体大小
    ax.tick_params(axis='both', labelsize=12)
    
    # 添加网格线
    ax.grid(True, alpha=0.3)
    
    # 优化图例
    legend = ax.legend(title='Group',
                      frameon=True,
                      fancybox=True,
                      framealpha=1,
                      shadow=True,
                      bbox_to_anchor=(1.05, 1),
                      loc='upper left',
                      fontsize=16)  # 原来是12，调整为1.3倍
    legend.get_title().set_fontsize(17)  # 原来是13，调整为1.3倍
    legend.get_title().set_fontweight('bold')
    
    # 设置背景为白色
    ax.set_facecolor('white')
    
    # 调整布局以适应图例
    plt.tight_layout()

def plot_tsne(data, batch_info, title, ax, perplexity=30, random_state=42):
    """
    使用t-SNE进行发表级别的降维可视化
    
    Parameters:
    -----------
    data: 表达矩阵 (基因 x 样本)
    batch_info: 包含批次信息的Series
    title: 图表标题
    ax: matplotlib的Axes对象
    """
    print(f"    正在运行t-SNE降维: {title}")
    
    # t-SNE降维
    tsne = TSNE(n_components=2, perplexity=perplexity, random_state=random_state, 
                n_iter=1000, verbose=0)
    embedding = tsne.fit_transform(data.T)
    
    # 创建用于绘图的DataFrame
    tsne_df = pd.DataFrame(data=embedding, columns=['tSNE1', 'tSNE2'], index=data.columns)
    tsne_df['Group'] = batch_info.values
    
    # Nature期刊常用配色
    nature_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f']
    
    # 专用于MYC/PVT1分组的配色
    myc_pvt1_colors = {
        'hi_hi': '#d62728',  # 红色 - 高表达
        'hi_lo': '#ff7f0e',  # 橙色
        'lo_hi': '#2ca02c',  # 绿色
        'lo_lo': '#1f77b4'   # 蓝色 - 低表达
    }
    
    unique_groups = tsne_df['Group'].unique()
    
    # 判断是否为MYC/PVT1分组
    is_myc_pvt1 = set(unique_groups).issubset({'hi_hi', 'hi_lo', 'lo_hi', 'lo_lo'})
    
    # 选择配色方案
    if is_myc_pvt1:
        # 如果是MYC/PVT1分组，使用专用配色
        colors = [myc_pvt1_colors.get(group, nature_colors[i]) for i, group in enumerate(unique_groups)]
    else:
        colors = [nature_colors[i % len(nature_colors)] for i in range(len(unique_groups))]

    # 批次名称映射
    batch_name_mapping = {
        'IMvigor210': 'Mariathasan et al.',
        'liu': 'Liu et al.',
        'ra': 'Riaz et al.',
        'hugo': 'Hugo et al.'
    }
    
    # 绘制散点图
    for i, group in enumerate(unique_groups):
        group_data = tsne_df[tsne_df['Group'] == group]
        # 使用映射后的批次名称
        display_name = batch_name_mapping.get(group, group)
        ax.scatter(group_data['tSNE1'], group_data['tSNE2'],
                  c=colors[i],
                  s=60,
                  alpha=0.8,
                  label=display_name,
                  edgecolors='white',
                  linewidth=0.5)
        
        # 只为cohort分组添加置信椭圆
        if not is_myc_pvt1 and len(group_data) >= 3:
            try:
                confidence_ellipse(
                    group_data['tSNE1'], group_data['tSNE2'],
                    ax, n_std=2.0,  # 2个标准差
                    edgecolor=colors[i],
                    facecolor=colors[i],
                    alpha=0.15,
                    linewidth=2
                )
            except Exception as e:
                print(f"跳过组 {group} 的椭圆绘制: {str(e)}")

    # 设置字体和样式
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['pdf.fonttype'] = 42
    
    # 设置标题和轴标签
    ax.set_title(title, fontsize=16, fontweight='bold', pad=20)
    ax.set_xlabel('t-SNE1', fontsize=14, fontweight='bold')
    ax.set_ylabel('t-SNE2', fontsize=14, fontweight='bold')
    
    # 显示完整的矩形边框
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_linewidth(1.5)
    
    # 设置刻度标签字体大小
    ax.tick_params(axis='both', labelsize=12)
    
    # 添加网格线
    ax.grid(True, alpha=0.3)
    
    # 优化图例
    legend = ax.legend(title='Group',
                      frameon=True,
                      fancybox=True,
                      framealpha=1,
                      shadow=True,
                      bbox_to_anchor=(1.05, 1),
                      loc='upper left',
                      fontsize=16)  # 原来是12，调整为1.3倍
    legend.get_title().set_fontsize(17)  # 原来是13，调整为1.3倍
    legend.get_title().set_fontweight('bold')
    
    # 设置背景为白色
    ax.set_facecolor('white')
    
    # 调整布局以适应图例
    plt.tight_layout()

def evaluate_clustering_metrics(data, batch_labels, title_prefix=""):
    """
    评估聚类质量的多种指标
    
    Parameters:
    -----------
    data : pandas.DataFrame
        表达矩阵 (基因 x 样本)，将被转置用于聚类
    batch_labels : pandas.Series
        批次标签
    title_prefix : str
        标题前缀
    
    Returns:
    --------
    dict : 包含各种聚类评估指标的字典
    """
    
    # 转置数据（样本 x 基因）用于聚类
    X = data.T
    
    # 获取批次标签的数值编码
    batch_labels_numeric = pd.Categorical(batch_labels).codes
    
    # 计算样本间的距离矩阵
    print(f"  {title_prefix}计算样本间距离矩阵...")
    
    # 使用相关系数距离（1-相关系数）
    # 计算样本间相关性
    correlation_matrix = X.T.corr()  # 样本间相关性
    distance_matrix = 1 - correlation_matrix
    
    # 确保距离矩阵非负且对称
    distance_matrix = np.abs(distance_matrix)
    distance_matrix = (distance_matrix + distance_matrix.T) / 2
    
    # 转换为numpy数组
    distance_matrix = distance_matrix.values
    np.fill_diagonal(distance_matrix, 0)
    
    # 层次聚类
    print(f"  {title_prefix}执行层次聚类...")
    linkage_matrix = linkage(squareform(distance_matrix), method='ward')
    
    # 确定聚类数量（使用批次数量）
    n_clusters = len(batch_labels.unique())
    cluster_labels = fcluster(linkage_matrix, n_clusters, criterion='maxclust')
    
    # 计算评估指标
    metrics = {}
    
    try:
        # 硅轮廓系数 (范围: -1 到 1, 越高越好)
        silhouette_avg = silhouette_score(distance_matrix, cluster_labels, metric='precomputed')
        metrics['silhouette_score'] = silhouette_avg
        
        # 调整兰德指数 (范围: -1 到 1, 越高越好)
        ari_score = adjusted_rand_score(batch_labels_numeric, cluster_labels)
        metrics['adjusted_rand_index'] = ari_score
        
        # 标准化互信息 (范围: 0 到 1, 越高越好)
        nmi_score = normalized_mutual_info_score(batch_labels_numeric, cluster_labels)
        metrics['normalized_mutual_info'] = nmi_score
        
        # 卡林斯基-哈拉巴兹指数 (越高越好)
        # 注意：这个指标需要原始特征空间，不能使用距离矩阵
        ch_score = calinski_harabasz_score(X, cluster_labels)
        metrics['calinski_harabasz_score'] = ch_score
        
        # 戴维斯-鲍尔丁指数 (越低越好)
        db_score = davies_bouldin_score(X, cluster_labels)
        metrics['davies_bouldin_score'] = db_score
        
        # 批次混合度评估
        batch_purity = calculate_batch_purity(cluster_labels, batch_labels_numeric)
        metrics['batch_purity'] = batch_purity
        
        # 批次分离度评估
        batch_separation = calculate_batch_separation(cluster_labels, batch_labels_numeric)
        metrics['batch_separation'] = batch_separation
        
    except Exception as e:
        print(f"    警告：计算某些指标时出错: {str(e)}")
    
    return metrics, linkage_matrix, cluster_labels

def calculate_batch_purity(cluster_labels, batch_labels):
    """
    计算批次纯度：衡量每个聚类中是否主要由一个批次组成
    返回值越接近1表示批次效应越强（校正效果差）
    """
    unique_clusters = np.unique(cluster_labels)
    total_samples = len(cluster_labels)
    weighted_purity = 0
    
    for cluster_id in unique_clusters:
        cluster_mask = cluster_labels == cluster_id
        cluster_batches = batch_labels[cluster_mask]
        
        if len(cluster_batches) > 0:
            # 计算该聚类中最主要批次的比例
            batch_counts = np.bincount(cluster_batches)
            max_batch_count = np.max(batch_counts)
            cluster_purity = max_batch_count / len(cluster_batches)
            
            # 按聚类大小加权
            cluster_weight = len(cluster_batches) / total_samples
            weighted_purity += cluster_purity * cluster_weight
    
    return weighted_purity

def calculate_batch_separation(cluster_labels, batch_labels):
    """
    计算批次分离度：衡量不同批次是否被分配到不同聚类
    返回值越接近1表示批次效应越强（校正效果差）
    """
    unique_batches = np.unique(batch_labels)
    total_separation = 0
    
    for batch_id in unique_batches:
        batch_mask = batch_labels == batch_id
        batch_clusters = cluster_labels[batch_mask]
        
        if len(batch_clusters) > 0:
            # 计算该批次样本的聚类多样性
            cluster_counts = np.bincount(batch_clusters)
            cluster_counts = cluster_counts[cluster_counts > 0]
            
            if len(cluster_counts) > 1:
                # 使用熵来衡量聚类多样性
                probabilities = cluster_counts / np.sum(cluster_counts)
                # 避免log(0)的情况
                probabilities = probabilities[probabilities > 0]
                entropy = -np.sum(probabilities * np.log2(probabilities))
                max_entropy = np.log2(len(cluster_counts))
                separation = entropy / max_entropy if max_entropy > 0 else 0
            else:
                separation = 0  # 所有样本都在同一个聚类中
            
            total_separation += separation
    
    return total_separation / len(unique_batches)

def plot_clustering_heatmap(data, batch_labels, linkage_matrix, cluster_labels, title, figsize=(12, 8)):
    """
    绘制聚类热图
    """
    try:
        # 准备数据
        X = data.T  # 样本 x 基因
        
        # 选择前1000个方差最大的基因进行可视化
        gene_vars = data.var(axis=1)
        n_genes = min(1000, len(gene_vars))
        top_genes = gene_vars.nlargest(n_genes).index
        data_subset = data.loc[top_genes]
        
        # 创建注释颜色
        batch_colors = {'IMvigor210': '#1f77b4', 'Liu et al., Nat Medicine 2019': '#ff7f0e', 'Riaz et al., Cell 2017': '#2ca02c', 'Hugo et al., Cell 2016': '#d62728'}
        
        # 创建颜色映射
        batch_color_map = {batch: color for batch, color in batch_colors.items() if batch in batch_labels.unique()}
        
        # 绘制热图
        g = sns.clustermap(data_subset, 
                          col_linkage=linkage_matrix,
                          col_cluster=True,
                          row_cluster=True,
                          col_colors=[batch_color_map.get(batch, '#gray') for batch in batch_labels],
                          cmap='RdBu_r',
                          center=0,
                          figsize=figsize,
                          dendrogram_ratio=0.15,
                          colors_ratio=0.02)
        
        g.fig.suptitle(title, fontsize=14, fontweight='bold')
        
        return g.fig
        
    except Exception as e:
        print(f"绘制聚类热图时出错: {str(e)}")
        # 返回一个简单的图形
        fig, ax = plt.subplots(figsize=figsize)
        ax.text(0.5, 0.5, f'聚类热图生成失败\n{str(e)}', 
                ha='center', va='center', transform=ax.transAxes)
        ax.set_title(title)
        return fig

def plot_clustering_evaluation(metrics_before, metrics_after, figsize=(14, 10)):
    """
    绘制聚类评估指标对比图
    """
    fig, axes = plt.subplots(2, 3, figsize=figsize)
    axes = axes.flatten()
    
    # 定义指标和其理想方向
    metrics_info = {
        'silhouette_score': ('Silhouette Score', 'higher_better'),
        'adjusted_rand_index': ('Adjusted Rand Index', 'higher_better'),
        'normalized_mutual_info': ('Normalized Mutual Info', 'higher_better'),
        'calinski_harabasz_score': ('Calinski-Harabasz Score', 'higher_better'),
        'davies_bouldin_score': ('Davies-Bouldin Score', 'lower_better'),
        'batch_purity': ('Batch Purity', 'lower_better')
    }
    
    for i, (metric_key, (metric_name, direction)) in enumerate(metrics_info.items()):
        if i >= len(axes):
            break
            
        ax = axes[i]
        
        # 获取数值
        before_val = metrics_before.get(metric_key, 0)
        after_val = metrics_after.get(metric_key, 0)
        
        # 绘制条形图
        categories = ['Before ComBat', 'After ComBat']
        values = [before_val, after_val]
        
        bars = ax.bar(categories, values, color=['#ff7f7f', '#7f7fff'], alpha=0.8)
        
        # 添加数值标签
        for bar, val in zip(bars, values):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + height*0.01,
                   f'{val:.3f}', ha='center', va='bottom', fontweight='bold')
        
        # 设置标题和标签
        ax.set_title(metric_name, fontsize=12, fontweight='bold')
        ax.set_ylabel('Score', fontsize=10)
        
        # 根据指标方向设置颜色
        if direction == 'higher_better':
            if after_val > before_val:
                ax.set_facecolor('#e8f5e8')  # 浅绿色背景表示改善
            else:
                ax.set_facecolor('#ffe8e8')  # 浅红色背景表示恶化
        else:  # lower_better
            if after_val < before_val:
                ax.set_facecolor('#e8f5e8')  # 浅绿色背景表示改善
            else:
                ax.set_facecolor('#ffe8e8')  # 浅红色背景表示恶化
        
        # 美化
        ax.grid(True, alpha=0.3)
        ax.set_axisbelow(True)
        
        # 设置y轴范围
        if metric_key in ['silhouette_score', 'adjusted_rand_index']:
            ax.set_ylim(-1, 1)
        elif metric_key in ['normalized_mutual_info', 'batch_purity']:
            ax.set_ylim(0, 1)
    
    plt.tight_layout()
    return fig


# ==============================================================================
# 1. 加载四个过滤后的表达矩阵文件
# ==============================================================================
# 定义文件路径
file_paths = {
    'IMvigor210': 'IMvigor210_filtered1.txt',
    'liu': 'liu_filtered1.txt', 
    'ra': 'ra_filtered1.txt',
    'hugo': 'hugo_filtered1.txt'
}

print("正在加载过滤后的表达矩阵文件...")
datasets = {}
all_samples = []
all_batches = []
sample_mapping = {}  # 存储原始样本名到前缀样本名的映射

# 加载每个文件
for batch_name, file_path in file_paths.items():
    try:
        # 读取表达矩阵（基因为行，样本为列）
        df = pd.read_csv(file_path, sep='\t', index_col=0)
        
        # 为样本名添加前缀以区分不同数据集
        original_samples = df.columns.tolist()
        prefixed_samples = [f"{batch_name}_{sample}" for sample in original_samples]
        df.columns = prefixed_samples
        
        # 记录原始样本名到前缀样本名的映射
        for orig, prefixed in zip(original_samples, prefixed_samples):
            sample_mapping[orig] = prefixed
        
        datasets[batch_name] = df
        
        # 收集样本信息
        all_samples.extend(prefixed_samples)
        all_batches.extend([batch_name] * len(prefixed_samples))
        
        print(f"  {batch_name}: {df.shape[0]}个基因, {df.shape[1]}个样本")
        print(f"    样本名示例: {prefixed_samples[0]} -> {prefixed_samples[-1]}")
        
    except FileNotFoundError:
        print(f"  警告：文件 {file_path} 不存在，跳过...")
    except Exception as e:
        print(f"  错误：加载 {file_path} 失败: {str(e)}")

# 检查是否有可用的数据集
if not datasets:
    print("错误：没有找到可用的数据文件！")
    exit(1)

print(f"\n成功加载 {len(datasets)} 个数据集")

# 找到所有数据集的共同基因
print("正在寻找共同基因...")
common_genes = None
for batch_name, df in datasets.items():
    if common_genes is None:
        common_genes = set(df.index)
    else:
        common_genes = common_genes.intersection(set(df.index))

common_genes = sorted(list(common_genes))
print(f"共同基因数量: {len(common_genes)}")

# 过滤每个数据集，只保留共同基因
print("正在过滤数据集...")
filtered_datasets = {}
for batch_name, df in datasets.items():
    filtered_df = df.loc[common_genes]
    filtered_datasets[batch_name] = filtered_df
    print(f"  {batch_name}: {filtered_df.shape[0]}个基因, {filtered_df.shape[1]}个样本")

# 合并所有表达矩阵
print("正在合并表达矩阵...")
expression_data = pd.concat(list(filtered_datasets.values()), axis=1)

# 创建批次信息
batch_info = pd.Series(all_batches, index=all_samples, name="Batch")

print("数据加载完成。")
print("合并后表达矩阵维度 (基因 x 样本):", expression_data.shape)
print("批次信息预览:")
print(batch_info.value_counts())
print(f"样本名添加前缀后总数: {len(all_samples)}")
print(f"样本名映射关系总数: {len(sample_mapping)}")
print("前缀样本名示例:")
for i, (orig, prefixed) in enumerate(list(sample_mapping.items())[:3]):
    print(f"  {orig} -> {prefixed}")
if len(sample_mapping) > 3:
    print(f"  ... 以及其他{len(sample_mapping)-3}个样本")

# 显示表达量基本统计信息
print("\n表达量基本统计:")
print(f"  最小值: {expression_data.min().min():.4f}")
print(f"  最大值: {expression_data.max().max():.4f}")
print(f"  平均值: {expression_data.mean().mean():.4f}")
print(f"  中位数: {expression_data.median().median():.4f}")
print("-" * 30)

# 确保表达矩阵的列名 (样本名) 和 batch_info 的索引完全一致
assert all(expression_data.columns == batch_info.index), "样本名不一致！"
print("样本名一致性检查通过。")

# 检查MYC和PVT1基因是否存在
print("\n关键基因检查:")
if 'MYC' in expression_data.index:
    print("  ✓ MYC基因存在")
    myc_stats = expression_data.loc['MYC'].describe()
    print(f"    MYC表达量统计: 平均值={myc_stats['mean']:.4f}, 中位数={myc_stats['50%']:.4f}, 标准差={myc_stats['std']:.4f}")
else:
    print("  ✗ MYC基因不存在")

if 'PVT1' in expression_data.index:
    print("  ✓ PVT1基因存在")
    pvt1_stats = expression_data.loc['PVT1'].describe()
    print(f"    PVT1表达量统计: 平均值={pvt1_stats['mean']:.4f}, 中位数={pvt1_stats['50%']:.4f}, 标准差={pvt1_stats['std']:.4f}")
else:
    print("  ✗ PVT1基因不存在")

# 加载MYC_PVT1注释文件
print("\n加载MYC_PVT1注释文件...")
try:
    annotation_df = pd.read_csv('MYC_PVT1_annotation.txt', sep='\t')
    print(f"  注释文件加载成功，包含{len(annotation_df)}个样本")
    
    # 检查每个样本是否在表达矩阵中
    valid_samples = []
    valid_statuses = []
    
    for idx, row in annotation_df.iterrows():
        sample = row['Sample']
        status = row['MYC_PVT1_Status']
        
        # 检查该样本是否在表达矩阵中
        if sample in expression_data.columns:
            valid_samples.append(sample)
            valid_statuses.append(status)
        else:
            print(f"    警告：样本 {sample} 不在表达矩阵中")
    
    # 创建hi_hi/lo_lo分组信息
    myc_pvt1_info = pd.Series(valid_statuses, index=valid_samples, name="MYC_PVT1_Status")
    
    print(f"  成功匹配样本数: {len(valid_samples)}")
    print(f"  分组统计:")
    print(f"    {myc_pvt1_info.value_counts().to_dict()}")
    
except FileNotFoundError:
    print("  警告：MYC_PVT1_annotation.txt文件不存在，请先运行数据处理脚本")
    myc_pvt1_info = None
except Exception as e:
    print(f"  错误：加载注释文件失败: {str(e)}")
    myc_pvt1_info = None


# ==============================================================================
# 2. 运行ComBat进行批次校正
# ==============================================================================
# ComBat函数需要表达数据（基因x样本）和批次信息
print("正在使用ComBat进行批次校正...")

# 检查是否有足够的批次进行ComBat分析
unique_batches = batch_info.unique()
if len(unique_batches) < 2:
    print("警告：只有一个批次，无法进行批次校正！")
    corrected_data = expression_data.copy()
else:
    print(f"检测到 {len(unique_batches)} 个批次: {list(unique_batches)}")
    
# 注意：ComBat处理的是数值型数据，如果数据中有非数值，需要提前处理
corrected_data = pycombat(expression_data, batch_info)
print("ComBat校正完成。")

print("校正后数据维度:", corrected_data.shape)
print("-" * 30)

# ==============================================================================
# 3. Cohort聚类分析评估批次效应校正效果
# ==============================================================================
print("正在进行Cohort聚类分析评估批次效应校正效果...")

# 评估校正前的聚类质量
print("\n=== 校正前聚类质量评估 ===")
metrics_before, linkage_before, clusters_before = evaluate_clustering_metrics(
    expression_data, batch_info, "校正前："
)

# 评估校正后的聚类质量
print("\n=== 校正后聚类质量评估 ===")
metrics_after, linkage_after, clusters_after = evaluate_clustering_metrics(
    corrected_data, batch_info, "校正后："
)

# 打印评估结果
print("\n" + "="*60)
print("COHORT聚类分析结果汇总")
print("="*60)

print(f"\n📊 聚类质量指标对比:")
print(f"{'指标':<25} {'校正前':<12} {'校正后':<12} {'变化':<10} {'效果'}")
print("-" * 70)

for metric_key, metric_name in [
    ('silhouette_score', 'Silhouette Score'),
    ('adjusted_rand_index', 'Adjusted Rand Index'),
    ('normalized_mutual_info', 'Normalized Mutual Info'),
    ('calinski_harabasz_score', 'Calinski-Harabasz'),
    ('davies_bouldin_score', 'Davies-Bouldin'),
    ('batch_purity', 'Batch Purity'),
    ('batch_separation', 'Batch Separation')
]:
    before_val = metrics_before.get(metric_key, 0)
    after_val = metrics_after.get(metric_key, 0)
    change = after_val - before_val
    
    # 判断效果（根据指标特性）
    if metric_key in ['davies_bouldin_score', 'batch_purity']:
        # 越低越好的指标
        effect = "✓ 改善" if change < 0 else "✗ 恶化"
    else:
        # 越高越好的指标
        effect = "✓ 改善" if change > 0 else "✗ 恶化"
    
    print(f"{metric_name:<25} {before_val:<12.3f} {after_val:<12.3f} {change:<+10.3f} {effect}")

# 绘制聚类评估指标对比图
print(f"\n📈 正在生成聚类评估指标对比图...")
fig_metrics = plot_clustering_evaluation(metrics_before, metrics_after)
fig_metrics.suptitle('Cohort Clustering Evaluation: Before vs After ComBat Correction', 
                     fontsize=16, fontweight='bold', y=0.98)
plt.show()

# 保存聚类评估图
fig_metrics.savefig('cohort_clustering_evaluation.pdf', dpi=300, bbox_inches='tight', facecolor='white')
print("聚类评估图已保存到: cohort_clustering_evaluation.pdf")

# 保存聚类结果
clustering_results = pd.DataFrame({
    'Sample': expression_data.columns,
    'Batch': batch_info.values,
    'Cluster_Before': clusters_before,
    'Cluster_After': clusters_after
})
clustering_results.to_csv('cohort_clustering_results.txt', sep='\t', index=False)
print("聚类结果已保存到: cohort_clustering_results.txt")

# ==============================================================================
# 4. 创建t-SNE分析图表
# ==============================================================================
print("\n🔄 开始t-SNE非线性降维分析...")

# t-SNE分析 - 去批次前
print("    正在运行t-SNE降维: Before ComBat - t-SNE by Batch")
fig_tsne1, ax_tsne1 = plt.subplots(figsize=(10, 8))
plot_tsne(expression_data, batch_info, 'Before ComBat - t-SNE by Batch', ax_tsne1)
plt.tight_layout()
plt.savefig('batch_effect_before_tsne.pdf', format='pdf', dpi=300, bbox_inches='tight', facecolor='white')
plt.close()

# t-SNE分析 - 去批次后
print("    正在运行t-SNE降维: After ComBat - t-SNE by Batch")
fig_tsne2, ax_tsne2 = plt.subplots(figsize=(10, 8))
plot_tsne(corrected_data, batch_info, 'After ComBat - t-SNE by Batch', ax_tsne2)
plt.tight_layout()
plt.savefig('batch_effect_after_tsne.pdf', format='pdf', dpi=300, bbox_inches='tight', facecolor='white')
plt.close()

# 保存合并图表
fig_batch_tsne, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))
plot_tsne(expression_data, batch_info, 'Before ComBat - t-SNE by Batch', ax1)
plot_tsne(corrected_data, batch_info, 'After ComBat - t-SNE by Batch', ax2)

# 添加总标题
fig_batch_tsne.suptitle('Batch Effect Correction - t-SNE Analysis by Batch', 
                        fontsize=20, fontweight='bold', y=0.98)
plt.tight_layout()
plt.savefig('batch_effect_tsne_by_batch.pdf', format='pdf', dpi=300, bbox_inches='tight', facecolor='white')
plt.close()

print("批次效应t-SNE分析图已保存为PDF格式")

# 如果有MYC/PVT1注释信息，也进行相应的t-SNE分析
if myc_pvt1_info is not None and len(myc_pvt1_info) > 0:
    print("\n🧬 按MYC/PVT1分组的t-SNE分析...")
    
    # 只使用有注释信息的样本
    annotated_samples = myc_pvt1_info.index
    filtered_expression = expression_data[annotated_samples]
    filtered_corrected = corrected_data[annotated_samples]
    
    if filtered_expression.shape[1] > 0:  # 检查是否有匹配的样本
        # t-SNE分析 - 去批次前
        fig_tsne3, ax_tsne3 = plt.subplots(figsize=(10, 8))
        plot_tsne(filtered_expression, myc_pvt1_info, 'Before ComBat - t-SNE by MYC/PVT1', ax_tsne3)
        plt.tight_layout()
        plt.savefig('myc_pvt1_before_tsne.pdf', format='pdf', dpi=300, bbox_inches='tight', facecolor='white')
        plt.close()
        
        # t-SNE分析 - 去批次后
        fig_tsne4, ax_tsne4 = plt.subplots(figsize=(10, 8))
        plot_tsne(filtered_corrected, myc_pvt1_info, 'After ComBat - t-SNE by MYC/PVT1', ax_tsne4)
        plt.tight_layout()
        plt.savefig('myc_pvt1_after_tsne.pdf', format='pdf', dpi=300, bbox_inches='tight', facecolor='white')
        plt.close()
        
        # 保存合并图表
        fig_myc_pvt1_tsne, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))
        plot_tsne(filtered_expression, myc_pvt1_info, 'Before ComBat - t-SNE by MYC/PVT1', ax1)
        plot_tsne(filtered_corrected, myc_pvt1_info, 'After ComBat - t-SNE by MYC/PVT1', ax2)
        
        # 添加总标题
        fig_myc_pvt1_tsne.suptitle('MYC/PVT1 Expression Groups - t-SNE Analysis', 
                                  fontsize=20, fontweight='bold', y=0.98)
        plt.tight_layout()
        plt.savefig('myc_pvt1_tsne_analysis.pdf', format='pdf', dpi=300, bbox_inches='tight', facecolor='white')
        plt.close()
        
        print("MYC/PVT1 t-SNE分析图已保存为PDF格式")
    else:
        print("警告：没有找到匹配的MYC/PVT1样本，跳过MYC/PVT1 t-SNE分析")
else:
    print("警告：MYC_PVT1_annotation.txt文件不存在或为空，跳过MYC/PVT1 t-SNE分析")

# ==============================================================================
# 5. 保存去批次后的结果
# ==============================================================================
print("正在保存去批次后的结果...")

# 保存去批次后的表达矩阵
corrected_data.to_csv('combined_expression_combat_corrected.txt', sep='\t')
print("去批次后的表达矩阵已保存到: combined_expression_combat_corrected.txt")

# 保存批次信息
batch_info.to_csv('batch_info.txt', sep='\t', header=True)
print("批次信息已保存到: batch_info.txt")

# 保存原始合并矩阵（用于对比）
expression_data.to_csv('combined_expression_before_combat.txt', sep='\t') 
print("原始合并矩阵已保存到: combined_expression_before_combat.txt")

# 如果有MYC_PVT1注释信息，也保存对应的子集
if myc_pvt1_info is not None:
    # 保存去批次后的MYC_PVT1样本子集
    annotated_samples = myc_pvt1_info.index
    filtered_corrected = corrected_data[annotated_samples]
    filtered_corrected.to_csv('myc_pvt1_expression_combat_corrected.txt', sep='\t')
    print("MYC_PVT1样本去批次后矩阵已保存到: myc_pvt1_expression_combat_corrected.txt")
    
    # 保存去批次前的MYC_PVT1样本子集
    filtered_expression = expression_data[annotated_samples]
    filtered_expression.to_csv('myc_pvt1_expression_before_combat.txt', sep='\t')
    print("MYC_PVT1样本去批次前矩阵已保存到: myc_pvt1_expression_before_combat.txt")

print("\n所有结果已保存完成！")
print("=" * 50)
print("\n输出文件说明:")
print("  📊 数据文件:")
print("    1. combined_expression_combat_corrected.txt - 去批次后的表达矩阵 (推荐用于后续分析)")
print("    2. combined_expression_before_combat.txt - 原始合并矩阵 (用于对比)")
print("    3. batch_info.txt - 批次信息文件")
print("    4. cohort_clustering_results.txt - Cohort聚类分析结果")
if myc_pvt1_info is not None:
    print("    5. myc_pvt1_expression_combat_corrected.txt - MYC_PVT1样本去批次后矩阵")
    print("    6. myc_pvt1_expression_before_combat.txt - MYC_PVT1样本去批次前矩阵")
print("\n  📈 图表文件:")
print("    1. cohort_clustering_evaluation.pdf - Cohort聚类评估图 (定量评估校正效果)")
print("    2. batch_effect_before_tsne.pdf - 批次效应校正前t-SNE图")
print("    3. batch_effect_after_tsne.pdf - 批次效应校正后t-SNE图")
print("    4. batch_effect_tsne_by_batch.pdf - 批次效应t-SNE对比图")
if myc_pvt1_info is not None:
    print("    5. myc_pvt1_before_tsne.pdf - MYC/PVT1校正前t-SNE图")
    print("    6. myc_pvt1_after_tsne.pdf - MYC/PVT1校正后t-SNE图")
    print("    7. myc_pvt1_tsne_analysis.pdf - MYC/PVT1 t-SNE对比图")
print("\n  🎨 图表特色:")
print("    • 使用Nature期刊配色方案")
print("    • 高分辨率300DPI输出")
print("    • 批次分组显示置信椭圆")
print("    • 清晰的去批次前后对比")
print("    • t-SNE非线性降维，强调局部聚类结构")