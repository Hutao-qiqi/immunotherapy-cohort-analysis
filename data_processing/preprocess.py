import pandas as pd
import numpy as np
import mygene
import os
import warnings
warnings.filterwarnings('ignore')

def convert_entrez_to_symbol(gene_ids):
    """
    使用mygene包将Entrez Gene ID转换为基因符号
    """
    mg = mygene.MyGeneInfo()
    
    # 过滤掉非数字ID
    valid_ids = []
    for gene_id in gene_ids:
        try:
            int(gene_id)
            valid_ids.append(gene_id)
        except:
            valid_ids.append(gene_id)
    
    # 批量查询
    result = mg.querymany(valid_ids, scopes='entrezgene', fields='symbol', species='human')
    
    # 创建ID到symbol的映射
    id_to_symbol = {}
    for item in result:
        if 'symbol' in item:
            id_to_symbol[str(item['query'])] = item['symbol']
        else:
            id_to_symbol[str(item['query'])] = str(item['query'])  # 保持原ID
    
    return id_to_symbol

def process_expression_matrix(file_path, file_type):
    """
    处理表达矩阵文件
    """
    print(f"正在处理 {file_path}...")
    
    if file_type == 'IMvigor210':
        # 基因符号为行名，样本为列名
        df = pd.read_csv(file_path, sep='\t', index_col=0)
        
    elif file_type == 'liu':
        # 基因符号为行名，样本为列名
        df = pd.read_csv(file_path, sep='\t', index_col=0)
        
    elif file_type in ['ra', 'hugo']:
        # 基因ID为行名，样本为列名，需要转换ID到symbol
        df = pd.read_csv(file_path, sep='\t', index_col=0)
        
        # 转换基因ID到symbol
        print(f"正在转换基因ID到symbol...")
        gene_ids = df.index.tolist()
        id_to_symbol = convert_entrez_to_symbol(gene_ids)
        
        # 创建新的索引
        new_index = [id_to_symbol.get(str(gene_id), str(gene_id)) for gene_id in gene_ids]
        df.index = new_index
        
        # 去重，如果有重复的基因符号，取平均值
        df = df.groupby(df.index).mean()
    
    return df

def analyze_myc_pvt1(df, dataset_name):
    """
    分析MYC和PVT1的表达量，标记样本状态
    """
    print(f"正在分析 {dataset_name} 中的MYC和PVT1表达...")
    
    # 检查是否有MYC和PVT1基因
    if 'MYC' not in df.index:
        print(f"警告：{dataset_name} 中未找到MYC基因")
        return None, None
    
    if 'PVT1' not in df.index:
        print(f"警告：{dataset_name} 中未找到PVT1基因")
        return None, None
    
    # 获取MYC和PVT1的表达量
    myc_expr = df.loc['MYC']
    pvt1_expr = df.loc['PVT1']
    
    # 检查数据类型，确保都是数值类型
    print(f"  检查数据类型...")
    
    # 检查MYC表达量
    non_numeric_myc = []
    for sample in myc_expr.index:
        value = myc_expr[sample]
        if not isinstance(value, (int, float, np.integer, np.floating)):
            non_numeric_myc.append(f"{sample}: {value} (type: {type(value).__name__})")
    
    # 检查PVT1表达量
    non_numeric_pvt1 = []
    for sample in pvt1_expr.index:
        value = pvt1_expr[sample]
        if not isinstance(value, (int, float, np.integer, np.floating)):
            non_numeric_pvt1.append(f"{sample}: {value} (type: {type(value).__name__})")
    
    # 报告非数值数据
    if non_numeric_myc:
        print(f"  警告：发现MYC非数值数据:")
        for item in non_numeric_myc[:5]:  # 只显示前5个
            print(f"    {item}")
        if len(non_numeric_myc) > 5:
            print(f"    ... 还有{len(non_numeric_myc)-5}个非数值数据")
    
    if non_numeric_pvt1:
        print(f"  警告：发现PVT1非数值数据:")
        for item in non_numeric_pvt1[:5]:  # 只显示前5个
            print(f"    {item}")
        if len(non_numeric_pvt1) > 5:
            print(f"    ... 还有{len(non_numeric_pvt1)-5}个非数值数据")
    
    # 转换为数值类型，非数值数据转为NaN
    myc_expr = pd.to_numeric(myc_expr, errors='coerce')
    pvt1_expr = pd.to_numeric(pvt1_expr, errors='coerce')
    
    # 检查是否有NaN值
    myc_nan_count = myc_expr.isna().sum()
    pvt1_nan_count = pvt1_expr.isna().sum()
    
    if myc_nan_count > 0:
        print(f"  警告：MYC表达量中有{myc_nan_count}个NaN值")
    if pvt1_nan_count > 0:
        print(f"  警告：PVT1表达量中有{pvt1_nan_count}个NaN值")
    
    # 删除含有NaN值的样本
    valid_samples = myc_expr.dropna().index.intersection(pvt1_expr.dropna().index)
    if len(valid_samples) < len(df.columns):
        removed_samples = len(df.columns) - len(valid_samples)
        print(f"  已删除{removed_samples}个含有NaN值的样本")
        myc_expr = myc_expr[valid_samples]
        pvt1_expr = pvt1_expr[valid_samples]
        df = df[valid_samples]  # 更新数据框
    
    print(f"  数据类型检查完成，有效样本数: {len(valid_samples)}")
    
    # 使用中位数作为阈值，>= 中位数为高表达，< 中位数为低表达
    myc_median = myc_expr.median()
    pvt1_median = pvt1_expr.median()
    
    print(f"  MYC中位数表达量: {myc_median:.4f}")
    print(f"  PVT1中位数表达量: {pvt1_median:.4f}")
    
    # 标记样本状态
    sample_status = {}
    hi_hi_samples = []
    hi_lo_samples = []
    lo_hi_samples = []
    lo_lo_samples = []
    
    for sample in df.columns:
        myc_val = myc_expr[sample]
        pvt1_val = pvt1_expr[sample]
        
        # 判断MYC表达水平（>= 中位数为高表达）
        if myc_val > myc_median:
            myc_level = 'high'
        else:
            myc_level = 'low'
        
        # 判断PVT1表达水平（>= 中位数为高表达）
        if pvt1_val > pvt1_median:
            pvt1_level = 'high'
        else:
            pvt1_level = 'low'
        
        # 分为四种状态
        if myc_level == 'high' and pvt1_level == 'high':
            sample_status[sample] = 'hi_hi'
            hi_hi_samples.append(sample)
        elif myc_level == 'high' and pvt1_level == 'low':
            sample_status[sample] = 'hi_lo'
            hi_lo_samples.append(sample)
        elif myc_level == 'low' and pvt1_level == 'high':
            sample_status[sample] = 'lo_hi'
            lo_hi_samples.append(sample)
        else:  # myc_level == 'low' and pvt1_level == 'low'
            sample_status[sample] = 'lo_lo'
            lo_lo_samples.append(sample)
    
    # 详细展示各个成分
    print(f"\n  {dataset_name} 详细样本分类:")
    print(f"    总样本数: {len(df.columns)}")
    
    # 统计单个基因的分布
    myc_high = sum(1 for sample in df.columns if myc_expr[sample] >= myc_median)
    myc_low = sum(1 for sample in df.columns if myc_expr[sample] < myc_median)
    
    pvt1_high = sum(1 for sample in df.columns if pvt1_expr[sample] >= pvt1_median)
    pvt1_low = sum(1 for sample in df.columns if pvt1_expr[sample] < pvt1_median)
    
    print(f"\n    单基因分布:")
    print(f"      MYC高表达 (>={myc_median:.4f}): {myc_high}个样本")
    print(f"      MYC低表达 (<{myc_median:.4f}): {myc_low}个样本")
    print(f"      PVT1高表达 (>={pvt1_median:.4f}): {pvt1_high}个样本")
    print(f"      PVT1低表达 (<{pvt1_median:.4f}): {pvt1_low}个样本")
    
    print(f"\n    最终分类:")
    print(f"      hi_hi (MYC高 & PVT1高): {len(hi_hi_samples)}个样本")
    print(f"      hi_lo (MYC高 & PVT1低): {len(hi_lo_samples)}个样本")
    print(f"      lo_hi (MYC低 & PVT1高): {len(lo_hi_samples)}个样本")
    print(f"      lo_lo (MYC低 & PVT1低): {len(lo_lo_samples)}个样本")
    
    # 显示各类样本列表
    sample_groups = [
        (hi_hi_samples, "hi_hi"),
        (hi_lo_samples, "hi_lo"),
        (lo_hi_samples, "lo_hi"),
        (lo_lo_samples, "lo_lo")
    ]
    
    for samples, group_name in sample_groups:
        if samples:
            print(f"\n  {group_name}样本列表:")
            for i, sample in enumerate(samples):
                if i < 5:  # 只显示前5个
                    print(f"    {sample}: MYC={myc_expr[sample]:.4f}, PVT1={pvt1_expr[sample]:.4f}")
                elif i == 5:
                    print(f"    ... (还有{len(samples)-5}个样本)")
                    break
    
    # 保留所有样本，不过滤
    all_samples = list(df.columns)
    all_status = {sample: sample_status[sample] for sample in all_samples}
    
    # 返回原始表达矩阵和样本状态
    print(f"\n  保留所有样本数: {len(all_samples)}")
    print(f"  基因数: {len(df.index)}")
    
    return df, all_status

def main():
    """
    主函数
    """
    # 文件路径和类型
    files = {
        'IMvigor_filtered.txt': 'IMvigor210',
        'liu_filtered.txt': 'liu',
        'ra_filtered.txt': 'ra',
        'hugo_filtered.txt': 'hugo'
    }
    
    all_sample_status = {}
    filtered_matrices = {}
    
    # 处理每个文件
    for file_path, file_type in files.items():
        if not os.path.exists(file_path):
            print(f"警告：文件 {file_path} 不存在，跳过...")
            continue
        
        try:
            # 处理表达矩阵
            df = process_expression_matrix(file_path, file_type)
            
            # 分析MYC和PVT1
            processed_df, sample_status = analyze_myc_pvt1(df, file_type)
            
            if processed_df is not None and sample_status is not None:
                # 保存处理后的矩阵（包含所有样本）
                filtered_matrices[file_type] = processed_df
                all_sample_status.update(sample_status)
                
                # 保存处理后的表达矩阵
                output_file = f"{file_type}_filtered1.txt"
                processed_df.to_csv(output_file, sep='\t')
                print(f"  已保存处理后的矩阵到: {output_file}")
            
        except Exception as e:
            print(f"处理 {file_path} 时出错: {str(e)}")
    
    # 生成注释文件
    if all_sample_status:
        print(f"\n生成注释文件...")
        annotation_df = pd.DataFrame(list(all_sample_status.items()), 
                                   columns=['Sample', 'MYC_PVT1_Status'])
        
        # 按状态排序
        annotation_df = annotation_df.sort_values('MYC_PVT1_Status')
        
        # 保存注释文件
        annotation_df.to_csv('MYC_PVT1_annotation.txt', sep='\t', index=False)
        print(f"已生成注释文件: MYC_PVT1_annotation.txt")
        
        # 打印统计信息
        status_counts = annotation_df['MYC_PVT1_Status'].value_counts()
        print(f"\n最终统计:")
        print(f"  总样本数: {len(annotation_df)}")
        print(f"  hi_hi样本数: {status_counts.get('hi_hi', 0)}")
        print(f"  hi_lo样本数: {status_counts.get('hi_lo', 0)}")
        print(f"  lo_hi样本数: {status_counts.get('lo_hi', 0)}")
        print(f"  lo_lo样本数: {status_counts.get('lo_lo', 0)}")
        
        # 按数据集分组统计
        datasets = {}
        for sample in annotation_df['Sample']:
            if sample.startswith('Patient'):
                dataset = 'liu'
            elif sample.startswith('Pt'):
                dataset = 'ra_or_hugo'
            elif sample.startswith('SAM'):
                dataset = 'IMvigor210'
            else:
                dataset = 'unknown'
            
            if dataset not in datasets:
                datasets[dataset] = []
            datasets[dataset].append(sample)
        
        print(f"\n各数据集样本数统计:")
        for dataset, samples in datasets.items():
            print(f"  {dataset}: {len(samples)}")
            if dataset == 'liu':
                # 详细统计liu数据集
                liu_samples = [s for s in annotation_df['Sample'] if s.startswith('Patient')]
                liu_anno = annotation_df[annotation_df['Sample'].isin(liu_samples)]
                liu_status_counts = liu_anno['MYC_PVT1_Status'].value_counts()
                print(f"    liu数据集详细统计:")
                for status in ['hi_hi', 'hi_lo', 'lo_hi', 'lo_lo']:
                    count = liu_status_counts.get(status, 0)
                    print(f"      {status}: {count}")
        
        # 显示前几行注释文件
        print(f"\n注释文件前10行:")
        print(annotation_df.head(10))
    
    print("\n数据处理完成！")

if __name__ == "__main__":
    main()
