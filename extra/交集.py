#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
样本名交集处理脚本 - 只处理liu数据集
处理liu表型数据和过滤数据文件，使样本名完全相同
"""

import pandas as pd
import numpy as np
import os

def get_samples_intersection(pheno_file, filtered_file):
    """
    获取表型数据和过滤数据的样本名交集
    """
    print(f"处理 liu 数据集...")
    
    # 读取表型数据
    pheno_df = pd.read_csv(pheno_file, sep='\t')
    pheno_samples = set(pheno_df.iloc[:, 0].tolist())
    print(f"表型数据样本数: {len(pheno_samples)}")
    
    # 读取过滤数据的第一行（样本名）
    with open(filtered_file, 'r') as f:
        header = f.readline().strip()
    
    # 解析样本名 - liu数据的第一列是空的，从第二列开始是样本名
    filtered_samples = set(header.split('\t')[1:])
    print(f"过滤数据样本数: {len(filtered_samples)}")
    
    # 计算交集
    intersection = pheno_samples.intersection(filtered_samples)
    print(f"交集样本数: {len(intersection)}")
    
    if len(intersection) == 0:
        print(f"警告：没有共同样本！")
        return None, None, None
    
    return pheno_samples, filtered_samples, intersection

def filter_and_save_pheno(pheno_file, intersection, output_file):
    """
    根据交集过滤表型数据并保存
    """
    pheno_df = pd.read_csv(pheno_file, sep='\t')
    
    # 过滤出交集中的样本
    filtered_pheno = pheno_df[pheno_df.iloc[:, 0].isin(intersection)]
    
    # 保存过滤后的文件
    filtered_pheno.to_csv(output_file, sep='\t', index=False)
    print(f"保存过滤后的表型数据到: {output_file}")

def filter_and_save_filtered(filtered_file, intersection, output_file):
    """
    根据交集过滤表达数据并保存
    """
    # 读取完整的过滤数据
    filtered_df = pd.read_csv(filtered_file, sep='\t')
    
    # 获取列名
    columns = filtered_df.columns.tolist()
    
    # liu数据第一列是空的，实际上是基因名
    gene_col = columns[0]
    sample_cols = [col for col in columns[1:] if col in intersection]
    
    # 创建新的DataFrame
    keep_columns = [gene_col] + sample_cols
    filtered_data = filtered_df[keep_columns]
    
    # 保存过滤后的文件
    filtered_data.to_csv(output_file, sep='\t', index=False)
    print(f"保存过滤后的表达数据到: {output_file}")

def main():
    """
    主函数：处理liu数据集
    """
    print("开始处理liu数据集的样本名交集...")
    
    # 定义文件路径
    pheno_file = 'liu_phenodata.txt'
    filtered_file = 'liu.txt'
    pheno_output = 'liu_phenodata_intersect.txt'
    filtered_output = 'liu_intersect.txt'
    
    try:
        # 检查文件是否存在
        if not os.path.exists(pheno_file):
            print(f"文件不存在: {pheno_file}")
            return
        if not os.path.exists(filtered_file):
            print(f"文件不存在: {filtered_file}")
            return
        
        # 获取交集
        pheno_samples, filtered_samples, intersection = get_samples_intersection(pheno_file, filtered_file)
        
        if intersection is None:
            print("处理失败：没有共同样本")
            return
        
        # 过滤并保存文件
        filter_and_save_pheno(pheno_file, intersection, pheno_output)
        filter_and_save_filtered(filtered_file, intersection, filtered_output)
        
        print(f"✓ liu 数据集处理完成")
        print(f"生成的文件：")
        print(f"  - {pheno_output}")
        print(f"  - {filtered_output}")
        
    except Exception as e:
        print(f"处理liu数据集时出错: {str(e)}")

if __name__ == "__main__":
    main()
