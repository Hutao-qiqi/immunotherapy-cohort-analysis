import pandas as pd
import numpy as np
import os

def main():
    """
    数据处理脚本：
    1. 对ra.txt只保留ra_phenodata.txt中出现的样本
    2. 对IMvigor210_TPM_annotated_by_SYMBOL.txt只保留IMvigor210_phenoData.txt中出现的样本
    3. 根据hugo_phenoData.txt的Biopsy Time列，只保留值为"pre-treatment"的hugo.txt样本，并去掉.baseline后缀
    4. 对liu.txt只保留liu_phenodata.txt中出现的样本
    """
    
    print("开始数据处理...")
    
    # 任务1: 处理ra.txt - 只保留ra_phenodata.txt中出现的样本
    print("1. 处理ra.txt文件...")
    try:
        # 读取ra_phenodata.txt获取有效样本列表
        ra_pheno_df = pd.read_csv('ra_phenodata.txt', sep='\t')
        valid_ra_samples = ra_pheno_df['Patient'].tolist()
        print(f"   ra_phenodata.txt中的样本数: {len(valid_ra_samples)}")
        
        # 读取ra.txt文件
        ra_df = pd.read_csv('ra.txt', sep='\t', index_col=0)
        print(f"   ra.txt原始数据形状: {ra_df.shape}")
        
        # 处理样本名称匹配问题：ra.txt的列名是"Pt1_Pre_AD101148-6"格式，需要提取Patient ID部分
        # 创建一个映射：从ra.txt的列名到Patient ID，同时过滤掉包含"On"的列
        ra_columns_map = {}
        for col in ra_df.columns:
            if '_' in col:
                # 从"Pt1_Pre_AD101148-6"中提取"Pt1"
                patient_id = col.split('_')[0]
                # 只保留ra_phenodata.txt中出现的样本，且排除包含"On"的列
                if patient_id in valid_ra_samples and 'On' not in col:
                    ra_columns_map[col] = patient_id
        
        # 只保留ra_phenodata.txt中出现的样本，且不包含"On"
        filtered_columns = list(ra_columns_map.keys())
        ra_df_filtered = ra_df[filtered_columns]
        print(f"   匹配的样本数: {len(filtered_columns)}")
        print(f"   过滤后数据形状: {ra_df_filtered.shape}")
        
        # 保存处理后的数据
        ra_df_filtered.to_csv('ra_filtered.txt', sep='\t')
        print("   ra.txt处理完成，保存为ra_filtered.txt")
        
    except Exception as e:
        print(f"   处理ra.txt时出错: {e}")
    
    # 任务2: 处理IMvigor数据 - 只保留IMvigor210_phenoData.txt中出现的样本
    print("\n2. 处理IMvigor数据...")
    try:
        # 读取IMvigor210_TPM_annotated_by_SYMBOL.txt表达矩阵
        imvigor_file = 'IMvigor210_TPM_annotated_by_SYMBOL.txt'
        print(f"   读取IMvigor表达矩阵: {imvigor_file}")
        
        # 读取表达矩阵
        imvigor_df = pd.read_csv(imvigor_file, sep='\t', index_col=0)
        print(f"   IMvigor表达矩阵原始形状: {imvigor_df.shape}")
        
        # 读取临床数据
        pheno_df = pd.read_csv('IMvigor210_phenoData.txt', sep='\t', index_col=0)
        print(f"   临床数据形状: {pheno_df.shape}")
        
        # 获取临床数据中的所有样本
        valid_samples = pheno_df.index.tolist()
        print(f"   临床数据中的样本数: {len(valid_samples)}")
        
        # 过滤表达矩阵，只保留临床数据中出现的样本
        common_samples = imvigor_df.columns.intersection(valid_samples)
        imvigor_df_filtered = imvigor_df[common_samples]
        print(f"   匹配的样本数: {len(common_samples)}")
        print(f"   过滤后IMvigor表达矩阵形状: {imvigor_df_filtered.shape}")
        
        # 保存过滤后的数据
        imvigor_df_filtered.to_csv('IMvigor_filtered.txt', sep='\t')
        print("   IMvigor数据处理完成，保存为IMvigor_filtered.txt")
        
    except Exception as e:
        print(f"   处理IMvigor数据时出错: {e}")
    
    # 任务3: 处理hugo.txt - 只保留pre-treatment样本
    print("\n3. 处理hugo.txt文件...")
    try:
        # 读取hugo临床数据
        hugo_pheno_df = pd.read_csv('hugo_phenoData.txt', sep='\t')
        print(f"   hugo临床数据形状: {hugo_pheno_df.shape}")
        
        # 找到Biopsy Time为"pre-treatment"的样本
        pre_treatment_samples = hugo_pheno_df[hugo_pheno_df['Biopsy Time'] == 'pre-treatment']['Patient ID'].tolist()
        print(f"   pre-treatment样本数: {len(pre_treatment_samples)}")
        
        # 读取hugo.txt表达数据
        hugo_df = pd.read_csv('hugo.txt', sep='\t', index_col=0)
        print(f"   hugo.txt原始数据形状: {hugo_df.shape}")
        
        # 处理样本名称匹配问题：hugo.txt的列名是"Pt1.baseline"格式，需要提取Patient ID部分
        # 创建一个映射：从hugo.txt的列名到Patient ID
        hugo_columns_map = {}
        for col in hugo_df.columns:
            if '.' in col:
                patient_id = col.split('.')[0]  # 提取Pt1、Pt2等
                if patient_id in pre_treatment_samples:
                    hugo_columns_map[col] = patient_id
        
        # 过滤样本，只保留pre-treatment样本
        filtered_columns = list(hugo_columns_map.keys())
        hugo_df_filtered = hugo_df[filtered_columns]
        
        # 去掉列名中的.baseline后缀
        new_column_names = {}
        for col in hugo_df_filtered.columns:
            if '.baseline' in col:
                new_name = col.replace('.baseline', '')
                new_column_names[col] = new_name
            else:
                new_column_names[col] = col
        
        hugo_df_filtered = hugo_df_filtered.rename(columns=new_column_names)
        print(f"   过滤后hugo.txt数据形状: {hugo_df_filtered.shape}")
        print(f"   匹配的样本: {len(filtered_columns)}个")
        print(f"   去掉.baseline后缀完成")
        
        # 保存过滤后的数据
        hugo_df_filtered.to_csv('hugo_filtered.txt', sep='\t')
        print("   hugo.txt处理完成，保存为hugo_filtered.txt")
        
    except Exception as e:
        print(f"   处理hugo.txt时出错: {e}")
    
    # 任务4: 处理liu.txt - 只保留liu_phenodata.txt中出现的样本
    print("\n4. 处理liu.txt文件...")
    try:
        # 读取liu_phenodata.txt获取有效样本列表
        # 尝试多种编码格式
        liu_pheno_df = None
        encodings = ['utf-8', 'utf-16', 'utf-16le', 'utf-16be', 'gbk', 'gb2312', 'latin-1', 'cp1252']
        
        for encoding in encodings:
            try:
                print(f"   尝试使用编码: {encoding}")
                with open("liu_phenodata.txt", "r", encoding=encoding) as f:
                    liu_pheno_df = pd.read_csv(f, sep='\t', index_col=0, low_memory=False)
                print(f"   成功使用编码: {encoding}")
                break
            except UnicodeDecodeError:
                continue
            except Exception as e:
                print(f"   编码 {encoding} 失败: {str(e)}")
                continue
        
        if liu_pheno_df is None:
            raise Exception("无法读取文件，尝试了所有常见编码格式")
            
        valid_liu_samples = liu_pheno_df.index.tolist()
        print(f"   liu_phenodata.txt中的样本数: {len(valid_liu_samples)}")
        
        # 读取liu.txt文件
        liu_df = pd.read_csv('liu.txt', sep='\t', index_col=0)
        print(f"   liu.txt原始数据形状: {liu_df.shape}")
        
        # 只保留liu_phenodata.txt中出现的样本
        common_samples = liu_df.columns.intersection(valid_liu_samples)
        liu_df_filtered = liu_df[common_samples]
        print(f"   匹配的样本数: {len(common_samples)}")
        print(f"   过滤后数据形状: {liu_df_filtered.shape}")
        
        # 保存处理后的数据
        liu_df_filtered.to_csv('liu_filtered.txt', sep='\t')
        print("   liu.txt处理完成，保存为liu_filtered.txt")
        
    except Exception as e:
        print(f"   处理liu.txt时出错: {e}")
    
    print("\n数据处理完成!")
    print("输出文件:")
    print("- ra_filtered.txt: 只保留ra_phenodata.txt中出现的样本")
    print("- IMvigor_filtered.txt: 只保留IMvigor210_phenoData.txt中出现的样本")
    print("- hugo_filtered.txt: 只保留pre-treatment样本的hugo.txt，已去掉.baseline后缀")
    print("- liu_filtered.txt: 只保留liu_phenodata.txt中出现的样本")

if __name__ == "__main__":
    main()
