import pandas as pd
import numpy as np

def load_expression_samples():
    """获取表达数据的样本ID列表"""
    # 读取表达数据的第一行（样本ID）
    with open('combined_expression_combat_corrected.txt', 'r') as f:
        header = f.readline().strip().split('\t')
    
    print(f"表达数据总样本数: {len(header)}")
    return header

def process_hugo_data():
    """处理hugo数据集的临床信息"""
    print("\n处理Hugo数据集...")
    df = pd.read_csv('hugo_phenoData.txt', sep='\t')
    
    # 提取生存信息
    # Overall Survival (days?) -> 转换为月份
    # Vital Status: Dead/Alive -> 转换为事件指示符
    
    survival_data = []
    for _, row in df.iterrows():
        patient_id = f"hugo_{row['Patient ID']}"
        
        # 生存时间（天数转换为月份）
        survival_time = pd.to_numeric(row['Overall Survival'], errors='coerce')
        if pd.notna(survival_time):
            survival_time_months = survival_time / 30.44  # 转换为月份
        else:
            survival_time_months = np.nan
        
        # 生存状态
        vital_status = row['Vital Status']
        if vital_status == 'Dead':
            event = 1
        elif vital_status == 'Alive':
            event = 0
        else:
            event = np.nan
        
        survival_data.append({
            'Sample_ID': patient_id,
            'OS_months': survival_time_months,
            'OS_event': event,
            'Dataset': 'Hugo'
        })
    
    hugo_df = pd.DataFrame(survival_data)
    print(f"Hugo数据集: {len(hugo_df)}个样本")
    print(f"有效生存数据: {hugo_df.dropna().shape[0]}个样本")
    
    return hugo_df

def process_imvigor210_data():
    """处理IMvigor210数据集的临床信息"""
    print("\n处理IMvigor210数据集...")
    df = pd.read_csv('IMvigor210_phenoData.txt', sep='\t')
    
    # 提取生存信息
    # os: 生存时间（月份）
    # censOS: 事件指示符 (0=censored/alive, 1=event/dead)
    
    survival_data = []
    for _, row in df.iterrows():
        patient_id = f"IMvigor210_{row.iloc[0]}"  # 第一列是样本ID，添加前缀
        
        # 生存时间（已经是月份）
        survival_time = pd.to_numeric(row['os'], errors='coerce')
        
        # 生存状态
        event = pd.to_numeric(row['censOS'], errors='coerce')
        
        survival_data.append({
            'Sample_ID': patient_id,
            'OS_months': survival_time,
            'OS_event': event,
            'Dataset': 'IMvigor210'
        })
    
    imvigor_df = pd.DataFrame(survival_data)
    print(f"IMvigor210数据集: {len(imvigor_df)}个样本")
    print(f"有效生存数据: {imvigor_df.dropna().shape[0]}个样本")
    
    return imvigor_df

def process_ra_data(expression_samples):
    """处理RA数据集的临床信息"""
    print("\n处理RA数据集...")
    df = pd.read_csv('ra_phenodata.txt', sep='\t')
    
    # 提取生存信息
    # Time to Death (weeks) -> 转换为月份
    # Dead/Alive (Dead = True) -> 转换为事件指示符
    
    survival_data = []
    for _, row in df.iterrows():
        # RA数据集的样本ID格式复杂，需要在表达数据中找到匹配的样本
        patient_base_id = row['Patient']  # 如 "Pt10"
        # 在表达数据中查找包含这个患者ID的样本（需要ra_前缀匹配）
        matched_samples = [sid for sid in expression_samples if sid.startswith('ra_') and patient_base_id in sid]
        
        if not matched_samples:
            continue  # 如果没有匹配的样本，跳过
        
        # 如果有多个匹配的样本，选择第一个
        patient_id = matched_samples[0]
        
        # 生存时间（周数转换为月份）
        survival_time = pd.to_numeric(row['Time to Death\r\n(weeks)'], errors='coerce')
        if pd.notna(survival_time):
            survival_time_months = survival_time / 4.34524  # 转换为月份
        else:
            survival_time_months = np.nan
        
        # 生存状态
        vital_status = row['Dead/Alive\r\n(Dead = True)']
        if vital_status == True or vital_status == 'TRUE':
            event = 1
        elif vital_status == False or vital_status == 'FALSE':
            event = 0
        else:
            event = np.nan
        
        survival_data.append({
            'Sample_ID': patient_id,
            'OS_months': survival_time_months,
            'OS_event': event,
            'Dataset': 'RA'
        })
    
    ra_df = pd.DataFrame(survival_data)
    print(f"RA数据集: {len(ra_df)}个样本")
    print(f"有效生存数据: {ra_df.dropna().shape[0]}个样本")
    
    return ra_df

def process_liu_data():
    """处理Liu数据集的临床信息"""
    print("\n处理Liu数据集...")
    df = pd.read_csv('liu_phenodata_filtered.txt', sep='\t')
    
    # 提取生存信息
    # OS: 总生存时间（天数？）-> 转换为月份
    # dead: 死亡状态 (1=dead, 0=alive)
    
    survival_data = []
    for _, row in df.iterrows():
        patient_id = f"liu_{row.iloc[0]}"  # 第一列是样本ID，添加前缀
        
        # 生存时间（天数转换为月份）
        survival_time = pd.to_numeric(row['OS'], errors='coerce')
        if pd.notna(survival_time):
            survival_time_months = survival_time / 30.44  # 转换为月份
        else:
            survival_time_months = np.nan
        
        # 生存状态
        event = pd.to_numeric(row['dead'], errors='coerce')
        
        survival_data.append({
            'Sample_ID': patient_id,
            'OS_months': survival_time_months,
            'OS_event': event,
            'Dataset': 'Liu'
        })
    
    liu_df = pd.DataFrame(survival_data)
    print(f"Liu数据集: {len(liu_df)}个样本")
    print(f"有效生存数据: {liu_df.dropna().shape[0]}个样本")
    
    return liu_df

def main():
    print("开始合并临床信息文件...")
    
    # 获取表达数据的样本列表
    expression_samples = load_expression_samples()
    
    # 处理各个数据集
    hugo_df = process_hugo_data()
    imvigor_df = process_imvigor210_data()
    ra_df = process_ra_data(expression_samples)
    liu_df = process_liu_data()
    
    # 合并所有数据
    print("\n合并所有数据集...")
    combined_df = pd.concat([hugo_df, imvigor_df, ra_df, liu_df], ignore_index=True)
    
    print(f"\n合并后总样本数: {len(combined_df)}")
    print(f"有效生存数据: {combined_df.dropna().shape[0]}个样本")
    
    # 与表达数据匹配
    print("\n与表达数据匹配...")
    matched_df = combined_df[combined_df['Sample_ID'].isin(expression_samples)]
    
    print(f"匹配的样本数: {len(matched_df)}")
    print(f"匹配的有效生存数据: {matched_df.dropna().shape[0]}个样本")
    
    # 按数据集统计匹配情况
    print("\n各数据集匹配情况:")
    for dataset in ['Hugo', 'IMvigor210', 'RA', 'Liu']:
        dataset_matched = matched_df[matched_df['Dataset'] == dataset]
        dataset_valid = dataset_matched.dropna()
        print(f"{dataset}: {len(dataset_matched)}个匹配样本, {len(dataset_valid)}个有效生存数据")
    
    # 保存合并的临床数据
    combined_df.to_csv('combined_clinical_data.txt', sep='\t', index=False)
    matched_df.to_csv('matched_clinical_data.txt', sep='\t', index=False)
    
    # 保存只包含有效生存数据的匹配样本
    valid_matched_df = matched_df.dropna()
    valid_matched_df.to_csv('updated_survival_data.txt', sep='\t', index=False)
    
    print(f"\n文件已保存:")
    print(f"- combined_clinical_data.txt: 所有合并的临床数据")
    print(f"- matched_clinical_data.txt: 与表达数据匹配的临床数据")
    print(f"- updated_survival_data.txt: 有效的生存分析数据（与表达数据严格对应）")
    
    # 显示最终结果摘要
    print(f"\n最终可用于生存分析的样本: {len(valid_matched_df)}个")
    
    return valid_matched_df

if __name__ == "__main__":
    result = main()
    print(result.head()) 