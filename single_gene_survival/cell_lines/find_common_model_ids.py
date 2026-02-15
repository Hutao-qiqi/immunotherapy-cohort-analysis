import pandas as pd
import os

# --- 配置 ---
file_crispr_effect = 'CRISPRGeneEffect.csv'
file_copy_number = 'OmicsCNGeneWGS.csv'
file_model_meta = 'Model.csv'
output_file = 'common_model_ids.txt'

def preview_csv_structure(filename, num_rows=5, num_cols=10):
    """打印CSV文件的行列名和前几行数据以供预览。"""
    print(f"\n{'='*20} 正在预览文件: {filename} {'='*20}")
    if not os.path.exists(filename):
        print(f"错误: 文件 '{filename}' 未找到。")
        print(f"{'='*20} 预览结束 {'='*20}\n")
        return
    try:
        # 始终将第一列作为索引进行预览，以便观察
        df_preview = pd.read_csv(filename, index_col=0, nrows=num_rows, low_memory=False)
        print("预览时，尝试将第一列作为行名（索引）。")
    except Exception:
        df_preview = pd.read_csv(filename, nrows=num_rows)
        print("预览时，按默认方式读取。")
    print(f"\n--- 前 {num_rows} 行数据预览 ---")
    print(df_preview)
    print(f"\n--- 前 {min(len(df_preview.columns), num_cols)} 个列名 ---")
    print(list(df_preview.columns[:num_cols]))
    print(f"\n--- 前 {min(len(df_preview.index), num_rows)} 个行名（索引） ---")
    print(list(df_preview.index[:num_rows]))
    print(f"\n{'='*20} 预览结束: {filename} {'='*20}\n")

# --- 更健壮的ID提取函数 ---

def get_model_ids_from_index(filename):
    """
    稳健的方法：从一个矩阵文件中读取行名（索引）作为ID。
    适用于 '模型 x 基因' 格式的文件。
    """
    if not os.path.exists(filename): return None
    try:
        # 明确告诉pandas第一列是索引，并只读取索引
        df = pd.read_csv(filename, index_col=0, usecols=[0])
        return set(df.index)
    except Exception as e:
        print(f"读取 '{filename}' 的索引时出错: {e}")
        return None

def get_model_ids_from_columns(filename):
    """
    稳健的方法：从一个矩阵文件中读取列名作为ID。
    适用于 '基因 x 模型' 格式的文件。
    """
    if not os.path.exists(filename): return None
    try:
        df = pd.read_csv(filename, nrows=0) # 只读表头
        # 假设第一列是基因名，其余是模型ID
        return set(df.columns[1:])
    except Exception as e:
        print(f"读取 '{filename}' 的列名时出错: {e}")
        return None

# --- 主程序 ---
if __name__ == "__main__":
    
    preview_csv_structure(file_crispr_effect)
    preview_csv_structure(file_copy_number)
    preview_csv_structure(file_model_meta)

    print("\n\n" + "*"*25 + " 文件预览完成，开始处理数据 " + "*"*25 + "\n")

    print("开始从文件中提取模型ID...")

    # 根据你的文件预览，我们现在可以确定正确的读取方式
    # CRISPRGeneEffect.csv: 行是模型, 列是基因 -> 从行名(index)读取
    crispr_ids = get_model_ids_from_index(file_crispr_effect)
    
    # OmicsCNGeneWGS.csv: 行是模型, 列是基因 -> 从行名(index)读取
    cnv_ids = get_model_ids_from_index(file_copy_number)
    
    # Model.csv: 第一列是模型ID，且它就是行名 -> 从行名(index)读取
    meta_ids = get_model_ids_from_index(file_model_meta)

    if crispr_ids is None or cnv_ids is None or meta_ids is None:
        print("\n由于文件读取错误，脚本终止。请检查文件路径和文件内容。")
    else:
        print("\n正在计算交集...")
        common_ids = crispr_ids & cnv_ids & meta_ids

        print("-" * 40)
        print(f"CRISPR 数据中的模型数量: {len(crispr_ids)}")
        print(f"拷贝数 (WGS) 数据中的模型数量: {len(cnv_ids)}")
        print(f"元数据 (Model.csv) 中的模型数量: {len(meta_ids)}")
        print("-" * 40)
        print(f"最终共同模型的数量: {len(common_ids)}")
        print("-" * 40)

        if common_ids:
            common_ids_list = sorted(list(common_ids))
            with open(output_file, 'w') as f:
                for model_id in common_ids_list:
                    f.write(model_id + '\n')
            print(f"成功！已将 {len(common_ids)} 个共同模型ID保存到文件: '{output_file}'")
        else:
            print("警告: 未找到任何共同的模型ID。")