library(GEOquery)

# 读取 SOFT 格式文件 (.soft.gz)
read_soft_file <- function(file_path) {
  # 直接读取压缩的 SOFT 文件
  gse <- getGEO(filename = file_path, GSEMatrix = TRUE)
  
  # 修正：正确处理多平台数据集
  if (is.list(gse)) {
    # 如果有多个平台，取第一个（通常只有一个）
    gse <- gse[[1]]
  }
  
  # 提取元数据 - 使用更可靠的方法
  metadata <- pData(gse)
  
  # 提取表达矩阵
  expression_matrix <- exprs(gse)
  
  # 提取特征数据（基因注释）
  feature_data <- fData(gse)
  
  return(list(
    metadata = metadata,
    expression = expression_matrix,
    features = feature_data
  ))
}

# 读取 MINiML 格式文件 (.xml.tgz)
read_miniml_file <- function(file_path) {
  # 直接读取压缩的 MINiML 文件
  gse <- getGEO(filename = file_path, GSEMatrix = TRUE)
  
  # 修正：正确处理多平台数据集
  if (is.list(gse)) {
    gse <- gse[[1]]
  }
  
  # 提取实验设计信息
  experiment_data <- experimentData(gse)
  
  # 提取样本处理协议
  protocols <- protocols(gse)
  
  # 提取样本元数据
  sample_metadata <- pData(gse)
  
  return(list(
    experiment = experiment_data@other,
    protocols = protocols,
    sample_metadata = sample_metadata
  ))
}

# 使用函数处理文件
# 处理 SOFT 文件
soft_data <- read_soft_file("GSE91061_family.soft.gz")

# 查看元数据
head(soft_data$metadata)

# 查看表达矩阵前几行
head(soft_data$expression[, 1:5])

# 查看基因注释
head(soft_data$features)

# 处理 MINiML 文件
miniml_data <- read_miniml_file("GSE91061_family.xml.tgz")

# 查看实验设计
str(miniml_data$experiment)

# 查看样本处理协议
miniml_data$protocols

# 查看样本元数据
head(miniml_data$sample_metadata)