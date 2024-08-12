# preprocess of seurat boj
# powered by LiquidMonsteR
library(Seurat)
library(stringdist)


init_seurat <- function(seurat_from_rds){
    if(seurat_from_rds){
      print("请指定rds文件的路径: ")
      rds_file_path <- readLines(con = "stdin", n = 1)
      seurat_obj <- readRDS(rds_file_path)
      seurat_obj <- NormalizeData(seurat_obj)
      seurat_obj <- FindVariableFeatures(seurat_obj)
      seurat_obj <- ScaleData(seurat_obj)
    }else{
      print("将使用 filtered_feature_bc_matrix.h5 文件生成Seurat对象")
      # 设置文件目录
      data_dir <- "./init_data/"

      # 获取所有 .h5 文件的路径
      h5_files <- list.files(data_dir, pattern = "\\.h5$", full.names = TRUE)

      # 检查是否有 .h5 文件
      if (length(h5_files) == 0) {
          print("数据缺失：没有找到任何 .h5 文件。")
      } else {
          print(paste("检测到", length(h5_files), "个 .h5 文件。"))
          
          # 遍历每个 .h5 文件并创建 Seurat 对象
          seurat_list <- list()
          for (i in seq_along(h5_files)) {
              # 读取10X Genomics数据
              data <- Read10X_h5(h5_files[i])
              
              # 创建Seurat对象
              seurat_temp <- CreateSeuratObject(counts = data, project = paste0("project", i))
              
              # 添加 group 元数据
              seurat_temp$group <- paste0("group", i)
              
              # 标准化表达矩阵
              seurat_temp <- NormalizeData(seurat_temp)

              # 将Seurat对象添加到列表中
              seurat_list[[i]] <- seurat_temp
          }
          
          # 合并所有Seurat对象
          seurat_obj <- merge(x = seurat_list[[1]], y = seurat_list[-1])
          
      }
    }
    # saveRDS(seurat_obj, "./init_data/seurat_obj.rds")

    return (check_and_preprocess_seurat(seurat_obj))
}


# 检查和预处理 Seurat 对象的函数
check_and_preprocess_seurat <- function(seurat_obj) {
  # 检查是否提供了 Seurat 对象
  if (missing(seurat_obj) || !inherits(seurat_obj, "Seurat")) {
    stop("请提供一个有效的 Seurat 对象。")
  }
  
  # 将默认 assay 设置为 RNA
  DefaultAssay(seurat_obj) <- "RNA"
  
  # 初始需要的列名
  required_cols <- c("cell_type", "x", "y", "orig.ident", "group", "nFeature_RNA", "nCount_RNA")
  
  # 模糊匹配函数，用于找到最接近的列名
  match_columns <- function(meta_data, required_cols) {
    matched_cols <- sapply(required_cols, function(req_col) {
      col_distances <- stringdist::stringdistmatrix(req_col, colnames(meta_data), method = "jw")
      closest_col <- colnames(meta_data)[which.min(col_distances)]
      # 如果最小距离大于某个阈值（如0.5），则认为没有匹配到
      if (min(col_distances) > 0.5) {
        return(NA)
      }
      return(closest_col)
    })
    names(matched_cols) <- required_cols
    return(matched_cols)
  }
  
  # 使用模糊匹配匹配列名
  matched_cols <- match_columns(seurat_obj@meta.data, required_cols)
  
  # 确定是否存在空间坐标
  use_spatial_coords <- "x" %in% matched_cols && "y" %in% matched_cols
  
  # 确定是否为多样本
  is_multiple <- "orig.ident" %in% matched_cols || "group" %in% matched_cols
  
  # 根据检测结果过滤所需的列名
  final_required_cols <- c()
  if (use_spatial_coords) {
    final_required_cols <- c(final_required_cols, "x", "y")
  }
  if (is_multiple) {
    final_required_cols <- c(final_required_cols, "orig.ident", "group")
  }
  
  # 将匹配到的列名重命名为标准名称
  for (req_col in final_required_cols) {
    col_name <- matched_cols[req_col]
    if (!req_col %in% colnames(seurat_obj@meta.data)) {
      seurat_obj@meta.data[[req_col]] <- seurat_obj@meta.data[[col_name]]
      seurat_obj@meta.data[[col_name]] <- NULL
    }
  }
  
  # 检查是否所有必需的列名都存在
  if (!all(final_required_cols %in% colnames(seurat_obj@meta.data))) {
    stop("Seurat 对象的 meta.data 中缺少一个或多个必需的列名: ", 
         paste(setdiff(final_required_cols, colnames(seurat_obj@meta.data)), collapse = ", "))
  }
  
  if (use_spatial_coords) {
    # 检查基因表达矩阵和空间位置矩阵的维度
    gene_expr_matrix <- seurat_obj@assays$RNA$counts
    spatial_matrix <- seurat_obj@meta.data[, c("x", "y")]
    
    if (ncol(gene_expr_matrix) != nrow(spatial_matrix)) {
      stop("基因表达矩阵中的细胞数量与空间位置矩阵中的细胞数量不匹配。")
    }
    
    if (!all(colnames(gene_expr_matrix) == rownames(spatial_matrix))) {
      stop("空间位置矩阵中的细胞顺序与基因表达矩阵中的顺序不匹配。")
    }
  }
  
  return(list(seurat_obj = seurat_obj, is_multiple = is_multiple, use_spatial_coords = use_spatial_coords))
}
