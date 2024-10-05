# preprocess of seurat obj
# powered by LiquidMonsteR
library(Seurat)
library(stringdist)
library(glmGamPoi)


init_seurat <- function(seurat_from_rds, use_spatial_coords) {

  options(future.globals.maxSize = 500 * 1024^3) 
  if (seurat_from_rds) {
    print("请指定rds文件的路径: ")
    rds_file_path <- readLines(con = "stdin", n = 1)
    seurat_obj <- readRDS(rds_file_path)
  } else {
    # 设置文件目录
    data_dir <- "./init_data/"

    if (use_spatial_coords) {
      print("将使用 filtered_feature_bc_matrix.h5 文件 + spatial 文件夹生成 Seurat 对象. ")
      # 获取所有文件夹路径
      sample_dirs <- list.dirs(data_dir, recursive = FALSE)

      # 初始化一个列表来存储所有样本的 Seurat 对象
      seurat_objects <- list()

      # 遍历每个样本目录，读取数据并处理
      for (i in seq_along(sample_dirs)) {
        sample_dir <- sample_dirs[i]
        
        # 检查 .h5 文件
        h5_file_path <- file.path(sample_dir, "filtered_feature_bc_matrix.h5")
        if (!file.exists(h5_file_path)) {
          stop(paste("缺失文件：在", sample_dir, "中没有找到 .h5 文件。"))
        }

        # 读取10X Genomics数据
        seurat_temp <- Load10X_Spatial(data.dir = sample_dir)
        seurat_temp$orig.ident <- paste0("Sample", i)
        
        # 为 seurat_temp 的 meta.data 增加一列，列名可以自定义为 "SampleID"
        seurat_temp <- AddMetaData(seurat_temp, metadata = paste0("Sample", i), col.name = "group")

        names(seurat_temp@images) <- seurat_temp$orig.ident[1]
        seurat_temp@images[names(seurat_temp@images)]$key <- paste0(names(seurat_temp@images), "_")
        
        # 获取空间坐标并添加到 meta.data
        coords_temp <- GetTissueCoordinates(seurat_temp)
        seurat_temp <- AddMetaData(seurat_temp, metadata = coords_temp[, "x"], col.name = "x")
        seurat_temp <- AddMetaData(seurat_temp, metadata = coords_temp[, "y"], col.name = "y")

        # 将处理后的 Seurat 对象添加到列表中
        seurat_objects[[i]] <- seurat_temp
      }

      # 合并所有 Seurat 对象
      seurat_obj <- merge(seurat_objects[[1]], y = seurat_objects[-1], add.cell.ids = paste0("Sample", seq_along(seurat_objects)))

      # 获取 Seurat 对象中的所有层数据
      all_layers <- seurat_obj[["Spatial"]]@layers
      # 初始化合并后的 counts 矩阵
      merged_counts <- NULL
      # 遍历所有层并合并 counts 数据
      for (layer in all_layers) {
        # 检查当前层是否为 counts 数据层
        if (inherits(layer, "Matrix")) {
          # 合并 counts 矩阵
          if (is.null(merged_counts)) {
            merged_counts <- layer
          } else {
            merged_counts <- cbind(merged_counts, layer)
          }
        }
      }
      # 更新到 Seurat 对象的 counts 中
      seurat_obj[["Spatial"]]$counts <- merged_counts
      # 确定 counts 层的索引
      counts_index <- which(names(seurat_obj[["Spatial"]]@layers) == "counts")
      # 设置 counts 为默认的 layer
      seurat_obj[["Spatial"]]@default <- counts_index

      # 检查结果
      print(seurat_obj)

    } else {
      print("将使用 filtered_feature_bc_matrix.h5 文件生成 Seurat 对象")
      
      # 获取所有 .h5 文件的路径
      h5_files <- list.files(data_dir, pattern = "\\.h5$", full.names = TRUE)

      if (length(h5_files) == 0) {
        stop("数据缺失：没有找到任何 .h5 文件。")
      } else {
        print(paste("检测到", length(h5_files), "个 .h5 文件。"))

        # 初始化空的 counts 矩阵和元数据
        combined_counts <- NULL
        combined_metadata <- data.frame()
        all_cell_names <- c()  # 存储所有样本的细胞名称，用于检测重复

        # 先遍历样本，收集所有细胞名称
        for (i in seq_along(h5_files)) {
          data <- Read10X_h5(h5_files[i])
          seurat_temp <- CreateSeuratObject(counts = data, project = paste0("project", i))
          cell_names <- Cells(seurat_temp)
          all_cell_names <- c(all_cell_names, cell_names)
        }

        # 检测是否存在重复的细胞名称
        if (any(duplicated(all_cell_names))) {
          add_orig_ident <- TRUE
        } else {
          add_orig_ident <- FALSE
        }

        # 重新遍历样本，进行数据合并
        for (i in seq_along(h5_files)) {
          data <- Read10X_h5(h5_files[i])
          seurat_temp <- CreateSeuratObject(counts = data, project = paste0("project", i))
          seurat_temp$group <- paste0("group", i)

          # 根据检测结果判断是否修改细胞名称
          if (add_orig_ident) {
            cell_names <- Cells(seurat_temp)
            new_cell_names <- paste0(seurat_temp$group, "_", cell_names)
            seurat_temp <- RenameCells(seurat_temp, new.names = new_cell_names)
          }

          counts_matrix <- GetAssayData(seurat_temp, slot = "counts")

          if (is.null(combined_counts)) {
            combined_counts <- counts_matrix
          } else {
            combined_counts <- cbind(combined_counts, counts_matrix)
          }

          combined_metadata <- rbind(combined_metadata, seurat_temp@meta.data)
        }

        seurat_obj <- CreateSeuratObject(counts = combined_counts, meta.data = combined_metadata)
      }
    }
  }

  return(check_and_preprocess_seurat(seurat_obj, use_spatial_coords))
}



# 检查和预处理 Seurat 对象的函数
check_and_preprocess_seurat <- function(seurat_obj, use_spatial_coords) {
    # 检查是否提供了 Seurat 对象
    if (missing(seurat_obj) || !inherits(seurat_obj, "Seurat")) {
        stop("请提供一个有效的 Seurat 对象. ")
    }


    active_assay <- seurat_obj@assays[[DefaultAssay(seurat_obj)]]
    # 检查 counts 数据的维度
    if (all(dim(active_assay$counts) == 0)) {
        stop("Seurat 对象中的 counts 数据为空（0x0）. ")
    }

    # 初始需要的列名
    required_cols <- c("cell_type", "x", "y", "orig.ident", "group")
    
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
    
    batch_key <- NULL

    # 检查 orig.ident 列是否存在且包含多个不同的元素
    if ("orig.ident" %in% matched_cols && length(unique(seurat_obj@meta.data$orig.ident)) > 1) {
        batch_key <- "orig.ident"
    }

    # 如果 orig.ident 不满足多样本条件，再检查 group 列
    if (is.null(batch_key) && "group" %in% matched_cols && length(unique(seurat_obj@meta.data$group)) > 1) {
        batch_key <- "group"
    }
    
    # 根据检测结果过滤所需的列名
    final_required_cols <- c()
    if (use_spatial_coords) {
        final_required_cols <- c(final_required_cols, "x", "y")
    }
    if (!is.null(batch_key)) {
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
        gene_expr_matrix <- active_assay$counts
        spatial_matrix <- seurat_obj@meta.data[, c("x", "y")]
        
        if (ncol(gene_expr_matrix) != nrow(spatial_matrix)) {
            stop("基因表达矩阵中的细胞数量与空间位置矩阵中的细胞数量不匹配. ")
        }
        
        if (!all(colnames(gene_expr_matrix) == rownames(spatial_matrix))) {
            stop("空间位置矩阵中的细胞顺序与基因表达矩阵中的顺序不匹配. ")
        }
    }

    if(use_spatial_coords){
        # 创建一个新的 RNA assay，使用 Spatial 的计数数据
        seurat_obj[["RNA"]] <- CreateAssayObject(counts = active_assay$counts)
    }

    DefaultAssay(seurat_obj) <- "RNA"

    return(list(seurat_obj = seurat_obj, batch_key = batch_key))
}
