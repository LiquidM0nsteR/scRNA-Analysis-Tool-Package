# batch_correction.R
# Powered by LiquidMonsteR
library(Seurat)
library(harmony)
library(sva)



# 主函数：调用上述不同的去批次方法
run_batch_correction <- function(seurat_obj, batch_correction_method) {
  
  print("正在进行批次校正. ")

  batch_key <- NULL
  if ("orig.ident" %in% colnames(seurat_obj@meta.data) && length(unique(seurat_obj@meta.data$orig.ident)) > 1) {
    batch_key <- "orig.ident"
  }

  # 如果 orig.ident 不满足多样本条件，再检查 group 列
  if (!is_multiple && "group" %in% colnames(seurat_obj@meta.data) && length(unique(seurat_obj@meta.data$group)) > 1) {
    batch_key <- "group"
  }

  
  # 根据 batch_correction_method 的值调用相应的去批次方法
  seurat_obj <- switch(
    batch_correction_method,
    `1` = run_seurat_integration(seurat_obj, batch_key),  # 1 对应 Seurat Integration
    `2` = run_harmony(seurat_obj, batch_key),             # 2 对应 Harmony
    `3` = run_combat(seurat_obj, batch_key),              # 3 对应 Combat
    `4` = run_cca(seurat_obj, batch_key),                 # 4 对应 CCA
  )
  
  print("批次校正完成. ")

  return(seurat_obj)
}


# 封装Seurat标准去批次流程
run_seurat_integration <- function(seurat_obj, batch_key) {
  seurat_list <- SplitObject(seurat_obj, split.by = batch_key)
  seurat_list <- lapply(seurat_list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    return(x)
  })
  
  anchors <- FindIntegrationAnchors(object.list = seurat_list, dims = 1:30)
  seurat_integrated <- IntegrateData(anchorset = anchors, dims = 1:30)
  
  DefaultAssay(seurat_integrated) <- "integrated"
  seurat_integrated <- ScaleData(seurat_integrated)
  
  # 不进行PCA降维，直接返回Seurat对象
  seurat_obj[["integrated"]] <- seurat_integrated@assays$integrated
  return(seurat_obj)
}

# 封装Harmony去批次方法
run_harmony <- function(seurat_obj, batch_key) {
  seurat_obj <- RunHarmony(seurat_obj, group.by.vars = batch_key)
  
  # 不进行PCA降维，直接返回Seurat对象
  seurat_obj@assays$harmony <- CreateAssayObject(counts = seurat_obj@reductions$harmony@cell.embeddings)
  return(seurat_obj)
}

# 封装Combat去批次方法
run_combat <- function(seurat_obj, batch_key) {
  # 获取表达矩阵
  expr_matrix <- as.matrix(GetAssayData(seurat_obj, slot = "data"))
  batch <- seurat_obj@meta.data[[batch_key]]
  
  # Combat批次校正
  combat_corrected <- ComBat(dat = expr_matrix, batch = batch)
  
  # 将校正后的数据放回Seurat对象
  seurat_obj[["combat"]] <- CreateAssayObject(counts = combat_corrected)
  return(seurat_obj)
}

# 封装CCA去批次方法
run_cca <- function(seurat_obj, batch_key) {
  seurat_list <- SplitObject(seurat_obj, split.by = batch_key)
  seurat_list <- lapply(seurat_list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
    x <- ScaleData(x)
    return(x)
  })
  
  seurat_obj <- RunCCA(object1 = seurat_list[[1]], object2 = seurat_list[[2]], num.cc = 30)
  # 不进行对齐和降维，直接返回Seurat对象
  seurat_obj@assays$cca <- CreateAssayObject(counts = seurat_obj@reductions$cca@cell.embeddings)
  return(seurat_obj)
}



