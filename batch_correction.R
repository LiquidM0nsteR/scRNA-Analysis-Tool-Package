# batch_correction.R
# Powered by LiquidMonsteR
library(Seurat)
library(harmony)
library(sva)



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

  seurat_obj[["integrated"]] <- seurat_integrated@assays$integrated

  print("去批次完成.")

  return(seurat_obj)
}



# 封装Harmony去批次方法
run_harmony <- function(seurat_obj, batch_key, reduction_method, num_components) {

  reduction_name <- switch(reduction_method,
                        `1` = "PCA",
                        `2` = "cica",
                        `3` = "scvi",
                        `4` = 'SpatialPCA',
                        `5` = 'SeuratSpatial',
                        stop("Invalid method"))

  seurat_obj <- RunHarmony(seurat_obj, batch_key, reduction_name)

  print("去批次完成.")

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

  DefaultAssay(seurat_obj) <- "combat"

  print("去批次完成.")

  return(seurat_obj)
}



