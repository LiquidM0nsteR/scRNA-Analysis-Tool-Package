# dimensionality_reduction.R
# powered by LiquidMonsteR
library(Seurat)
library(ica)  # For c-ICA
library(reticulate)  # For Python interaction


perform_dimensionality_reduction <- function(seurat_obj, reduction_method = 1, batch_corre_method = NULL) {
  # 运行PCA
  seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)

  # 使用ElbowPlot寻找最佳的低维维度
  p1 <- ElbowPlot(seurat_obj) +
    ggtitle("Elbow Plot for PCA Components") +
    theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))


  # 保存Elbow Plot为PDF文件
  pdf("./report2.pdf")
  print(p1)
  dev.off()

  # 提取PCA结果
  pca_result <- seurat_obj[["pca"]]

  # 计算解释方差比例
  explained_variance <- cumsum(pca_result@stdev^2) / sum(pca_result@stdev^2)
  num_components <- which(explained_variance >= 0.90)[1]

  print(paste("至少解释90%方差的主成分数量:", num_components, "，设置为降维结果的维度"))

  # 根据方法选择降维算法
  if (reduction_method == 1) {
    seurat_obj <- run_pca(seurat_obj, num_components, batch_corre_method)
  } else if (reduction_method == 2) {
    seurat_obj <- run_cica(seurat_obj, num_components, batch_corre_method)
  } else if (reduction_method == 3) {
    seurat_obj <- run_scvi(seurat_obj, num_components)
  } else {
    stop("Invalid reduction method. Choose a number between 1 and 3.")
  }
  
  # 标准化降维结果
  seurat_obj <- standardize_reduction(seurat_obj, reduction_method)
  
  print("降维完成.")

  # 返回结果
  return(list(seurat_obj = seurat_obj, num_components = num_components))
}



run_pca <- function(seurat_obj, num_components, batch_corre_method) {
  if (is.null(batch_corre_method)) {
    # 单样本情况，直接使用原始数据进行PCA降维
    DefaultAssay(seurat_obj) <- "RNA"  # 假设原始数据存储在 "RNA" assay 中
  } else {
    # 多样本情况，根据batch_corre_method选择去批次后的数据
    assay_name <- switch(
      batch_corre_method,
      `1` = "integrated",  # Seurat Integration 的结果
      `2` = "harmony",     # Harmony 的结果
      `3` = "combat",      # Combat 的结果
      `4` = "cca",         # CCA 的结果
    )

    # 检查是否存在对应的去批次校正结果
    if (!assay_name %in% names(seurat_obj@assays)) {
      stop(paste("The assay", assay_name, "is not found in the Seurat object. Make sure you have run batch correction."))
    }

    # 设置默认assay为指定的批次校正结果
    DefaultAssay(seurat_obj) <- assay_name
  }
  
  # 使用设定的默认assay进行PCA降维
  seurat_obj <- RunPCA(seurat_obj, npcs = num_components, verbose = FALSE)
  embeddings <- seurat_obj[["pca"]]@cell.embeddings
  colnames(embeddings) <- paste0("reduction_", 1:ncol(embeddings))
  seurat_obj[["pca"]] <- CreateDimReducObject(embeddings = embeddings, key = "pca_")
  return(seurat_obj)
}


run_cica <- function(seurat_obj, num_components, batch_corre_method = NULL) {

  if (is.null(batch_corre_method)) {
    # 单样本情况，直接使用原始数据进行CICA降维
    DefaultAssay(seurat_obj) <- "RNA"  # 假设原始数据存储在 "RNA" assay 中
  } else {
    # 多样本情况，根据batch_corre_method选择去批次后的数据
    assay_name <- switch(
      batch_corre_method,
      `1` = "integrated",  # Seurat Integration 的结果
      `2` = "harmony",     # Harmony 的结果
      `3` = "combat",      # Combat 的结果
      `4` = "cca",         # CCA 的结果
    )

    # 检查是否存在对应的去批次校正结果
    if (!assay_name %in% names(seurat_obj@assays)) {
      stop(paste("The assay", assay_name, "is not found in the Seurat object. Make sure you have run batch correction."))
    }

    # 设置默认assay为指定的批次校正结果
    DefaultAssay(seurat_obj) <- assay_name
  }
  
  # 提取基因表达矩阵并运行CICA降维

  gene_expr_matrix <- GetAssayData(seurat_obj, slot = "data")
  print(paste0("基因表达矩阵的大小为： ", nrow(gene_expr_matrix), "*", ncol(gene_expr_matrix)))
  
  cica_result <- icafast(t(gene_expr_matrix), nc = num_components) 

  # 创建降维对象并存储在Seurat对象中
  embeddings <- cica_result$S
  if (nrow(embeddings) < ncol(embeddings)) {
    embeddings <- t(embeddings)
  }
  print(paste0("降维矩阵的大小为：", nrow(embeddings), "*", ncol(embeddings)))
  rownames(embeddings) <- colnames(gene_expr_matrix)
  colnames(embeddings) <- paste0("cica_", 1:ncol(embeddings))
  seurat_obj[["cica"]] <- CreateDimReducObject(embeddings = embeddings, key = "cica_")
  
  return(seurat_obj)
}


run_scvi <- function(seurat_obj, num_components) {

  use_python(Sys.which("python"), required = TRUE)

  features <- VariableFeatures(seurat_obj)
  # Export Seurat object data to a file
  expr_data <- t(seurat_obj@assays$RNA$counts[features, ])
  meta_data <- seurat_obj@meta.data

  write.csv(expr_data, file = "./cache/expr_data.csv", row.names = TRUE)
  write.csv(meta_data, file = "./cache/meta_data.csv", row.names = TRUE)

  print(paste0("基因表达矩阵的大小为： ", nrow(expr_data), "*", ncol(expr_data)))

  py_run_string(sprintf("

import os
import sys
import scvi
import scanpy as sc
import pandas as pd
import torch
from joblib import parallel_backend


# Set up multi-core CPU computation
num_cores = len(os.sched_getaffinity(0))
torch.set_num_threads(max(1, num_cores - 2))  # Adjust to the actual number of available CPU cores

# Load the expression data and metadata
expr_data = pd.read_csv('./cache/expr_data.csv', index_col=0)
meta_data = pd.read_csv('./cache/meta_data.csv', index_col=0)

# Create an AnnData object
adata = sc.AnnData(X=expr_data.values, obs=meta_data)

# Set up the data for scVI
if 'orig.ident' in adata.obs and adata.obs['orig.ident'].nunique() > 1:
    scvi.model.SCVI.setup_anndata(adata, batch_key='orig.ident')
elif 'group' in adata.obs and adata.obs['group'].nunique() > 1:
    scvi.model.SCVI.setup_anndata(adata, batch_key='group')
else:
    scvi.model.SCVI.setup_anndata(adata)


# Create and train the scVI model
vae = scvi.model.SCVI(adata, n_latent = %d)

# Use parallel computing to speed up training
with parallel_backend('threading', n_jobs = max(1, num_cores - 2)):  # Adjust to the actual number of available CPU cores
    vae.train(max_epochs = 400)

# Get the latent representation
latent = vae.get_latent_representation()

# Save the latent representation
pd.DataFrame(latent, index=adata.obs.index).to_csv('./cache/latent.csv')

", num_components))
  

  # 提取编码数据
  latent_representation <- as.matrix(read.csv("./cache/latent.csv", row.names = 1))

  if (nrow(latent_representation) < ncol(latent_representation)) {
    latent_representation <- t(latent_representation)
  }

  print(paste0("降维矩阵的大小为：", nrow(latent_representation), "*", ncol(latent_representation)))

  rownames(latent_representation) <- rownames(expr_data)
  colnames(latent_representation) <- paste0("scvi_", 1:ncol(latent_representation))
  
  seurat_obj[["scvi"]] <- CreateDimReducObject(embeddings = latent_representation, key = "scvi_")
  return(seurat_obj)
}




standardize_reduction <- function(seurat_obj, reduction_method) {
  reduction_name <- switch(reduction_method,
                           `1` = "pca",
                           `2` = "cica",
                           `3` = "scvi",
                           `4` = "sca")
  reduction <- seurat_obj[[reduction_name]]
  
  # 获取降维结果的细胞名和维度名
  cell_names <- colnames(reduction@cell.embeddings)
  dim_names <- rownames(reduction@cell.embeddings)
  
  # 检查并设置行名（维度名）和列名（细胞名）
  if (is.null(dim_names)) {
    num_dims <- ncol(reduction@cell.embeddings)
    dim_names <- paste0(reduction_name, "_", 1:num_dims)
    rownames(reduction@cell.embeddings) <- dim_names
  }
  
  if (any(colnames(reduction@cell.embeddings) != cell_names)) {
    reduction@cell.embeddings <- reduction@cell.embeddings[, cell_names, drop = FALSE]
  }
  
  # 统一命名reduction和key
  reduction@key <- paste0(reduction_name, "_")
  seurat_obj[[reduction_name]] <- reduction
  
  return(seurat_obj)
}

