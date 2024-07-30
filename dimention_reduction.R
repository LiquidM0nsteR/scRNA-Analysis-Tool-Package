# dimensionality_reduction.R
# powered by LiquidMonsteR
library(Seurat)
library(ica)  # For c-ICA
library(reticulate)  # For Python interaction


perform_dimensionality_reduction <- function(seurat_obj, method) {
  # 预处理数据
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj)
  
  if (method == 1) {
    seurat_obj <- run_pca(seurat_obj)
  } else if (method == 2) {
    seurat_obj <- run_cica(seurat_obj)
  } else if (method == 3) {
    seurat_obj <- run_scvi(seurat_obj)
  } else if (method == 4) {
    seurat_obj <- run_sca(seurat_obj)
  } else {
    stop("Invalid reduction method. Choose a number between 1 and 5.")
  }
  
  seurat_obj <- standardize_reduction(seurat_obj, method)
  return(seurat_obj)
}

run_pca <- function(seurat_obj) {
  seurat_obj <- RunPCA(seurat_obj, npcs = 30)
  embeddings <- seurat_obj[["pca"]]@cell.embeddings
  colnames(embeddings) <- paste0("reduction_", 1:ncol(embeddings))
  seurat_obj[["pca"]] <- CreateDimReducObject(embeddings = embeddings, key = "pca_")
  return(seurat_obj)
}

run_cica <- function(seurat_obj) {
  DefaultAssay(seurat_obj) <- "RNA"  # 确保设置默认的Assay
  features <- VariableFeatures(seurat_obj)
  gene_expr_matrix <- GetAssayData(seurat_obj, slot = "data")[features, ]
  print(paste0("size of gene_expr_matrix is ", nrow(gene_expr_matrix), "*", ncol(gene_expr_matrix)))
  cica_result <- icafast(t(gene_expr_matrix), nc = 30)  # 假设我们选择提取30个独立成分
  embeddings <- cica_result$S
  if (nrow(embeddings) < ncol(embeddings)) {
    embeddings <- t(embeddings)
  }
  print(paste0("size of embeddings is ", nrow(embeddings), "*", ncol(embeddings)))
  rownames(embeddings) <- colnames(gene_expr_matrix)
  colnames(embeddings) <- paste0("cica_", 1:ncol(embeddings))
  seurat_obj[["cica"]] <- CreateDimReducObject(embeddings = embeddings, key = "cica_")
  return(seurat_obj)
}


run_scvi <- function(seurat_obj) {
  use_condaenv("rkit", required = TRUE)
  features <- VariableFeatures(seurat_obj)
  # Export Seurat object data to a file
  expr_data <- t(seurat_obj@assays$RNA$count[features, ])
  meta_data <- seurat_obj@meta.data
  write.csv(expr_data, file = "./cache/expr_data.csv", row.names = TRUE)
  write.csv(meta_data, file = "./cache/meta_data.csv", row.names = TRUE)

  print(paste0("size of gene_expr_matrix is ", nrow(expr_data), "*", ncol(expr_data)))

  py_run_string("

import os

import sys
print(sys.executable)

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
if 'stage' in adata.obs:
    scvi.model.SCVI.setup_anndata(adata, batch_key='stage')
else:
    scvi.model.SCVI.setup_anndata(adata)

# Create and train the scVI model
vae = scvi.model.SCVI(adata, n_latent=24)

# Use parallel computing to speed up training
with parallel_backend('threading', n_jobs = max(1, num_cores - 2)):  # Adjust to the actual number of available CPU cores
    vae.train()

# Get the latent representation
latent = vae.get_latent_representation()

# Save the latent representation
pd.DataFrame(latent, index=adata.obs.index).to_csv('./cache/latent.csv')

")
  

  # 提取编码数据
  latent_representation <- as.matrix(read.csv("./cache/latent.csv", row.names = 1))

  if (nrow(latent_representation) < ncol(latent_representation)) {
    latent_representation <- t(latent_representation)
  }

  print(paste0("size of latent is ", nrow(latent_representation), "*", ncol(latent_representation)))

  rownames(latent_representation) <- rownames(expr_data)
  colnames(latent_representation) <- paste0("scvi_", 1:ncol(latent_representation))
  
  seurat_obj[["scvi"]] <- CreateDimReducObject(embeddings = latent_representation, key = "scvi_")
  return(seurat_obj)
}



run_sca <- function(seurat_obj) {
  # 确保当前工作目录下有一个缓存文件夹
  if (!dir.exists("./cache")) {
    dir.create("./cache")
  }

  # 获取Variable Features
  features <- VariableFeatures(seurat_obj)
  gene_expr_matrix <- GetAssayData(seurat_obj, slot = "data")[features, ]

  # 打印基因表达矩阵的大小
  print(paste0("size of gene_expr_matrix is ", nrow(gene_expr_matrix), "*", ncol(gene_expr_matrix)))

  # 将数据写入临时文件
  expr_data <- t(gene_expr_matrix)
  meta_data <- seurat_obj@meta.data
  write.csv(expr_data, file = "./cache/expr_data.csv", row.names = TRUE)
  write.csv(meta_data, file = "./cache/meta_data.csv", row.names = TRUE)

  # 运行Python代码
  py_run_string("
import scanpy as sc
from shannonca.dimred import reduce_scanpy
import pandas as pd

# 读取数据
expr_data = pd.read_csv('./cache/expr_data.csv', index_col=0)
meta_data = pd.read_csv('./cache/meta_data.csv', index_col=0)

# 创建 AnnData 对象
adata = sc.AnnData(X=expr_data.values, obs=meta_data)

# 执行 SCA 降维
reduce_scanpy(adata, keep_scores=True, keep_loadings=True, keep_all_iters=True, layer=None, key_added='sca', iters=1, n_comps=30)

# 保存降维结果
adata.obsm['X_sca'].to_csv('./cache/sca_embeddings.csv')
")

  # 读取降维结果
  sca_embeddings <- read.csv("./cache/sca_embeddings.csv", row.names = 1)
  if (nrow(sca_embeddings) < ncol(sca_embeddings)) {
    sca_embeddings <- t(sca_embeddings)
  }
  rownames(sca_embeddings) <- colnames(gene_expr_matrix)
  colnames(sca_embeddings) <- paste0("sca_", 1:ncol(sca_embeddings))

  # 打印降维结果的大小
  print(paste0("size of embeddings is ", nrow(sca_embeddings), "*", ncol(sca_embeddings)))

  # 将降维结果保存到Seurat对象中
  seurat_obj[["sca"]] <- CreateDimReducObject(embeddings = sca_embeddings, key = "sca_")
  return(seurat_obj)
}



standardize_reduction <- function(seurat_obj, method) {
  reduction_name <- switch(method,
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