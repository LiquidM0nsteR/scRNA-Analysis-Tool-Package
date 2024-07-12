# dimensionality_reduction.R
# powered by LiquidMonsteR
library(Seurat)
library(ica)  # For c-ICA
library(autoencoder)  # For Autoencoder
library(reticulate)  # For Python interaction
library(scater)  # For SCA


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
    seurat_obj <- run_autoencoder(seurat_obj)
  } else if (method == 4) {
    seurat_obj <- run_scvi(seurat_obj)
  } else if (method == 5) {
    seurat_obj <- run_sca(seurat_obj)
  } else {
    stop("Invalid reduction method. Choose a number between 1 and 5.")
  }
  
  seurat_obj <- standardize_reduction(seurat_obj)
  return(seurat_obj)
}

run_pca <- function(seurat_obj) {
  seurat_obj <- RunPCA(seurat_obj, npcs = 30)
  embeddings <- seurat_obj[["pca"]]@cell.embeddings
  colnames(embeddings) <- paste0("reduction_", 1:ncol(embeddings))
  seurat_obj[["reduction"]] <- CreateDimReducObject(embeddings = embeddings, key = "reduction_")
  return(seurat_obj)
}

run_cica <- function(seurat_obj) {
  gene_expr_matrix <- GetAssayData(seurat_obj, slot = "data")
  cica_result <- icafast(gene_expr_matrix, nc = 30)  # 假设我们选择提取30个独立成分
  embeddings <- cica_result$S
  if (nrow(embeddings) > ncol(embeddings)) {
    embeddings <- t(embeddings)
  }
  colnames(embeddings) <- colnames(gene_expr_matrix)
  rownames(embeddings) <- paste0("reduction_", 1:nrow(embeddings))
  seurat_obj[["reduction"]] <- CreateDimReducObject(embeddings = embeddings, key = "reduction_")
  return(seurat_obj)
}

run_autoencoder <- function(seurat_obj) {
  gene_expr_matrix <- GetAssayData(seurat_obj, slot = "data")
  ae_model <- autoencoder.train(gene_expr_matrix, hidden.structure = c(128, 64, 32), learning.rate = 0.01)
  encoded_data <- autoencoder.encode(ae_model, gene_expr_matrix)
  if (nrow(encoded_data) > ncol(encoded_data)) {
    encoded_data <- t(encoded_data)
  }
  colnames(encoded_data) <- colnames(gene_expr_matrix)
  rownames(encoded_data) <- paste0("reduction_", 1:nrow(encoded_data))
  seurat_obj[["reduction"]] <- CreateDimReducObject(embeddings = encoded_data, key = "reduction_")
  return(seurat_obj)
}

run_scvi <- function(seurat_obj) {
  # 安装并加载scvi-tools包
  py_install("scvi-tools")
  # 将数据转换为anndata对象
  gene_expr_matrix <- GetAssayData(seurat_obj, slot = "data")
  anndata <- import("anndata", convert = FALSE)
  scvi <- import("scvi", convert = FALSE)
  adata <- anndata$AnnData(X = t(gene_expr_matrix))
  
  # 训练scVI模型
  scvi$model$SCVI$setup_anndata(adata)
  vae <- scvi$model$SCVI(adata)
  vae$train()
  
  # 提取编码数据
  latent_representation <- vae$get_latent_representation()
  if (nrow(latent_representation) > ncol(latent_representation)) {
    latent_representation <- t(latent_representation)
  }
  rownames(latent_representation) <- colnames(gene_expr_matrix)
  colnames(latent_representation) <- paste0("reduction_", 1:ncol(latent_representation))
  
  seurat_obj[["reduction"]] <- CreateDimReducObject(embeddings = latent_representation, key = "reduction_")
  return(seurat_obj)
}

run_sca <- function(seurat_obj) {
  gene_expr_matrix <- GetAssayData(seurat_obj, slot = "data")
  sca_result <- sca(gene_expr_matrix, k = 30)  # 假设我们选择提取30个Surprisal成分
  embeddings <- sca_result$components
  if (nrow(embeddings) > ncol(embeddings)) {
    embeddings <- t(embeddings)
  }
  colnames(embeddings) <- colnames(gene_expr_matrix)
  rownames(embeddings) <- paste0("reduction_", 1:nrow(embeddings))
  seurat_obj[["reduction"]] <- CreateDimReducObject(embeddings = embeddings, key = "reduction_")
  return(seurat_obj)
}

standardize_reduction <- function(seurat_obj) {
  reduction <- seurat_obj[["reduction"]]
  
  # 获取降维结果的细胞名和维度名
  cell_names <- colnames(reduction@cell.embeddings)
  dim_names <- rownames(reduction@cell.embeddings)
  
  # 检查并设置行名（维度名）和列名（细胞名）
  if (is.null(dim_names)) {
    num_dims <- ncol(reduction@cell.embeddings)
    dim_names <- paste0("reduction_", 1:num_dims)
    rownames(reduction@cell.embeddings) <- dim_names
  }
  
  if (any(colnames(reduction@cell.embeddings) != cell_names)) {
    reduction@cell.embeddings <- reduction@cell.embeddings[, cell_names, drop = FALSE]
  }
  
  # 统一命名reduction和key
  reduction@key <- "reduction_"
  seurat_obj[["reduction"]] <- reduction
  
  return(seurat_obj)
}
