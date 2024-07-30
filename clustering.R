# clustering.R
# powered by LiquidMonsteR
library(Seurat)
library(cluster)  # For K-medoids
library(mclust)  # For Model-Based Clustering
library(igraph)  # For Walktrap
library(scran) # For Walktrap


perform_clustering <- function(seurat_obj, method, reduction_method) {
  # 映射 method 到降维方法的名称
  reduction_name <- switch(as.character(reduction_method),
                           `1` = "pca",
                           `2` = "cica",
                           `3` = "scvi",
                           `4` = "sca",
                           stop("Invalid method"))
  
  if (method == 1) {
    seurat_obj <- run_findclusters(seurat_obj, reduction_name)
  } else if (method == 2) {
    seurat_obj <- run_walktrap(seurat_obj, reduction_name)
  } else if (method == 3) {
    seurat_obj <- run_kmedoids(seurat_obj, reduction_name)
  } else if (method == 4) {
    seurat_obj <- run_model_based_clustering(seurat_obj, reduction_name)
  } else {
    stop("Invalid clustering method. Choose a number between 1 and 4.")
  }
  
  return(seurat_obj)
}

run_findclusters <- function(seurat_obj, reduction_method) {
  seurat_obj <- FindNeighbors(seurat_obj, reduction = reduction_method, dims = 1:20)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
  return(seurat_obj)
}

run_walktrap <- function(seurat_obj, reduction_method) {
  seurat_obj <- FindNeighbors(seurat_obj, reduction = reduction_method, dims = 1:20)
  # 提取邻居图
  graph_matrix <- as.matrix(seurat_obj@graphs[[1]])
  
  # 将邻接矩阵转换为 igraph 对象
  graph <- igraph::graph_from_adjacency_matrix(graph_matrix, mode = "undirected", weighted = TRUE)
  
  clusters <- cluster_walktrap(graph)
  seurat_obj$seurat_clusters <- factor(clusters$membership)
  return(seurat_obj)
}

run_kmedoids <- function(seurat_obj, reduction_method) {
  # 获取cell_type种类的数量
  num_clusters <- length(unique(seurat_obj@meta.data$cell_type))
  
  # 使用给定的降维方法获取嵌入
  data <- Embeddings(seurat_obj, reduction = reduction_method)
  
  # 运行K-medoids聚类
  kmedoids_result <- pam(data, k = num_clusters)  # 使用cell_type种类的数量作为簇的数量
  seurat_obj$seurat_clusters <- factor(kmedoids_result$clustering)
  
  return(seurat_obj)
}


run_model_based_clustering <- function(seurat_obj, reduction_method) {
  data <- Embeddings(seurat_obj, reduction = reduction_method)
  mclust_result <- Mclust(data)
  seurat_obj$seurat_clusters <- factor(mclust_result$classification)
  return(seurat_obj)
}
