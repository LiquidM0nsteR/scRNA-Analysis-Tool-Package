# clustering.R
# powered by LiquidMonsteR
library(Seurat)
library(cluster)  # For K-medoids
library(mclust)  # For Model-Based Clustering
library(igraph)  # For Walktrap

perform_clustering <- function(seurat_obj, method) {
  if (method == 1) {
    seurat_obj <- run_findclusters(seurat_obj)
  } else if (method == 2) {
    seurat_obj <- run_walktrap(seurat_obj)
  } else if (method == 3) {
    seurat_obj <- run_kmedoids(seurat_obj)
  } else if (method == 4) {
    seurat_obj <- run_model_based_clustering(seurat_obj)
  } else {
    stop("Invalid clustering method. Choose a number between 1 and 4.")
  }
  
  return(seurat_obj)
}

run_findclusters <- function(seurat_obj) {
  seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", dims = 1:10)
  seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
  return(seurat_obj)
}

run_walktrap <- function(seurat_obj) {
  graph <- buildSNNGraph(seurat_obj, use.dimred = "pca", k = 20)
  clusters <- cluster_walktrap(graph)
  seurat_obj$seurat_clusters <- factor(clusters$membership)
  return(seurat_obj)
}

run_kmedoids <- function(seurat_obj) {
  data <- Embeddings(seurat_obj, reduction = "pca")
  kmedoids_result <- pam(data, k = 5)  # 假设选择5个簇
  seurat_obj$seurat_clusters <- factor(kmedoids_result$clustering)
  return(seurat_obj)
}


run_model_based_clustering <- function(seurat_obj) {
  data <- Embeddings(seurat_obj, reduction = "pca")
  mclust_result <- Mclust(data)
  seurat_obj$seurat_clusters <- factor(mclust_result$classification)
  return(seurat_obj)
}

# powered by LiquidMonsteR