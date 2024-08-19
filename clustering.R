# clustering.R
# powered by LiquidMonsteR
library(Seurat)
library(cluster)  # For K-medoids
library(mclust)  # For Model-Based Clustering
library(igraph)  # For Walktrap
library(scran) # For Walktrap


perform_clustering <- function(seurat_obj, method, reduction_method, num_components, is_multiple) {
  # 映射 method 到降维方法的名称
  reduction_name <- switch(reduction_method,
                           `1` = "pca",
                           `2` = "cica",
                           `3` = "scvi",
                           stop("Invalid method"))
  
  if (method == 1) {
    seurat_obj <- run_findclusters(seurat_obj, reduction_name)
  } else if (method == 2) {
    seurat_obj <- run_walktrap(seurat_obj, reduction_name)
  } else if (method == 3) {
    seurat_obj <- run_model_based_clustering(seurat_obj, reduction_name)
  } else {
    stop("Invalid clustering method. Choose a number between 1 and 4.")
  }

  seurat_obj <- plot_umap(seurat_obj, reduction_method, method, num_components, is_multiple)
  
  return(seurat_obj)
}

run_findclusters <- function(seurat_obj, reduction_method) {

  neighbor_graph_name <- "my_custom_nn"
  snn_graph_name <- "my_custom_snn"
  # 使用指定的邻居图名称运行 FindNeighbors
  seurat_obj <- FindNeighbors(
    seurat_obj, 
    reduction = reduction_method, 
    dims = 1:10, 
    k.param = 10, 
    graph.name = c(neighbor_graph_name, snn_graph_name)
  )

  # 使用 FindClusters 调用指定的 SNN 图
  seurat_obj <- FindClusters(
    seurat_obj, 
    graph.name = snn_graph_name, 
    resolution = 0.3
  )
  return(seurat_obj)
}

run_walktrap <- function(seurat_obj, reduction_method) {
  seurat_obj <- FindNeighbors(seurat_obj, reduction = reduction_method)
  # 提取邻居图
  graph_matrix <- as.matrix(seurat_obj@graphs[[1]])
  
  # 将邻接矩阵转换为 igraph 对象
  graph <- igraph::graph_from_adjacency_matrix(graph_matrix, mode = "undirected", weighted = TRUE)
  
  clusters <- cluster_walktrap(graph)
  seurat_obj$seurat_clusters <- factor(clusters$membership)
  return(seurat_obj)
}


run_model_based_clustering <- function(seurat_obj, reduction_method) {
  data <- Embeddings(seurat_obj, reduction = reduction_method)
  mclust_result <- Mclust(data)
  seurat_obj$seurat_clusters <- factor(mclust_result$classification)
  return(seurat_obj)
}



plot_umap <- function(seurat_obj, reduction_method, clustering_method, num_components, is_multiple){

  reduction_name <- switch(as.character(reduction_method),
                      `1` = "pca",
                      `2` = "cica",
                      `3` = "scvi",
                      stop("Invalid method"))

  seurat_obj <- Seurat::RunUMAP(seurat_obj, reduction = reduction_name, dims = 1:num_components)
  seurat_obj <- Seurat::RunTSNE(seurat_obj, reduction = reduction_name, dims = 1:num_components)


  reduction_methods <- c("PCA", "c-ICA", "scVI")
  clustering_methods <- c("Louvain", "Walktrap", "Mclust")
  reduction_name <- reduction_methods[reduction_method]
  clustering_name <- clustering_methods[clustering_method]

  # 绘制UMAP图像，颜色为聚类结果，并将图例设置为两列
  p1 <- DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters") +
    ggtitle(paste("UMAP -", reduction_name, "+", clustering_name)) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20)) +
    guides(color = guide_legend(ncol = 2, override.aes = list(size = 4)))


  if(is_multiple){
    # 绘制UMAP图像，颜色为不同批次，并将图例设置为两列
    p2 <- DimPlot(seurat_obj, reduction = "umap", group.by = "orig.ident") +
      ggtitle(paste("UMAP -", reduction_name, "+", "Colored by batch")) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20)) +
      guides(color = guide_legend(override.aes = list(size = 4)))
  }


  # 聚类图转换成pdf
  pdf("./report3.pdf")
    if(is_multiple){print(p2)}
    print(p1)
  dev.off()

  print("聚类完成.")

  # saveRDS(seurat_obj, file = "./init_data/seurat_obj.rds")

  return(seurat_obj)
}