# plot_and_ari.R
# powered by LiquidMonsteR

library(Seurat)
library(ggplot2)
library(clustree)  # For calculating ARI
library(mclust)  # For ARI
library(gridExtra)
library(grid)

plot_umap <- function(seurat_obj, reduction_method, clustering_method) {

  reduction_name <- switch(as.character(reduction_method),
                          `1` = "pca",
                          `2` = "cica",
                          `3` = "scvi",
                          `4` = "sca",
                          stop("Invalid method"))

  seurat_obj <- RunUMAP(seurat_obj, reduction = reduction_name, dims = 1:10)

  # 获取降维和聚类方法的名称
  reduction_methods <- c("PCA", "c-ICA", "Autoencoder", "scVI", "SCA")
  clustering_methods <- c("Louvain", "Walktrap", "K-medoids", "Mclust")
  
  reduction_name <- reduction_methods[reduction_method]
  clustering_name <- clustering_methods[clustering_method]
  
  # 绘制UMAP图像，颜色为原始细胞类型
  p1 <- DimPlot(seurat_obj, reduction = "umap", group.by = "cell_type") +
    ggtitle(paste("UMAP -", reduction_name, "- Manual Annotation")) +
    theme(plot.title = element_text(hjust = 0.5))
  
  # 绘制UMAP图像，颜色为聚类结果
  p2 <- DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters") +
    ggtitle(paste("UMAP -", reduction_name, "+", clustering_name)) +
    theme(plot.title = element_text(hjust = 0.5))
  
  # 计算ARI指标
  true_labels <- seurat_obj$cell_type
  predicted_labels <- seurat_obj$seurat_clusters
  ari <- adjustedRandIndex(true_labels, predicted_labels)
  
  # 在图像上标注ARI指标
  p2 <- p2 + annotate("text", x = Inf, y = Inf, label = paste("ARI:", round(ari, 3)), hjust = 1.1, vjust = 1.1, size = 5, color = "red")
  
  return(list(p1, p2))
}
