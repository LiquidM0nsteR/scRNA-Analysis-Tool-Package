# clustering.R
# powered by LiquidMonsteR
library(Seurat)
library(cluster)  # For K-medoids
library(mclust)  # For Model-Based Clustering
library(igraph)  # For Walktrap
library(scran) # For Walktrap


perform_clustering <- function(seurat_obj, clustering_method, reduction_method, num_components, use_spatial_coords, batch_key) {

    print("正在进行聚类. ")
    # 映射 method 到降维方法的名称
    reduction_name <- switch(reduction_method,
                            `1` = "PCA",
                            `2` = "cica",
                            `3` = "scvi",
                            `4` = "SpatialPCA",
                            `5` = "SeuratSpatial",
                            stop("不存在的降维方法. "))
    
    # 如果包含Harmony降维结果，就设置为harmony
    orig_reduction_name <- NULL
    if ("harmony" %in% names(seurat_obj@reductions)) {
        print("使用 harmony 结果进行聚类.")
        orig_reduction_name <- reduction_name
        reduction_name <- "harmony"
    }
    
    # 根据聚类方法，执行相应的聚类函数
    if (clustering_method == 1) {
      seurat_obj <- run_findclusters(seurat_obj, reduction_name)
    } else if (clustering_method == 2) {
      seurat_obj <- run_walktrap(seurat_obj, reduction_name)
    } else if (clustering_method == 3) {
      seurat_obj <- run_model_based_clustering(seurat_obj, reduction_name)
    } else {
      stop("Invalid clustering method. Choose a number between 1 and 3.")
    }

    # 根据聚类结果画umap图，同时也画出不同批次的umap图，观察去批次的效果
    if(!is.null(orig_reduction_name)){reduction_name <- orig_reduction_name}
    seurat_obj <- plot_umap(seurat_obj, reduction_name, clustering_method, num_components, use_spatial_coords, batch_key)
    
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



plot_umap <- function(seurat_obj, reduction_name, clustering_method, num_components, use_spatial_coords, batch_key){

    if(reduction_name == "SpatialPCA"){num_components = min(20, num_components)} # "SpatialPCA"结果默认20个维度
    seurat_obj <- Seurat::RunUMAP(seurat_obj, reduction = reduction_name, dims = 1:num_components, verbose = FALSE)
    seurat_obj <- Seurat::RunTSNE(seurat_obj, reduction = reduction_name, dims = 1:num_components, verbose = FALSE)


    clustering_methods <- c("Louvain", "Walktrap", "Mclust")
    clustering_name <- clustering_methods[clustering_method]


    # 绘制UMAP图像，颜色为聚类结果，并将图例设置为两列
    p1 <- DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters") +
      ggtitle(paste("UMAP -", reduction_name, "+", clustering_name)) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20)) +
      guides(color = guide_legend(ncol = 2, override.aes = list(size = 4)))

    # 绘制tSNE图像，颜色为聚类结果，并将图例设置为两列
    p2 <- DimPlot(seurat_obj, reduction = "tsne", group.by = "seurat_clusters") +
      ggtitle(paste("tSNE -", reduction_name, "+", clustering_name)) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20)) +
      guides(color = guide_legend(ncol = 2, override.aes = list(size = 4)))


    if(!is.null(batch_key)){
        # 绘制UMAP图像，颜色为不同批次，并将图例设置为两列
        p3 <- DimPlot(seurat_obj, reduction = "umap", group.by = "orig.ident") +
          ggtitle("UMAP Colored by corrected batch") +
          theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20)) +
          guides(color = guide_legend(override.aes = list(size = 4)))
    }
    

    if (use_spatial_coords) {
        # 画出细胞的空间分布图，x和y轴为细胞的空间坐标
        spatial_data <- seurat_obj@meta.data
        
        if(is.null(batch_key)){
            # 绘制空间分布图
            p4 <- ggplot(spatial_data, aes(x = x, y = y, color = factor(seurat_clusters))) +
              geom_point(size = 1.5, alpha = 0.8) +  # 绘制点图
              scale_color_discrete(name = "Seurat Clusters") +  # 设置颜色图例
              theme_minimal() +  # 使用简洁的主题
              labs(title = "Spatial Distribution of Cells",
                  x = "X Coordinate",
                  y = "Y Coordinate") +  # 设置标题和坐标轴标签
              theme(legend.position = "right", 
                    plot.title = element_text(hjust = 0.5, face = "bold", size = 20)) +
              guides(color = guide_legend(override.aes = list(size = 4)))

        }else{
            # 创建一个空列表来存放每个批次的图
            plot_list <- list()
            # 获取所有批次
            batch_list <- unique(spatial_data[[batch_key]])

            # 循环遍历每个批次，生成并存储图表到 plot_list 中
            for (batch in batch_list) {
              # 筛选出当前批次的数据
              batch_data <- subset(spatial_data, spatial_data[[batch_key]] == batch)
              
              # 绘制空间分布图
              p4 <- ggplot(batch_data, aes(x = x, y = y, color = factor(seurat_clusters))) +
                geom_point(size = 1.5, alpha = 0.8) +  # 绘制点图
                scale_color_discrete(name = "Seurat Clusters") +  # 设置颜色图例
                theme_minimal() +  # 使用简洁的主题
                labs(title = paste("Spatial Distribution of Cells", "\nSample:", batch),
                    x = "X Coordinate",
                    y = "Y Coordinate") +  # 设置标题和坐标轴标签
                theme(legend.position = "right", 
                      plot.title = element_text(hjust = 0.5, face = "bold", size = 20)) +
                guides(color = guide_legend(override.aes = list(size = 4)))
              
              # 将图表存储到 plot_list 中，使用 batch 作为图表名称
              plot_list[[paste("Batch", batch)]] <- p4
            }
        }

    }


    # 聚类图转换成pdf
    pdf("./report3.pdf", width = 10, height = 8 + length(unique(seurat_obj@meta.data[, "seurat_clusters"])) / 16)
        if(!is.null(batch_key)){print(p3)}
        print(p1)
        print(p2)
        if(use_spatial_coords){
            if(is.null(batch_key)){
                print(p4)
            }else{for(plot in plot_list){print(plot)}}
        }
    dev.off()

    print("聚类完成. ")

    return(seurat_obj)
}