# clustering.R
# powered by LiquidMonsteR
library(Seurat)
library(cluster)  # For K-medoids
library(mclust)  # For Model-Based Clustering
library(igraph)  # For Walktrap
library(scran) # For Walktrap
library(monocle3) # For monocle3-clustering


perform_clustering <- function(seurat_obj, clustering_method, reduction_method, num_components, use_spatial_coords, batch_key) {

    print("正在进行聚类. ")
    # 映射 method 到降维方法的名称
    reduction_name <- switch(reduction_method,
                            `1` = "PCA",
                            `2` = "cica",
                            `3` = "scvi",
                            `4` = "SpatialPCA",
                            `5` = "SeuratSpatial",
                            `6` = "SEDR",
                            stop("不存在的降维方法编号. "))

    if(reduction_name == "SpatialPCA" || reduction_name == "SEDR"){ # 多样本 SpatialPCA 结果固定20个维度
        num_components <- min(ncol(Embeddings(seurat_obj, reduction = reduction_name)), num_components)
    }

    # 根据聚类方法编号，执行相应的聚类函数
    seurat_obj <- switch(
        clustering_method,
        `1` = run_findclusters(seurat_obj, reduction_name, num_components),
        `2` = run_walktrap(seurat_obj, reduction_name),
        `3` = run_model_based_clustering(seurat_obj, reduction_name),
        `4` = run_monocle3_clustering(seurat_obj, reduction_name),
        stop("无效的聚类方法编号，请重选1到6之间的数字. ")
    )


    # 根据聚类结果画umap图，同时也画出不同批次的umap图，观察去批次的效果
    seurat_obj <- plot_umap(seurat_obj, reduction_name, clustering_method, num_components, use_spatial_coords, batch_key)
    
    return(seurat_obj)
}



run_findclusters <- function(seurat_obj, reduction_name, num_components) {

    neighbor_graph_name <- "my_custom_nn"
    snn_graph_name <- "my_custom_snn"
    # 使用指定的邻居图名称运行 FindNeighbors
        seurat_obj <- FindNeighbors(
        seurat_obj, 
        reduction = reduction_name, 
        dims = 1:num_components,  # 使用全部的低维矩阵
        k.param = 10,             # 设置邻居数
        graph.name = c(neighbor_graph_name, snn_graph_name)
    )

    # 使用 FindClusters 调用指定的 SNN 图
    seurat_obj <- FindClusters(
        seurat_obj, 
        graph.name = snn_graph_name, 
        resolution = 0.2  # 代表聚类的分辨率，值越高聚类数就越多
    )
    return(seurat_obj)
}



run_walktrap <- function(seurat_obj, reduction_name) {

    seurat_obj <- FindNeighbors(seurat_obj, reduction = reduction_name)
    # 提取邻居图
    graph_matrix <- as.matrix(seurat_obj@graphs[[1]])
    
    # 将邻接矩阵转换为 igraph 对象
    graph <- igraph::graph_from_adjacency_matrix(graph_matrix, mode = "undirected", weighted = TRUE)
    
    clusters <- cluster_walktrap(graph)
    seurat_obj$seurat_clusters <- factor(clusters$membership)
    return(seurat_obj)
}



run_model_based_clustering <- function(seurat_obj, reduction_name) {

    # 获取指定降维方法的嵌入数据
    data <- Embeddings(seurat_obj, reduction = reduction_name)
    
    # 使用Mclust进行模型驱动的聚类
    mclust_result <- Mclust(data)
    
    # 处理未分配的细胞
    if (any(is.na(mclust_result$classification))) {
        # 找到未分配的细胞的索引
        unclassified_indices <- which(is.na(mclust_result$classification))
        
        for (i in unclassified_indices) {
            # 找到该细胞所属的最近的聚类
            nearest_cluster <- which.max(mclust_result$z[i, ])
            mclust_result$classification[i] <- nearest_cluster
        }
    }
    
    # 将聚类结果添加到Seurat对象
    seurat_obj$seurat_clusters <- factor(mclust_result$classification)
    
    return(seurat_obj)
}



run_monocle3_clustering <- function(seurat_obj, reduction_name){
    
    # 将 Seurat 对象转换为 Monocle3 的 CellDataSet 对象
    cds <- new_cell_data_set(seurat_obj@assays$RNA$counts, cell_metadata = seurat_obj@meta.data)
    
    # 将 Seurat 对象中的指定降维方法嵌入数据添加到 Monocle3 对象中
    # 这里由于monocle3的降维方法名称限制，使用PCA进行伪装
    reducedDims(cds)[["PCA"]] <- Embeddings(seurat_obj, reduction_name)
    
    # 使用 Monocle3 的 cluster_cells 函数进行聚类, 其中k为细胞邻居数
    cds <- cluster_cells(cds, reduction_method = "PCA", k = 5)
    
    # 将 Monocle3 的聚类结果提取并添加到 Seurat 对象中
    seurat_obj$seurat_clusters <- factor(cds@clusters[["PCA"]]$clusters)
    
    return(seurat_obj)
}






plot_umap <- function(seurat_obj, reduction_name, clustering_method, num_components, use_spatial_coords, batch_key){

    seurat_obj <- Seurat::RunUMAP(seurat_obj, reduction = reduction_name, dims = 1:num_components, verbose = FALSE)
    seurat_obj <- Seurat::RunTSNE(seurat_obj, reduction = reduction_name, dims = 1:num_components, verbose = FALSE)


    clustering_methods <- c("Louvain", "Walktrap", "Mclust", "Monocle3_clustering")
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