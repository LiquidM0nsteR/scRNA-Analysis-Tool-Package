# dimensionality_reduction.R
# powered by LiquidMonsteR
library(Seurat)
library(ica)  # For c-ICA
library(SpatialPCA) # For SpatialPCA
library(reticulate)  # For Python interaction
library(doParallel) # for multicores
library(future)
library(future.apply)
library(foreach)
library(parallel)
source("./batch_correction.R")

perform_dimensionality_reduction <- function(seurat_obj, reduction_method, batch_correction_method, batch_key) {

    print("正在进行降维. ")

    seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
    seurat_obj <- RunUMAP(seurat_obj, dims = 1:30 ,verbose = FALSE)
  
    # 如果是多样本，先执行去批次操作（Harmony方法除外）
    if(!is.null(batch_key)){ 
        if(!(reduction_method %in% c(3, 4))){ # scvi 和 SpatialPCA 自带去批次功能，所以跳过这一段
            # 根据 batch_correction_method 的值调用相应的去批次方法
            seurat_obj <- switch(
              batch_correction_method,
              `1` = run_seurat_integration(seurat_obj, batch_key),  # 1 对应 Seurat Integration            
              `2` = run_combat(seurat_obj, batch_key),              # 2 对应 Combat
              `3` = seurat_obj,                                     # 3 对应 Harmony，等降维之后再执行
              `4` = seurat_obj                                      # 4 对应 不去批次，直接返回原Seurat对象
            )
        }
    }

    # 使用ElbowPlot寻找最佳的低维维度
    p1 <- ElbowPlot(seurat_obj) +
        ggtitle("Elbow Plot for PCA Components") +
        theme(plot.title = element_text(size = 20, face = "bold", hjust = 0.5))
    
    if(!is.null(batch_key)){
        # 绘制UMAP图像，颜色为不同批次，并将图例设置为两列
        p2 <- DimPlot(seurat_obj, reduction = "umap", group.by = "orig.ident") +
            ggtitle("UMAP Colored by initial batch") +
            theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20)) +
            guides(color = guide_legend(override.aes = list(size = 4)))
    }

    # 保存Elbow Plot为PDF文件
    pdf("./report2.pdf", width = 10)
    print(p1)
    if(!is.null(batch_key)){print(p2)}
    dev.off()

    # 提取PCA结果
    pca_result <- seurat_obj[["pca"]]

    # 计算解释方差比例
    explained_variance <- cumsum(pca_result@stdev^2) / sum(pca_result@stdev^2)
    num_components <- which(explained_variance >= 0.90)[1]

    print(paste("至少解释90%方差的主成分数量:", num_components, "，设置为降维结果的维度. "))


    # 根据方法选择降维算法
    seurat_obj <- switch(reduction_method,
                        `1` = run_pca(seurat_obj, num_components),
                        `2` = run_cica(seurat_obj, num_components),
                        `3` = run_scvi(seurat_obj, num_components),
                        `4` = run_SpatialPCA(seurat_obj, batch_key, num_components),
                        `5` = run_Seurat_Spatial(seurat_obj, batch_key, num_components),
                        stop("Invalid reduction method. Choose a number between 1 and 5.")
    )


    # 如果选用了Harmony去批次，那么就降维之后再去批次。
    if(!is.null(batch_key) && batch_correction_method == 3 && !(reduction_method %in% c(3, 4))){
        print("开始进行 Harmony 去批次.")
        seurat_obj <- run_harmony(seurat_obj, batch_key, reduction_method, num_components)
    }
    

    print("降维完成. ")

    # 返回结果
    return(list(seurat_obj = seurat_obj, num_components = num_components))
}



run_pca <- function(seurat_obj, num_components) {
    
    # 检查当前assay中是否有 VariableFeatures
    if (length(VariableFeatures(seurat_obj)) == 0) {
        seurat_obj <- FindVariableFeatures(seurat_obj, verbose = FALSE)
    }

    # 获取当前assay的名字
    current_assay <- DefaultAssay(seurat_obj)
    # 检查当前assay中是否有scale.data
    if (is.null(seurat_obj[[current_assay]]$scale.data) || nrow(seurat_obj[[current_assay]]$scale.data) == 0) {
        seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
    }


    # 使用设定的默认assay进行PCA降维
    seurat_obj <- RunPCA(seurat_obj, npcs = num_components, verbose = FALSE)
    embeddings <- seurat_obj[["pca"]]@cell.embeddings
    colnames(embeddings) <- paste0("PCA_", 1:ncol(embeddings))
    rownames(embeddings) <- colnames(seurat_obj@assays$RNA$data)
    seurat_obj[["PCA"]] <- CreateDimReducObject(embeddings = embeddings, key = "PCA_")
    return(seurat_obj)
}



run_cica <- function(seurat_obj, num_components) {
  
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

    print(paste0("降维矩阵的大小为：", nrow(latent_representation), "*", ncol(latent_representation)))

    rownames(latent_representation) <- rownames(expr_data)
    colnames(latent_representation) <- paste0("scvi_", 1:ncol(latent_representation))
    
    seurat_obj[["scvi"]] <- CreateDimReducObject(embeddings = latent_representation, key = "scvi_")
    return(seurat_obj)
}



run_SpatialPCA <- function(seurat_obj, batch_key, num_components) {

    num_cores <- max(1, parallel::detectCores() - 2)

    if (!is.null(batch_key)) {
        # 处理多样本情况
        
        # 提取基因表达矩阵和位置信息列表，按batch_key区分
        sample_ids <- unique(seurat_obj@meta.data[[batch_key]])
        
        count_list <- lapply(sample_ids, function(sample_id) {
            sample_cells <- which(seurat_obj@meta.data[[batch_key]] == sample_id)
            as.matrix(GetAssayData(seurat_obj, slot = "counts")[, sample_cells])
        })
        
        location_list <- lapply(sample_ids, function(sample_id) {
            sample_cells <- which(seurat_obj@meta.data[[batch_key]] == sample_id)
            as.matrix(seurat_obj@meta.data[sample_cells, c("x", "y")])
        })
        
        # 调用多样本SpatialPCA函数
        spatial_pca_result <- SpatialPCA_Multiple_Sample(
            count_list = count_list,        # 各个批次的基因表达数据
            location_list = location_list,  # 各个批次的细胞位置
            gene.type = "spatial",          # 需要分析的基因类型设置为空间高可变基因
            sparkversion = "sparkx",        # 识别空间高可变基因的算法版本
            numCores_spark = num_cores,     # spark核心数设置为当前可用核心数-2
            gene.number = 2000,             # 保留用于分析的基因数量
            customGenelist = NULL,          # 指定需要分析的基因
            min.loctions = 30,              # 过滤在小于min.loctions个位置上表达的基因，也就是大部分位置不表达的基因
            min.features = 30,              # 过滤只有小于min.features个基因表达的细胞，也就是大部分基因不表达的细胞
            bandwidth_common = 0.1          # 体现空间位置相关性的重要程度，越大空间相关性越重要
        )
        
        # 提取并存储降维结果
        embeddings <- t(do.call(cbind, spatial_pca_result$SpatialPC_list))

        # 使用 gsub 去掉行名中的 "Sample1_"
        rownames(embeddings) <- gsub("^Sample[0-9]+_", "", rownames(embeddings))
        colnames(embeddings) <- paste0("SpatialPCA_", 1:ncol(embeddings))
        
        seurat_obj[["SpatialPCA"]] <- CreateDimReducObject(embeddings = embeddings, key = "SpatialPCA_")
        
    } else {
        # 处理单样本情况
        
        # 提取基因表达矩阵和位置信息
        count_matrix <- as.matrix(GetAssayData(seurat_obj, slot = "counts"))
        location_matrix <- as.matrix(seurat_obj@meta.data[, c("x", "y")])
        
        # 调用单样本SpatialPCA函数
        spatial_pca_result <- CreateSpatialPCAObject(
            counts = count_matrix,
            location = location_matrix,
            gene.type = "spatial",
            sparkversion = "sparkx",
            numCores_spark = num_cores,
            gene.number = 2000,
            customGenelist = NULL,
            min.loctions = 50,
            min.features = 200
        )
        
        # 估算空间主成分
        spatial_pca_result <- SpatialPCA_buildKernel(spatial_pca_result, kerneltype = "gaussian")
        spatial_pca_result <- SpatialPCA_EstimateLoading(spatial_pca_result, fast = TRUE, SpatialPCnum = num_components)
        spatial_pca_result <- SpatialPCA_SpatialPCs(spatial_pca_result, fast = TRUE)
        
        # 提取并存储降维结果
        embeddings <- t(spatial_pca_result@SpatialPCs)
        
        colnames(embeddings) <- paste0("SpatialPCA_", 1:ncol(embeddings))
        
        seurat_obj[["SpatialPCA"]] <- CreateDimReducObject(embeddings = embeddings, key = "SpatialPCA_")
    }
    
    return(seurat_obj)
}



run_Seurat_Spatial <- function(seurat_obj, batch_key, num_components) {

    if (!is.null(batch_key)) {
        # 处理多样本情况
        # 将 Seurat 对象按照 batch_key 列拆分为多个子集
        seurat_list <- SplitObject(seurat_obj, split.by = batch_key)
        
        # 初始化列表来存储每个批次的高可变基因
        all_variable_features <- list()

        # 设置并行计算
        plan(multisession, workers = max(1, parallel::detectCores() - 2))
        # 并行处理每个批次的数据
        all_variable_features <- future_lapply(seurat_list, function(batch_data) {
            
            # 1. SCTransform 标准化处理
            batch_data <- SCTransform(batch_data, verbose = FALSE)
            
            # 查找高变基因
            batch_data <- FindVariableFeatures(batch_data, assay = "SCT", slot = "data", verbose = FALSE)
            
            # 获取与当前批次相关联的 image 数据
            sample_name <- unique(batch_data$group)
            batch_data@images <- seurat_obj@images[sample_name]
            
            # 2. 寻找空间高可变基因
            batch_data <- FindSpatiallyVariableFeatures(
                object = batch_data,
                assay = "SCT",
                selection.method = "moransi",
                features = VariableFeatures(batch_data),  # 选择每个样本的高可变基因
                nfeatures = 2000
            )
            
            # 返回空间高可变基因
            return(VariableFeatures(batch_data))
        })

        # 合并所有批次的高可变基因，并去重
        combined_variable_features <- unique(unlist(all_variable_features))

        # 3. 缩放数据以适合 PCA 等降维方法
        seurat_obj <- ScaleData(seurat_obj, features = combined_variable_features, verbose = FALSE)

        # 4. 使用这些合并的空间高可变基因进行 PCA 降维
        seurat_obj <- RunPCA(seurat_obj, features = combined_variable_features, npcs = num_components, verbose = FALSE)

        # 5. 将 PCA 结果存储在 Seurat 对象中
        embeddings <- Embeddings(seurat_obj, "pca")
        colnames(embeddings) <- paste0("SeuratSpatial_", 1:ncol(embeddings))
        seurat_obj[["SeuratSpatial"]] <- CreateDimReducObject(embeddings = embeddings, key = "SeuratSpatial_")

    } else {
        # 处理单样本情况
        # 1. SCTransform 标准化处理
        seurat_obj <- SCTransform(seurat_obj, verbose = FALSE)
        seurat_obj <- FindVariableFeatures(seurat_obj, assay = "SCT", slot = "data", verbose = FALSE)
        
        # 2. 寻找空间高可变基因
        seurat_obj <- FindSpatiallyVariableFeatures(
            object = seurat_obj,
            assay = "SCT",
            selection.method = "moransi",
            features = VariableFeatures(seurat_obj),  # 使用经过SCTransform标准化后的高可变基因
            nfeatures = 2000
        )

        # 3. 缩放数据以适合 PCA 等降维方法
        seurat_obj <- ScaleData(seurat_obj, features = VariableFeatures(seurat_obj), verbose = FALSE)
        
        # 4. 使用这些空间高可变基因进行 PCA 降维
        seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(seurat_obj), npcs = num_components, verbose = FALSE)
        
        # 5. 将 PCA 结果存储在 Seurat 对象中
        embeddings <- Embeddings(seurat_obj, "pca")
        colnames(embeddings) <- paste0("SeuratSpatial_", 1:ncol(embeddings))
        seurat_obj[["SeuratSpatial"]] <- CreateDimReducObject(embeddings = embeddings, key = "SeuratSpatial_")
    }
    return(seurat_obj)
}





