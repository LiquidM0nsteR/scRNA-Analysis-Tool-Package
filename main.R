# Powered by LiquidMonsteR
# From Hainan institute of WhHan University of Technology, Sanya, Hainan.
# Firstly created in July-7-2024


#source("./requirements.R")
source("./QC.R")
source("./preprocess.R")
# source("./batch_correction.R")
source("./dimension_reduction.R")
source("./clustering.R")
source("./downstream_analysis.R")
options(warn = -1)
cat("Powered by LiquidMonsteR.\n\n ")
cat("请将包含原始数据的 .rds 文件，或者 .h5 文件，或者带有 spatial + .h5 的文件夹，放入 init_data 文件夹中. \n\n\n")
Sys.sleep(3.6)

####################################### 参数设置面板 ####################################################

# 参数1: 是否从RDS文件导入Seurat对象
cat("是否从RDS文件导入Seurat对象？请输入 TRUE 或 FALSE：\n")
seurat_from_rds <- as.logical(readLines(con = "stdin", n = 1))


# 参数2：是否使用空间坐标
cat("是否为空间转录组数据？请输入 TRUE 或 FALSE：\n")
use_spatial_coords <- as.logical(readLines(con = "stdin", n = 1))
if(use_spatial_coords){
    if(seurat_from_rds){
        cat("请确保rds的meta.data中有 x 和 y 两列，即空间坐标. ")
    }
    else{
        cat("请确保存放表达矩阵的h5文件名称为 filtered_feature_bc_matrix.h5 \n")
    }
    Sys.sleep(3)
}


# 参数3: 细胞类型注释参考数据集
cat("
请指定细胞类型注释的参考数据集编号（可输入多个编号，以空格分隔）：
1: HumanPrimaryCellAtlasData（人类）
2: MouseRNAseqData（老鼠）
3: BlueprintEncodeData（免疫细胞）
4: DatabaseImmuneCellExpressionData（免疫细胞）
5: ImmGenData（小鼠免疫细胞）
6: MonacoImmuneData（人免疫细胞）
7: NovershternHematopoieticData（人造血干细胞和祖细胞）
8: 自定义参考集\n")
# 获取用户输入的注释集编号
anno_ref <- as.numeric(strsplit(readLines(con = "stdin", n = 1), "\\s+")[[1]])


# 参数4: 物种选择
cat("请问物种选择是什么？请输入数字：1: 人类  2: 老鼠\n")
species <- as.numeric(readLines(con = "stdin", n = 1))


# 参数5: 去批次方法
cat("请指定去批次方法编号：\n")
cat("1: Seurat标准流程\n2: ComBat\n3: Harmony\n4: 不进行去批次(单样本)\n")
batch_correction_method <- as.numeric(readLines(con = "stdin", n = 1))

# 去批次方法参考文献
# 1: Seurat标准流程
#    Stuart T, Butler A, et al. (2019). Comprehensive Integration of Single-Cell Data. Cell, 177(7), 1888-1902.e21. https://doi.org/10.1016/j.cell.2019.05.031
# 2: Harmony
#    Korsunsky I, et al. (2019). Fast, sensitive and accurate integration of single-cell data with Harmony. Nat Methods, 16, 1289–1296. https://doi.org/10.1038/s41592-019-0619-0
# 3: ComBat
#    Johnson WE, Li C, Rabinovic A. (2007). Adjusting batch effects in microarray expression data using empirical Bayes methods. Biostatistics, 8(1), 118-127. https://doi.org/10.1093/biostatistics/kxj037


# 参数6: 降维方法
cat("请指定降维方法编号：")
cat("
1: PCA（Seurat标准流程）
2: c-ICA
3: scVI（自带去批次功能）
4: SpatialPCA（仅适用于空间转录组数据，自带去批次功能）
5: SeuratSpatial（Seurat标准流程，仅适用于空间转录组数据）\n")
reduction_method <- as.numeric(readLines(con = "stdin", n = 1))
# 检查use_spatial_coords的值
if (!use_spatial_coords && (reduction_method == 4 || reduction_method == 5)) {
    stop("错误：当前未使用空间坐标，不能选择4或5作为降维方法。请重新选择. ")
}

# 降维方法参考文献
# 1: PCA (主成分分析)
#    Jolliffe IT. (2002). Principal Component Analysis. Springer Series in Statistics. Springer-Verlag, New York. https://doi.org/10.1007/b98835
# 2: c-ICA (独立成分分析)
#    Teschendorff AE, et al. (2007). Elucidating the altered transcriptional programs in breast cancer using independent component analysis. PLOS Computational Biology, 3(8), e161. https://doi.org/10.1371/journal.pcbi.0030161
# 3: scVI (变分自动编码器)
#    Lopez R, Regier J, et al. (2018). Deep generative modeling for single-cell transcriptomics. Nat Methods, 15(12), 1053-1058. https://doi.org/10.1038/s41592-018-0229-2
# 4: SpatialPCA (空间主成分分析)
#    Shang, L., & Zhou, X. (2022). Spatially Aware Dimension Reduction for Spatial Transcriptomics. Nature Communications. https://doi.org/10.1038/s41467-022-34879-1
# 5: SeuratSpatial（Seurat标准流程）
#    Satija, R., Farrell, J.A., Gennert, D., Schier, A.F., & Regev, A. (2015). Spatial reconstruction of single-cell gene expression data. Nature Biotechnology, 33(5), 495-502. https://doi.org/10.1038/nbt.3192


# 参数7: 聚类方法
cat("请指定聚类方法编号：")
cat("
1: FindClusters (Louvain)
2: Walktrap
3: Model-Based Clustering (Mclust)\n")
clustering_method <- as.numeric(readLines(con = "stdin", n = 1))

# 聚类方法参考文献
# 1: FindClusters (Louvain算法)
#    Blondel VD, et al. (2008). Fast unfolding of communities in large networks. J Stat Mech, 2008(10), P10008. https://doi.org/10.1088/1742-5468/2008/10/P10008
# 2: Walktrap
#    Pons P, Latapy M. (2005). Computing communities in large networks using random walks. J Graph Algorithms Appl, 10(2), 191-218. https://doi.org/10.7155/jgaa.00124
# 3: Model-Based Clustering (Mclust)
#    Fraley C, Raftery AE. (2002). Model-Based Clustering, Discriminant Analysis, and Density Estimation. J Am Stat Assoc, 97(458), 611-631. https://doi.org/10.1198/016214502760047131


# 检查是否选择了自定义参考集路径
if (8 %in% anno_ref) {
    cat("请指定自定义参考集 count_matrix.txt.gz 文件的绝对路径: \n")
    count_matrix_path <- readLines(con = "stdin", n = 1)

    cat("请指定自定义参考集 cell_annotation.txt.gz 文件的绝对路径: \n")
    cell_annotation_path <- readLines(con = "stdin", n = 1)
} else {
    count_matrix_path <- NULL
    cell_annotation_path <- NULL
}



########################################### 流程执行 ##########################################################

# 初始化Seurat对象
result1 <- init_seurat(seurat_from_rds, use_spatial_coords)
    seurat_obj <- result1$seurat_obj
    batch_key <- result1$batch_key

# 运行质量控制
seurat_obj <- run_QC(seurat_obj, seurat_from_rds, batch_key)

# 执行降维
result2 <- perform_dimensionality_reduction(seurat_obj, reduction_method, batch_correction_method, batch_key)
    seurat_obj <- result2$seurat_obj
    num_components <- result2$num_components

# 执行聚类
seurat_obj <- perform_clustering(seurat_obj, clustering_method, reduction_method, num_components, use_spatial_coords, batch_key)

# 下游分析
downstream_analysis(seurat_obj, anno_ref, count_matrix_path, cell_annotation_path, species)
