# Powered by LiquidMonsteR
# From Hainan institute of WhHan University of Technology, Sanya, Hainan.
# Firstly created in July-7-2024


source("./requirements.R")
source("./QC.R")
source("./preprocess.R")
source("./dimention_reduction.R")
source("./clustering.R")
source("./downstream_analysis.R")
options(warn = -1)


####################################### 参数设置面板 ####################################################


# 如果已经有准备好的seurat对象存储在rds文件中，请放置在init_data文件夹下，设置seurat_from_rds为TRUE
# 如果数据来自 filtered_feature_bc_matrix.h5，则将seurat_from_rds置为FALSE
seurat_from_rds <- FALSE 

# 指定细胞自动标注的参考数据集：
# 1. HumanPrimaryCellAtlasData（人类）； 2. MouseRNAseqData（老鼠）
# 3. BlueprintEncodeData（免疫细胞）; 4. DatabaseImmuneCellExpressionData（免疫细胞）
# 5. ImmGenData（小鼠免疫细胞）
# 6. MonacoImmuneData（人免疫细胞）; 7. NovershternHematopoieticData（人造血干细胞和祖细胞）
# 8. 其它的注释，稍后将会要求指定count_matrix_path, cell_annotation_path。
#     注意，若要搜索对应物种的测序数据，可搜索以下网站：
#     https://www.ncbi.nlm.nih.gov/geo/
print("请指定细胞类型注释的参考数据集，可输入多个数字，以空格分隔：")
anno_ref <- readLines(con = "stdin", n = 1)

species <- 1 # 1为人类，2为老鼠（目前仅支持人和老鼠）

# 降维方法编号提示
# 1: PCA   // 主成分分析
# 2: c-ICA // 独立成分分析
# 3: scvi  // 变分自动编码器
reduction_method <- 3  # 选择降维方法，1到5之间的数字

# 聚类方法编号提示
# 1: FindClusters (Louvain)
# 2: Walktrap
# 3: Model-Based Clustering (Mclust)
clustering_method <- 3 # 选择聚类方法，1到4的数字


########################################### 流程执行 ##########################################################

# 将输入的字符串分割为单个元素，并转换为数值向量
anno_ref <- as.numeric(strsplit(anno_ref, "\\s+")[[1]])

if(8 %in% anno_ref){
    print("请指定自定义参考集 count_matrix.txt.gz 文件的绝对路径: ")
    count_matrix_path <- readLines(con = "stdin", n = 1)

    print("请指定自定义参考集 cell_annotation.txt.gz 文件的绝对路径: ")
    cell_annotation_path <- readLines(con = "stdin", n = 1)
}

result1 <- init_seurat(seurat_from_rds)
    seurat_obj <- result1$seurat_obj
    is_multiple <- result1$is_multiple
    use_spatial_coords <- result1$use_spatial_coords

seurat_obj <- run_QC(seurat_obj, seurat_from_rds)

result2 <- perform_dimensionality_reduction(seurat_obj, reduction_method)
    seurat_obj <- result2$seurat_obj
    num_components <- result2$num_components

seurat_obj <- perform_clustering(seurat_obj, clustering_method, reduction_method, num_components)

downstream_analysis(seurat_obj, anno_ref, count_matrix_path, cell_annotation_path, species)

