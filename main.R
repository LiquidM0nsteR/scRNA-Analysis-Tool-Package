# Powered by LiquidMonsteR
# From Hainan institute of WhHan University of Technology, Sanya, Hainan.
# Firstly created in July-7-2024

source("./SingleR.R")
source("./requirements.R")
source("./check_seurat.R")
source("./dimention_reduction.R")
source("./clustering.R")
source("./plot_and_ari.R")


####################################### 参数设置面板 ####################################################

# 如果已经有准备好的seurat对象存储在rds文件中，请放置在init_data文件夹下，设置seurat_from_rds为TRUE
# 如果数据来自barcodes.tsv.gz、features.tsv.gz和matrix.mtx.gz，则将seurat_from_rds置为FALSE
seurat_from_rds <- FALSE  

# 是否为多样本
is_multi_sample <- FALSE  

# 是否使用空间位置信息进行降维
use_spatial_coords <- FALSE

# 降维方法编号提示
# 1: PCA   // 主成分分析
# 2: c-ICA // 独立成分分析
# 3: scvi  // 变分自动编码器
# 4: SCA   // 意外成分分析, 注意：截至2024/7/29，该包在pypi上仍处于维护状态，暂时不可以使用
reduction_method <- 1  # 选择降维方法，1到5之间的数字


# 聚类方法编号提示
# 1: FindClusters (Louvain)
# 2: Walktrap
# 3: K-medoids
# 4: Model-Based Clustering (Mclust)
clustering_method <- 1 # 选择聚类方法，1到4的数字


######################################################################################################


if(seurat_from_rds){
  rds_file_path <- readline(prompt = "请指定rds文件的路径: ")
  seurat_obj <- readRDS(rds_file_path)
}else{
  print("将使用barcodes.tsv.gz、features.tsv.gz 和 matrix.mtx.gz生成Seurat对象")
  # run_singleR()
  seurat_obj <- readRDS("./init_data/seurat_obj.rds")
}

seurat_obj <- check_and_preprocess_seurat(seurat_obj, is_multi_sample, use_spatial_coords)
seurat_obj <- perform_dimensionality_reduction(seurat_obj, reduction_method)
seurat_obj <- standardize_reduction(seurat_obj, reduction_method)
seurat_obj <- perform_clustering(seurat_obj, clustering_method, reduction_method)
plots <- plot_umap(seurat_obj, reduction_method, clustering_method)

# 保存图像到PDF文件
pdf("./report.pdf")
for (plot in plots) {
  print(plot)
}
dev.off()

print("Successfully finished. ")
