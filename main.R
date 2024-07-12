source("./requirements.R")
source("./check_seurat.R")
source("./dimention_reduction.R")
source("./clustering.R")
source("./plot_and_ari.R")
# powered by LiquidMonsteR
# 参数设置面板

is_multi_sample <- FALSE  # 如果是多样本，设置为TRUE，否则设置为FALSE

seurat_obj_path <- "/share/home/bgi_clicx/test_flower/obj_rice_flower_merge_filter.rds" # 设置seurat对象的存储路径

# 降维方法编号提示
# 1: PCA
# 2: c-ICA
# 3: Autoencoder
# 4: AVAE
# 5: SCA
reduction_method <- 1  # 选择降维方法，1到5之间的数字

# 聚类方法编号提示
# 1: FindClusters (Louvain)
# 2: Walktrap
# 3: K-medoids
# 4: Model-Based Clustering (Mclust)
clustering_method <- 1 # 选择聚类方法，1到4的数字




seurat_obj = readRDS(seurat_obj_path)
seurat_obj <- check_and_preprocess_seurat(seurat_obj, is_multi_sample)
seurat_obj <- perform_dimensionality_reduction(seurat_obj, reduction_method)
seurat_obj <- standardize_reduction(seurat_obj)
seurat_obj <- perform_clustering(seurat_obj, clustering_method)
plots <- plot_umap(seurat_obj, reduction_method, clustering_method)

# 保存图像到PDF文件
pdf("./report.pdf")
for (plot in plots) {
  print(plot)
}
dev.off()
