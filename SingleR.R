# 加载SingleR和Seurat包
library(SingleR)
library(Seurat)
library(celldex)
library(Matrix)


run_singleR <- function(){

    print("正在检测init_data文件夹下的必要数据：barcodes.tsv.gz、features.tsv.gz 和 matrix.mtx.gz.")

    # 设置文件路径和文件名
    data_dir <- "./init_data/"  
    required_files <- c("barcodes.tsv.gz", "features.tsv.gz", "matrix.mtx.gz")

    # 检查文件是否存在
    missing_files <- required_files[!file.exists(file.path(data_dir, required_files))]

    # 输出检查结果
    if (length(missing_files) > 0) {
        cat("缺失文件: ", paste(missing_files, collapse = ", "), "\n")
    } else {
        cat("数据检查完毕.\n")
    }

    # 设定文件路径
    barcode_file <- "./init_data/barcodes.tsv.gz"
    features_file <- "./init_data/features.tsv.gz"
    matrix_file <- "./init_data/matrix.mtx.gz"

    # 读取barcodes, features, and matrix
    barcodes <- readLines(barcode_file)
    features <- read.delim(features_file, header = FALSE)$V2 # 使用基因标志(Gene Symbols)作为行名
    matrix <- readMM(matrix_file)

    # 创建稀疏矩阵并指定行名和列名
    rownames(matrix) <- features
    colnames(matrix) <- barcodes

    # 除去重复的基因
    matrix <- matrix[unique(features), ]


    # 创建Seurat对象
    seurat_obj <- CreateSeuratObject(counts = matrix)
    seurat_obj <- NormalizeData(seurat_obj)

    # 加载HumanPrimaryCellAtlasData参考数据集
    ref <- HumanPrimaryCellAtlasData()

    # 运行SingleR进行细胞类型标记
    singleR_results <- SingleR(test = seurat_obj@assays$RNA$data, ref = ref, labels = ref$label.fine)

    # 将SingleR结果添加到Seurat对象的meta.data中
    seurat_obj$cell_type <- singleR_results$labels


    saveRDS(seurat_obj, file = "./init_data/seurat_obj.rds")
}
