# 加载SingleR和Seurat包
library(SingleR)
library(Seurat)
library(celldex)
library(cowplot)
library(Matrix)
library(scales)
library(patchwork)
library(grid)
library(extrafont)

run_singleR <- function(){

    # 读取数据
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


    print(paste("基因数量为：", nrow(matrix)))
    print(paste("细胞过滤前的细胞数为：", ncol(matrix)))


    ##################################### 质量控制(QC) ##############################################


    # 计算线粒体基因比例
    seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-") 


    # 绘制细胞基因测量数量小提琴图
    p1 <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3) +
        theme(plot.title = element_text(hjust = 0.5, size = 16))
    # 设置总标题
    p1 <- wrap_plots(p1) + 
        plot_annotation(
            title = "Vlolin plot of Quality Control index",
            theme = theme(
            plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
            )
        )


    # 读取Seurat自带的细胞周期基因集，用来看细胞周期的分布
    s.genes <-cc.genes$s.genes
    g2m.genes<-cc.genes$g2m.genes
    # 计算细胞周期得分
    seurat_obj <- CellCycleScoring(seurat_obj, s.features = s.genes, g2m.features = g2m.genes)
    # 提取S.Score和G2M.Score
    scores <- seurat_obj@meta.data[, c("S.Score", "G2M.Score")]


    # 绘制S.Score的直方图
    p2 <- ggplot(scores, aes(x = S.Score)) +
        geom_histogram(binwidth = 0.05, fill = "blue", color = "black", alpha = 0.5) +
        ggtitle("Histogram of S.Score") +
        xlab("S.Score") +
        ylab("Cell Count") +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold", size = 20), # 调整标题字体大小并居中
            plot.title.position = "plot"
        )


    # 绘制G2M.Score的直方图
    p3 <- ggplot(scores, aes(x = G2M.Score)) +
        geom_histogram(binwidth = 0.05, fill = "red", color = "black", alpha = 0.5) +
        ggtitle("Histogram of G2M.Score") +
        xlab("G2M.Score") +
        ylab("Cell Count") +
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold", size = 20), # 调整标题字体大小并居中
            plot.title.position = "plot"
        )

    # 保存图像到PDF文件
    pdf("./report1.pdf")
    print(p1)
    print(p2)
    print(p3)
    dev.off()


    # 设置质控条件，筛除检测到的基因数量在200以下的细胞，以及线粒体基因占总基因30%以上的细胞
    seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & percent.mt < 30)

    print("细胞筛选完毕.")
    print(paste("细胞过滤后的细胞数为：", nrow(seurat_obj@meta.data)))
    
    ##################################### SingleR细胞注释 ##############################################


    # 加载HumanPrimaryCellAtlasData参考数据集
    ref <- HumanPrimaryCellAtlasData()

    # 运行SingleR进行细胞类型标记
    singleR_results_main <- SingleR(test = seurat_obj@assays$RNA$data, ref = ref, labels = ref$label.main)

    singleR_results_detail <- SingleR(test = seurat_obj@assays$RNA$data, ref = ref, labels = ref$label.fine)

    # 将SingleR结果添加到Seurat对象的meta.data中
    seurat_obj$cell_type <- singleR_results_main$labels
    seurat_obj$sub_cell_type <-  singleR_results_detail$labels


    saveRDS(seurat_obj, file = "./init_data/seurat_obj.rds")
    
    print("SingleR注释完成. ")
}
