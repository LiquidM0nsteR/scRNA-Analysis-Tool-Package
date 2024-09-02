# 加载SingleR和Seurat包
library(cowplot)
library(Matrix)
library(scales)
library(patchwork)
library(grid)
library(dplyr)
library(ggplot2)
library(doParallel)
library(foreach)
library(parallel)




run_QC <- function(seurat_obj, seurat_from_rds, batch_key){
    
    if(!seurat_from_rds){

        print(paste("基因数量为：", nrow(seurat_obj@assays$RNA$counts)))
        print(paste("细胞过滤前的细胞数为：", ncol(seurat_obj@assays$RNA$counts)))


        # 计算线粒体基因比例
        seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-") 

        # 设置质控条件，筛除检测到的基因数量在200以下的细胞，以及线粒体基因占总基因30%以上的细胞
        seurat_obj <- subset(seurat_obj, subset = nFeature_RNA > 200 & percent.mt < 30)

        print("细胞筛选完毕.")
        print(paste("细胞过滤后的细胞数为：", nrow(seurat_obj@meta.data)))
        
        # 读取Seurat自带的细胞周期基因集，用来看细胞周期的分布
        # 从 Seurat 加载新版的 cc.genes
        cc.genes.updated.2019 <- Seurat::cc.genes.updated.2019
        s.genes <- cc.genes.updated.2019$s.genes
        g2m.genes <- cc.genes.updated.2019$g2m.genes
        

        # # 计算细胞周期得分
        seurat_obj <- NormalizeData(seurat_obj)
        seurat_obj <- CellCycleScoring(seurat_obj, s.features = s.genes, g2m.features = g2m.genes)


        # 消除与细胞周期相关基因的影响
        # 设置并行计算
        print("正在消除与细胞周期相关基因的影响.")

        library(future)
        num_cores <- max(1, parallel::detectCores() - 2)
        plan("multicore", workers = num_cores)  
        options(future.globals.maxSize = 224 * 1024^3)  
        # 加速ScaleData过程
        seurat_obj <- ScaleData(seurat_obj, vars.to.regress = c("S.Score", "G2M.Score", batch_key), 
            features = rownames(seurat_obj), do.scale = TRUE, do.center = TRUE)  
        plan("sequential") 

        print("消除完成.")
       
    }

    seurat_obj <- FindVariableFeatures(seurat_obj, verbose = FALSE)

    QC_plot(seurat_obj, seurat_from_rds)

    return(seurat_obj)
}


QC_plot <- function(seurat_obj, seurat_from_rds){
    
    if(seurat_from_rds){
        # 计算线粒体基因比例
        seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-") 

        # 计算细胞周期分数
        s.genes <-cc.genes$s.genes
        g2m.genes<-cc.genes$g2m.genes
        seurat_obj <- CellCycleScoring(seurat_obj, s.features = s.genes, g2m.features = g2m.genes)
    }
    

    # 绘制细胞基因测量数量小提琴图
    p1 <- VlnPlot(seurat_obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, group.by = 'orig.ident') +
        theme(plot.title = element_text(hjust = 0.5, size = 16))
    # 设置总标题
    p1 <- wrap_plots(p1) + 
        plot_annotation(
            title = "Vlolin plot of Quality Control index",
            theme = theme(
            plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
            )
        )


    # 画出细胞周期得分的小提琴图
    p2 <- VlnPlot(seurat_obj ,features = c('S.Score','G2M.Score'), group.by = 'orig.ident')+ 
        theme(plot.title = element_text(hjust = 0.5, size = 16))
    p2 <- wrap_plots(p2) + 
        plot_annotation(
            title = "Violin plot of Cell Cycle Score",
            theme = theme(
                plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
            )
        )
    
    # 保存图像到PDF文件
    pdf("./report1.pdf", width = 10)
    print(p1)
    print(p2)
    dev.off()
}

