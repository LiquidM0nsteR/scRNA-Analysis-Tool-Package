# downstream_analysis.R
# Powered by LiquidMonsteR.
library(Seurat)
library(future)
library(monocle3)
library(BiocParallel)
library(parallel)
library(reticulate)
library(SingleR)
library(dplyr)
library(ggplot2)
library(SummarizedExperiment)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(pdftools)


downstream_analysis <- function(seurat_obj, anno_ref, count_matrix_path = NULL, cell_annotation_path =NULL, species){

    # seurat_obj = readRDS("./init_data/seurat_obj.rds")

    Idents(seurat_obj) <- seurat_obj$seurat_clusters

    ####################################### 差异表达分析 ###############################################

    print("进行差异表达分析. ")

    # 配置并行计算
    num_cores <- max(1, parallel::detectCores() - 2)
    plan("multicore", workers = num_cores) # 根据实际CPU核心数量设置workers数量

    # 使用FindAllMarkers进行标记基因计算
    all_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

    # 恢复为单线程
    plan("sequential")

    write.csv(all_markers, "./cache/all_markers.csv")


    # 1. 计算每个细胞类型的细胞数并排序
    cell_type_counts <- seurat_obj@meta.data %>%
      group_by(seurat_clusters) %>%
      summarize(count = n()) %>%
      arrange(desc(count))

    # 2. 计算累积和并选择能覆盖95%细胞的细胞类型
    cell_type_counts <- cell_type_counts %>%
      mutate(cumulative_sum = cumsum(count),
            cumulative_perc = cumulative_sum / sum(count))

    selected_cell_types <- cell_type_counts %>%
      filter(cumulative_perc <= 0.95) %>%
      pull(seurat_clusters)

    # 如果覆盖95%的细胞数量的细胞类型不够多，可以放宽条件
    if (length(selected_cell_types) <= 7) {
      selected_cell_types <- cell_type_counts %>%
        filter(cumulative_perc <= 0.98) %>%
        pull(seurat_clusters)
    }

    # 选择符合条件的细胞
    selected_cells <- seurat_obj@meta.data %>%
      filter(seurat_clusters %in% selected_cell_types) %>%
      rownames()

    # 创建一个新的 Seurat 对象，只包含选定的细胞
    seurat_obj_selected <- subset(seurat_obj, cells = selected_cells)


    # 使用slice_max选择每个群集前5个显著的标记基因
    top_markers <- all_markers %>%
      group_by(cluster) %>%
      slice_max(order_by = avg_log2FC, n = 5) %>%
      ungroup()

    selected_genes <- unique(top_markers$gene)


    # 5. 创建热图
    p1 <- DoHeatmap(seurat_obj, features = selected_genes) + NoLegend() +
      ggtitle("Top 5 Marker Genes for Each Cluster") +
      theme(plot.title = element_text(size = 20, face = "bold",hjust = 0.5), plot.margin = margin(5, 50, 5, 5))

    p2 <- DotPlot(seurat_obj, features = top_markers$gene[1:10], group.by = "seurat_clusters") + 
      RotatedAxis()+
      ggtitle("Top 5 Marker Genes Expression Across Clusters")+
      theme(plot.title = element_text(size = 20, face = "bold",hjust = 0.5))

    print("表达差异分析完成. ")


    ######################################## SingleR注释 ############################################

    # 对聚类结果用 singleR 标注细胞类型
    print("开始进行 SingleR 细胞类型注释")

    # 加载和准备参考数据集
    ref_list <- list()
    label_list <- list()

    # 检查并加载 celldex 提供的参考数据集
    for (ref_num in anno_ref) {
      ref <- switch(as.character(ref_num),
        `1` = celldex::HumanPrimaryCellAtlasData(),
        `2` = celldex::MouseRNAseqData(),
        `3` = celldex::BlueprintEncodeData(),
        `4` = celldex::DatabaseImmuneCellExpressionData(),
        `5` = celldex::ImmGenData(),
        `6` = celldex::MonacoImmuneData(),
        `7` = celldex::NovershternHematopoieticData(),
        `8` = {
          print("使用指定的自定义参考数据集.")
          ref_data <- read.table(count_matrix_path, header = TRUE, row.names = 1)
          ref_label <- read.table(cell_annotation_path, header = TRUE, row.names = 1, sep = "\t")
          
          custom_ref <- SummarizedExperiment(
            assays = list(counts = as.matrix(ref_data)),
            colData = DataFrame(cell_type = as.vector(ref_label$Cell_subtype))
          )
          
          # 并行计算 log1p
          parallel_log1p <- function(counts_matrix) {
            print("正在计算表达矩阵的log1p.")
            use_condaenv("rkit", required = TRUE)
            torch <- import("torch")
            
            has_gpu <- torch$cuda$is_available()
            
            if (has_gpu) {
              print("检测到GPU, 使用GPU进行计算.")
              counts_tensor <- torch$tensor(counts_matrix, device = "cuda")
              logcounts_tensor <- torch$log1p(counts_tensor)
              logcounts_matrix <- as.matrix(logcounts_tensor$cpu()$numpy())
            } else {
              print("未检测到GPU, 使用CPU并行计算.")
              counts_tensor <- torch$tensor(counts_matrix, device = "cpu")
              logcounts_tensor <- torch$log1p(counts_tensor)
              logcounts_matrix <- as.matrix(logcounts_tensor$numpy())
            }
            
            return(logcounts_matrix)
          }
          
          logcounts_matrix <- parallel_log1p(assay(custom_ref, "counts"))
          assay(custom_ref, "logcounts", withDimnames = FALSE) <- logcounts_matrix
          print("log1p 计算完成")
          
          ref_list[[length(ref_list) + 1]] <- custom_ref
          label_list[[length(label_list) + 1]] <- ref_label$Cell_subtype  # 使用自定义的标签
        },
        NULL  # 如果 ref_num 不匹配上述情况，则返回 NULL
      )
      
      if (!is.null(ref)) {
        ref_list[[length(ref_list) + 1]] <- ref
        
        # 添加相应的标签到 label_list
        label_list[[length(label_list) + 1]] <- if ("label.fine" %in% colnames(colData(ref))) {
          ref$label.fine
        } else {
          ref$label.main  # 或使用其他合适的标签列
        }
      }
    }

    if(length(ref_list) == 1){
        ref_list = ref_list[[1]]
        label_list = label_list[[1]]
    }

    # 获取聚类信息
    clusters <- Idents(seurat_obj)

    # 使用 SingleR 进行注释，确保 ref 和 labels 列表长度一致
    singleR_results <- SingleR(
      test = seurat_obj@assays$RNA$data, 
      ref = ref_list,  # 确保 ref_list 为一个列表
      labels = label_list,  # 确保 label_list 与 ref_list 长度一致
      clusters = clusters,
      de.method = "wilcox"
    )

    # 把SingleR的注释结果放到cell_type中
    seurat_obj$cell_type <- singleR_results$labels[match(clusters, rownames(singleR_results))]

    # 绘制UMAP图，根据cell_type分组
    p3 <- DimPlot(seurat_obj, reduction = "umap", group.by = "cell_type") + 
      ggtitle("Cell Type Annotation on UMAP")+
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20))
    
    # 画出 SingleR 注释准确度的得分图。
    #p4 <- plotScoreHeatmap(singleR_results)

    # 注释完成之后，看一下细胞类型的组成：
    # 计算每个样本内细胞类型的比例
    cell_type_proportions <- seurat_obj@meta.data %>%
        group_by(orig.ident, cell_type) %>%
        summarise(count = n()) %>%
        mutate(proportion = count / sum(count))

    # 可视化每个样本内细胞类型的比例
    p5 <- ggplot(cell_type_proportions, aes(x = orig.ident, y = proportion, fill = cell_type)) +
        geom_bar(stat = "identity", position = "fill") +
        theme_minimal() +
        labs(x = "Sample", y = "Proportion", fill = "Cell Type") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1))+ 
        plot_annotation(
            title = "Proportion of each cell type",
            theme = theme(
            plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
            )
        )

    print("SingleR 注释完成.")


    ###################################### Monocle 轨迹推断 ###########################################
    
    print("开始进行 Monocle 轨迹推断.")

    data <- seurat_obj@assays$RNA$counts
    gene_annotation <- data.frame(gene_short_name = rownames(data))
    rownames(gene_annotation) <- rownames(data)

    cds <- new_cell_data_set(data, cell_metadata = seurat_obj@meta.data, gene_metadata = gene_annotation)
    reducedDims(cds)$UMAP <- Embeddings(seurat_obj, "umap")
    cds <- cluster_cells(cds)
    cds <- learn_graph(cds)
    p6 <- plot_cells(cds, color_cells_by = "cell_type") + 
      ggtitle("Cell Trajectory Inference Map")+
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20))

    print("轨迹推断完成.")


    ######################################## 功能富集分析 ###################################################
    
    print("开始进行功能富集分析.") 

    # 功能富集分析
    DE_results <- all_markers
    DEG <- DE_results$gene[DE_results$p_val_adj < 0.05]

    enrich_results <- enrichGO(gene = DEG, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", 
                              pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, 
                              minGSSize = 10, maxGSSize = 500, readable = TRUE) 

    # 可视化富集分析结果
    p7 <- dotplot(enrich_results, showCategory = 10) + 
      ggtitle("GO Enrichment Analysis")+
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20))

    print("功能富集分析完成.")


    ######################################## 信号通路分析 ###################################################

    print("开始进行信号通路分析. ")

    # 把基因名称从 symbol 转为 enterzid
    DEG <- bitr(DEG, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

    if(species == 1){
      organism = "human"
    }else if(species == 2){
      organism = "mouse"
    }

    # 信号通路分析
    pathway_results <- enrichPathway(gene = DEG$ENTREZID, organism = organism, 
                                    pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, 
                                    minGSSize = 10, maxGSSize = 500, readable = TRUE) 

    # 可视化信号通路分析结果
    p8 <- barplot(pathway_results, showCategory = 10) + 
      ggtitle("Pathway Enrichment Analysis") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20))

    print("信号通路分析完成. ")

    saveRDS(seurat_obj, "./init_data/seurat_obj.rds")


    ###################################### 保存图于pdf ############################################


    # 将图保存于 pdf 中
    pdf("./report4.pdf", width = 10, height = 8)
      print(p1)
      print(p2)
    dev.off()

    pdf("./report5.pdf", width = 12, height = 8)
      print(p3)
    dev.off()

    pdf("./report6.pdf", width = 12, height = 14)
      # print(p4)
      plotScoreHeatmap(singleR_results)
    dev.off()

    file.remove("./Rplot.pdf")
    
    pdf("./report7.pdf")
      print(p5)
      print(p6)
      print(p7)
      print(p8)
    dev.off()

   
    # 合并所有的pdf
    # 获取当前目录下所有以 "report" 开头的 PDF 文件
    pdf_report_files <- list.files(pattern = "\\.pdf$")

    # 排除文件名为 "report.pdf" 的文件
    pdf_report_files <- pdf_report_files[pdf_report_files != "report.pdf"]

    # 合并这些PDF文件
    if (length(pdf_report_files) > 0) {
      pdf_combine(pdf_report_files, "report.pdf")
      print("PDF文件已合并到report.pdf\n")
      
      # 删除所有以 "report" 开头的 PDF 文件（除了最终的 "report.pdf"）
      file.remove(pdf_report_files)
      
    } else {
      print("当前目录下没有找到以'report'开头的PDF文件。\n")
    }

    print("所有模块运行完成")
}












# done by xzy at 2024/8/11, 11:20 am.