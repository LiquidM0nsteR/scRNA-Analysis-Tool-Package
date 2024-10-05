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
library(ggraph)
library(SummarizedExperiment)
library(clusterProfiler)
library(ReactomePA)
library(org.Hs.eg.db)
library(pdftools)
library(CellChat)
library(GENIE3)
library(doRNG)



# 下游分析函数，包括差异表达分析、细胞注释等多个流程
downstream_analysis <- function(seurat_obj, anno_ref, count_matrix_path = NULL, cell_annotation_path =NULL, species, batch_key, use_spatial_coords){

    DefaultAssay(seurat_obj) <- "RNA"
    Idents(seurat_obj) <- seurat_obj$seurat_clusters


    ####################################### 差异表达分析 ###############################################


    print("进行差异表达分析. ")

    # 配置并行计算
    num_cores <- max(1, parallel::detectCores() - 2)
    plan("multicore", workers = num_cores) # 根据实际CPU核心数量设置workers数量

    # 使用FindAllMarkers进行标记基因计算
    all_markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)

    # 恢复为单线程
    plan("sequential")

    write.csv(all_markers, "./cache/all_markers.csv")

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
        ggtitle("Top 10 Marker Genes Expression Across Clusters")+
        theme(plot.title = element_text(size = 20, face = "bold",hjust = 0.5))

    print("表达差异分析完成. ")


    ######################################## SingleR注释 ############################################


    # 对聚类结果用 singleR 标注细胞类型
    print("开始进行 SingleR 细胞类型注释")

    mapped_ref_list <- list()
    for (ref_num in anno_ref) {
        mapped_ref_num <- switch(as.character(species),
            `1` = switch(as.character(ref_num),
                `1` = 1,  # HumanPrimaryCellAtlasData
                `2` = 3,  # BlueprintEncodeData
                `3` = 4,  # DatabaseImmuneCellExpressionData
                `4` = 6,  # MonacoImmuneData
                `5` = 7,  # NovershternHematopoieticData
                `6` = 8   # 自定义参考集
            ),
            `2` = switch(as.character(ref_num),
                `1` = 2,  # MouseRNAseqData
                `2` = 5,  # ImmGenData
                `3` = 8   # 自定义参考集
            ),
            `3` = 8,  # 自定义参考集
            NULL  # 如果 species 不匹配，则返回 NULL
        )
        # 添加映射后的参考编号到列表
        if (!is.null(mapped_ref_num)) {
            mapped_ref_list[[length(mapped_ref_list) + 1]] <- mapped_ref_num
        } else {
            stop("无效的参考集编号。")
        }
    }

    # 加载和准备参考数据集
    ref_list <- list()
    label_list_main <- list()
    label_list_fine <- list()

    # 检查并加载 celldex 提供的参考数据集
    for (ref_num in mapped_ref_list) {
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
                ref_matrix_anno = read_file(count_matrix_path, cell_annotation_path)
                ref_data <- ref_matrix_anno$ref_data
                ref_label <- ref_matrix_anno$ref_label
                
                custom_ref <- SummarizedExperiment(
                  assays = list(counts = as.matrix(ref_data)),
                  colData = DataFrame(cell_type = as.vector(ref_label))
                )
                

                # 并行计算 log1p
                parallel_log1p <- function(counts_matrix) {

                    print("正在计算表达矩阵的log1p.")
                    use_python(Sys.which("python"), required = TRUE)
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
                print("log1p 计算完成. ")
                
                ref_list[[length(ref_list) + 1]] <- custom_ref
                label_list_main[[length(label_list_main) + 1]] <- ref_label  # 使用自定义的标签
                label_list_fine[[length(label_list_fine) + 1]] <- ref_label  # 使用自定义的标签
            },
            NULL  # 如果 ref_num 不匹配上述情况，则返回 NULL
        )
        

        if (ref_num < 8 && !is.null(ref)) {
            ref_list[[length(ref_list) + 1]] <- ref
            # 添加相应的标签到 label_list
            label_list_main[[length(label_list_main) + 1]] <- ref$label.main
            label_list_fine[[length(label_list_fine) + 1]] <- ref$label.fine
        }
    }

    if(length(ref_list) == 1){
        ref_list = ref_list[[1]]
        label_list_main = label_list_main[[1]]
        label_list_fine = label_list_fine[[1]]
    }

    # 获取聚类信息
    clusters <- Idents(seurat_obj)

    # 使用每个群集前10个显著的标记基因进行SingleR标注
    top_markers <- all_markers %>%
        group_by(cluster) %>%
        slice_max(order_by = avg_log2FC, n = 10) %>%
        ungroup()

    selected_genes <- unique(top_markers$gene)

    print(paste("使用 ", length(selected_genes), " 个 Marker 基因进行细胞注释. "))

    # 设置并行参数，使用多核心（保留两个核心给系统）
    nCores <- max(detectCores() - 2, 1)
    BPPARAM <- MulticoreParam(workers = nCores)

    # 使用 SingleR 进行细胞类型注释，并使用并行计算
    # 第一步：使用main labels进行粗粒度注释
    singleR_main <- SingleR(
        test = seurat_obj@assays$RNA$data[selected_genes, ],  
        ref = ref_list,  
        labels = label_list_main,  # 注释集为main标签
        clusters = clusters,  # 按聚类结果注释
        method = "cluster",  
        de.method = "wilcox",  
        BPPARAM = BPPARAM  # 并行参数
    )

    # 提取 main 注释结果（按簇的注释结果）
    main_labels <- singleR_main$labels  # 名称为 cluster 的编号
    names(main_labels) <- levels(clusters)

    # 根据注释结果对 clusters 进行分组
    same_type_clusters_list <- split(names(main_labels), main_labels)

    # 找到所有注释相同的不同 clusters
    multi_cluster_types <- names(same_type_clusters_list)[sapply(same_type_clusters_list, length) > 1]

    # 汇总所有具有相同注释的 clusters
    all_same_type_clusters <- unlist(same_type_clusters_list[multi_cluster_types])

    # 获取这些 clusters 对应的所有细胞索引
    cell_indices <- which(as.character(clusters) %in% all_same_type_clusters)

    print(paste("聚类不同但注释相同的 clusters 包含", length(cell_indices), "个细胞，将进行细粒度注释. "))

    if (length(cell_indices) > 0) {
        # 第二步：对这些细胞进行细化注释
        singleR_fine <- SingleR(
            test = seurat_obj@assays$RNA@data[selected_genes, cell_indices],  # 只对这些细胞的数据进行注释
            ref = ref_list,  
            labels = label_list_fine,  # 注释集为 fine 标签
            clusters = as.character(clusters[cell_indices]),  # 使用相应细胞的 clusters
            method = "cluster",  
            de.method = "wilcox",  
            BPPARAM = BPPARAM   # 并行参数
        )
    }

    # 使用 singleR_main 的注释作为初始结果
    seurat_obj$cell_type <- singleR_main$labels[match(clusters, rownames(singleR_main))]

    # 覆盖 singleR_fine 有注释的部分
    fine_match <- match(clusters, rownames(singleR_fine))
    seurat_obj$cell_type[!is.na(fine_match)] <- singleR_fine$labels[fine_match[!is.na(fine_match)]]

    # 绘制UMAP图，根据cell_type分组
    p3 <- DimPlot(seurat_obj, reduction = "umap", group.by = "cell_type") + 
        ggtitle("Cell Type Annotation on UMAP")+
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20))
    

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
    p6 <- monocle3::plot_cells(cds, color_cells_by = "cell_type", 
                    cell_size = 0.5,  # 调整细胞节点大小
                    graph_label_size = 6,  # 调整轨迹图的标签大小
                    group_label_size = 6  # 调整细胞类型标签的字体大小
                    ) + 
        ggtitle("Cell Trajectory Inference Map") +
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20))

    print("轨迹推断完成.")


    ######################################## CellChat 细胞通讯分析 ###############################################
    
    if(species != 3){ # cellchat目前只有人和老鼠的数据集，其他物种请自定义 CellChatDB

        print("开始进行细胞通讯分析.") 

        # 使用 species 参数确定数据库
        CellChatDB <- if (species == 1) CellChatDB.human else CellChatDB.mouse

        # Step 1: 创建 CellChat 对象
        if (!is.null(batch_key)) {
            # 根据 batch_key 将 Seurat 对象分割为多个对象
            seurat_list <- SplitObject(seurat_obj, split.by = batch_key)

            # 为每个分割的 Seurat 对象创建 CellChat 对象并处理数据
            cellchat_list <- lapply(seurat_list, function(x) {
                x@meta.data$cell_type <- factor(x@meta.data$cell_type)  # 转换为因子
                x@meta.data$cell_type <- droplevels(x@meta.data$cell_type)  # 去除未使用的因子水平
                
                # 创建 CellChat 对象，考虑 use_spatial_coords 参数
                if (use_spatial_coords) {
                    cellchat_sample <- createCellChat(object = x, group.by = "cell_type", coordinates = as.matrix(x@meta.data[, c("x", "y")]))
                } else {
                    cellchat_sample <- createCellChat(object = x, group.by = "cell_type")
                }

                cellchat_sample@DB <- CellChatDB
                cellchat_sample <- subsetData(cellchat_sample)
                return(cellchat_sample)
            })

            # 合并所有已处理的 CellChat 对象
            cellchat <- mergeCellChat(cellchat_list, add.names = names(cellchat_list))
            # 一般idents在 mergeCellChat 后为一个 list，最后一个 idents 就是合并之后的idents
            cellchat@idents <- cellchat@idents[[length(cellchat@idents)]]

        } else {
            # 如果没有 batch_key，仅根据 use_spatial_coords 处理单个 Seurat 对象
            if (use_spatial_coords) {
                cellchat <- createCellChat(object = seurat_obj, group.by = "cell_type", coordinates = as.matrix(seurat_obj@meta.data[, c("x", "y")]))
            } else {
                cellchat <- createCellChat(object = seurat_obj, group.by = "cell_type")
            }
            cellchat@DB <- CellChatDB
            cellchat <- subsetData(cellchat)
        }
        
        # 除去重复的idents
        cellchat@idents <- droplevels(cellchat@idents)

        # Step 2: 鉴定过表达的基因和配体-受体
        cellchat <- identifyOverExpressedGenes(cellchat)
        cellchat <- identifyOverExpressedInteractions(cellchat)

        # Step 3: 计算通讯概率和信号通路
        cellchat <- computeCommunProb(cellchat)
        cellchat <- filterCommunication(cellchat, min.cells = 10)
        cellchat <- computeCommunProbPathway(cellchat)
        cellchat <- aggregateNet(cellchat)

        # Step 4: 可视化细胞通讯网络
        groupSize <- as.numeric(table(seurat_obj@meta.data$cell_type))
        
        # heatmap 可视化
        p7 <- netVisual_heatmap(cellchat, measure = "count", 
                                font.size = 10, font.size.title = 20, 
                                title.name = "Cell Communication Heatmap")

        print("细胞通讯分析完成.")
    }


    ####################################### Genie3 基因调控网络分析 ###############################################
    
    
    print("开始进行基因调控网络分析. ")

    # 准备 top 5% Marker 基因的基因表达矩阵
    expr_matrix <- as.matrix(seurat_obj@assays$RNA$counts[selected_genes, ])

    # 使用并行化计算加速 GENIE3 推断
    nCores <- max(detectCores() - 2, 1)  # 自动检测可用核心数，并使用所有核心减去2个
    registerDoParallel(cores = nCores)
    
    # 使用 GENIE3 推断基因调控网络，并指定并行核心数
    weight_matrix <- GENIE3(expr_matrix, nCores = nCores)

    # 完成后注销并行环境
    stopImplicitCluster()

    # 将结果转换为 igraph 对象，方便网络图展示
    g <- graph_from_adjacency_matrix(weight_matrix, weighted = TRUE, mode = "directed")

    # 获取前10个权重最高的边
    E(g)$weight <- E(g)$weight
    top_edges <- head(sort(E(g)$weight, decreasing = TRUE), 10)

    # 过滤网络图中权重较低的边，只保留前10个
    g_top <- subgraph.edges(g, which(E(g)$weight %in% top_edges))

    # 使用 ggraph 和 igraph 绘制网络图
    p8 <- ggraph(g_top, layout = 'fr', niter = 1200, area = vcount(g_top)^1.2) +  
        geom_edge_link(aes(edge_alpha = weight, edge_width = weight), show.legend = TRUE, alpha = 0.3) +  
        geom_node_point(size = 2, color = "steelblue") +  
        geom_node_text(aes(label = name), repel = TRUE, size = 5, color = "black", nudge_y = 0.2, nudge_x = 0.2) +  
        theme_void() +  
        ggtitle("Top 10 Gene Regulatory Network Interactions") +  
        theme(
            plot.title = element_text(hjust = 0.5, face = "bold", size = 20),  
            legend.position = "right",  # 将图例移到右侧
            legend.text = element_text(size = 14),  
            legend.title = element_text(size = 14),  
            legend.key.size = unit(3, "lines")  # 调整图例的大小
        )

    print("基因调控网络分析完成.")


    ######################################## 功能富集分析 ###################################################
    

    if(species != 3){ # 功能富集分析目前只支持人和老鼠，其他物种请自定义 OrgDb 参数

        print("开始进行功能富集分析.") 

        # 功能富集分析
        DE_results <- all_markers
        DEG <- DE_results$gene[DE_results$p_val_adj < 0.05]

        enrich_results <- enrichGO(gene = DEG, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP", 
                                pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, 
                                minGSSize = 10, maxGSSize = 500, readable = TRUE) 

        # 可视化富集分析结果
        p9 <- dotplot(enrich_results, showCategory = 10) + 
        ggtitle("GO Enrichment Analysis")+
        theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20))

        print("功能富集分析完成.")
    }

    ######################################## 信号通路分析 ###################################################


    if(species != 3){ # 信号通路分析目前只支持人和老鼠，其他物种请自定义 OrgDb 参数

        print("开始进行信号通路分析. ")

        if(species == 1){
            organism = "human"
        }else if(species == 2){
            organism = "mouse"
        }

        # 把基因名称从 symbol 转为 enterzid
        DEG <- bitr(DEG, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)

        # 信号通路分析
        pathway_results <- enrichPathway(gene = DEG$ENTREZID, organism = organism, 
            pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.05, 
            minGSSize = 10, maxGSSize = 500, readable = TRUE) 

        # 可视化信号通路分析结果
        p10 <- barplot(pathway_results, showCategory = 10) + 
            ggtitle("Pathway Enrichment Analysis") +
            theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 20))

        print("信号通路分析完成. ")
    }

    saveRDS(seurat_obj, "./init_data/seurat_obj.rds")


    ###################################### 保存图于pdf ############################################
    

    plot_index = max(0.9, length(unique(clusters)) / 16)

    # 将图保存于 pdf 中
    pdf("./report4.pdf", width = 12 * plot_index , height = 10 * plot_index)
        print(p1)
        print(p2)
    dev.off()

    pdf("./report5.pdf", width = 12 * plot_index, height = 8 * plot_index)
        print(p3)
    dev.off()

    pdf("./report6.pdf", width = 12 * plot_index * (1 + 0.3 * (length(anno_ref) %/% 3)), height = 14 * plot_index * length(anno_ref))
        plotScoreHeatmap(results = singleR_main, max.labels = 25)
        file.remove("./Rplot.pdf")
        if(!is.null(singleR_fine)){
            plotScoreHeatmap(singleR_fine, max.labels = 25)
            file.remove("./Rplot.pdf")
        }
    dev.off()

    
    pdf("./report7.pdf", width = 12 * plot_index , height = 10 * plot_index)
        print(p5)
        print(p6)
        print(p7)
        print(p8)
        print(p9)
        print(p10)
    dev.off()

   
    # 合并所有的pdf
    # 获取当前目录下所有以 "report" 开头的 PDF 文件
    pdf_report_files <- list.files(pattern = "\\.pdf$")

    # 排除文件名为 "report.pdf" 的文件
    pdf_report_files <- pdf_report_files[pdf_report_files != "report.pdf"]

    # 合并这些PDF文件
    if (length(pdf_report_files) > 0) {
        pdf_combine(pdf_report_files, "report.pdf")
        print("PDF文件已合并到report.pdf \n")
        
        # 删除所有以 "report" 开头的 PDF 文件（除了最终的 "report.pdf"）
        file.remove(pdf_report_files)
      
    } else {
        print("当前目录下没有找到以'report'开头的PDF文件.\n")
    }

    print("所有模块运行完成.")
}



# 读取文件函数，接受两个文件路径，分别处理表达矩阵和细胞注释文件
read_file <- function(count_matrix_path, cell_annotation_path) {
  
    # 定义读取表达矩阵的逻辑
    count_matrix_ext <- tools::file_ext(count_matrix_path)
    ref_data <- switch(count_matrix_ext,
        "csv" = read.csv(count_matrix_path, header = TRUE, row.names = 1),
        "txt" = read.table(count_matrix_path, header = TRUE, sep = "\t", row.names = 1),
        "gz" = read.table(count_matrix_path, header = TRUE, sep = "\t", row.names = 1),
        "h5" = {
            if (!requireNamespace("rhdf5", quietly = TRUE)) {
                stop("请安装 rhdf5 包来读取 h5 文件.")
            }
            rhdf5::h5read(count_matrix_path, "/")
        },
        stop("不支持的表达矩阵文件类型: ", count_matrix_ext)
    )
    
    # 定义读取细胞注释的逻辑
    annotation_ext <- tools::file_ext(cell_annotation_path)
    ref_label <- switch(annotation_ext,
        "csv" = read.csv(cell_annotation_path, header = TRUE),
        "txt" = read.table(cell_annotation_path, header = TRUE, sep = "\t"),
        "gz" = read.table(cell_annotation_path, header = TRUE, sep = "\t"),
        stop("不支持的细胞注释文件类型: ", annotation_ext)
    )
    
    # 如果细胞注释文件是行向量，将其转置为列向量
    if (ncol(ref_label) == 1) {
        ref_label <- ref_label[[1]]  # 如果是单列，提取为列向量
    } else if (nrow(ref_label) == 1) {
        # 如果是单行，转置为列向量
        ref_label <- t(ref_label)[, 1]
    }
    
    # 返回包含表达矩阵和细胞注释的列表
    return(list(ref_data = ref_data, ref_label = ref_label))
}











# done by xzy at 2024/8/11, 11:20 am.