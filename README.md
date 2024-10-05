 
# LiquidMonsteR: Single-Cell Transcriptomics Analysis Pipeline

### Powered by LiquidMonsteR  
**From Hainan Institute of WuHan University of Technology, Sanya, Hainan**  
**First Created on July 7, 2024**

## Overview

LiquidMonsteR is a comprehensive R-based pipeline designed for analyzing single-cell transcriptomics data, including spatial transcriptomics. It provides functions for quality control, dimensionality reduction, clustering, and downstream analysis, ensuring a robust workflow tailored to various experimental conditions.

## Prerequisites

Before running the script, make sure you have installed all necessary R packages as listed in the `requirements.R` file.

 

## Input Data

Place your raw data files in the `init_data` directory. Supported file types include:

- `.rds` files containing Seurat objects.
- `.h5` files with expression matrices.
- Folders with spatial data folder and `.h5`file. **Ensure correct file naming**: Make sure your .h5 file is named ```filtered_feature_bc_matrix.h5```. If the file is not correctly named, the script will fail.



## Usage

1.**Navigate to the directory**: Change to the directory containing the scripts.

2.**Run the pipeline**: Execute the R script using Rscript.

- ```main.R```

You will be prompted to provide the following inputs:

1. **Import Seurat Object from RDS File**: Indicate whether to import a Seurat object from an RDS file (`TRUE` or `FALSE`).
2. **Spatial Transcriptomics Data**: Specify if the data includes spatial coordinates (`TRUE` or `FALSE`).
3. **Cell Type Annotation Reference Dataset**: Choose from available datasets or provide a custom reference dataset.
4. **Species Selection**: Choose between human (`1`) or mouse (`2`).
5. **Batch Correction Method**: Select a batch correction method from the provided options.
6. **Dimensionality Reduction Method**: Choose a dimensionality reduction method suitable for your data.
7. **Clustering Method**: Specify the clustering method to use.

## Workflow

The script executes the following steps:

1. **Initialization**: Loads the data into a Seurat object and processes spatial information if applicable.
2. **Quality Control**: Performs quality control checks on the data.
3. **Dimensionality Reduction**: Applies the selected dimensionality reduction method.
4. **Clustering**: Performs clustering based on the reduced dimensions.
5. **Downstream Analysis**: Conducts downstream analyses, including cell type annotation using reference datasets.

## References

The methods used in this pipeline are based on well-established scientific literature:

- **Batch Correction**: Seurat, Harmony, and ComBat methods.
- **Dimensionality Reduction**: PCA, c-ICA, scVI, SpatialPCA, and SeuratSpatial.
- **Clustering**: Louvain, Walktrap, and Model-Based Clustering.

  ###  Batch Correction Methods References

1. **Seurat Standard Workflow**  
   Stuart, T., Butler, A., et al. (2019). Comprehensive Integration of Single-Cell Data. *Cell*, 177(7), 1888-1902.e21. [Link to article](https://doi.org/10.1016/j.cell.2019.05.031) 
(GitHub: [Seurat](https://github.com/satijalab/seurat))


2. **Harmony**  
   Korsunsky, I., et al. (2019). Fast, sensitive and accurate integration of single-cell data with Harmony. *Nature Methods*, 16, 1289–1296. [Link to article](https://doi.org/10.1038/s41592-019-0619-0) (GitHub: [Harmony](https://github.com/immunogenomics/harmony))

3. **ComBat**  
   Johnson, W. E., Li, C., & Rabinovic, A. (2007). Adjusting batch effects in microarray expression data using empirical Bayes methods. *Biostatistics*, 8(1), 118-127. [Link to article](https://doi.org/10.1093/biostatistics/kxj037)  (ComBat itself does not have a standalone GitHub repository, but it is implemented as part of the `sva` R package.   GitHub: [sva-devel](https://github.com/zhangyuqing/sva-devel))

  ### Dimensionality Reduction Methods References

1. **PCA (Principal Component Analysis)**  
   Jolliffe, I. T. (2002). *Principal Component Analysis*. Springer Series in Statistics. Springer-Verlag, New York. [Link to book](https://doi.org/10.1007/b98835) (PCA is a classic method without a dedicated GitHub repository but is widely used in R and Python libraries.  Example: Implemented in [Seurat](https://github.com/satijalab/seurat))

2. **c-ICA (Independent Component Analysis)**  
   Teschendorff, A. E., et al. (2007). Elucidating the altered transcriptional programs in breast cancer using independent component analysis. *PLOS Computational Biology*, 3(8), e161. [Link to article](https://doi.org/10.1371/journal.pcbi.0030161) (c-ICA is typically part of Independent Component Analysis and does not have a dedicated GitHub repository, but similar functionality is provided by the `fastICA` package in R. GitHub: [fastICA](https://github.com/cran/fastICA)
)

3. **scVI (Variational Autoencoder)**  
   Lopez, R., Regier, J., et al. (2018). Deep generative modeling for single-cell transcriptomics. *Nature Methods*, 15(12), 1053-1058. [Link to article](https://doi.org/10.1038/s41592-018-0229-2) (GitHub: [scVI](https://github.com/YosefLab/scVI))

4. **SpatialPCA (Spatial Principal Component Analysis)**  
   Shang, L., & Zhou, X. (2022). Spatially Aware Dimension Reduction for Spatial Transcriptomics. *Nature Communications*. [Link to article](https://doi.org/10.1038/s41467-022-34879-1) (GitHub: [SpatialPCA](https://github.com/shangll123/SpatialPCA))

5. **SeuratSpatial (Seurat Standard Workflow)**  
   Satija, R., Farrell, J. A., Gennert, D., Schier, A. F., & Regev, A. (2015). Spatial reconstruction of single-cell gene expression data. *Nature Biotechnology*, 33(5), 495-502. [Link to article](https://doi.org/10.1038/nbt.3192) (GitHub: [Seurat](https://github.com/satijalab/seurat))

  ### Clustering Methods References 

1. **FindClusters (Louvain Algorithm)**  
   Blondel, V. D., Guillaume, J. L., Lambiotte, R., & Lefebvre, E. (2008). Fast unfolding of communities in large networks. *Journal of Statistical Mechanics: Theory and Experiment*, 2008(10), P10008. [Link to article](https://doi.org/10.1088/1742-5468/2008/10/P10008) (Implemented in [Seurat](https://github.com/satijalab/seurat))

2. **Walktrap**  
   Pons, P., & Latapy, M. (2005). Computing communities in large networks using random walks. *Journal of Graph Theory Algorithms and Applications*, 10(2), 191-218. [Link to article](https://doi.org/10.7155/jgaa.00124) (GitHub: [igraph](https://github.com/igraph/rigraph))

3. **Model-Based Clustering (Mclust)**  
   Fraley, C., & Raftery, A. E. (2002). Model-Based Clustering, Discriminant Analysis, and Density Estimation. *Journal of the American Statistical Association*, 97(458), 611-631. [Link to article](https://doi.org/10.1198/016214502760047131) (GitHub: [mclust](https://github.com/cran/mclust))

  ### Additional Tools

1. **Monocle3 (Trajectory Inference)**  
   GitHub: [Monocle3](https://github.com/cole-trapnell-lab/monocle3)

2. **SingleR (Cell Type Annotation)**  
   GitHub: [SingleR](https://github.com/dviraran/SingleR)

3. **clusterProfiler and ReactomePA (Functional Enrichment Analysis)**  
   - [clusterProfiler](https://github.com/YuLab-SMU/clusterProfiler)  
   - [ReactomePA](https://github.com/YuLab-SMU/ReactomePA)

## Required Resources of Running demo

- **Data Type**: Large-scale data
- **CPU**: 32 cores or more
- **Memory**: 256 GB or more
