# preprocess of seurat boj
# powered by LiquidMonsteR
library(Seurat)
library(stringdist)

# Function to check and preprocess Seurat object
check_and_preprocess_seurat <- function(seurat_obj, is_multi_sample, use_spatial_coords) {
  # Check if Seurat object is provided
  if (missing(seurat_obj) || !inherits(seurat_obj, "Seurat")) {
    stop("Please provide a valid Seurat object.")
  }
  
  # Set default assay to RNA if not already set
  DefaultAssay(seurat_obj) <- "RNA"
  
  # Required meta.data columns
  if(use_spatial_coords){
    required_cols <- c("cell_type", "x", "y")
  }else{
    required_cols <- c("cell_type")
  }

  
  # Additional columns required for multi-sample data
  if (is_multi_sample) {
    required_cols <- c(required_cols, "sample")
  }
  
  # Fuzzy matching function to find the closest match for required columns
  match_columns <- function(meta_data, required_cols) {
    matched_cols <- sapply(required_cols, function(req_col) {
      col_distances <- stringdist::stringdistmatrix(req_col, colnames(meta_data), method = "jw")
      closest_col <- colnames(meta_data)[which.min(col_distances)]
      return(closest_col)
    })
    names(matched_cols) <- required_cols
    return(matched_cols)
  }
  
  # Match required columns using fuzzy matching
  matched_cols <- match_columns(seurat_obj@meta.data, required_cols)
  
  # Rename matched columns to standard names
  for (req_col in required_cols) {
    col_name <- matched_cols[req_col]
    if (!req_col %in% colnames(seurat_obj@meta.data)) {
      seurat_obj@meta.data[[req_col]] <- seurat_obj@meta.data[[col_name]]
      seurat_obj@meta.data[[col_name]] <- NULL
    }
  }
  
  # Check if all required columns are present after matching
  if (!all(required_cols %in% colnames(seurat_obj@meta.data))) {
    stop("The Seurat object is missing one or more required columns in meta.data: ", 
         paste(setdiff(required_cols, colnames(seurat_obj@meta.data)), collapse = ", "))
  }
  
  if(use_spatial_coords){
    # Check the dimensions of the gene expression matrix and spatial position matrix
    gene_expr_matrix <- seurat_obj@assays$RNA$counts
    spatial_matrix <- seurat_obj@meta.data[, c("x", "y")]
    
    if (ncol(gene_expr_matrix) != nrow(spatial_matrix)) {
      stop("The number of cells in the gene expression matrix does not match the number of cells in the spatial position matrix.")
    }
    
    if (!all(colnames(gene_expr_matrix) == rownames(spatial_matrix))) {
      stop("The order of cells in the spatial position matrix does not match the order in the gene expression matrix.")
    }
  }
  return(seurat_obj)
}