# requirements.R
# Powered by LiquidMonsteR
# List of required packages
required_packages <- c(
  "Seurat",
  "ggplot2",
  "dplyr",
  "FNN",
  "stringdist",
  "igraph",
  "cluster",
  "mclust",
  "clustree",
  "gridExtra",
  "grid",
  "ica",
  "autoencoder",
  "scater",
  "reticulate"
)


# Function to set CRAN mirror
set_cran_mirror <- function(primary = "https://mirrors.huaweicloud.com/CRAN/", secondary = "https://cran.csiro.au/") {
  options(repos = c(CRAN = primary))
  message("Using primary CRAN mirror: ", primary)
  tryCatch({
    available.packages()
  }, error = function(e) {
    message("Primary CRAN mirror failed, switching to secondary: ", secondary)
    options(repos = c(CRAN = secondary))
    available.packages()
  })
}

# Function to set Bioconductor mirror
set_bioc_mirror <- function() {
  options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
}

# Function to install Bioconductor if not already installed
install_bioc <- function() {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager", repos = getOption("repos"))
  }
  library(BiocManager)
}

# Function to install a package using conda
install_conda <- function(pkg) {
  tryCatch({
    system(paste("conda install -y r-", pkg, sep = ""))
  }, error = function(e) {
    message("Failed to install package using conda: ", pkg)
  })
}


# Function to install and load packages
install_and_load <- function(packages) {
  for (pkg in packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      tryCatch({
        install.packages(pkg, dependencies = TRUE)
      }, error = function(e) {
        # If install.packages fails, use BiocManager to install
        tryCatch({
          BiocManager::install(pkg)
        }, error = function(e) {
          # If BiocManager install fails, use conda to install
          install_conda(pkg)
        })
      })
      tryCatch({
        library(pkg, character.only = TRUE)
      }, error = function(e) {
        message("Failed to load package: ", pkg)
      })
    } else {
      library(pkg, character.only = TRUE)
    }
  }
}

# Set CRAN and Bioconductor mirrors
set_cran_mirror()
set_bioc_mirror()

# Install Bioconductor and then required packages
install_bioc()
install_and_load(required_packages)

# LquidMonsteR, all rights reserved.