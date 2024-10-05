# requirements.R
# Powered by LiquidMonsteR

# 必需的R包列表
# 如果还需要安装别的R包，直接在下面添加，然后命令行运行Rscript requirements.R即可
required_packages <- c(
  "Seurat",
  "SingleR",
  "ica",
  "reticulate",
  "mclust",
  "gridExtra",
  "stringdist",
  "future",
  "scran",
  "celldex",
  "Matrix",
  "cowplot",
  "scales",
  "devtools",
  "clusterProfiler",
  "ReactomePA",
  "org.Hs.eg.db",
  "doParallel",
  "foreach",
  "hdf5r",
  "rhdf5",
  "Rfast2",
  "FNN",
  "pdftools",
  "harmony",
  "sva",
  "glmGamPoi",
  "GENIE3",
  "hgu133plus2.db",
  "ComplexHeatmap",
  "NMF",
  "doRNG",
  "ggraph"
) 


required_pip_packages <- list(
  list(name = "torch", version = "2.3.1"),
  list(name = "scvi-tools", version = "1.1.5"),
  list(name = "scVelo", version = "0.3.2")
)


##################################################################################################


# 安装和加载R包
install_and_load <- function(packages) {
  conda_channels <- c("conda-forge", "anaconda", "bioconda", "r")
  cran_mirrors <- c("https://mirrors.tuna.tsinghua.edu.cn/CRAN/",
                    "https://mirrors.huaweicloud.com/CRAN/",
                    "https://cloud.r-project.org/")

  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      installed <- FALSE
      # 尝试使用 conda 安装
      for (channel in conda_channels) {
        message("尝试使用 ", channel, " 安装包: ", pkg)
        system(paste("conda install -c", channel, "-y r-", pkg, sep = " "))
        tryCatch({
          suppressMessages(library(pkg, character.only = TRUE))
          installed <- TRUE
          break
        }, error = function(e) {
          message("通过 ", channel, " 安装包失败: ", pkg)
        })
      }

      # 如果 conda 安装失败，尝试使用 CRAN 安装
      if (!installed) {
        for (cran_mirror in cran_mirrors) {
          message("尝试使用 CRAN 镜像 ", cran_mirror, " 安装包: ", pkg)
          options(repos = c(CRAN = cran_mirror))
          tryCatch({
            install.packages(pkg, dependencies = TRUE)
            suppressMessages(library(pkg, character.only = TRUE))
            installed <- TRUE
            break
          }, error = function(e) {
            message("使用 CRAN 镜像 ", cran_mirror, " 安装包失败: ", pkg)
          })
        }
      }

      # 如果 CRAN 安装失败，尝试使用 Bioconductor 安装
      if (!installed) {
        message("尝试使用 Bioconductor 安装包: ", pkg)
        tryCatch({
          BiocManager::install(pkg, dependencies = TRUE, update = TRUE)
          suppressMessages(library(pkg, character.only = TRUE))
          installed <- TRUE
        }, error = function(e) {
          message("使用 Bioconductor 安装包失败: ", pkg)
        })

        # 尝试使用备用的 Bioconductor 镜像
        if (!installed) {
          for (mirror in c(6, 7, 4, 1)) {
            message("尝试备用镜像: ", mirror)
            chooseBioCmirror(ind = mirror) 
            tryCatch({
              BiocManager::install(pkg, dependencies = TRUE)
              suppressMessages(library(pkg, character.only = TRUE))
              installed <- TRUE
              break
            }, error = function(e) {
              message("备用镜像安装包失败: ", mirror)
            })
          }
        }
      }

      # 如果所有方法均失败，提示用户
      if (!installed) {
        stop("无法安装和加载包: ", pkg)
      }
    } else {
      suppressMessages(library(pkg, character.only = TRUE))
    }
  }
}



# 检查并安装 pip 包
check_and_install_pip_package <- function(pkg_info) {
  installed <- system(paste("pip show", pkg_info$name), intern = TRUE)
  if (length(installed) == 0) {
    message("安装 pip 包: ", pkg_info$name)
    system(paste("pip3 install", pkg_info$name))
  } else {
    installed_version <- sub("Version: ", "", grep("Version:", installed, value = TRUE))
    if (!is.null(pkg_info$version) && installed_version != pkg_info$version) {
      message("更新 pip 包: ", pkg_info$name)
      system(paste("pip3 install", pkg_info$name))
    } else {
      message("pip 包已安装且版本符合要求: ", pkg_info$name)
    }
  }
}


####################################################################################################


# 设置镜像
options(repos = c(CRAN = "https://mirrors.huaweicloud.com/CRAN/"))
Sys.setenv(PIP_INDEX_URL = "https://pypi.tuna.tsinghua.edu.cn/simple")

# 安装bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", dependencies = TRUE)
}
library(BiocManager)


if (Sys.which("gsl-config") == "" && system("conda install -c conda-forge gsl -y", intern = FALSE) != 0) {
    stop("conda 安装 gsl 出错.")
}

# 检查并安装 pip 包
for (pkg in required_pip_packages) {
  check_and_install_pip_package(pkg)
}

# 设置bioconductor镜像
chooseBioCmirror(ind = 5) 
# 设置超时时间,单位为秒
options(timeout = 2000)
# 安装其他的R包
install_and_load(required_packages)



###########################################################################################################



# 安装monocle3等，并自动升级所有依赖包
if (!requireNamespace("monocle3", quietly = TRUE)){
    devtools::install_github('cole-trapnell-lab/monocle3', upgrade = "always", dependencies =TRUE)
}
if (!requireNamespace("presto", quietly = TRUE)){
    devtools::install_github('immunogenomics/presto', upgrade = "always", dependencies =TRUE)
}
if (!requireNamespace("SpatialPCA", quietly = TRUE)){
    devtools::install_github("shangll123/SpatialPCA", upgrade = "always", dependencies =TRUE)
}
if (!requireNamespace("circlize", quietly = TRUE)){
    devtools::install_github("jokergoo/circlize", upgrade = "always", dependencies =TRUE)
}
if (!requireNamespace("CellChat", quietly = TRUE)){
    devtools::install_github("jinworks/CellChat", upgrade = "always", dependencies =TRUE)
}
if (length(system("pip show SEDR", intern = TRUE)) == 0) {
  system("pip install git+https://github.com/JinmiaoChenLab/SEDR.git")
  sedr_installed <- system("pip show SEDR", intern = TRUE)
}

if (length(system("pip show SEDR", intern = TRUE)) == 0) {stop("SEDR 安装失败. ")}
suppressMessages(library(SpatialPCA))
suppressMessages(library(monocle3))
suppressMessages(library(presto))
suppressMessages(library(circlize))
suppressMessages(library(CellChat))