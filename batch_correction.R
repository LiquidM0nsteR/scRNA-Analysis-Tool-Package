# batch_correction.R
# Powered by LiquidMonsteR
library(Seurat)
library(harmony)
library(sva)
library(reticulate)



# 封装Seurat标准去批次流程
run_seurat_integrating <- function(seurat_obj, batch_key) {
    seurat_list <- SplitObject(seurat_obj, split.by = batch_key)
    seurat_list <- lapply(seurat_list, FUN = function(x) {
        x <- NormalizeData(x)
        x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
        return(x)
    })
    
    anchors <- FindIntegrationAnchors(object.list = seurat_list, dims = 1:30)
    seurat_integrated <- IntegrateData(anchorset = anchors, dims = 1:30)
    
    DefaultAssay(seurat_integrated) <- "integrated"
    seurat_integrated <- ScaleData(seurat_integrated)

    seurat_obj[["integrated"]] <- seurat_integrated@assays$integrated

    print("去批次完成.")

    return(seurat_obj)
}



# 封装Harmony去批次方法
run_harmony <- function(seurat_obj, batch_key, reduction_name) {

    seurat_obj <- RunHarmony(seurat_obj, batch_key, reduction_name)

    # 将 harmony 的 embeddings 传递给指定的 reduction_name
    harmony_embeddings <- Embeddings(seurat_obj[["harmony"]])

    # 更改 embeddings 矩阵的列名，确保它与 reduction_name 一致
    colnames(harmony_embeddings) <- paste0(reduction_name, "_", 1:ncol(harmony_embeddings))

    # 替换 reduction_name 的 embeddings 为 harmony 的 embeddings
    seurat_obj[[reduction_name]]@cell.embeddings <- harmony_embeddings

    # 修改 key 以匹配新的 reduction_name
    seurat_obj[[reduction_name]]@key <- paste0(reduction_name, "_")

    print("去批次完成.")

    return(seurat_obj)
}



# 封装Combat去批次方法
run_combat <- function(seurat_obj, batch_key) {
  # 获取表达矩阵
  expr_matrix <- as.matrix(GetAssayData(seurat_obj, slot = "data"))
  batch <- seurat_obj@meta.data[[batch_key]]
  
  # Combat批次校正
  combat_corrected <- ComBat(dat = expr_matrix, batch = batch)
  
  # 将校正后的数据放回Seurat对象
  seurat_obj[["combat"]] <- CreateAssayObject(counts = combat_corrected)

  DefaultAssay(seurat_obj) <- "combat"

  print("去批次完成.")

  return(seurat_obj)
}



run_deepMNN <- function(seurat_obj, batch_key, reduction_name) {

    # 定义网络参数
    n_blocks = 2          # 指定 ResNet 中的残差块
    learning_rate = 1e-3  # 初始学习率
    max_epochs = 200      # 最大迭代次数
    batch_size = 32       # 每次梯度更新时使用的样本数
    decay_step_size = 20  # 学习率每隔多少个 epoch 衰减
    lr_decay_factor = 0.8 # 学习率衰减因子
    min_lr = 1e-6         # 最小学习率

    # 1. 提取降维后的矩阵和批次信息
    embedding_matrix <- as.matrix(Embeddings(seurat_obj, reduction = reduction_name))
    batch <- seurat_obj@meta.data[[batch_key]]

    # 设置 Python 环境
    use_python(Sys.which("python"), required = TRUE)

    # 使用 reticulate 执行 Python 代码
    py_run_string("
import torch
import numpy as np
from sklearn.preprocessing import normalize

# 定义 ResNet block
class ResnetBlock(torch.nn.Module):
    def __init__(self, dim):
        super(ResnetBlock, self).__init__()
        self.block = self.build_resnet_block(dim)
        
    def build_resnet_block(self, dim):
        block = [torch.nn.Linear(dim, 2*dim),
                 torch.nn.BatchNorm1d(2*dim),
                 torch.nn.PReLU()]
        block += [torch.nn.Linear(2*dim, dim),
                  torch.nn.BatchNorm1d(dim),
                  torch.nn.PReLU()]
        return torch.nn.Sequential(*block)
    
    def forward(self, x):
        out = x + self.block(x)
        return out

# 定义网络
class Net(torch.nn.Module):
    def __init__(self, input_dim, n_blocks, device):
        super(Net, self).__init__()
        model = []
        for i in range(n_blocks):
            model += [ResnetBlock(input_dim)]
        self.model = torch.nn.Sequential(*model)
        self.model.to(device=device)
        
    def forward(self, input):
        out = self.model(input)
        return out

# 学习率更新函数
def update_lr(optim, epoch, init_lr, min_lr=1e-6, decay_step_size=20, lr_decay_factor=0.8):
    exponent = int(np.floor((epoch + 1) / decay_step_size))
    lr = init_lr * np.power(lr_decay_factor, exponent)
    if lr < min_lr:
        optim.param_groups[0]['lr'] = min_lr
    else:
        optim.param_groups[0]['lr'] = lr
    print('Learning rate = %.7f' % optim.param_groups[0]['lr'])

# 主批次校正逻辑
def run_deepMNN_py(embedding_matrix, batch, n_blocks, learning_rate = 1e-6, max_epochs = 200, batch_size = 32, decay_step_size = 20, lr_decay_factor = 0.8, min_lr = 1e-6):
    
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    # 基于 batch 信息构造 indices 来分割 embedding_matrix
    datasets = []
    unique_batches = np.unique(batch)  # 获取批次的唯一值
    for b in unique_batches:
        indices = [i for i, batch_val in enumerate(batch) if batch_val == b]
        if len(indices) == 0:
            continue  # 跳过空的批次
        datasets.append(embedding_matrix[indices, :])
    
    if len(datasets) == 0:
        raise ValueError('所有批次均为空，无法继续执行批次校正.')

    input_dim = datasets[0].shape[1]
    net = Net(input_dim, int(n_blocks), device).to(device)
    optimizer = torch.optim.Adam(net.parameters(), lr=learning_rate)

    losses = []
    net.train(True)
    
    embedding_tensor = torch.Tensor(embedding_matrix).to(device)
    data_loader = torch.utils.data.DataLoader(list(range(embedding_tensor.size(0))), batch_size=int(batch_size), shuffle=True)
    
    for epoch in range(int(max_epochs)):
        batch_losses = []
        for batch_indices in data_loader:
            batch_data = embedding_tensor[batch_indices, :]
            
            # 计算损失：直接比较校正后的嵌入矩阵与原始的嵌入矩阵
            corrected_batch = net(batch_data)
            loss = torch.nn.functional.mse_loss(corrected_batch, batch_data)
            
            optimizer.zero_grad()
            loss.backward()
            optimizer.step()
            
            batch_losses.append(loss.item())

        epoch_loss = np.mean(batch_losses)
        losses.append(epoch_loss)
        print(f'Epoch {epoch + 1}/{max_epochs} - Loss: {epoch_loss:.4f}')

        # 每个 epoch 结束后更新学习率
        update_lr(optimizer, epoch, learning_rate, min_lr, decay_step_size, lr_decay_factor)

    net.eval()
    corrected_embeddings = net(embedding_tensor).cpu().detach().numpy()
    
    return corrected_embeddings
")

    # 将整个 embedding_matrix、batch 和 cells 传递给 Python 环境
    py$embedding_matrix <- embedding_matrix
    py$batch <- batch

    # 调用 Python 函数进行批次校正
    py$corrected_embeddings <- py$run_deepMNN_py(py$embedding_matrix, py$batch, n_blocks, learning_rate, max_epochs, batch_size, decay_step_size, lr_decay_factor, min_lr)

    # 将校正后的嵌入矩阵放回 Seurat 对象
    corrected_embeddings <- py$corrected_embeddings
    colnames(corrected_embeddings) <- paste0(reduction_name, "_", 1:ncol(corrected_embeddings))
    rownames(corrected_embeddings) <- rownames(embedding_matrix)
    seurat_obj[[reduction_name]]@cell.embeddings <- corrected_embeddings

    # 修改 key 以匹配新的 reduction_name
    seurat_obj[[reduction_name]]@key <- paste0(reduction_name, "_")

    print("去批次完成.")

    return(seurat_obj)
}
