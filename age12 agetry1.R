# =====================================================
# STEP 01：检查6只小鼠的原始数据
# =====================================================

# 1. 原始数据路径
raw_dir <- "E:/新建文件夹/expressions"

# 2. 新项目保存位置
out_dir <- paste0(
  "D:/R data/Age_vs_AgeTrp_scRNA/",
  "New_6Mice_3vs3_FromRaw"
)

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# 3. 本次使用的6只小鼠
samples <- c(
  "Age1",
  "Age2",
  "Age_Trp1",
  "Age_Trp2",
  "Age_Trp3",
  "Age_Trp4"
)

# 4. 实验分组
groups <- c(
  "Age",
  "Age",
  "Age",
  "Age_Trp",
  "Age_Trp",
  "Age_Trp"
)

# 5. 检查目录
cat("原始目录存在：", dir.exists(raw_dir), "\n")

if (!dir.exists(raw_dir)) {
  stop("原始数据目录不存在，请检查路径。")
}

# 6. 显示原始目录内容
cat("\n===== 原始目录文件 =====\n")

print(
  list.files(
    raw_dir,
    recursive = FALSE
  )
)

# 7. 检查每只鼠是否有同名目录
check_table <- data.frame(
  Sample = samples,
  Group = groups,
  SampleDirectory = file.path(raw_dir, samples),
  DirectoryExists = dir.exists(
    file.path(raw_dir, samples)
  )
)

cat("\n===== 六只小鼠检查 =====\n")

print(check_table)

# 8. 保存分组信息
write.csv(
  check_table,
  file.path(
    out_dir,
    "01_Sample_Information.csv"
  ),
  row.names = FALSE
)

# 9. 保存原始文件清单
all_files <- list.files(
  raw_dir,
  recursive = TRUE,
  full.names = TRUE
)

writeLines(
  all_files,
  file.path(
    out_dir,
    "01_Raw_File_List.txt"
  )
)

cat("\n第一步检查完成。\n")



# =====================================================
# STEP 02：检查原始10x表达矩阵
# 不进行过滤、不修改数据
# =====================================================

library(Seurat)

raw_dir <- "E:/新建文件夹/expressions"

out_dir <- paste0(
  "D:/R data/Age_vs_AgeTrp_scRNA/",
  "New_6Mice_3vs3_FromRaw"
)

samples <- c(
  "Age1",
  "Age2",
  "Age_Trp1",
  "Age_Trp2",
  "Age_Trp3",
  "Age_Trp4"
)

groups <- c(
  "Age", "Age", "Age",
  "Age_Trp", "Age_Trp", "Age_Trp"
)

results <- list()

for (i in seq_along(samples)) {
  
  sample_name <- samples[i]
  
  sample_dir <- file.path(
    raw_dir,
    sample_name
  )
  
  cat("\n====================================\n")
  cat("样本：", sample_name, "\n")
  cat("====================================\n")
  
  # 查看文件夹中的文件
  sample_files <- list.files(
    sample_dir,
    recursive = TRUE,
    full.names = TRUE
  )
  
  print(basename(sample_files))
  
  # 检查是否包含标准10x文件
  matrix_files <- sample_files[
    grepl(
      "^matrix\\.mtx(\\.gz)?$",
      basename(sample_files),
      ignore.case = TRUE
    )
  ]
  
  if (length(matrix_files) != 1) {
    stop(
      sample_name,
      "：无法唯一定位matrix.mtx文件，",
      "请把该样本的文件列表发给我。"
    )
  }
  
  matrix_dir <- dirname(matrix_files[1])
  
  # 读取10x矩阵
  counts <- Read10X(
    data.dir = matrix_dir
  )
  
  if (is.list(counts)) {
    if ("Gene Expression" %in% names(counts)) {
      counts <- counts[["Gene Expression"]]
    } else {
      stop(
        sample_name,
        "包含多个数据类型，需要确认RNA矩阵。"
      )
    }
  }
  
  results[[sample_name]] <- data.frame(
    Sample = sample_name,
    Group = groups[i],
    Genes = nrow(counts),
    Cells = ncol(counts),
    TotalUMI = sum(counts)
  )
  
  cat("基因数：", nrow(counts), "\n")
  cat("细胞数：", ncol(counts), "\n")
  cat("UMI总数：", sum(counts), "\n")
  
  rm(counts)
  gc(verbose = FALSE)
}

# 合并统计结果
results_df <- do.call(
  rbind,
  results
)

rownames(results_df) <- NULL

cat("\n========== 六只小鼠原始数据汇总 ==========\n")

print(results_df)

# 保存结果
write.csv(
  results_df,
  file.path(
    out_dir,
    "02_Raw_Matrix_Summary.csv"
  ),
  row.names = FALSE
)

cat("\n第二步完成！\n")


# =====================================================
# STEP 03：建立6只小鼠的Seurat对象并进行QC评估
# 不过滤细胞
# =====================================================

library(Seurat)
library(ggplot2)
library(patchwork)

# -----------------------------------------------------
# 1. 路径
# -----------------------------------------------------

raw_dir <- "E:/新建文件夹/expressions"

out_dir <- paste0(
  "D:/R data/Age_vs_AgeTrp_scRNA/",
  "New_6Mice_3vs3_FromRaw"
)

qc_dir <- file.path(out_dir, "03_QC_BeforeFiltering")

dir.create(
  qc_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# -----------------------------------------------------
# 2. 分组
# -----------------------------------------------------

samples <- c(
  "Age1",
  "Age2",
  "Age_Trp1",
  "Age_Trp2",
  "Age_Trp3",
  "Age_Trp4"
)

groups <- c(
  "Age",
  "Age",
  "Age",
  "Age_Trp",
  "Age_Trp",
  "Age_Trp"
)

# -----------------------------------------------------
# 3. 逐样本建立Seurat对象
# -----------------------------------------------------

seurat_list <- list()

for (i in seq_along(samples)) {
  
  sample_name <- samples[i]
  
  cat("\n正在处理：", sample_name, "\n")
  
  sample_dir <- file.path(raw_dir, sample_name)
  
  counts <- Read10X(
    data.dir = sample_dir
  )
  
  if (is.list(counts)) {
    counts <- counts[["Gene Expression"]]
  }
  
  # 不设置最低基因数过滤
  seu <- CreateSeuratObject(
    counts = counts,
    project = sample_name,
    min.cells = 0,
    min.features = 0
  )
  
  seu$sample <- sample_name
  seu$Group <- groups[i]
  
  # 小鼠线粒体基因通常以mt-开头
  mt_genes <- grep(
    "^mt-",
    rownames(seu),
    value = TRUE,
    ignore.case = TRUE
  )
  
  cat("识别到线粒体基因：", length(mt_genes), "\n")
  
  if (length(mt_genes) == 0) {
    stop(
      sample_name,
      "：未识别到线粒体基因，",
      "需要检查原始基因名称。"
    )
  }
  
  seu[["percent.mt"]] <- PercentageFeatureSet(
    seu,
    features = mt_genes
  )
  
  seurat_list[[sample_name]] <- seu
  
  rm(counts, seu)
  gc(verbose = FALSE)
}

# -----------------------------------------------------
# 4. 汇总所有细胞的QC信息
# -----------------------------------------------------

qc_data <- do.call(
  rbind,
  lapply(
    names(seurat_list),
    function(s) {
      
      x <- seurat_list[[s]]
      
      data.frame(
        Sample = s,
        Group = x$Group,
        nFeature_RNA = x$nFeature_RNA,
        nCount_RNA = x$nCount_RNA,
        percent.mt = x$percent.mt
      )
    }
  )
)

qc_data$Sample <- factor(
  qc_data$Sample,
  levels = samples
)

# -----------------------------------------------------
# 5. 汇总每只鼠QC统计
# -----------------------------------------------------

qc_summary <- do.call(
  rbind,
  lapply(
    samples,
    function(s) {
      
      x <- qc_data[
        qc_data$Sample == s,
      ]
      
      data.frame(
        Sample = s,
        Group = unique(x$Group),
        Cells = nrow(x),
        
        Median_nFeature = median(x$nFeature_RNA),
        Mean_nFeature = mean(x$nFeature_RNA),
        
        Median_nCount = median(x$nCount_RNA),
        Mean_nCount = mean(x$nCount_RNA),
        
        Median_MT = median(x$percent.mt),
        Mean_MT = mean(x$percent.mt),
        
        MT_Above10 = sum(x$percent.mt > 10),
        MT_Above15 = sum(x$percent.mt > 15),
        MT_Above20 = sum(x$percent.mt > 20)
      )
    }
  )
)

rownames(qc_summary) <- NULL

cat("\n========== QC汇总 ==========\n")

print(qc_summary)

write.csv(
  qc_summary,
  file.path(qc_dir, "03_QC_Summary.csv"),
  row.names = FALSE
)

write.csv(
  qc_data,
  file.path(qc_dir, "03_QC_AllCells.csv"),
  row.names = FALSE
)

# -----------------------------------------------------
# 6. 小提琴图：每个样本分别显示
# -----------------------------------------------------

p1 <- ggplot(
  qc_data,
  aes(x = Sample, y = nFeature_RNA, fill = Group)
) +
  geom_violin(trim = TRUE) +
  geom_boxplot(
    width = 0.12,
    fill = "white",
    outlier.shape = NA
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1)
  ) +
  labs(title = "Detected genes per cell", x = NULL)

p2 <- ggplot(
  qc_data,
  aes(x = Sample, y = nCount_RNA, fill = Group)
) +
  geom_violin(trim = TRUE) +
  geom_boxplot(
    width = 0.12,
    fill = "white",
    outlier.shape = NA
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1)
  ) +
  labs(title = "UMI counts per cell", x = NULL)

p3 <- ggplot(
  qc_data,
  aes(x = Sample, y = percent.mt, fill = Group)
) +
  geom_violin(trim = TRUE) +
  geom_boxplot(
    width = 0.12,
    fill = "white",
    outlier.shape = NA
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1)
  ) +
  labs(title = "Mitochondrial percentage", x = NULL)

p_qc <- (p1 | p2 | p3) +
  plot_layout(guides = "collect")

ggsave(
  file.path(qc_dir, "03_QC_Violin.png"),
  p_qc,
  width = 18,
  height = 5.5,
  dpi = 500,
  bg = "white"
)

ggsave(
  file.path(qc_dir, "03_QC_Violin.pdf"),
  p_qc,
  width = 18,
  height = 5.5
)

# -----------------------------------------------------
# 7. 散点图：UMI数、基因数和线粒体比例
# -----------------------------------------------------

p_scatter1 <- ggplot(
  qc_data,
  aes(x = nCount_RNA, y = nFeature_RNA)
) +
  geom_point(alpha = 0.15, size = 0.4) +
  facet_wrap(~ Sample, scales = "free") +
  theme_classic() +
  labs(
    title = "UMI counts vs detected genes"
  )

p_scatter2 <- ggplot(
  qc_data,
  aes(x = nCount_RNA, y = percent.mt)
) +
  geom_point(alpha = 0.15, size = 0.4) +
  facet_wrap(~ Sample, scales = "free") +
  theme_classic() +
  labs(
    title = "UMI counts vs mitochondrial percentage"
  )

ggsave(
  file.path(qc_dir, "03_QC_Scatter_Features.png"),
  p_scatter1,
  width = 12,
  height = 8,
  dpi = 500,
  bg = "white"
)

ggsave(
  file.path(qc_dir, "03_QC_Scatter_Mito.png"),
  p_scatter2,
  width = 12,
  height = 8,
  dpi = 500,
  bg = "white"
)

# -----------------------------------------------------
# 8. 保存尚未过滤的Seurat对象
# -----------------------------------------------------

saveRDS(
  seurat_list,
  file.path(out_dir, "03_Seurat_BeforeQC.rds")
)

cat("\n第三步完成：已生成QC结果，尚未过滤细胞。\n")


# =====================================================
# STEP 04：六只小鼠QC过滤
# 不做双细胞检测、不合并对象
# =====================================================

library(Seurat)
library(ggplot2)

# 1. 路径
out_dir <- paste0(
  "D:/R data/Age_vs_AgeTrp_scRNA/",
  "New_6Mice_3vs3_FromRaw"
)

qc_dir <- file.path(out_dir, "04_QC_Filtering")

dir.create(
  qc_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# 2. 读取第三步保存的对象
seurat_list <- readRDS(
  file.path(out_dir, "03_Seurat_BeforeQC.rds")
)

samples <- c(
  "Age1", "Age2", "Age_Trp1",
  "Age_Trp2", "Age_Trp3", "Age_Trp4"
)

groups <- c(
  "Age", "Age", "Age",
  "Age_Trp", "Age_Trp", "Age_Trp"
)

# 3. 设置过滤阈值
min_features <- 300
max_features <- 8000
min_counts <- 500
max_mt <- 10

filtered_list <- list()
qc_results <- list()

# 4. 逐样本过滤
for (i in seq_along(samples)) {
  
  s <- samples[i]
  seu <- seurat_list[[s]]
  
  md <- seu@meta.data
  
  keep <- (
    md$nFeature_RNA >= min_features &
      md$nFeature_RNA <= max_features &
      md$nCount_RNA >= min_counts &
      md$percent.mt < max_mt
  )
  
  keep[is.na(keep)] <- FALSE
  
  kept_cells <- rownames(md)[keep]
  
  if (length(kept_cells) == 0) {
    stop(s, "：过滤后没有剩余细胞。")
  }
  
  filtered_list[[s]] <- subset(
    seu,
    cells = kept_cells
  )
  
  qc_results[[s]] <- data.frame(
    Sample = s,
    Group = groups[i],
    Before = ncol(seu),
    After = length(kept_cells),
    Removed = ncol(seu) - length(kept_cells),
    RetentionPct = round(
      length(kept_cells) / ncol(seu) * 100,
      2
    ),
    LowFeatures = sum(md$nFeature_RNA < min_features),
    HighFeatures = sum(md$nFeature_RNA > max_features),
    LowCounts = sum(md$nCount_RNA < min_counts),
    HighMT = sum(md$percent.mt >= max_mt)
  )
  
  cat(
    s,
    ": ",
    ncol(seu),
    " -> ",
    length(kept_cells),
    " cells\n",
    sep = ""
  )
}

# 5. 汇总过滤结果
qc_table <- do.call(
  rbind,
  qc_results
)

rownames(qc_table) <- NULL

cat("\n========== QC过滤结果 ==========\n")

print(qc_table)

# 注意：各项不合格数量可能相互重叠，
# 因此不能直接相加作为Removed。

write.csv(
  qc_table,
  file.path(qc_dir, "04_QC_Filtering_Summary.csv"),
  row.names = FALSE
)

# 6. 绘制过滤前后细胞数量图
plot_data <- rbind(
  data.frame(
    Sample = qc_table$Sample,
    Group = qc_table$Group,
    Stage = "Before",
    Cells = qc_table$Before
  ),
  data.frame(
    Sample = qc_table$Sample,
    Group = qc_table$Group,
    Stage = "After",
    Cells = qc_table$After
  )
)

plot_data$Sample <- factor(
  plot_data$Sample,
  levels = samples
)

p <- ggplot(
  plot_data,
  aes(
    x = Sample,
    y = Cells,
    fill = Stage
  )
) +
  geom_col(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  geom_text(
    aes(label = Cells),
    position = position_dodge(width = 0.8),
    vjust = -0.3,
    size = 3
  ) +
  scale_fill_manual(
    values = c(
      Before = "#B7B7B7",
      After = "#3674AB"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.12))
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    )
  ) +
  labs(
    title = "Cell retention after QC",
    x = NULL,
    y = "Number of cells"
  )

ggsave(
  file.path(qc_dir, "04_Cell_Retention.png"),
  p,
  width = 9,
  height = 6,
  dpi = 500,
  bg = "white"
)

ggsave(
  file.path(qc_dir, "04_Cell_Retention.pdf"),
  p,
  width = 9,
  height = 6
)

# 7. 保存过滤后的对象
saveRDS(
  filtered_list,
  file.path(out_dir, "04_Seurat_AfterQC.rds")
)

cat("\n第四步完成！\n")
cat("过滤后总细胞数：",
    sum(qc_table$After), "\n")



# =====================================================
# STEP 05：逐样本检测并去除Doublets
# Age n=3 vs Age_Trp n=3
# =====================================================

# 1. 安装与加载
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("scDblFinder", quietly = TRUE)) {
  BiocManager::install(
    "scDblFinder",
    ask = FALSE,
    update = FALSE
  )
}

library(Seurat)
library(SingleCellExperiment)
library(scDblFinder)
library(Matrix)
library(ggplot2)

set.seed(123)

# 2. 路径
out_dir <- paste0(
  "D:/R data/Age_vs_AgeTrp_scRNA/",
  "New_6Mice_3vs3_FromRaw"
)

doublet_dir <- file.path(
  out_dir,
  "05_Doublet_Detection"
)

dir.create(
  doublet_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# 3. 读取上一阶段QC后的对象
seurat_list <- readRDS(
  file.path(out_dir, "04_Seurat_AfterQC.rds")
)

samples <- c(
  "Age1",
  "Age2",
  "Age_Trp1",
  "Age_Trp2",
  "Age_Trp3",
  "Age_Trp4"
)

groups <- c(
  "Age", "Age", "Age",
  "Age_Trp", "Age_Trp", "Age_Trp"
)

# 4. 初始化
singlet_list <- list()
doublet_results <- list()
cell_classifications <- list()

# 5. 每只鼠分别检测
for (i in seq_along(samples)) {
  
  s <- samples[i]
  
  cat("\n=================================\n")
  cat("正在检测：", s, "\n")
  cat("=================================\n")
  
  seu <- seurat_list[[s]]
  
  DefaultAssay(seu) <- "RNA"
  
  if (inherits(seu[["RNA"]], "Assay5")) {
    seu <- JoinLayers(seu, assay = "RNA")
  }
  
  counts <- GetAssayData(
    seu,
    assay = "RNA",
    layer = "counts"
  )
  
  # 转换为SingleCellExperiment
  sce <- SingleCellExperiment(
    assays = list(counts = counts)
  )
  
  # 检测双细胞
  sce <- scDblFinder(
    sce,
    samples = NULL,
    clusters = FALSE,
    BPPARAM = BiocParallel::SerialParam(),
    verbose = FALSE
  )
  
  # 取出每个细胞的评分和分类
  df <- data.frame(
    Cell = colnames(sce),
    Sample = s,
    Group = groups[i],
    DoubletScore = sce$scDblFinder.score,
    Classification = as.character(
      sce$scDblFinder.class
    )
  )
  
  cell_classifications[[s]] <- df
  
  # 只保留Singlet
  singlet_cells <- df$Cell[
    df$Classification == "singlet"
  ]
  
  singlet_list[[s]] <- subset(
    seu,
    cells = singlet_cells
  )
  
  # 每只鼠统计
  n_before <- ncol(seu)
  n_singlet <- length(singlet_cells)
  n_doublet <- sum(df$Classification == "doublet")
  
  doublet_results[[s]] <- data.frame(
    Sample = s,
    Group = groups[i],
    Before = n_before,
    Singlets = n_singlet,
    Doublets = n_doublet,
    DoubletPct = round(
      n_doublet / n_before * 100,
      2
    ),
    RetentionPct = round(
      n_singlet / n_before * 100,
      2
    )
  )
  
  cat(
    "过滤前：", n_before,
    "  Singlets：", n_singlet,
    "  Doublets：", n_doublet,
    "\n"
  )
  
  rm(sce, counts, seu)
  gc(verbose = FALSE)
}

# 6. 汇总
doublet_table <- do.call(
  rbind,
  doublet_results
)

rownames(doublet_table) <- NULL

cat("\n========== Doublet检测汇总 ==========\n")

print(doublet_table)

write.csv(
  doublet_table,
  file.path(
    doublet_dir,
    "05_Doublet_Summary.csv"
  ),
  row.names = FALSE
)

all_classifications <- do.call(
  rbind,
  cell_classifications
)

write.csv(
  all_classifications,
  file.path(
    doublet_dir,
    "05_Cell_Doublet_Scores.csv"
  ),
  row.names = FALSE
)

# 7. 双细胞比例柱状图
doublet_table$Sample <- factor(
  doublet_table$Sample,
  levels = samples
)

p <- ggplot(
  doublet_table,
  aes(
    x = Sample,
    y = DoubletPct,
    fill = Group
  )
) +
  geom_col(width = 0.65) +
  geom_text(
    aes(label = paste0(DoubletPct, "%")),
    vjust = -0.4,
    size = 4
  ) +
  scale_fill_manual(
    values = c(
      Age = "#E66A61",
      Age_Trp = "#16A7B1"
    )
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.15))
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    ),
    legend.position = "none"
  ) +
  labs(
    title = "Doublet detection across six samples",
    x = NULL,
    y = "Predicted doublets (%)"
  )

ggsave(
  file.path(doublet_dir, "05_Doublet_Rate.png"),
  p,
  width = 8,
  height = 5.5,
  dpi = 500,
  bg = "white"
)

ggsave(
  file.path(doublet_dir, "05_Doublet_Rate.pdf"),
  p,
  width = 8,
  height = 5.5
)

# 8. 保存去除Doublet的对象
saveRDS(
  singlet_list,
  file.path(
    out_dir,
    "05_Seurat_Singlets.rds"
  )
)

cat("\n第五步完成！\n")
cat(
  "六只小鼠剩余Singlets：",
  sum(doublet_table$Singlets),
  "\n"
)


# =====================================================
# STEP 06：合并、标准化、高变基因、PCA
# 不进行Harmony或聚类
# =====================================================

library(Seurat)
library(ggplot2)
library(patchwork)

set.seed(123)

# 1. 路径
out_dir <- paste0(
  "D:/R data/Age_vs_AgeTrp_scRNA/",
  "New_6Mice_3vs3_FromRaw"
)

pca_dir <- file.path(out_dir, "06_PCA")

dir.create(
  pca_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# 2. 读取去双细胞后的对象
seurat_list <- readRDS(
  file.path(out_dir, "05_Seurat_Singlets.rds")
)

samples <- c(
  "Age1",
  "Age2",
  "Age_Trp1",
  "Age_Trp2",
  "Age_Trp3",
  "Age_Trp4"
)

# 确认样本完整
stopifnot(all(samples %in% names(seurat_list)))

seurat_list <- seurat_list[samples]

# 3. 合并对象
combined <- merge(
  x = seurat_list[[1]],
  y = seurat_list[-1],
  add.cell.ids = samples,
  project = "Ovary_6Mice_3vs3"
)

DefaultAssay(combined) <- "RNA"

# Seurat v5合并RNA层
if (inherits(combined[["RNA"]], "Assay5")) {
  combined <- JoinLayers(
    combined,
    assay = "RNA"
  )
}

# 4. 明确实验分组
age_samples <- c(
  "Age1",
  "Age2",
  "Age_Trp1"
)

trp_samples <- c(
  "Age_Trp2",
  "Age_Trp3",
  "Age_Trp4"
)

combined$sample <- as.character(combined$sample)

combined$Group6 <- factor(
  ifelse(
    combined$sample %in% age_samples,
    "Age",
    "Age_Trp"
  ),
  levels = c("Age", "Age_Trp")
)

cat("\n===== 六只小鼠分组 =====\n")
print(table(combined$sample, combined$Group6))

cat("\n合并细胞总数：", ncol(combined), "\n")

stopifnot(
  ncol(combined) == 37304,
  all(combined$sample %in% samples)
)

# 5. RNA标准化
combined <- NormalizeData(
  combined,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = FALSE
)

# 6. 寻找高变基因
combined <- FindVariableFeatures(
  combined,
  selection.method = "vst",
  nfeatures = 3000,
  verbose = FALSE
)

cat(
  "高变基因数：",
  length(VariableFeatures(combined)),
  "\n"
)

# 7. ScaleData
combined <- ScaleData(
  combined,
  features = VariableFeatures(combined),
  verbose = FALSE
)

# 8. PCA
combined <- RunPCA(
  combined,
  features = VariableFeatures(combined),
  npcs = 40,
  seed.use = 123,
  verbose = FALSE
)

# 9. PCA：按样本着色
p_sample <- DimPlot(
  combined,
  reduction = "pca",
  group.by = "sample",
  raster = TRUE
) +
  ggtitle("PCA - Six samples")

# 10. PCA：按组别着色
p_group <- DimPlot(
  combined,
  reduction = "pca",
  group.by = "Group6",
  cols = c(
    "Age" = "#E66A61",
    "Age_Trp" = "#16A7B1"
  ),
  raster = TRUE
) +
  ggtitle("PCA - Age vs Age_Trp")

p_combined <- p_sample | p_group

ggsave(
  file.path(pca_dir, "06_PCA_Sample_Group.png"),
  p_combined,
  width = 15,
  height = 6,
  dpi = 400,
  bg = "white"
)

ggsave(
  file.path(pca_dir, "06_PCA_Sample_Group.pdf"),
  p_combined,
  width = 15,
  height = 6
)

# 11. ElbowPlot
p_elbow <- ElbowPlot(
  combined,
  ndims = 40
) +
  ggtitle("PCA Elbow Plot")

ggsave(
  file.path(pca_dir, "06_PCA_Elbow.png"),
  p_elbow,
  width = 8,
  height = 5.5,
  dpi = 400,
  bg = "white"
)

# 12. 前30个PC解释方差
pca_stdev <- Stdev(combined, reduction = "pca")

variance_pct <- 100 * pca_stdev^2 / sum(pca_stdev^2)

pca_summary <- data.frame(
  PC = seq_along(variance_pct),
  VariancePct = variance_pct,
  CumulativePct = cumsum(variance_pct)
)

cat("\n===== PCA方差贡献（前30个） =====\n")
print(head(pca_summary, 30))

write.csv(
  pca_summary,
  file.path(pca_dir, "06_PCA_Variance.csv"),
  row.names = FALSE
)

# 13. 保存PCA对象
saveRDS(
  combined,
  file.path(out_dir, "06_Seurat_PCA.rds")
)

cat("\n第六步完成！\n")
cat("细胞总数：", ncol(combined), "\n")
cat("结果目录：", pca_dir, "\n")



# =====================================================
# STEP 07: Harmony integration
# Six mice, Age n=3 vs Age_Trp n=3
# =====================================================

library(Seurat)
library(harmony)
library(ggplot2)
library(patchwork)

set.seed(123)

# 1. 路径
out_dir <- paste0(
  "D:/R data/Age_vs_AgeTrp_scRNA/",
  "New_6Mice_3vs3_FromRaw"
)

harmony_dir <- file.path(out_dir, "07_Harmony")

dir.create(
  harmony_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# 2. 读取PCA对象
combined <- readRDS(
  file.path(out_dir, "06_Seurat_PCA.rds")
)

# 3. 再次检查分组
cat("\n===== 样本分组 =====\n")

print(
  table(combined$sample, combined$Group6)
)

# 4. Harmony整合
# 只使用前30个PC
# 仅将sample作为校正变量，不使用Group6

combined <- harmony::RunHarmony(
  object = combined,
  group.by.vars = "sample",
  reduction.use = "pca",
  dims.use = 1:30,
  reduction.save = "harmony",
  project.dim = FALSE,
  verbose = TRUE
)

# 5. 检查Harmony结果
cat("\n===== Harmony结果 =====\n")

print(Reductions(combined))

cat(
  "Harmony维度：",
  ncol(Embeddings(combined, "harmony")),
  "\n"
)

# 6. PCA与Harmony按样本着色
p_before <- DimPlot(
  combined,
  reduction = "pca",
  group.by = "sample",
  raster = TRUE
) +
  ggtitle("Before Harmony")

p_after <- DimPlot(
  combined,
  reduction = "harmony",
  group.by = "sample",
  raster = TRUE
) +
  ggtitle("After Harmony")

p_compare <- p_before | p_after

ggsave(
  file.path(
    harmony_dir,
    "07_Before_After_Harmony_Samples.png"
  ),
  p_compare,
  width = 15,
  height = 6,
  dpi = 400,
  bg = "white"
)

ggsave(
  file.path(
    harmony_dir,
    "07_Before_After_Harmony_Samples.pdf"
  ),
  p_compare,
  width = 15,
  height = 6
)

# 7. Harmony按组别着色
p_group <- DimPlot(
  combined,
  reduction = "harmony",
  group.by = "Group6",
  cols = c(
    Age = "#E66A61",
    Age_Trp = "#16A7B1"
  ),
  raster = TRUE
) +
  ggtitle("Harmony - Age vs Age_Trp")

ggsave(
  file.path(
    harmony_dir,
    "07_Harmony_ByGroup.png"
  ),
  p_group,
  width = 8,
  height = 6,
  dpi = 400,
  bg = "white"
)

# 8. Harmony按样本着色
p_sample <- DimPlot(
  combined,
  reduction = "harmony",
  group.by = "sample",
  raster = TRUE
) +
  ggtitle("Harmony - Six Samples")

ggsave(
  file.path(
    harmony_dir,
    "07_Harmony_BySample.png"
  ),
  p_sample,
  width = 9,
  height = 6,
  dpi = 400,
  bg = "white"
)

# 9. 保存对象
saveRDS(
  combined,
  file.path(
    out_dir,
    "07_Seurat_Harmony.rds"
  )
)

cat("\n第七步Harmony整合完成！\n")
cat("细胞数量：", ncol(combined), "\n")




# =====================================================
# STEP 08
# UMAP + t-SNE + Clustering
#
# Age: Age1, Age2, Age_Trp1
# Age_Trp: Age_Trp2, Age_Trp3, Age_Trp4
# =====================================================

library(Seurat)
library(ggplot2)
library(patchwork)

set.seed(123)

# 1. 文件路径
out_dir <- paste0(
  "D:/R data/Age_vs_AgeTrp_scRNA/",
  "New_6Mice_3vs3_FromRaw"
)

step_dir <- file.path(out_dir, "08_Clustering")

dir.create(
  step_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# 2. 读取Harmony对象
combined <- readRDS(
  file.path(out_dir, "07_Seurat_Harmony.rds")
)

# 3. 检查Harmony
stopifnot(
  "harmony" %in% Reductions(combined),
  ncol(Embeddings(combined, "harmony")) >= 30
)

cat("细胞总数：", ncol(combined), "\n")

print(table(
  combined$sample,
  combined$Group6
))

# 4. 构建邻居图
combined <- FindNeighbors(
  combined,
  reduction = "harmony",
  dims = 1:30,
  verbose = FALSE
)

# 5. 聚类
combined <- FindClusters(
  combined,
  resolution = 0.5,
  random.seed = 123,
  verbose = FALSE
)

# 6. UMAP
combined <- RunUMAP(
  combined,
  reduction = "harmony",
  dims = 1:30,
  n.neighbors = 30,
  min.dist = 0.3,
  seed.use = 123,
  verbose = FALSE
)

# 7. t-SNE
combined <- RunTSNE(
  combined,
  reduction = "harmony",
  dims = 1:30,
  seed.use = 123,
  check_duplicates = FALSE,
  verbose = FALSE
)

# 8. 绘图函数
save_fig <- function(
    plot,
    name,
    width = 9,
    height = 7
) {
  
  ggsave(
    file.path(step_dir, paste0(name, ".png")),
    plot,
    width = width,
    height = height,
    dpi = 500,
    bg = "white"
  )
  
  ggsave(
    file.path(step_dir, paste0(name, ".pdf")),
    plot,
    width = width,
    height = height,
    bg = "white"
  )
}

# =====================================================
# 9. UMAP - 聚类编号
# =====================================================

p_cluster <- DimPlot(
  combined,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  raster = FALSE
) +
  ggtitle("UMAP - Cell clusters")

save_fig(
  p_cluster,
  "08_UMAP_Clusters"
)

# =====================================================
# 10. UMAP - 组别
# =====================================================

p_group <- DimPlot(
  combined,
  reduction = "umap",
  group.by = "Group6",
  cols = c(
    "Age" = "#E66A61",
    "Age_Trp" = "#16A7B1"
  ),
  raster = FALSE
) +
  ggtitle("UMAP - Age vs Age_Trp")

save_fig(
  p_group,
  "08_UMAP_Groups"
)

# =====================================================
# 11. UMAP - 六只小鼠
# =====================================================

p_sample <- DimPlot(
  combined,
  reduction = "umap",
  group.by = "sample",
  raster = FALSE
) +
  ggtitle("UMAP - Six samples")

save_fig(
  p_sample,
  "08_UMAP_Samples",
  width = 10,
  height = 7
)

# =====================================================
# 12. UMAP - 按组分别显示
# =====================================================

p_split <- DimPlot(
  combined,
  reduction = "umap",
  group.by = "seurat_clusters",
  split.by = "Group6",
  raster = FALSE
)

save_fig(
  p_split,
  "08_UMAP_SplitGroups",
  width = 15,
  height = 7
)

# =====================================================
# 13. t-SNE - 聚类
# =====================================================

p_tsne <- DimPlot(
  combined,
  reduction = "tsne",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  raster = FALSE
) +
  ggtitle("t-SNE - Cell clusters")

save_fig(
  p_tsne,
  "08_TSNE_Clusters"
)

# =====================================================
# 14. t-SNE - 组别
# =====================================================

p_tsne_group <- DimPlot(
  combined,
  reduction = "tsne",
  group.by = "Group6",
  cols = c(
    "Age" = "#E66A61",
    "Age_Trp" = "#16A7B1"
  ),
  raster = FALSE
)

save_fig(
  p_tsne_group,
  "08_TSNE_Groups"
)

# =====================================================
# 15. 各Cluster细胞数
# =====================================================

cluster_table <- as.data.frame(
  table(
    Sample = combined$sample,
    Group = combined$Group6,
    Cluster = combined$seurat_clusters
  )
)

colnames(cluster_table)[4] <- "Cells"

write.csv(
  cluster_table,
  file.path(
    step_dir,
    "08_Cluster_Counts_BySample.csv"
  ),
  row.names = FALSE
)

cat("\n========== Cluster统计 ==========\n")

print(table(combined$seurat_clusters))

cat(
  "\nCluster总数：",
  length(unique(combined$seurat_clusters)),
  "\n"
)

# 16. 保存对象
saveRDS(
  combined,
  file.path(
    out_dir,
    "08_Seurat_Clustered.rds"
  )
)

cat("\n第八步完成！\n")
cat("保存目录：", step_dir, "\n")


# ============================================================
# STEP 09
# 六只小鼠卵巢单细胞：最终细胞大类注释
#
# Age: Age1, Age2, Age_Trp1
# Age_Trp: Age_Trp2, Age_Trp3, Age_Trp4
#
# 删除Cluster: 10, 15, 23
# Cluster 12: Cycling cells
#
# 不重新降维，不重新聚类
# 所有图片保存为 PDF + SVG
# ============================================================

library(Seurat)
library(ggplot2)
library(patchwork)

if (!requireNamespace("svglite", quietly = TRUE)) {
  install.packages("svglite")
}

set.seed(123)

# ============================================================
# 1. 路径
# ============================================================

out_dir <- paste0(
  "D:/R data/Age_vs_AgeTrp_scRNA/",
  "New_6Mice_3vs3_FromRaw"
)

step_dir <- file.path(
  out_dir,
  "09_Final_CellType_Annotation"
)

dir.create(
  step_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# 2. 读取第八步的聚类对象
# ============================================================

combined <- readRDS(
  file.path(out_dir, "08_Seurat_Clustered.rds")
)

cat("原始细胞数：", ncol(combined), "\n")

cat("\n原始Cluster数量：\n")
print(table(combined$seurat_clusters))

# ============================================================
# 3. 删除Cluster 10、15、23
# ============================================================

Idents(combined) <- "seurat_clusters"

remove_clusters <- c("10", "15", "23")

remove_cells <- WhichCells(
  combined,
  idents = remove_clusters
)

cat(
  "\n准备删除的细胞数：",
  length(remove_cells),
  "\n"
)

# 删除指定细胞
combined <- subset(
  combined,
  idents = remove_clusters,
  invert = TRUE
)

# 去掉未使用的Cluster因子水平
combined$seurat_clusters <- droplevels(
  factor(as.character(combined$seurat_clusters))
)

cat("\n删除后细胞数：", ncol(combined), "\n")

cat("\n保留的Cluster：\n")
print(table(combined$seurat_clusters))

# ============================================================
# 4. 细胞大类注释
# ============================================================

cluster_annotation <- c(
  
  "0" = "Stroma",
  "1" = "Granulosa",
  "2" = "T cells",
  "3" = "Stroma",
  "4" = "NK cells",
  "5" = "Steroidogenic cells",
  "6" = "Endothelium",
  "7" = "Granulosa",
  "8" = "Epithelium",
  "9" = "Phagocytes",
  
  "11" = "Stroma",
  "12" = "Cycling cells",
  "13" = "Phagocytes",
  "14" = "T cells",
  "16" = "Stroma",
  "17" = "Endothelium",
  "18" = "Steroidogenic cells",
  "19" = "Epithelium",
  "20" = "Pericytes",
  "21" = "B cells",
  "22" = "Neutrophils"
)

# 将Cluster编号映射为细胞类型
cluster_ids <- as.character(
  combined$seurat_clusters
)

combined$Major_CellType <- unname(
  cluster_annotation[cluster_ids]
)

if (anyNA(combined$Major_CellType)) {
  stop("存在未注释Cluster，请检查映射表。")
}

# 指定细胞类型显示顺序
celltype_order <- c(
  "Phagocytes",
  "Granulosa",
  "Stroma",
  "Steroidogenic cells",
  "Endothelium",
  "Epithelium",
  "T cells",
  "B cells",
  "NK cells",
  "Pericytes",
  "Neutrophils",
  "Cycling cells"
)

combined$Major_CellType <- factor(
  combined$Major_CellType,
  levels = celltype_order
)

# ============================================================
# 5. 配色
# ============================================================

celltype_colors <- c(
  "Phagocytes" = "#E64B35",
  "Granulosa" = "#D89000",
  "Stroma" = "#4DBBD5",
  "Steroidogenic cells" = "#00A087",
  "Endothelium" = "#3C5488",
  "Epithelium" = "#7E62A3",
  "T cells" = "#8491B4",
  "B cells" = "#B09C85",
  "NK cells" = "#91D1C2",
  "Pericytes" = "#F39B7F",
  "Neutrophils" = "#C77CFF",
  "Cycling cells" = "#FF69B4"
)

group_colors <- c(
  "Age" = "#E66A61",
  "Age_Trp" = "#16A7B1"
)

samples <- c(
  "Age1",
  "Age2",
  "Age_Trp1",
  "Age_Trp2",
  "Age_Trp3",
  "Age_Trp4"
)

# ============================================================
# 6. 确认分组
# ============================================================

age_samples <- c(
  "Age1",
  "Age2",
  "Age_Trp1"
)

combined$sample <- as.character(combined$sample)

stopifnot(
  all(combined$sample %in% samples)
)

combined$Group6 <- factor(
  ifelse(
    combined$sample %in% age_samples,
    "Age",
    "Age_Trp"
  ),
  levels = c("Age", "Age_Trp")
)

cat("\n===== 六只小鼠分组 =====\n")

print(
  table(
    combined$sample,
    combined$Group6
  )
)

# ============================================================
# 7. 统一图片保存函数：PDF + SVG
# ============================================================

save_fig <- function(
    p,
    name,
    width = 10,
    height = 7
) {
  
  ggsave(
    filename = file.path(
      step_dir,
      paste0(name, ".pdf")
    ),
    plot = p,
    device = grDevices::cairo_pdf,
    width = width,
    height = height,
    bg = "white"
  )
  
  ggsave(
    filename = file.path(
      step_dir,
      paste0(name, ".svg")
    ),
    plot = p,
    device = svglite::svglite,
    width = width,
    height = height,
    bg = "white"
  )
}

# ============================================================
# 8. UMAP：12个细胞大类
# ============================================================

p_celltype <- DimPlot(
  combined,
  reduction = "umap",
  group.by = "Major_CellType",
  cols = celltype_colors,
  label = TRUE,
  repel = TRUE,
  raster = FALSE
) +
  ggtitle("Major ovarian cell types") +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    )
  )

save_fig(
  p_celltype,
  "09_UMAP_Major_CellTypes",
  width = 12,
  height = 8
)

# ============================================================
# 9. UMAP：Age和Age_Trp分组
# ============================================================

p_group <- DimPlot(
  combined,
  reduction = "umap",
  group.by = "Group6",
  cols = group_colors,
  raster = FALSE
) +
  ggtitle("Age vs Age_Trp") +
  theme_classic(base_size = 13)

save_fig(
  p_group,
  "09_UMAP_Groups",
  width = 10,
  height = 7
)

# ============================================================
# 10. UMAP：两组分别显示细胞类型
# ============================================================

p_split <- DimPlot(
  combined,
  reduction = "umap",
  group.by = "Major_CellType",
  split.by = "Group6",
  cols = celltype_colors,
  raster = FALSE
) +
  ggtitle("Major cell types by group")

save_fig(
  p_split,
  "09_UMAP_CellTypes_SplitGroups",
  width = 17,
  height = 8
)

# ============================================================
# 11. UMAP：原Cluster编号
# ============================================================

p_cluster <- DimPlot(
  combined,
  reduction = "umap",
  group.by = "seurat_clusters",
  label = TRUE,
  repel = TRUE,
  raster = FALSE
) +
  ggtitle("Retained clusters")

save_fig(
  p_cluster,
  "09_UMAP_Retained_Clusters",
  width = 11,
  height = 8
)

# ============================================================
# 12. t-SNE：细胞大类
# ============================================================

p_tsne <- DimPlot(
  combined,
  reduction = "tsne",
  group.by = "Major_CellType",
  cols = celltype_colors,
  label = TRUE,
  repel = TRUE,
  raster = FALSE
) +
  ggtitle("Major ovarian cell types - t-SNE")

save_fig(
  p_tsne,
  "09_TSNE_Major_CellTypes",
  width = 12,
  height = 8
)

# ============================================================
# 13. 导出Cluster与细胞类型对应关系
# ============================================================

annotation_table <- data.frame(
  Cluster = names(cluster_annotation),
  CellType = unname(cluster_annotation)
)

write.csv(
  annotation_table,
  file.path(
    step_dir,
    "09_Cluster_Annotation.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 14. 导出每只小鼠的细胞数量
# ============================================================

cell_count <- as.data.frame(
  table(
    Sample = factor(
      combined$sample,
      levels = samples
    ),
    CellType = factor(
      combined$Major_CellType,
      levels = celltype_order
    )
  )
)

colnames(cell_count)[3] <- "Cells"

write.csv(
  cell_count,
  file.path(
    step_dir,
    "09_CellCounts_BySample_CellType.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 15. 查看细胞大类数量
# ============================================================

cat("\n===== 细胞大类数量 =====\n")

print(
  table(combined$Major_CellType)
)

cat("\n===== 保留的Cluster数量 =====\n")

print(
  length(unique(
    as.character(combined$seurat_clusters)
  ))
)

# ============================================================
# 16. 保存最终注释对象
# ============================================================

Idents(combined) <- "Major_CellType"

saveRDS(
  combined,
  file.path(
    out_dir,
    "09_Seurat_Final_Annotated.rds"
  )
)

cat("\n===================================\n")
cat("第九步完成！\n")
cat("删除Cluster：10、15、23\n")
cat("Cluster 12：Cycling cells\n")
cat("剩余细胞数：", ncol(combined), "\n")
cat("细胞大类数量：",
    length(unique(combined$Major_CellType)), "\n")
cat("所有图片均为PDF + SVG\n")
cat("===================================\n")

# =====================================================
# STEP 09 补充：t-SNE按组展示细胞大类
# 与09_UMAP_CellTypes_SplitGroups保持一致
# 输出PDF + SVG
# =====================================================

library(Seurat)
library(ggplot2)
library(svglite)

out_dir <- paste0(
  "D:/R data/Age_vs_AgeTrp_scRNA/",
  "New_6Mice_3vs3_FromRaw"
)

step_dir <- file.path(
  out_dir,
  "09_Final_CellType_Annotation"
)

dir.create(step_dir, recursive = TRUE, showWarnings = FALSE)

# 读取已注释对象，不重新计算
combined <- readRDS(
  file.path(out_dir, "09_Seurat_Final_Annotated.rds")
)

# 与UMAP使用相同的颜色
celltype_colors <- c(
  "Phagocytes" = "#E64B35",
  "Granulosa" = "#D89000",
  "Stroma" = "#4DBBD5",
  "Steroidogenic cells" = "#00A087",
  "Endothelium" = "#3C5488",
  "Epithelium" = "#7E62A3",
  "T cells" = "#8491B4",
  "B cells" = "#B09C85",
  "NK cells" = "#91D1C2",
  "Pericytes" = "#F39B7F",
  "Neutrophils" = "#C77CFF",
  "Cycling cells" = "#FF69B4"
)

# t-SNE：按组分面，细胞大类着色
p_tsne_split <- DimPlot(
  combined,
  reduction = "tsne",
  group.by = "Major_CellType",
  split.by = "Group6",
  cols = celltype_colors,
  raster = FALSE,
  pt.size = 0.25
) +
  labs(title = "Major ovarian cell types by group") +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    legend.title = element_blank(),
    strip.text = element_text(
      face = "bold",
      size = 13
    )
  )

# PDF
ggsave(
  file.path(
    step_dir,
    "09_TSNE_CellTypes_SplitGroups.pdf"
  ),
  plot = p_tsne_split,
  device = grDevices::cairo_pdf,
  width = 17,
  height = 8,
  bg = "white"
)

# SVG
ggsave(
  file.path(
    step_dir,
    "09_TSNE_CellTypes_SplitGroups.svg"
  ),
  plot = p_tsne_split,
  device = svglite::svglite,
  width = 17,
  height = 8,
  bg = "white"
)

print(p_tsne_split)

cat("t-SNE分组图已保存为PDF和SVG。\n")


# ============================================================
# STEP 10
# 六只小鼠细胞组成与占比分析
#
# Age: Age1, Age2, Age_Trp1
# Age_Trp: Age_Trp2, Age_Trp3, Age_Trp4
#
# 输入：09_Seurat_Final_Annotated.rds
#
# 输出：
# 1. 各样本细胞数量和比例
# 2. 六样本堆积柱状图
# 3. 组平均堆积柱状图
# 4. 细胞类型比例比较图
# 5. Granulosa比例比较图
# 6. 组间统计
#
# 图片：PDF + SVG
# ============================================================

library(Seurat)
library(ggplot2)
library(dplyr)
library(tidyr)

if (!requireNamespace("svglite", quietly = TRUE)) {
  install.packages("svglite")
}

library(svglite)

set.seed(123)

# ============================================================
# 1. 文件路径
# ============================================================

out_dir <- paste0(
  "D:/R data/Age_vs_AgeTrp_scRNA/",
  "New_6Mice_3vs3_FromRaw"
)

step_dir <- file.path(
  out_dir,
  "10_Cell_Proportion"
)

dir.create(
  step_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ============================================================
# 2. 读取最终注释对象
# ============================================================

combined <- readRDS(
  file.path(
    out_dir,
    "09_Seurat_Final_Annotated.rds"
  )
)

stopifnot(
  all(
    c("sample", "Major_CellType") %in%
      colnames(combined@meta.data)
  )
)

# ============================================================
# 3. 正确的3 vs 3分组
# ============================================================

age_samples <- c(
  "Age1",
  "Age2",
  "Age_Trp1"
)

trp_samples <- c(
  "Age_Trp2",
  "Age_Trp3",
  "Age_Trp4"
)

samples_all <- c(
  age_samples,
  trp_samples
)

group_colors <- c(
  "Age" = "#E66A61",
  "Age_Trp" = "#16A7B1"
)

# ============================================================
# 4. 细胞大类及配色
# 与第九步完全一致
# ============================================================

celltype_order <- c(
  "Phagocytes",
  "Granulosa",
  "Stroma",
  "Steroidogenic cells",
  "Endothelium",
  "Epithelium",
  "T cells",
  "B cells",
  "NK cells",
  "Pericytes",
  "Neutrophils",
  "Cycling cells"
)

celltype_colors <- c(
  "Phagocytes" = "#E64B35",
  "Granulosa" = "#D89000",
  "Stroma" = "#4DBBD5",
  "Steroidogenic cells" = "#00A087",
  "Endothelium" = "#3C5488",
  "Epithelium" = "#7E62A3",
  "T cells" = "#8491B4",
  "B cells" = "#B09C85",
  "NK cells" = "#91D1C2",
  "Pericytes" = "#F39B7F",
  "Neutrophils" = "#C77CFF",
  "Cycling cells" = "#FF69B4"
)

# ============================================================
# 5. 建立细胞元数据
# ============================================================

meta <- combined@meta.data %>%
  mutate(
    Sample = as.character(sample),
    CellType = as.character(Major_CellType)
  ) %>%
  filter(Sample %in% samples_all)

if (nrow(meta) != ncol(combined)) {
  stop("存在不属于指定6只鼠的细胞。")
}

if (any(!meta$CellType %in% celltype_order)) {
  stop("发现未定义的细胞类型，请检查第九步注释。")
}

meta <- meta %>%
  mutate(
    Group = ifelse(
      Sample %in% age_samples,
      "Age",
      "Age_Trp"
    )
  )

# ============================================================
# 6. 统计每只鼠的各细胞类型数量
# ============================================================

cell_count <- meta %>%
  count(
    Sample,
    Group,
    CellType,
    name = "Cells"
  ) %>%
  complete(
    Sample = samples_all,
    CellType = celltype_order,
    fill = list(Cells = 0)
  ) %>%
  mutate(
    Group = ifelse(
      Sample %in% age_samples,
      "Age",
      "Age_Trp"
    )
  ) %>%
  group_by(Sample) %>%
  mutate(
    TotalCells = sum(Cells),
    Percentage = Cells / TotalCells * 100
  ) %>%
  ungroup()

cell_count$Sample <- factor(
  cell_count$Sample,
  levels = samples_all
)

cell_count$Group <- factor(
  cell_count$Group,
  levels = c("Age", "Age_Trp")
)

cell_count$CellType <- factor(
  cell_count$CellType,
  levels = celltype_order
)

cat("\n===== 各样本细胞组成 =====\n")

print(
  cell_count %>%
    arrange(Sample, CellType)
)

write.csv(
  cell_count,
  file.path(
    step_dir,
    "10_CellCounts_And_Percentages.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 7. 输出各样本比例宽表
# ============================================================

percentage_wide <- cell_count %>%
  select(
    Sample,
    Group,
    CellType,
    Percentage
  ) %>%
  pivot_wider(
    names_from = CellType,
    values_from = Percentage
  )

write.csv(
  percentage_wide,
  file.path(
    step_dir,
    "10_CellPercentage_Wide.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 8. 计算组平均比例
#
# 重要：
# 先计算每只鼠比例，
# 然后对3只鼠求平均，
# 避免细胞数多的小鼠权重过高
# ============================================================

group_mean <- cell_count %>%
  group_by(Group, CellType) %>%
  summarise(
    MeanPercentage = mean(Percentage),
    SD = sd(Percentage),
    SEM = sd(Percentage) / sqrt(n()),
    .groups = "drop"
  )

write.csv(
  group_mean,
  file.path(
    step_dir,
    "10_GroupMean_CellPercentage.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 9. 统一图片保存函数
# ============================================================

save_fig <- function(
    p,
    filename,
    width = 10,
    height = 7
) {
  
  ggsave(
    file.path(
      step_dir,
      paste0(filename, ".pdf")
    ),
    plot = p,
    device = grDevices::cairo_pdf,
    width = width,
    height = height,
    bg = "white"
  )
  
  ggsave(
    file.path(
      step_dir,
      paste0(filename, ".svg")
    ),
    plot = p,
    device = svglite::svglite,
    width = width,
    height = height,
    bg = "white"
  )
}

# ============================================================
# 10. 六只小鼠的细胞比例堆积柱状图
# 每只鼠一根柱子
# ============================================================

p_sample <- ggplot(
  cell_count,
  aes(
    x = Sample,
    y = Percentage,
    fill = CellType
  )
) +
  geom_col(
    width = 0.78,
    color = "white",
    linewidth = 0.15
  ) +
  scale_fill_manual(
    values = celltype_colors,
    drop = FALSE
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    expand = c(0, 0)
  ) +
  labs(
    title = "Cell type composition across six mice",
    x = NULL,
    y = "Percentage of cells (%)",
    fill = "Cell type"
  ) +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    ),
    legend.position = "right"
  )

save_fig(
  p_sample,
  "10_StackedBar_6Mice",
  width = 12,
  height = 7
)

# ============================================================
# 11. 两组平均细胞比例堆积图
# 每组一根柱子
# ============================================================

p_group <- ggplot(
  group_mean,
  aes(
    x = Group,
    y = MeanPercentage,
    fill = CellType
  )
) +
  geom_col(
    width = 0.65,
    color = "white",
    linewidth = 0.2
  ) +
  scale_fill_manual(
    values = celltype_colors,
    drop = FALSE
  ) +
  scale_y_continuous(
    limits = c(0, 100),
    expand = c(0, 0)
  ) +
  theme_classic(base_size = 13) +
  labs(
    title = "Mean cell type composition",
    subtitle = "Age (n=3) vs Age_Trp (n=3)",
    x = NULL,
    y = "Mean cell percentage (%)",
    fill = "Cell type"
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      hjust = 0.5
    )
  )

save_fig(
  p_group,
  "10_StackedBar_GroupMean",
  width = 9,
  height = 7
)

# ============================================================
# 12. 各细胞大类的组间比例比较
# 每个点代表一只小鼠
# ============================================================

p_compare <- ggplot(
  cell_count,
  aes(
    x = Group,
    y = Percentage,
    color = Group
  )
) +
  geom_boxplot(
    width = 0.5,
    outlier.shape = NA
  ) +
  geom_point(
    size = 2.3,
    position = position_jitter(
      width = 0.08,
      height = 0
    )
  ) +
  facet_wrap(
    ~ CellType,
    scales = "free_y",
    ncol = 4
  ) +
  scale_color_manual(
    values = group_colors
  ) +
  theme_classic(base_size = 12) +
  labs(
    title = "Cell type proportion comparison",
    subtitle = "Each dot represents one mouse",
    x = NULL,
    y = "Cell percentage (%)"
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      hjust = 0.5
    ),
    legend.position = "none",
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    ),
    strip.text = element_text(
      face = "bold",
      size = 10
    )
  )

save_fig(
  p_compare,
  "10_CellType_Proportion_Comparison",
  width = 15,
  height = 11
)

# ============================================================
# 13. 重点：Granulosa细胞比例比较
# ============================================================

gc_prop <- cell_count %>%
  filter(CellType == "Granulosa")

p_gc <- ggplot(
  gc_prop,
  aes(
    x = Group,
    y = Percentage,
    color = Group
  )
) +
  geom_boxplot(
    width = 0.45,
    outlier.shape = NA,
    linewidth = 0.8
  ) +
  geom_point(
    size = 3.5,
    position = position_jitter(
      width = 0.07,
      height = 0
    )
  ) +
  ggrepel::geom_text_repel(
    aes(label = Sample),
    color = "black",
    size = 3.4,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = group_colors
  ) +
  theme_classic(base_size = 14) +
  labs(
    title = "Granulosa cell proportion",
    subtitle = "Age (n=3) vs Age_Trp (n=3)",
    x = NULL,
    y = "Granulosa cells (%)"
  ) +
  theme(
    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),
    plot.subtitle = element_text(
      hjust = 0.5
    ),
    legend.position = "none"
  )

save_fig(
  p_gc,
  "10_Granulosa_Proportion",
  width = 7,
  height = 6
)

# ============================================================
# 14. 组间统计
# 每只小鼠为独立重复
#
# n=3 vs n=3，使用Welch t检验作探索性比较
# 同时保留Wilcoxon P值
# ============================================================

stat_results <- lapply(
  celltype_order,
  function(ct) {
    
    x <- cell_count %>%
      filter(CellType == ct)
    
    age_values <- x$Percentage[
      x$Group == "Age"
    ]
    
    trp_values <- x$Percentage[
      x$Group == "Age_Trp"
    ]
    
    p_t <- tryCatch(
      t.test(
        trp_values,
        age_values
      )$p.value,
      error = function(e) NA_real_
    )
    
    p_w <- tryCatch(
      wilcox.test(
        trp_values,
        age_values,
        exact = FALSE
      )$p.value,
      error = function(e) NA_real_
    )
    
    data.frame(
      CellType = ct,
      AgeMean = mean(age_values),
      AgeTrpMean = mean(trp_values),
      Difference = mean(trp_values) - mean(age_values),
      Welch_P = p_t,
      Wilcox_P = p_w
    )
  }
)

stat_results <- do.call(
  rbind,
  stat_results
)

stat_results$Welch_FDR <- p.adjust(
  stat_results$Welch_P,
  method = "BH"
)

write.csv(
  stat_results,
  file.path(
    step_dir,
    "10_CellType_Proportion_Statistics.csv"
  ),
  row.names = FALSE
)

# ============================================================
# 15. 查看结果
# ============================================================

cat("\n====================================\n")
cat("第十步：细胞比例分析完成！\n")
cat("====================================\n")

cat("\n各组细胞平均比例：\n")
print(group_mean)

cat("\n颗粒细胞占比：\n")
print(gc_prop)

cat("\n组间统计结果：\n")
print(stat_results)

cat("\n所有图片均为PDF + SVG\n")
cat("结果目录：", step_dir, "\n")



# ============================================================
# STEP 11
# Granulosa细胞提取与Pseudobulk表达矩阵构建
#
# Granulosa = Cluster 1 + Cluster 7
#
# Age: Age1, Age2, Age_Trp1
# Age_Trp: Age_Trp2, Age_Trp3, Age_Trp4
#
# 不重新进行QC/Harmony/聚类
# 生成各样本的原始UMI汇总矩阵
# 输出图片：PDF + SVG
# ============================================================

library(Seurat)
library(Matrix)
library(ggplot2)

if (!requireNamespace("svglite", quietly = TRUE)) {
  install.packages("svglite")
}

# ------------------------------------------------------------
# 1. 路径
# ------------------------------------------------------------

out_dir <- paste0(
  "D:/R data/Age_vs_AgeTrp_scRNA/",
  "New_6Mice_3vs3_FromRaw"
)

step_dir <- file.path(
  out_dir,
  "11_Granulosa_DEG"
)

dir.create(
  step_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# ------------------------------------------------------------
# 2. 读取最终注释对象
# ------------------------------------------------------------

combined <- readRDS(
  file.path(
    out_dir,
    "09_Seurat_Final_Annotated.rds"
  )
)

DefaultAssay(combined) <- "RNA"

# ------------------------------------------------------------
# 3. 提取Cluster 1和7
# ------------------------------------------------------------

Idents(combined) <- "seurat_clusters"

granulosa <- subset(
  combined,
  idents = c("1", "7")
)

cat("\nGranulosa细胞总数：",
    ncol(granulosa), "\n")

cat("\n各Cluster细胞数量：\n")
print(table(granulosa$seurat_clusters))

# ------------------------------------------------------------
# 4. 确认六只鼠及分组
# ------------------------------------------------------------

samples <- c(
  "Age1",
  "Age2",
  "Age_Trp1",
  "Age_Trp2",
  "Age_Trp3",
  "Age_Trp4"
)

groups <- c(
  "Age", "Age", "Age",
  "Age_Trp", "Age_Trp", "Age_Trp"
)

sample_info <- data.frame(
  Sample = samples,
  Group = groups,
  stringsAsFactors = FALSE
)

granulosa$sample <- as.character(
  granulosa$sample
)

stopifnot(
  all(granulosa$sample %in% samples)
)

granulosa$Group6 <- factor(
  ifelse(
    granulosa$sample %in% samples[1:3],
    "Age",
    "Age_Trp"
  ),
  levels = c("Age", "Age_Trp")
)

sample_info$CellCount <- as.integer(
  table(
    factor(
      granulosa$sample,
      levels = samples
    )
  )
)

cat("\n===== 各样本Granulosa细胞数量 =====\n")
print(sample_info)

if (any(sample_info$CellCount == 0)) {
  stop("某只小鼠没有Granulosa细胞，无法完成6样本Pseudobulk分析。")
}

write.csv(
  sample_info,
  file.path(
    step_dir,
    "11_Granulosa_SampleInfo.csv"
  ),
  row.names = FALSE
)

# ------------------------------------------------------------
# 5. 原始UMI counts
# ------------------------------------------------------------

# Seurat v5多层RNA数据先合并
if (inherits(granulosa[["RNA"]], "Assay5")) {
  granulosa <- JoinLayers(
    granulosa,
    assay = "RNA"
  )
}

counts <- GetAssayData(
  granulosa,
  assay = "RNA",
  layer = "counts"
)

cat("\n原始UMI矩阵维度：\n")
print(dim(counts))

# ------------------------------------------------------------
# 6. 逐只鼠汇总原始UMI
#
# Pseudobulk:
# 每只鼠的同一基因UMI计数相加
#
# 最终得到：
# Genes × 6 mice
# ------------------------------------------------------------

pb_counts <- sapply(
  samples,
  function(s) {
    
    cell_ids <- colnames(granulosa)[
      granulosa$sample == s
    ]
    
    Matrix::rowSums(
      counts[, cell_ids, drop = FALSE]
    )
  }
)

rownames(pb_counts) <- rownames(counts)
colnames(pb_counts) <- samples

storage.mode(pb_counts) <- "integer"

cat("\n===== Pseudobulk矩阵 =====\n")
print(dim(pb_counts))

cat("\n每只鼠总UMI：\n")
print(colSums(pb_counts))

# ------------------------------------------------------------
# 7. 保存矩阵及样本信息
# ------------------------------------------------------------

saveRDS(
  granulosa,
  file.path(
    step_dir,
    "11_Granulosa_Seurat.rds"
  )
)

saveRDS(
  pb_counts,
  file.path(
    step_dir,
    "11_Granulosa_Pseudobulk_Counts.rds"
  )
)

write.csv(
  pb_counts,
  file.path(
    step_dir,
    "11_Granulosa_Pseudobulk_Counts.csv"
  )
)

# ------------------------------------------------------------
# 8. 绘制Granulosa UMAP
# 保留原UMAP坐标
# ------------------------------------------------------------

p_umap <- DimPlot(
  granulosa,
  reduction = "umap",
  group.by = "Group6",
  cols = c(
    "Age" = "#E66A61",
    "Age_Trp" = "#16A7B1"
  ),
  raster = FALSE
) +
  ggtitle("Granulosa cells: Age vs Age_Trp") +
  theme_classic(base_size = 13)

ggsave(
  file.path(
    step_dir,
    "11_Granulosa_UMAP.pdf"
  ),
  p_umap,
  device = grDevices::cairo_pdf,
  width = 9,
  height = 7,
  bg = "white"
)

ggsave(
  file.path(
    step_dir,
    "11_Granulosa_UMAP.svg"
  ),
  p_umap,
  device = svglite::svglite,
  width = 9,
  height = 7,
  bg = "white"
)

# ------------------------------------------------------------
# 9. 输出各样本细胞数量图
# ------------------------------------------------------------

sample_info$Sample <- factor(
  sample_info$Sample,
  levels = samples
)

sample_info$Group <- factor(
  sample_info$Group,
  levels = c("Age", "Age_Trp")
)

p_count <- ggplot(
  sample_info,
  aes(
    x = Sample,
    y = CellCount,
    fill = Group
  )
) +
  geom_col(
    width = 0.7
  ) +
  scale_fill_manual(
    values = c(
      "Age" = "#E66A61",
      "Age_Trp" = "#16A7B1"
    )
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(
      angle = 35,
      hjust = 1
    )
  ) +
  labs(
    title = "Granulosa cell counts",
    x = NULL,
    y = "Number of cells"
  )

ggsave(
  file.path(
    step_dir,
    "11_Granulosa_CellCounts.pdf"
  ),
  p_count,
  device = grDevices::cairo_pdf,
  width = 9,
  height = 6,
  bg = "white"
)

ggsave(
  file.path(
    step_dir,
    "11_Granulosa_CellCounts.svg"
  ),
  p_count,
  device = svglite::svglite,
  width = 9,
  height = 6,
  bg = "white"
)

# ------------------------------------------------------------
# 10. 完成
# ------------------------------------------------------------

cat("\n====================================\n")
cat("STEP 11 完成\n")
cat("Granulosa = Cluster 1 + Cluster 7\n")
cat("总细胞数：", ncol(granulosa), "\n")
cat("Pseudobulk样本数：", ncol(pb_counts), "\n")
cat("结果目录：", step_dir, "\n")
cat("====================================\n")



# ==========================================================
# STEP 12：Granulosa Pseudobulk edgeR差异表达分析
# Age_Trp vs Age
# n = 3 vs 3
# 输出 PDF + SVG
# ==========================================================

# 1. 安装并加载R包
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

if (!requireNamespace("edgeR", quietly = TRUE)) {
  BiocManager::install(
    "edgeR",
    ask = FALSE,
    update = FALSE
  )
}

if (!requireNamespace("svglite", quietly = TRUE)) {
  install.packages("svglite")
}

library(edgeR)
library(ggplot2)

# 2. 路径
out_dir <- paste0(
  "D:/R data/Age_vs_AgeTrp_scRNA/",
  "New_6Mice_3vs3_FromRaw"
)

step_dir <- file.path(out_dir, "12_Granulosa_edgeR")

dir.create(
  step_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

# 3. 读取Pseudobulk矩阵
pb <- readRDS(
  file.path(
    out_dir,
    "11_Granulosa_DEG",
    "11_Granulosa_Pseudobulk_Counts.rds"
  )
)

pb <- as.matrix(pb)

# 4. 确定分组和样本顺序
samples <- c(
  "Age1", "Age2", "Age_Trp1",
  "Age_Trp2", "Age_Trp3", "Age_Trp4"
)

stopifnot(all(samples %in% colnames(pb)))

pb <- pb[, samples, drop = FALSE]

stopifnot(
  ncol(pb) == 6,
  all(is.finite(pb)),
  all(pb >= 0),
  all(pb == round(pb))
)

group <- factor(
  c(rep("Age", 3), rep("Age_Trp", 3)),
  levels = c("Age", "Age_Trp")
)

design <- model.matrix(~ Group, data = data.frame(
  Group = group
))

rownames(design) <- samples

cat("\n分组设计矩阵：\n")
print(design)

# 5. 构建edgeR对象
y <- DGEList(
  counts = pb,
  group = group
)

# 6. 过滤低表达基因
keep <- filterByExpr(
  y,
  design = design
)

cat("\n过滤前基因数：", nrow(y), "\n")
cat("过滤后基因数：", sum(keep), "\n")

if (sum(keep) < 100) {
  stop("过滤后基因过少，请检查Pseudobulk矩阵。")
}

y <- y[keep, , keep.lib.sizes = FALSE]

# 7. TMM标准化
y <- calcNormFactors(
  y,
  method = "TMM"
)

# 8. 估计离散度
y <- estimateDisp(
  y,
  design = design
)

# 9. 准似然模型
fit <- glmQLFit(
  y,
  design = design
)

# 10. Age_Trp vs Age
qlf <- glmQLFTest(
  fit,
  coef = "GroupAge_Trp"
)

# 11. 全部检测基因的差异分析结果
deg <- topTags(
  qlf,
  n = nrow(y)
)$table

deg$Gene <- rownames(deg)

deg <- deg[
  ,
  c(
    "Gene",
    setdiff(colnames(deg), "Gene")
  )
]

# 12. 差异表达状态
logfc_cutoff <- log2(1.5)

deg$Status <- "Not significant"

deg$Status[
  deg$FDR < 0.05 &
    deg$logFC >= logfc_cutoff
] <- "Up"

deg$Status[
  deg$FDR < 0.05 &
    deg$logFC <= -logfc_cutoff
] <- "Down"

# 13. 保存结果
write.csv(
  deg,
  file.path(
    step_dir,
    "12_All_Genes_edgeR.csv"
  ),
  row.names = FALSE
)

sig <- deg[deg$FDR < 0.05, ]

write.csv(
  sig,
  file.path(
    step_dir,
    "12_FDR05_DEGs.csv"
  ),
  row.names = FALSE
)

sig_fc <- deg[
  deg$Status %in% c("Up", "Down"),
]

write.csv(
  sig_fc,
  file.path(
    step_dir,
    "12_FDR05_FC1.5_DEGs.csv"
  ),
  row.names = FALSE
)

up <- deg[deg$Status == "Up", ]

down <- deg[deg$Status == "Down", ]

write.csv(
  up,
  file.path(step_dir, "12_Upregulated.csv"),
  row.names = FALSE
)

write.csv(
  down,
  file.path(step_dir, "12_Downregulated.csv"),
  row.names = FALSE
)

nominal <- deg[
  deg$PValue < 0.05 &
    abs(deg$logFC) >= logfc_cutoff,
]

write.csv(
  nominal,
  file.path(
    step_dir,
    "12_NominalP05_FC1.5_Exploratory.csv"
  ),
  row.names = FALSE
)

# 14. 重点基因
focus_genes <- c(
  "Col1a1", "Nat10", "Esr2",
  "Fshr", "Amh", "Inha", "Foxl2",
  "Mki67", "Top2a"
)

focus <- deg[
  tolower(deg$Gene) %in% tolower(focus_genes),
]

write.csv(
  focus,
  file.path(
    step_dir,
    "12_Focus_Genes.csv"
  ),
  row.names = FALSE
)

# 15. 保存标准化表达矩阵
logcpm <- cpm(
  y,
  log = TRUE,
  prior.count = 2
)

saveRDS(
  logcpm,
  file.path(
    step_dir,
    "12_Normalized_logCPM.rds"
  )
)

# 16. 保存样本测序计数和标准化信息
sample_qc <- data.frame(
  Sample = samples,
  Group = as.character(group),
  LibrarySize = y$samples$lib.size,
  TMMFactor = y$samples$norm.factors,
  EffectiveLibrarySize =
    y$samples$lib.size * y$samples$norm.factors
)

write.csv(
  sample_qc,
  file.path(
    step_dir,
    "12_Sample_Library_QC.csv"
  ),
  row.names = FALSE
)

# 17. 保存图片函数
save_fig <- function(p, filename, width, height) {
  
  ggsave(
    file.path(step_dir, paste0(filename, ".pdf")),
    plot = p,
    device = grDevices::cairo_pdf,
    width = width,
    height = height,
    bg = "white"
  )
  
  ggsave(
    file.path(step_dir, paste0(filename, ".svg")),
    plot = p,
    device = svglite::svglite,
    width = width,
    height = height,
    bg = "white"
  )
}

# 18. Pseudobulk PCA
pca <- prcomp(
  t(logcpm),
  center = TRUE,
  scale. = FALSE
)

variance <- 100 *
  pca$sdev^2 / sum(pca$sdev^2)

pca_df <- data.frame(
  Sample = samples,
  Group = group,
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2]
)

write.csv(
  pca_df,
  file.path(step_dir, "12_PCA_Coordinates.csv"),
  row.names = FALSE
)

group_colors <- c(
  Age = "#E66A61",
  Age_Trp = "#16A7B1"
)

p_pca <- ggplot(
  pca_df,
  aes(x = PC1, y = PC2, color = Group)
) +
  geom_point(size = 3.5) +
  geom_text(
    aes(label = Sample),
    vjust = -0.9,
    size = 3.5,
    show.legend = FALSE
  ) +
  scale_color_manual(values = group_colors) +
  labs(
    title = "Granulosa pseudobulk PCA",
    x = paste0("PC1 (", round(variance[1], 1), "%)"),
    y = paste0("PC2 (", round(variance[2], 1), "%)")
  ) +
  theme_classic(base_size = 13)

save_fig(
  p_pca,
  "12_Pseudobulk_PCA",
  9, 7
)

# 19. 火山图
deg$negLog10FDR <- -log10(
  pmax(deg$FDR, .Machine$double.xmin)
)

volcano_colors <- c(
  "Up" = "#E64B35",
  "Down" = "#3C78B5",
  "Not significant" = "#BDBDBD"
)

label_df <- sig_fc[
  order(sig_fc$FDR),
  ,
  drop = FALSE
]

label_df <- head(label_df, 12)

p_volcano <- ggplot(
  deg,
  aes(
    x = logFC,
    y = negLog10FDR,
    color = Status
  )
) +
  geom_point(
    size = 0.9,
    alpha = 0.7
  ) +
  geom_vline(
    xintercept = c(-logfc_cutoff, logfc_cutoff),
    linetype = "dashed",
    color = "grey50"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dashed",
    color = "grey50"
  ) +
  geom_text(
    data = label_df,
    aes(label = Gene),
    size = 3,
    vjust = -0.6,
    check_overlap = TRUE,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = volcano_colors
  ) +
  labs(
    title = "Granulosa: Age_Trp vs Age",
    subtitle = "edgeR pseudobulk | FDR < 0.05 | FC >= 1.5",
    x = "log2 Fold Change (Age_Trp / Age)",
    y = "-log10(FDR)",
    color = "DEG status"
  ) +
  theme_classic(base_size = 13)

save_fig(
  p_volcano,
  "12_Granulosa_Volcano",
  10, 8
)

# 20. 统计结果
summary_table <- data.frame(
  Category = c(
    "Tested genes",
    "FDR < 0.05",
    "FDR < 0.05 and FC >= 1.5",
    "Upregulated",
    "Downregulated",
    "Nominal P < 0.05 and FC >= 1.5"
  ),
  Count = c(
    nrow(deg),
    nrow(sig),
    nrow(sig_fc),
    nrow(up),
    nrow(down),
    nrow(nominal)
  )
)

write.csv(
  summary_table,
  file.path(
    step_dir,
    "12_DEG_Summary.csv"
  ),
  row.names = FALSE
)

# 21. 保存edgeR模型
saveRDS(
  list(
    y = y,
    design = design,
    fit = fit,
    qlf = qlf
  ),
  file.path(
    step_dir,
    "12_edgeR_Analysis.rds"
  )
)

cat("\n========== DEG结果 ==========\n")
print(summary_table)

cat("\n========== 重点基因 ==========\n")
print(focus)

cat("\n========== 样本测序量 ==========\n")
print(sample_qc)

cat("\nSTEP 12 完成！\n")