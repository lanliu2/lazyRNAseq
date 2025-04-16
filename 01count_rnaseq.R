#import library
library(ggplot2)
library(DESeq2)
library(clusterProfiler)
library(org.Dr.eg.db)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(vegan)
library(ape)             # 用于 pcoa()
library(RColorBrewer)
library(Rtsne)
library(ggforce)


### 1.环境配置 ###
setwd("D:/desktop/工作/冉晨瑞/RNA_202504/新的注释文件/02质控样本/")  #设置工作路径

### 2.Ensembl_id转换为gene_name ###
## 2.1读入上游生成的count.txt文件,更改列名 ##
input_data <- read.csv("gene_count_matrix.csv",
                       header=TRUE, row.names=1) #将基因名作为行名


sample_names <- colnames(input_data)
# > sample_names
# [1] "NG1"  "NG10" "NG11" "NG12" "NG2"  "NG3"  "NI1"  "NI10" "NI11" "NI12" "NI2"  "NI3"  "PG1"  "PG10"
# [15] "PG11" "PG12" "PG2"  "PG3"  "PI1"  "PI10" "PI11" "PI12" "PI2"  "PI3" 


#重命名样本名称
#查看基因是否重复 事实证明stringtie的结果没有重复基因
symbols <- rownames(input_data)
table(duplicated(symbols))
# FALSE 
# 31338 

## 清除基因名的重复显示问题 ##
row.names(input_data) <- sub("\\|.*", "", row.names(input_data))

### 3.初步过滤低表达基因（筛选标准不唯一、依情况而定）###
## 3.1筛选出至少在重复样本数量内的表达量counts大于1的行（基因）##
keep_feature <- rowSums(input_data>1) >= 2
table(keep_feature)  #查看筛选情况，FALSE为低表达基因数（行数），TURE为要保留基因数
# keep_feature
# FALSE  TRUE 
# 4923 26415 
counts_filt <- input_data[keep_feature, ] #替换counts为筛选后的基因矩阵（保留较高表达量的基因）
#### 4.保存数据 ####
counts_raw=input_data#这里重新命名方便后续分析调用
counts=counts_filt
write.table(counts,"counts.txt",sep="\t")


#### DESeq差异分析 ####

# 设置分组
sample_info <- data.frame(
  sample = c("WT-Vi-G-1", "WT-G-1", "WT-G-2", "WT-G-3", "WT-Vi-G-2", "WT-Vi-G-3", "WT-Vi-I-1", "WT-I-1", "WT-I-2", "WT-I-3", "WT-Vi-I-2", "WT-Vi-I-3",
             "TG-Vi-G-1", "TG-G-1", "TG-G-2", "TG-G-3", "TG-Vi-G-2", "TG-Vi-G-3", "TG-Vi-I-1", "TG-I-1", "TG-I-2", "TG-I-3", "TG-Vi-I-2", "TG-Vi-I-3"),
  
  # sample = c("NG1", "NG10", "NG11", "NG12", "NG2", "NG3", "NI1", "NI10", "NI11", "NI12", "NI2", "NI3",
  #            "PG1", "PG10", "PG11", "PG12", "PG2", "PG3", "PI1", "PI10", "PI11", "PI12", "PI2", "PI3"),
  
  genotype = c(rep("WT",12), rep("TG", 12)), 
  treatment = c("Soak24h", "Naive", "Naive", "Naive", "Soak24h", "Soak24h",  # NG系列
                "Soak24h", "Naive", "Naive", "Naive", "Soak24h", "Soak24h",  # NI系列
                "Soak24h", "Naive", "Naive", "Naive", "Soak24h", "Soak24h",  # PG系列
                "Soak24h", "Naive", "Naive", "Naive", "Soak24h", "Soak24h"), # PI系列
  tissue = c(rep("Gills", 6), rep("Intestine", 6),rep("Gills", 6), rep("Intestine", 6))  # Gills 和 Intestine
)

rownames(sample_info) <- sample_info$sample
write.csv(sample_info, "sample_info.csv")

# 将 genotype, treatment, tissue 列转换为因子类型
sample_info$genotype <- factor(sample_info$genotype)
sample_info$treatment <- factor(sample_info$treatment)
sample_info$tissue <- factor(sample_info$tissue)

# 确保 counts 列名与 sample_info$sample 一致
colnames(counts) <- sample_info$sample

# 使用DESeq2包建立DESeqDataSet对象
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = sample_info,  # 使用 sample_info 作为 colData
                              design = ~ genotype + treatment + tissue)


###### PCA 聚类 ###########
# 1. 使用 VST 或 rlog 转换数据（推荐 VST，效率更高）
vsd <- vst(dds, blind = FALSE)  # blind=FALSE 保留设计信息用于样本分组可视化

# 2. 绘制 PCA 图（ggplot2）
pcaData <- plotPCA(vsd, intgroup = c("genotype", "treatment", "tissue"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# 添加样本名称标签
p <- ggplot(pcaData, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = treatment, shape = genotype), size = 2) +
  geom_text_repel(aes(label = row.names(pcaData)), size = 2, max.overlaps = 90) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA of RNA-seq Samples") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))






# p_pcoa <- ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2)) +
#   geom_point(aes(color = treatment, shape = genotype), size = 2) +
#   geom_text_repel(aes(label = Sample), size = 2.5, max.overlaps = 100) +
#   stat_ellipse(aes(group = group, color = treatment), level = 0.90, linetype = "dashed", size = 0.8, alpha = 0.4) +
#   xlab(paste0("PCoA1: ", percentVar[1], "%")) +
#   ylab(paste0("PCoA2: ", percentVar[2], "%")) +
#   ggtitle("PCoA of RNA-seq Samples") +
#   theme_bw() +
#   theme(plot.title = element_text(hjust = 0.5))


p

# 保存 PCA 图
ggsave("PCA_plot.pdf", plot = p, width = 8, height = 6)

################pcoa########################
# # 1. 构建 DESeqDataSet 并做 VST 转换
# dds <- DESeqDataSetFromMatrix(countData = counts,
#                               colData = sample_info,
#                               design = ~ genotype + treatment + tissue)
# vsd <- vst(dds, blind = FALSE)
# 1. 获取表达矩阵
expr_matrix <- assay(vsd)

# 2. 计算样本间距离矩阵（欧式距离）
dist_matrix <- vegdist(t(expr_matrix), method = "euclidean")

# 3. 执行 PCoA 分析
pcoa_res <- ape::pcoa(dist_matrix)

# 4. 整理 PCoA 输出数据
pcoa_df <- as.data.frame(pcoa_res$vectors[, 1:2])
pcoa_df$Sample <- rownames(pcoa_df)
pcoa_df <- merge(pcoa_df, as.data.frame(colData(vsd)), by.x = "Sample", by.y = "row.names")

# 5. 添加组合分组变量（四组）
pcoa_df$group <- paste0(pcoa_df$tissue, "_", pcoa_df$treatment)

# 6. 可视化 PCoA + 置信椭圆
percentVar <- round(100 * pcoa_res$values$Relative_eig[1:2], 2)

p_pcoa <- ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2)) +
  geom_point(aes(color = treatment, shape = genotype), size = 2) +
  geom_text_repel(aes(label = Sample), size = 2.5, max.overlaps = 100) +
  stat_ellipse(aes(group = group, color = treatment), level = 0.90, linetype = "dashed", size = 0.8, alpha = 0.4) +
  xlab(paste0("PCoA1: ", percentVar[1], "%")) +
  ylab(paste0("PCoA2: ", percentVar[2], "%")) +
  ggtitle("PCoA of RNA-seq Samples") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
p_pcoa

# 7. 保存图像
ggsave("PCoA_plot_with_90CI.pdf", plot = p_pcoa, width = 8, height = 6, dpi = 300)

######################################

##########################热图绘制####################
# 计算皮尔森相关性矩阵
cor_matrix <- cor(assay(vsd), method = "pearson")
cor_dist <- as.dist(1 - cor_matrix)

# 格式化数字标签
formatted_matrix <- matrix(sprintf("%.1f", cor_matrix),
                           nrow = nrow(cor_matrix),
                           dimnames = dimnames(cor_matrix))

# 分组信息（样本注释），行名应为样本名
annotation_df <- sample_info[, c("genotype", "treatment", "tissue")]

# 设置分组颜色（可选）
ann_colors <- list(
  genotype = c(WT = "#d6bdfb", TG = "#9dfceb"),
  treatment = c(Naive = "tomato", Soak24h = "skyblue"),
  tissue = c(Gills = "lightgreen", Intestine = "yellow")
)

# 自定义颜色渐变：white -> #d6bdfb
# my_colors <- colorRampPalette(c("#eee1ff","#250fff"))(100)

# 绘图
pheatmap(cor_matrix,
         # color = my_colors,                      # ✅ 自定义颜色
         clustering_distance_rows = cor_dist,
         clustering_distance_cols = cor_dist,
         display_numbers = formatted_matrix,
         number_color = "black",
         annotation_col = annotation_df,
         annotation_row = annotation_df,
         annotation_colors = ann_colors,
         main = "Pearson Correlation Heatmap (Custom Colors)",
         fontsize = 10)


# 保存为 PDF
pdf("Pearson_Correlation_Heatmap.pdf", width = 9, height = 8)
pheatmap(cor_matrix,
         # color = my_colors,                      # ✅ 自定义颜色
         clustering_distance_rows = cor_dist,
         clustering_distance_cols = cor_dist,
         display_numbers = formatted_matrix,
         number_color = "black",
         annotation_col = annotation_df,
         annotation_row = annotation_df,
         annotation_colors = ann_colors,
         main = "Pearson Correlation Heatmap (Custom Colors)",
         fontsize = 10)
dev.off()

write.csv(pcaData, file = "PCA_data.csv", row.names = FALSE)
write.csv(pcoa_df, file = "PCoA_data.csv", row.names = FALSE)
write.csv(cor_matrix, file = "Correlation_matrix.csv")
write.csv(as.matrix(cor_dist), file = "Correlation_distance_matrix.csv")



##################全基因热图#####################

######################################################

######################批量绘制pca和pcoa小组####################################
run_conditional_pca <- function(dds,
                                fixed_factors = list(genotype = "WT", tissue = "Gills"),
                                variable_factor = "treatment",
                                output_dir = "PCA_conditional") {
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # ---- 1. 子集筛选 ----
  sel_samples <- rep(TRUE, ncol(dds))
  for (fac in names(fixed_factors)) {
    sel_samples <- sel_samples & (colData(dds)[[fac]] %in% fixed_factors[[fac]])
  }
  dds_sub <- dds[, sel_samples]
  
  # 去除 factor 多余水平
  for (col in colnames(colData(dds_sub))) {
    if (is.factor(colData(dds_sub)[[col]])) {
      colData(dds_sub)[[col]] <- droplevels(colData(dds_sub)[[col]])
    }
  }
  
  # ---- 2. 检查分组是否足够 ----
  if (length(unique(colData(dds_sub)[[variable_factor]])) < 2) {
    warning(paste0("Not enough levels in ", variable_factor, " after filtering. Skip."))
    return(NULL)
  }
  
  # ---- 3. vst变换 ----
  design(dds_sub) <- as.formula(paste("~", variable_factor))
  vsd <- vst(dds_sub, blind = FALSE)
  mat <- assay(vsd)
  
  # -------------------------------------
  # ✅ PCA 图
  # -------------------------------------
  pcaData <- plotPCA(vsd, intgroup = variable_factor, returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  pcaData$name <- rownames(pcaData)
  
  p_pca <- ggplot(pcaData, aes(x = PC1, y = PC2, color = .data[[variable_factor]], label = name)) +
    geom_point(size = 4) +
    stat_ellipse(aes(group = .data[[variable_factor]]), level = 0.90, type = "norm", linetype = "dashed", linewidth = 0.8) +
    geom_text_repel(size = 3, max.overlaps = 100) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    ggtitle(paste0("PCA | ", 
                   paste(paste0(names(fixed_factors), "=", fixed_factors), collapse = ", "),
                   " | by ", variable_factor)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  
  # -------------------------------------
  # ✅ PCoA 图（使用 Bray-Curtis 距离）
  # -------------------------------------
  dist_mat <- vegdist(t(mat), method = "bray")
  pcoa_coords <- cmdscale(dist_mat, k = 2, eig = TRUE)
  pcoaData <- as.data.frame(pcoa_coords$points)
  colnames(pcoaData) <- c("PCoA1", "PCoA2")
  pcoaData[[variable_factor]] <- colData(dds_sub)[[variable_factor]]
  pcoaData$name <- rownames(pcoaData)
  
  eigs <- round(100 * pcoa_coords$eig[1:2] / sum(pcoa_coords$eig[pcoa_coords$eig > 0]))
  
  p_pcoa <- ggplot(pcoaData, aes(x = PCoA1, y = PCoA2, color = .data[[variable_factor]], label = name)) +
    geom_point(size = 4) +
    stat_ellipse(aes(group = .data[[variable_factor]]), level = 0.90, type = "norm", linetype = "dashed", linewidth = 0.8) +
    geom_text_repel(size = 3, max.overlaps = 100) +
    xlab(paste0("PCoA1: ", eigs[1], "% variance")) +
    ylab(paste0("PCoA2: ", eigs[2], "% variance")) +
    ggtitle(paste0("PCoA | ", 
                   paste(paste0(names(fixed_factors), "=", fixed_factors), collapse = ", "),
                   " | by ", variable_factor)) +
    theme_bw() +
    theme(plot.title = element_text(hjust = 0.5))
  
  # ---- 4. 输出文件 ----
  suffix <- paste(paste0(names(fixed_factors), "_", fixed_factors), collapse = "_")
  
  file_pca <- file.path(output_dir, paste0("PCA_", suffix, "_by_", variable_factor, ".pdf"))
  file_pcoa <- file.path(output_dir, paste0("PCoA_", suffix, "_by_", variable_factor, ".pdf"))
  
  ggsave(file_pca, plot = p_pca, width = 8, height = 6)
  ggsave(file_pcoa, plot = p_pcoa, width = 8, height = 6)
}



for (gen in c("WT", "TG")) {
  for (tis in c("Gills", "Intestine")) {
    run_conditional_pca(dds,
                        fixed_factors = list(genotype = gen, tissue = tis),
                        variable_factor = "treatment")
  }
}


for (tre in c("Naive", "Soak24h")) {
  for (tis in c("Gills", "Intestine")) {
    run_conditional_pca(dds,
                        fixed_factors = list(tissue = tis,treatment = tre),
                        variable_factor = "genotype")
  }
}

for (tre in c("Naive", "Soak24h")) {
  for (gen in c("WT", "TG")) {
    run_conditional_pca(dds,
                        fixed_factors = list(treatment = tre, genotype = gen),
                        variable_factor = "tissue")
  }
}




###################################

# 运行DESeq2差异分析
dds <- DESeq(dds)

# 获取差异表达结果
res <- results(dds)

resOrdered <- res[order(res$padj), ]

#创建一个专门用于存放差异表达输出结果的目录
dir.create("./diffout", showWarnings = FALSE, recursive = TRUE)

# 定义一个用于导出差异表达结果并添加分类的函数
export_differential_results <- function(dds, contrast, output_file, significance_level = 0.05, log2fc_threshold = 1) {
  # 获取差异表达结果
  res <- results(dds, contrast = contrast)
  
  # 添加类别列 (up, down, normal)
  res$category <- "normal"  # 默认类别为 "normal"
  
  # 根据 log2 fold change 和 p-value 分类
  res$category[res$log2FoldChange > log2fc_threshold & res$padj < significance_level] <- "up"
  res$category[res$log2FoldChange < -log2fc_threshold & res$padj < significance_level] <- "down"
  
  # 对结果按调整后的 p 值排序
  resOrdered <- res[order(res$padj), ]
  
  # 将结果保存为 CSV 文件
  write.csv(as.data.frame(resOrdered), file = output_file)
}

# 创建目录用于存储输出结果（如果不存在）
dir.create("./diffout", showWarnings = FALSE, recursive = TRUE)

# 运行差异表达分析并导出不同分组的结果

# 比较 Tg vs wt
export_differential_results(dds, c("genotype", "TG", "WT"), "./diffout/genotype_comparison_Tg_vs_wt.csv")

# 比较 Soak24h vs Naive
export_differential_results(dds, c("treatment", "Soak24h", "Naive"), "./diffout/treatment_comparison_Soak24h_vs_Naive.csv")
  
# 比较 Gills vs Intestine
export_differential_results(dds, c("tissue", "Gills", "Intestine"), "./diffout/tissue_comparison_Gills_vs_Intestine.csv")



### 划分组织比较 #####

# 提取肠组织（Intestine）样本
intestine_samples <- sample_info[sample_info$tissue == "Intestine", ]
intestine_dds <- dds[, rownames(intestine_samples)]  # 筛选肠组织样本对应的 dds 数据

gills_samples <- sample_info[sample_info$tissue == "Gills",]
gills_dds <- dds[, rownames(gills_samples)]  # 筛选腮组织样本对应的 dds 数据

export_differential_results(intestine_dds, c("treatment", "Soak24h", "Naive"),
                            "./diffout/intestine_comparison_Soak24h_vs_Naive.csv")
export_differential_results(gills_dds, c("treatment", "Soak24h", "Naive"),
                            "./diffout/gills_comparison_Soak24h_vs_Naive.csv")


######################小组划分#######################################
run_conditional_deseq <- function(dds,
                                  fixed_factors = list(genotype = "WT", tissue = "Gills"),
                                  variable_factor = "treatment",
                                  group1 = "Soak24h",
                                  group2 = "Naive",
                                  output_dir = "diffout",
                                  significance_level = 0.05,
                                  log2fc_threshold = 1) {

  
  # 创建输出目录
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # ---- 1. 子集筛选 ----
  sel_samples <- rep(TRUE, ncol(dds))
  for (fac in names(fixed_factors)) {
    sel_samples <- sel_samples & (colData(dds)[[fac]] %in% fixed_factors[[fac]])
  }
  dds_sub <- dds[, sel_samples]
  
  # 去除未使用水平
  for (col in colnames(colData(dds_sub))) {
    if (is.factor(colData(dds_sub)[[col]])) {
      colData(dds_sub)[[col]] <- droplevels(colData(dds_sub)[[col]])
    }
  }
  
  # ---- 2. 确保比较组存在 ----
  var_levels <- levels(colData(dds_sub)[[variable_factor]])
  if (!all(c(group1, group2) %in% var_levels)) {
    warning(paste("Skip comparison:", group1, "vs", group2, "in", variable_factor, "— not present in this subset"))
    return(NULL)
  }
  
  # ---- 3. 设定 design 并运行差异分析 ----
  design(dds_sub) <- as.formula(paste("~", variable_factor))
  dds_sub <- DESeq(dds_sub)
  
  # ---- 4. 获取差异结果 ----
  res <- results(dds_sub, contrast = c(variable_factor, group1, group2))
  
  # ---- 5. 添加分类 ----
  res$category <- "normal"
  res$category[res$log2FoldChange > log2fc_threshold & res$padj < significance_level] <- "up"
  res$category[res$log2FoldChange < -log2fc_threshold & res$padj < significance_level] <- "down"
  resOrdered <- res[order(res$padj), ]
  
  # ---- 6. 导出结果 ----
  suffix <- paste(paste0(names(fixed_factors), "_", fixed_factors), collapse = "_")
  file_name <- paste0("DEG_", suffix, "_", variable_factor, "_", group1, "_vs_", group2, ".csv")
  out_path <- file.path(output_dir, file_name)
  
  write.csv(as.data.frame(resOrdered), file = out_path)
  message(paste("✓ Exported:", out_path))
}

for (tis in c("Gills", "Intestine")) {
  for (gen in c("WT", "TG")) {
    run_conditional_deseq(dds,
                          fixed_factors = list(tissue = tis,genotype = gen),
                          variable_factor = "treatment",
                          group1 = "Soak24h",
                          group2 = "Naive")
  }
}

for (tis in c("Gills", "Intestine")) {
  for (tre in c("Soak24h", "Naive")) {
    run_conditional_deseq(dds,
                          fixed_factors = list(treatment = tre, tissue = tis),
                          variable_factor = "genotype",
                          group1 = "TG",
                          group2 = "WT")
  }
}

for (tre in c("Soak24h", "Naive")) {
  for (gen in c("WT", "TG")) {
    run_conditional_deseq(dds,
                          fixed_factors = list(genotype = gen, treatment = tre),
                          variable_factor = "tissue",
                          group1 = "Intestine",
                          group2 = "Gills")
  }
}
