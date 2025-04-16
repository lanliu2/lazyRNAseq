#import library
library(ggplot2)
library(DESeq2)
library(clusterProfiler)
library(org.Dr.eg.db)
library(pheatmap)
library(RColorBrewer)
library(ggrepel)
library(ggplotify)
library(grid)
library(ComplexHeatmap)
library(circlize)


### 1.环境配置 ###
setwd("D:/desktop/工作/冉晨瑞/RNA_202504/新的注释文件/02质控样本/")  #设置工作路径
# 创建保存目录
filename <- "Intestine_12samples"
# filename <- "Gills_12samples"
dir.create(paste0("./",filename), showWarnings = FALSE, recursive = TRUE)
setwd(file.path(getwd(),filename))

### 2.Ensembl_id转换为gene_name ###
## 2.1读入上游生成的count.txt文件,更改列名 ##
input_data <- read.csv("../gene_count_matrix.csv",
                       header=TRUE, row.names=1) #将基因名作为行名
table(duplicated(rownames(input_data)))
# FALSE 
# 53709

sample_names <- colnames(input_data)
# > sample_names
# [1] "NG1"  "NG10" "NG11" "NG12" "NG2"  "NG3"  "NI1"  "NI10" "NI11" "NI12" "NI2"  "NI3"  "PG1"  "PG10"
# [15] "PG11" "PG12" "PG2"  "PG3"  "PI1"  "PI10" "PI11" "PI12" "PI2"  "PI3" 


### 保存只需要的样本 ##########
# 指定你想保留的样本列名
# selected_samples <- c("NG10", "NG11","NG12", "PG10", "PG11","PG12")
# selected_samples <- c("NG1","NG2", "NG3", "NG10", "NG11","NG12",
#                       "PG1", "PG2", "PG3", "PG10", "PG11","PG12")
selected_samples <- c("NI1","NI2", "NI3", "NI10", "NI11","NI12",
                      "PI1", "PI2", "PI3", "PI10", "PI11","PI12")

# 从 input_data 中筛选这些列
subset_data <- input_data[, selected_samples]
colnames(subset_data) <- c("WT-Vi-I-1","WT-Vi-I-2","WT-Vi-I-3","WT-I-1", "WT-I-2","WT-I-3",
                          "TG-Vi-I-1","TG-Vi-I-2","TG-Vi-I-3","TG-I-1", "TG-I-2","TG-I-3")
# 将筛选后的数据保存为新 CSV 文件
write.csv(subset_data, file = paste0(filename,".csv"))


#重命名样本名称
#查看基因是否重复 事实证明stringtie的结果没有重复基因
symbols <- rownames(subset_data)
table(duplicated(symbols))
# FALSE 
# 31338 

## 清除基因名的重复显示问题 ##
row.names(subset_data) <- sub("\\|.*", "", row.names(subset_data))

### 3.初步过滤低表达基因（筛选标准不唯一、依情况而定）###
## 3.1筛选出至少在重复样本数量内的表达量counts大于1的行（基因）##
keep_feature <- rowSums(subset_data>1) >= 2
table(keep_feature)  #查看筛选情况，FALSE为低表达基因数（行数），TURE为要保留基因数
# keep_feature
# FALSE  TRUE 
# 4923 26415 
counts_filt <- subset_data[keep_feature, ] #替换counts为筛选后的基因矩阵（保留较高表达量的基因）
#### 4.保存数据 ####
counts_raw=subset_data#这里重新命名方便后续分析调用
counts=counts_filt


###############构建分组信息###################
sample_info <- data.frame(
  row.names = colnames(counts),
  tissue = "Intestine",
  # tissue = "Gills",
  treatment = c(rep("Soak24h", 3), rep("Naive", 3), rep("Soak24h", 3), rep("Naive", 3)),
  genotype = c(rep("WT", 6), rep("TG",6)),
  group = c(rep("WT", 6), rep("TG",6))
)
# sample_info <- data.frame(
#   row.names = colnames(counts),
#   tissue = "Gills",
#   treatment = "naive",
#   genotype = c("WT", "WT", "WT", "TG", "TG", "TG"),
#   group = c("WT", "WT", "WT", "TG", "TG", "TG")
# )

################6.构建 DESeq2 数据集#############
dds <- DESeqDataSetFromMatrix(
  countData = round(counts),
  colData = sample_info,
  design = ~ group
)

##############PCA分析##############
vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = c("treatment", "genotype"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

p1 <- ggplot(pcaData, aes(x = PC1, y = PC2, color = treatment, shape=genotype, label = name)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA: WT vs TG in Gills") +
  theme_bw()

ggsave("01pca.pdf",p1,dpi=800)

p2 <- ggplot(pcaData, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = treatment, shape = genotype), size = 3.5) +
  geom_text_repel(aes(label = name), size = 3) +
  stat_ellipse(aes(group = treatment, color = treatment), level = 0.90, linetype = "dashed", size = 1) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA: WT vs TG in Intestine") +
  theme_bw()

# 保存图像
ggsave("PCA_with_90CI.pdf", p2, dpi = 800)


###########8. 聚类分析（热图）##################
### ========= 样本相关性热图 (皮尔森) ========= ###
cor_matrix <- cor(assay(vsd), method = "pearson")
cor_dist <- as.dist(1 - cor_matrix)
formatted_matrix <- matrix(sprintf("%.2f", cor_matrix),
                           nrow = nrow(cor_matrix),
                           dimnames = dimnames(cor_matrix))

annotation_df <- data.frame(Treatment = pcaData$treatment, Genotype= pcaData$genotype)
rownames(annotation_df) <- rownames(pcaData)

ann_colors <- list(
  Treatment = c("Naive" = "#619cff", "Soak24h" = "#f8766d"),
  Genotype = c("WT" = "green", "TG" = "pink")
)

# pdf("02Pearson_Correlation_Heatmap.pdf", width = 9, height = 8)
p3 <- pheatmap(cor_matrix,
         clustering_distance_rows = cor_dist,
         clustering_distance_cols = cor_dist,
         display_numbers = formatted_matrix,
         number_color = "black",
         annotation_col = annotation_df,
         annotation_row = annotation_df,
         annotation_colors = ann_colors,
         main = "Pearson Correlation Heatmap",
         fontsize = 10)
# dev.off()
pp3 <- as.ggplot(p3)
ggsave("02Pearson_Correlation_Heatmap.pdf", pp3, width = 9, height = 8, dpi=800)



############差异表达分析#############
dds <- DESeq(dds)
res <- results(dds, contrast = c("group", "TG", "WT"))  # TG vs WT

# 去除 NA
res <- na.omit(res)
res_df <- as.data.frame(res)

# 添加 gene symbol 列
res_df$gene <- rownames(res_df)

# 添加差异表达类别标签
res_df$category <- ifelse(res_df$padj < 0.05 & res_df$log2FoldChange > 1, "up",
                          ifelse(res_df$padj < 0.05 & res_df$log2FoldChange < -1, "down", "normal"))

table(res_df$category)

# 调整列顺序（可选）
res_df <- res_df[, c("gene", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "category")]

# 保存为CSV
write.csv(res_df, file = paste0("DEG_", filename, ".csv"), row.names = FALSE)

##############







