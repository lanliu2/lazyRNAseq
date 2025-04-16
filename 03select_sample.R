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
filename <- "tissue_Gills_treatment_naive_genotype_WT_vs_TG_counts"
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

# 给新样本命名
filename <- "tissue_Gills_treatment_naive_genotype_WT_vs_TG_counts"

### 保存只需要的样本 ##########
# 指定你想保留的样本列名
selected_samples <- c("NG10", "NG11", "PG10", "PG12")

# 从 input_data 中筛选这些列
subset_data <- input_data[, selected_samples]

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
  tissue = "Gills",
  treatment = "naive",
  genotype = c("WT", "WT", "TG", "TG"),
  group = c("WT", "WT", "TG", "TG")
)

################6.构建 DESeq2 数据集#############
dds <- DESeqDataSetFromMatrix(
  countData = round(counts),
  colData = sample_info,
  design = ~ group
)

##############PCA分析##############
vsd <- vst(dds, blind = FALSE)
pcaData <- plotPCA(vsd, intgroup = "group", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

p1 <- ggplot(pcaData, aes(x = PC1, y = PC2, color = group, label = name)) +
  geom_point(size = 3) +
  geom_text_repel(size = 3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA: WT vs TG in Gills") +
  theme_bw()

ggsave("01pca.pdf",p1,dpi=800)

###########8. 聚类分析（热图）##################
sampleDists <- dist(t(assay(vsd)))  # 样本之间距离矩阵
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(vsd)
colnames(sampleDistMatrix) <- colnames(vsd)

# 设置颜色
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)

# 画热图
p2 <- Heatmap(sampleDistMatrix,
              name = "Distance",  # 图例标题
              col = colorRamp2(c(min(sampleDistMatrix), max(sampleDistMatrix)), c("#084594","white" )),
              clustering_distance_rows = sampleDists,
              clustering_distance_columns = sampleDists,
              column_names_side = "bottom",
              row_names_side = "left",
              column_title = "Sample Distance Heatmap")

# 保存热图

# 转换为 ggplot 对象
p_gg2 <- as.ggplot(p2)

# 保存（现在可以用 ggsave）
ggsave("02Sample_Distance_Heatmap_ggplot.pdf", p_gg2, width = 9, height = 8)



############差异表达分析#############
dds <- DESeq(dds)
res <- results(dds, contrast = c("group", "TG", "WT"))  # TG vs WT

# 去除 NA
res <- na.omit(res)
res_df <- as.data.frame(res)

# 添加 gene symbol 列
res_df$gene <- rownames(res_df)

# 添加差异表达类别标签
res_df$category <- ifelse(res_df$padj < 0.1 & res_df$log2FoldChange > 1, "up",
                          ifelse(res_df$padj < 0.1 & res_df$log2FoldChange < -1, "down", "normal"))

table(res_df$category)

# 调整列顺序（可选）
res_df <- res_df[, c("gene", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "category")]

# 保存为CSV
write.csv(res_df, file = paste0("DEG_", filename, ".csv"), row.names = FALSE)

##############







