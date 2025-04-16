#导入包
library(ggplot2)
library(DESeq2)
library(clusterProfiler)
library(org.Dr.eg.db)
library(dplyr)
library(enrichplot)

# 设置工作路径
setwd("D:/desktop/工作/冉晨瑞/RNA_202504/新的注释文件/02质控样本/")
filename <- "treatment_Naive_tissue_Gills_genotype_TG_vs_WT"
dir.create(paste0("./",filename), showWarnings = FALSE, recursive = TRUE)
# tissue_Gills_genotype_WT_treatment_Soak24h_vs_Naive
# tissue_Intestine_genotype_WT_treatment_Soak24h_vs_Naive
# tissue_Gills_genotype_TG_treatment_Soak24h_vs_Naive
# tissue_Intestine_genotype_TG_treatment_Soak24h_vs_Naive

# treatment_Soak24h_tissue_Intestine_genotype_TG_vs_WT
# treatment_Soak24h_tissue_Gills_genotype_TG_vs_WT
# treatment_Naive_tissue_Intestine_genotype_TG_vs_WT
# treatment_Naive_tissue_Gills_genotype_TG_vs_WT




setwd(file.path(getwd(),filename))
#导入差异表达数据
file <- paste0("DEG_",filename,".csv")
plottitle <- gsub(".csv","",file)
data <- read.csv(paste0("../diffout/",file), header = TRUE, row.names = 1)
#gene=30550

##################################差异基因火山图######################
dir.create("./volcano", showWarnings = FALSE, recursive = TRUE)
plot_volcano <- function(res, title) {
  # 设置阈值
  significance_level <- 0.05
  log2fc_threshold <- 1
  
  # 替换NA值，以便绘图
  res$padj[is.na(res$padj)] <- 1
  res$log2FoldChange[is.na(res$log2FoldChange)] <- 0
  # 统计上调和下调的基因个数
  upregulated <- sum(res$padj < significance_level & res$log2FoldChange > log2fc_threshold, na.rm = TRUE)
  downregulated <- sum(res$padj < significance_level & res$log2FoldChange < -log2fc_threshold, na.rm = TRUE)
  # 开始绘图
  plot(res$log2FoldChange, -log10(res$padj),
       pch = 20, cex = 0.6,
       main = title,
       xlab = "Log2(FC)",
       ylab = "-Log10P-value",
       # xlim = c(-15, 15), ylim = c(0, 25))
       xlim = c(-15, 15), ylim = c(0, 40))
       # xaxp = c(-10, 10, 5), yaxp = c(0, 100, 5))
  
  # 根据阈值上色
  col_points <- ifelse(res$padj < significance_level & res$log2FoldChange < -log2fc_threshold, '#0577e2',
                       ifelse(res$padj < significance_level & res$log2FoldChange > log2fc_threshold, '#ea0f54', 'black'))
  
  points(res$log2FoldChange, -log10(res$padj), pch = 20, col = col_points)
  
  # 添加垂直虚线
  abline(v = 0, col = "black", lty = 2, lwd = 2)
  # 在图的左上角和右上角添加上调和下调的基因个数
  # text(-12, 20, paste0(downregulated), col = '#0577e2', cex = 1.6)
  # text(12, 20, paste0(upregulated), col = '#ea0f54', cex = 1.6)
  text(-10, 35, paste0(downregulated), col = '#0577e2', cex = 1.6)
  text(10, 35, paste0(upregulated), col = '#ea0f54', cex = 1.6)
}

# plot_volcano(data,plottitle)
volcano_plot <- plot_volcano(data,plottitle)
volcano_plot
# 保存pdf
pdf(file = paste0("./volcano/", plottitle, ".pdf"), width = 7, height = 7)  # 设置画布大小
plot_volcano(data, plottitle)  # 生成图
dev.off()  # 关闭设备，完成保存


#########################只用草鱼基因富集################################################
# #准备基因列表
# genes_up <- rownames(data[data$category == "up", ])
# genes_down <- rownames(data[data$category == "down", ])

########################草鱼基因对到斑马鱼基因富集######################################
# 准备基因列表
gene_up <- rownames(data[data$category == "up", ])
gene_down <- rownames(data[data$category == "down", ])

#替换基因名
#将草鱼基因名替换为斑马鱼基因名
# homolog_gene <- read.csv("./Grasscarp_Zebrafish_Human_Homolog_gene.csv")
# # 定义函数 Trans_gene
# Trans_gene <- function(genelist, homolog_df) {
#   tempdf <- as.data.frame(genelist)
#   colnames(tempdf) <- "X"
#   tempmerge <- merge(tempdf, homolog_df, by.x = "X", by.y = "Grasscarp_gene", all.x = TRUE)
#   tempmerge <- tempmerge %>%
#     filter(Zebrafish_gene != "#N/A")
#   newlist <- tempmerge$Zebrafish_gene
#   mylist <- unique(newlist)
#   list_name <- deparse(substitute(genelist))
#   write.csv(tempmerge,file=paste0("./GOgenes/",list_name,"_",file),row.names = FALSE)
#   return(mylist)
# }
# genes_up <- Trans_gene(gene_up,homolog_gene)
# genes_down <- Trans_gene(gene_down,homolog_gene)
# genes_diff <- c(genes_up,genes_down)

#########################草鱼基因对到斑马鱼基因################################


# GO 富集分析
ego_up <- enrichGO(gene = gene_up, 
                   OrgDb = org.Dr.eg.db, 
                   keyType = "SYMBOL", 
                   ont = "ALL", 
                   pAdjustMethod = "BH", 
                   qvalueCutoff = 0.05, 
                   readable = TRUE)

ego_down <- enrichGO(gene = gene_down, 
                     OrgDb = org.Dr.eg.db, 
                     keyType = "SYMBOL", 
                     ont = "ALL", 
                     pAdjustMethod = "BH", 
                     qvalueCutoff = 0.05, 
                     readable = TRUE)

gene_diff <- c(gene_down,gene_up)
ego_diff <- enrichGO(gene = gene_diff, 
                     OrgDb = org.Dr.eg.db, 
                     keyType = "SYMBOL", 
                     ont = "ALL", 
                     pAdjustMethod = "BH", 
                     qvalueCutoff = 0.05, 
                     readable = TRUE)

#排序保存富集结果
summary(ego_down)
dir.create("./GOgenes", showWarnings = FALSE, recursive = TRUE)
ego_down_df <- as.data.frame(ego_down)
write.csv(ego_down_df,file=paste0("./GOgenes/","GOdown_",file))
ego_up_df <- as.data.frame(ego_up)
write.csv(ego_up_df,file=paste0("./GOgenes/","GOup_",file))
ego_diff_df <- as.data.frame(ego_diff)
write.csv(ego_diff_df,file=paste0("./GOgenes/","GOdiff_",file))
##############绘制基因通路图#########
####上调基因的GO富集结果条形图
bp_up <- barplot(ego_up, showCategory=20, title=paste0("GOup_",gsub(".csv","",file)))
bp_up
ggsave(paste0("./GOgenes/GOup_",gsub(".csv",".pdf",file)),bp_up,height = 15, width = 10,dpi=800)
# 下调基因的GO富集结果条形图
bp_down <- barplot(ego_down, showCategory=20, title=paste0("GOdown_",gsub(".csv","",file)))
bp_down
ggsave(paste0("./GOgenes/GOdown_",gsub(".csv",".pdf",file)),bp_down,height = 15, width = 10,dpi=800)
#所有差异基因富集结果条形图
# 下调基因的GO富集结果条形图
bp_diff <- barplot(ego_diff, showCategory=20, title=paste0("GOdiff_",gsub(".csv","",file)))
bp_diff
ggsave(paste0("./GOgenes/GOdiff_",gsub(".csv",".pdf",file)),bp_diff,height = 15, width = 10,dpi=800)

###### ego_diff ##########
dotplot <- ggplot(ego_diff_df, aes(x = Count, y = reorder(Description, Count), color = p.adjust, size = Count)) +
  geom_point() +  # 使用geom_point()来绘制气泡图
  scale_color_gradient(low = "blue", high = "red") +  # 设置p.adjust的颜色梯度
  scale_size_continuous(range = c(1, 8)) +  # 控制气泡大小范围
  theme_bw() +
  labs(title = paste0("GOdiff",plottitle), x = "Count", y = NULL) +
  theme(axis.text.x = element_text(hjust = 1))  # 让x轴标签倾斜，便于显示
dir.create("./GOgenes", showWarnings = FALSE, recursive = TRUE)
ggsave(dotplot,file=paste0("./GOgenes/","GO_DotDiff_",plottitle,".pdf"),
       height = 15, width=10, dpi=800)

barplot <- ggplot(ego_diff_df, aes(x = Count, y = reorder(Description, Count),fill = p.adjust)) +
  geom_bar(stat = "identity") +  # 绘制条形图
  scale_fill_gradient(low = "blue", high = "red") +  # 设置p.adjust的颜色梯度
  theme_bw() +
  labs(title = paste0("GOdiff",plottitle), y = NULL, x = "Count") +
  theme(axis.text.x = element_text(angle = 45)) 

ggsave(barplot,file=paste0("./GOgenes/","GO_BarDiff_",plottitle,".pdf"),
       height = 15, width=10, dpi=800)

##################################################################

###### ego_up ##########
dotplot <- ggplot(ego_up_df, aes(x = Count, y = reorder(Description, Count), color = p.adjust, size = Count)) +
  geom_point() +  # 使用geom_point()来绘制气泡图
  scale_color_gradient(low = "blue", high = "red") +  # 设置p.adjust的颜色梯度
  scale_size_continuous(range = c(1, 8)) +  # 控制气泡大小范围
  theme_bw() +
  labs(title = paste0("GOup",plottitle), x = "Count", y = NULL) +
  theme(axis.text.x = element_text(hjust = 1))  # 让x轴标签倾斜，便于显示
dir.create("./GOgenes", showWarnings = FALSE, recursive = TRUE)
ggsave(dotplot,file=paste0("./GOgenes/","GO_Dotup_",plottitle,".pdf"),
       height = 15, width=10, dpi=800)

barplot <- ggplot(ego_up_df, aes(x = Count, y = reorder(Description, Count),fill = p.adjust)) +
  geom_bar(stat = "identity") +  # 绘制条形图
  scale_fill_gradient(low = "blue", high = "red") +  # 设置p.adjust的颜色梯度
  theme_bw() +
  labs(title = paste0("GOup",plottitle), y = NULL, x = "Count") +
  theme(axis.text.x = element_text(angle = 45)) 

ggsave(barplot,file=paste0("./GOgenes/","GO_Barup_",plottitle,".pdf"),
       height = 15, width=10, dpi=800)

##################################################################

###### ego_up ##########
dotplot <- ggplot(ego_down_df, aes(x = Count, y = reorder(Description, Count), color = p.adjust, size = Count)) +
  geom_point() +  # 使用geom_point()来绘制气泡图
  scale_color_gradient(low = "blue", high = "red") +  # 设置p.adjust的颜色梯度
  scale_size_continuous(range = c(1, 8)) +  # 控制气泡大小范围
  theme_bw() +
  labs(title = paste0("GOdown_",plottitle), x = "Count", y = NULL) +
  theme(axis.text.x = element_text(hjust = 1))  # 让x轴标签倾斜，便于显示
dir.create("./GOgenes", showWarnings = FALSE, recursive = TRUE)
ggsave(dotplot,file=paste0("./GOgenes/","GO_Dotdown__",plottitle,".pdf"),
       height = 15, width=10, dpi=800)

barplot <- ggplot(ego_down_df, aes(x = Count, y = reorder(Description, Count),fill = p.adjust)) +
  geom_bar(stat = "identity") +  # 绘制条形图
  scale_fill_gradient(low = "blue", high = "red") +  # 设置p.adjust的颜色梯度
  theme_bw() +
  labs(title = paste0("GOdown_",plottitle), y = NULL, x = "Count") +
  theme(axis.text.x = element_text(angle = 45)) 

ggsave(barplot,file=paste0("./GOgenes/","GO_Bardown__",plottitle,".pdf"),
       height = 15, width=10, dpi=800)

##################################################################

#############################KEGG富集#################
# # 将基因符号转换为 ENTREZ ID
ups <- bitr(gene_up, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Dr.eg.db)
downs <- bitr(gene_down, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Dr.eg.db)
diffs <- rbind(ups,downs)


# KEGG富集分析（上调基因），直接使用基因符号
kegg_up <- enrichKEGG(gene = ups$ENTREZID,
                      organism = 'dre',    # 斑马鱼为 "dre",草鱼为cide
                      keyType = "ncbi-geneid",  # 直接使用基因符号
                      pvalueCutoff = 0.05)

kegg_down <- enrichKEGG(gene = downs$ENTREZID,
                        organism = 'dre',
                        keyType = "ncbi-geneid",
                        pvalueCutoff = 0.05)

# KEGG富集分析（所有差异基因），直接使用基因符号
kegg_diffs <- enrichKEGG(gene = diffs$ENTREZID,
                        organism = 'dre',
                        keyType = "ncbi-geneid",
                        pvalueCutoff = 0.05)



# 保存结果为 csv 文件
dir.create("./kegg/",showWarnings = FALSE, recursive = TRUE)

# 导出 KEGG 上调基因富集结果
write.csv(kegg_up@result, file = "./kegg/kegg_up_all_results.csv", row.names = FALSE)
write.csv(kegg_down@result, file = "./kegg/kegg_down_all_results.csv", row.names = FALSE)
write.csv(kegg_diffs@result, file = "./kegg/kegg_diffs_all_results.csv", row.names = FALSE)


write.csv(kegg_up, paste0("./kegg/kegg_up_",file), row.names = FALSE)
write.csv(kegg_up, paste0("./kegg/kegg_down_",file), row.names = FALSE)
write.csv(kegg_diffs, paste0("./kegg/kegg_diffs_",file), row.names = FALSE)

# 绘制 KEGG 富集分析结果图
# 下调基因 KEGG 富集气泡图
dotdiffs <- dotplot(kegg_diffs, showCategory = 20) + ggtitle(paste0("KEGGdiffs_",gsub(".csv","",file)))
ggsave(paste0("./kegg/diffs_",gsub(".csv",".pdf",file)),dotdiffs)
# ggsave(paste0("./kegg/diffs_",gsub(".csv",".pdf",file)),dotdiffs,height = 10, width = 10,dpi=800)
# 上调基因 KEGG 富集气泡图
dotup <- dotplot(kegg_up, showCategory = 20) + ggtitle(paste0("KEGGup_",gsub(".csv","",file)))
ggsave(paste0("./kegg/up_",gsub(".csv",".pdf",file)),dotup)
# 下调基因 KEGG 富集气泡图
dotdown <- dotplot(kegg_down, showCategory = 20) + ggtitle(paste0("KEGGdown_",gsub(".csv","",file)))
ggsave(paste0("./kegg/down_",gsub(".csv",".pdf",file)),dotdown)


