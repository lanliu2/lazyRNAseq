# lazyRNAseq

## 简介

这是一个用于RNA-seq下游分析的脚本，暂未封装。输入数据是count文件，并且要求数据读入后将基因名（symbol）所在的列设置为`row.names`

然后进一步对数据进行过滤（这里只设置了过滤要求是基因的count数要大于1，并且至少要在2个样本中表达）

设置样本信息后（样本名、具体的分组）后就进行pca、pcoa、聚类热图、差异基因表达分析，进一步则绘制差异基因火山图，GO通路富集， Kegg等。

## pca分组结果
![](./figures/pca.png)
## pca聚类热图结果
![](./figures/heatmap.png)
## pcoa结果
### pcoa结果
![](./figures/pcoa.png)
### 90%置信区间
![](./figures/pcoa90IC.png)
### 95%置信区间
![](./figures/pcoa95IC.png)
## 差异基因火山图
![](./figures/volcano.png)
