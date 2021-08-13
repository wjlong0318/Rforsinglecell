library(slingshot)
library(SingleCellExperiment)
library(mclust, quietly = TRUE)
library(Seurat)
library(scater)
library(uwot)
library(tradeSeq)

#devtools::install_github('satijalab/seurat-data')

# BiocManager::install(c('scater', 'scran', 'uwot'))
# BiocManager::install("slingshot")
# BiocManager::install("tradeSeq")

setwd("singlecell/P-MATRIX")
print("load library complete...")
print(Sys.time())
setwd("singlecell/P-MATRIX")


print("load RDS complete...")
pbmc<- readRDS("sample.combined_T10_others_singleR_ref3_main_fine.rds")
scdata1=subset(pbmc,Treatment!=5)
scdata2=subset(scdata1,Response!="Na")

#scdata=subset(scdata2,singleR.main2=="Bcell")
scdata=subset(scdata2,b_is_cell=="true")

pbmc_sce=as.SingleCellExperiment(scdata, assay = 'RNA')
# pbmc_sce=scater::addPerCellQC(pbmc_sce)
# pbmc_sce=scater::addPerFeatureQC(pbmc_sce)
print("col and row...")
colData(pbmc_sce)
rowData(pbmc_sce)
pbmc_sce

pbmc_sce=scran::computeSumFactors(pbmc_sce)
#pbmc_sce<- scater::normalize(pbmc_sce)


geneFilter <- apply(assays(pbmc_sce)$counts,1,function(x){
  sum(x >= 3) >= 10
})
pbmc_sce<- pbmc_sce[geneFilter, ]

FQnorm <- function(counts){
  rk <- apply(counts,2,rank,ties.method='min')
  counts.sort <- apply(counts,2,sort)
  refdist <- apply(counts.sort,1,median)
  norm <- apply(rk,2,function(r){ refdist[r] })
  rownames(norm) <- rownames(counts)
  return(norm)
}

assays(pbmc_sce)$norm <- FQnorm(assays(pbmc_sce)$counts)

saveRDS(mycds, file="T10_Tcell_slingshot_sce.rds")

pca <- prcomp(t(log1p(assays(pbmc_sce)$norm)), scale. = FALSE)
rd1 <- pca$x[,1:2]
head(rd1)
#        PC1       PC2
#c1 -15.55731 -6.218455
#c2 -14.48979 -7.060021
#c3 -16.10326 -7.630399
#c4 -15.59147 -7.443861
#c5 -15.86607 -6.180724
#c6 -15.78163 -6.733366

plot(rd1, col = rgb(0,0,0,.5), pch=16, asp = 1)

library(uwot)
## Loading required package: Matrix
## 
## Attaching package: 'Matrix'
## The following object is masked from 'package:S4Vectors':
## 
##     expand
# 进行UMAP降维
rd2 <- umap(t(log1p(assays(pbmc_sce)$norm)))
colnames(rd2) <- c('UMAP1', 'UMAP2')
head(rd2)
#        UMAP1     UMAP2
#[1,] 6.357223 -4.358757
#[2,] 6.442034 -4.216842
#[3,] 6.478745 -4.421883
#[4,] 6.712623 -4.487934
#[5,] 6.497382 -4.276127
#[6,] 6.302136 -4.531643

plot(rd2, col = rgb(0,0,0,.5), pch=16, asp = 1)

reducedDims(pbmc_sce) <- pbmc_scepleList(PCA = rd1, UMAP = rd2)
saveRDS(mycds, file="T10_Tcell_slingshot_sce_PCA_UMAP.rds")
pbmc_sce

library(mclust, quietly = TRUE)
## Package 'mclust' version 5.4.6
## Type 'citation("mclust")' for citing this R package in publications.
## 
## Attaching package: 'mclust'
## The following object is masked from 'package:mgcv':
## 
##     mvn
# 使用Mclust函数进行细胞聚类
cl1 <- Mclust(rd1)$classification
head(cl1)
#c1 c2 c3 c4 c5 c6 
# 1  1  1  1  1  1 
colData(pbmc_sce)$GMM <- cl1

saveRDS(mycds, file="T10_Tcell_slingshot_sce_PCA_UMAP_ Mclust.rds")
pbmc_sce
#class: SingleCellExperiment 
#dim: 734 300 
#metadata(0):
#assays(2): counts norm
#rownames(734): G2 G3 ... G749 G750
#rowData names(0):
#colnames(300): c1 c2 ... c299 c300
#colData names(1): GMM
#reducedDimNames(2): PCA UMAP
#spikeNames(0):

# 可视化聚类分群的结果
library(RColorBrewer)
# plot(rd1, col = brewer.pal(9,"Set1")[cl1], pch=16, asp = 1)
# 
# cl2 <- kmeans(rd1, centers = 4)$cluster
# head(cl2)
# #c1 c2 c3 c4 c5 c6 
# # 2  2  2  2  2  2 
# colData(pbmc_sce)$kmeans <- cl2
# 
# # 可视化聚类分群的结果
# #plot(rd1, col = brewer.pal(9,"Set1")[cl2], pch=16, asp = 1)

pbmc_sce <- slingshot(pbmc_sce, clusterLabels = 'GMM', reducedDim = 'PCA')
## Using full covariance matrix
pbmc_sce

summary(pbmc_sce$slingPseudotime_1)
##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
##   0.000   8.633  21.118  21.415  34.367  43.186
pdf("T10_Tcell_plot_slingPseudotime.pdf")
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(pbmc_sce$slingPseudotime_1, breaks=100)]
plot(reducedDims(pbmc_sce)$PCA, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(pbmc_sce), lwd=2, col='black')
dev.off()

pdf("T10_Tcell_plot_sling_GMM.pdf")
plot(reducedDims(pbmc_sce)$PCA, col = brewer.pal(9,'Set1')[pbmc_sce$GMM], pch=16, asp = 1)
lines(SlingshotDataSet(pbmc_sce), lwd=2, type = 'lineages', col = 'black')
dev.off()

# BiocManager::install("tradeSeq")
# library(tradeSeq)
# 
# # fit negative binomial GAM
# pbmc_sce <- fitGAM(pbmc_sce)
# 
# # test for dynamic expression
# ATres <- associationTest(pbmc_sce)
# topgenes <- rownames(ATres[order(ATres$pvalue), ])[1:250]
# pst.ord <- order(pbmc_sce$slingPseudotime_1, na.last = NA)
# heatdata <- assays(pbmc_sce)$counts[topgenes, pst.ord]
# heatclus <- pbmc_sce$GMM[pst.ord]
# 
# heatmap(log1p(heatdata), Colv = NA,
#         ColSideColors = brewer.pal(9,"Set1")[heatclus])
# 
# ---------
# ## Using full covariance matrix
# 
# lin1
# ## class: SlingshotDataSet 
# ## 
# ##  Samples Dimensions
# ##      140          2
# ## 
# ## lineages: 2 
# ## Lineage1: 1  2  3  5  
# ## Lineage2: 1  2  3  4  
# ## 
# ## curves: 0
# plot(rd, col = brewer.pal(9,"Set1")[cl], asp = 1, pch = 16)
# lines(lin1, lwd = 3, col = 'black')
