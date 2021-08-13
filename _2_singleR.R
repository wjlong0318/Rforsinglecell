# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("celldex")
library(Seurat)
library(celldex)
library(SingleR)
print("load library complete...")

setwd("singlecell/P-MATRIX")
pbmc<- readRDS("IntegrateSample_tsne_umap_T10.rds")
print("load data complete...")

pbmc_for_SingleR <- GetAssayData(pbmc, slot="data")
print("read pbmc_for_SingleR...")
head(pbmc_for_SingleR)

hpca.se <- celldex::HumanPrimaryCellAtlasData()
Blue.se=celldex::BlueprintEncodeData() 
Immune.se=celldex::DatabaseImmuneCellExpressionData()
print("read hpca.se...")
#head(hpca.se)

clusters=pbmc@meta.data$seurat_clusters
print("read clusters...")
head(clusters)

pred.hesc <- SingleR(test = pbmc_for_SingleR, ref = list(BP=Blue.se, HPCA=hpca.se,Im=Immune.se), labels =list(Blue.se$label.main, hpca.se$label.main,Immune.se$label.main),
                     method = "cluster", clusters = clusters, 
                     assay.type.test = "logcounts", assay.type.ref = "logcounts")
print("[result]pred.hesc...")
pred.hesc$labels
#table(pred.hesc$labels,clusters)

celltype = data.frame(ClusterID=rownames(pred.hesc), celltype=pred.hesc$labels, stringsAsFactors = F) 
print("[result]celltype...")
celltype
pbmc@meta.data$singleR=celltype[match(clusters,celltype$ClusterID),'celltype']
write.csv(pbmc@meta.data,file="pbmc_meta.data.ref3.csv")

# dice=celldex::DatabaseImmuneCellExpressionData()
# 
# pred.dice <- SingleR(test = pbmc_for_SingleR, ref = dice, labels = dice$label.main,
#                      method = "cluster", clusters = clusters, 
#                      assay.type.test = "logcounts", assay.type.ref = "logcounts")
# print("[result]pred.dice...")
# pred.dice$labels
# 
# celltype_dice = data.frame(ClusterID=rownames(pred.dice), celltype_dice=pred.dice$labels, stringsAsFactors = F) 
# print("[result]dice celltype...")
# celltype_dice
# pbmc@meta.data$singleR_dice=celltype_dice[match(clusters,celltype_dice$ClusterID),'celltype_dice']



#print("[result]singleR...")
#pbmc@meta.data$singleR
print("start write pdf..." )
saveRDS(pbmc, file="sample.combined_T10_singleR_ref3.rds")
pdf("umap_singleR_ref3_dimplot.pdf")
DimPlot(pbmc, reduction = "umap", group.by = "singleR")
dev.off()
pdf("tsne_singleR_ref3_dimplot.pdf")
DimPlot(pbmc, reduction = "tsne", group.by = "singleR")
dev.off()

pdf("umap_Response_singleR_ref3_dimplot.pdf")
DimPlot(pbmc, reduction = "umap", group.by = "singleR",split.by = 'Response')
dev.off()
pdf("tsne_Response_singleR_ref3_dimplot.pdf")
DimPlot(pbmc, reduction = "tsne", group.by = "singleR",split.by = 'Response')
dev.off()


