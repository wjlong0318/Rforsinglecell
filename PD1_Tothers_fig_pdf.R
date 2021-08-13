library(Seurat)
library(cowplot)
library(Matrix)
library(dplyr)
library(ggplot2)
library(reshape2)
setwd("singlecell/P-MATRIX")

print("start...\n")

# 
# nCoV.integrated<- readRDS("IntegrateSample_Tothers.rds")
# 
# 
# sample_info = as.data.frame(colnames(nCoV.integrated))
# colnames(sample_info) = c('ID')
# rownames(sample_info) = sample_info$ID
# sample_info$sample = nCoV.integrated@meta.data$orig.ident
# samples=read.csv("../meta.csv",header=T,sep=",")
# sample_info = dplyr::left_join(sample_info,samples)
# rownames(sample_info) = sample_info$ID
# nCoV.integrated = AddMetaData(object = nCoV.integrated, metadata = sample_info)
# 
# DefaultAssay(nCoV.integrated) <- "RNA"
# nCoV.integrated[['percent.mito']] <- PercentageFeatureSet(nCoV.integrated, pattern = "^MT-")
# nCoV.integrated <- NormalizeData(object = nCoV.integrated, normalization.method = "LogNormalize", scale.factor = 1e4)
# nCoV.integrated <- FindVariableFeatures(object = nCoV.integrated, selection.method = "vst", nfeatures = 2000,verbose = FALSE)
# nCoV.integrated <- ScaleData(nCoV.integrated, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito"))
# nCoV.integrated <- RunPCA(nCoV.integrated, verbose = FALSE,npcs = 100)
# nCoV.integrated <- ProjectDim(object = nCoV.integrated)
# nCoV.integrated <- FindNeighbors(object = nCoV.integrated, dims = 1:50)
# nCoV.integrated <- FindClusters(object = nCoV.integrated, resolution = 1.2)
# nCoV.integrated <- RunTSNE(object = nCoV.integrated, dims = 1:50)
# nCoV.integrated  <- RunUMAP(nCoV.integrated  , dims = 1:50)
# saveRDS(nCoV.integrated, file="IntegrateSample_tsne_umap_Tothers.rds")
# print("SaveData IntegrateSample_tsne_umap_Tothers.rds complete...\n")
# pdf("umap_Tothers_dimplot.pdf")
# DimPlot(nCoV.integrated , reduction = "umap", label = TRUE)
# dev.off()
# 
# pdf("tsne_Tothers_dimplot.pdf")
# DimPlot(nCoV.integrated , reduction = "tsne", label = TRUE)
# dev.off()
# 
# pdf("tchain_tsne_Tothers_dimplot.pdf")
# DimPlot(nCoV.integrated , reduction = "tsne", group.by = "t_chain")
# dev.off()
# 
# pdf("bchain_tsne_Tothers_dimplot.pdf")
# DimPlot(nCoV.integrated , reduction = "tsne", group.by = "b_chain")
# dev.off()
# 
# pdf("tchain_umap_Tothers_dimplot.pdf")
# DimPlot(nCoV.integrated , reduction = "umap", group.by = "t_chain")
# dev.off()
# 
# pdf("bchain_umap_Tothers_dimplot.pdf")
# DimPlot(nCoV.integrated , reduction = "umap", group.by = "b_chain")
# dev.off()
# 
# pdf("umap_Response_Tothers_dimplot.pdf")
# DimPlot(nCoV.integrated , reduction = "umap", label = TRUE,split.by = 'Response')
# dev.off()
# 
# pdf("tsne_Response_Tothers_dimplot.pdf")
# DimPlot(nCoV.integrated , reduction = "tsne", label = TRUE,split.by = 'Response')
# dev.off()
# 
# pdf("tchain_tsne_Response_Tothers_dimplot.pdf")
# DimPlot(nCoV.integrated , reduction = "tsne", group.by = "t_chain",split.by = 'Response')
# dev.off()
# 
# pdf("bchain_tsne_Response_Tothers_dimplot.pdf")
# DimPlot(nCoV.integrated , reduction = "tsne", group.by = "b_chain",split.by = 'Response')
# dev.off()
# 
# pdf("tchain_umap_Response_Tothers_dimplot.pdf")
# DimPlot(nCoV.integrated , reduction = "umap", group.by = "t_chain",split.by = 'Response')
# dev.off()
# 
# pdf("bchain_umap_Response_Tothers_dimplot.pdf")
# DimPlot(nCoV.integrated , reduction = "umap", group.by = "b_chain",split.by = 'Response')
# dev.off()

print("start singleR...")
library(celldex)
library(SingleR)
pbmc=readRDS("IntegrateSample_tsne_umap_Tothers.rds")

pbmc_for_SingleR <- GetAssayData(pbmc, slot="data")
print("read pbmc_for_SingleR...")

hpca.se <- celldex::HumanPrimaryCellAtlasData()
Blue.se=celldex::BlueprintEncodeData() 
Immune.se=celldex::DatabaseImmuneCellExpressionData()
print("read hpca.se,blue.se,Immune.se...")


clusters=pbmc@meta.data$seurat_clusters
print("read clusters...")
head(clusters)

pred.hesc <- SingleR(test = pbmc_for_SingleR, ref = list(BP=Blue.se, HPCA=hpca.se,Im=Immune.se), labels =list(Blue.se$label.fine, hpca.se$label.fine,Immune.se$label.fine),
                     method = "cluster", clusters = clusters, 
                     assay.type.test = "logcounts", assay.type.ref = "logcounts")
print("[result]pred.hesc...")
pred.hesc$labels

celltype = data.frame(ClusterID=rownames(pred.hesc), celltype=pred.hesc$labels, stringsAsFactors = F) 
print("[result]celltype...")
celltype
pbmc@meta.data$singleR=celltype[match(clusters,celltype$ClusterID),'celltype']
#write.csv(pbmc@meta.data,file="pbmc_Tothers_meta.data.ref3_fine.csv")


pred.main <- SingleR(test = pbmc_for_SingleR, ref = list(BP=Blue.se, HPCA=hpca.se,Im=Immune.se), labels =list(Blue.se$label.main, hpca.se$label.main,Immune.se$label.main),
                     method = "cluster", clusters = clusters, 
                     assay.type.test = "logcounts", assay.type.ref = "logcounts")
print("[result]pred.main...")
pred.main$labels

celltype.main = data.frame(ClusterID=rownames(pred.main), celltype.main=pred.main$labels, stringsAsFactors = F) 
print("[result]celltype.main...")
celltype.main
pbmc@meta.data$singleR.main=celltype.main[match(clusters,celltype.main$ClusterID),'celltype.main']
print("[result]write metadata fine and main...")

tem=pbmc@meta.data$singleR
tem[which(tem=="Monocytes")] <-"Monocyte"
pbmc@meta.data$singleR2=tem

tem2=pbmc@meta.data$singleR.main
tem2[which(tem2=="Monocytes")] <-"Monocyte"
pbmc@meta.data$singleR.main2=tem2

write.csv(pbmc@meta.data,file="pbmc_Tothers_meta.data.ref3_main_fine.csv")
saveRDS(pbmc, file="sample.combined_Tothers_singleR_ref3_main_fine.rds")

print("start write pdf..." )
pdf("umap_Tothers_singleR2_ref3_dimplot.pdf")
DimPlot(pbmc, reduction = "umap", group.by = "singleR2")
dev.off()
pdf("tsne_Tothers_singleR2_ref3_dimplot.pdf")
DimPlot(pbmc, reduction = "tsne", group.by = "singleR2")
dev.off()

pdf("umap_Tothers_Response_singleR2_ref3_dimplot.pdf")
DimPlot(pbmc, reduction = "umap", group.by = "singleR2",split.by = 'Response')
dev.off()
pdf("tsne_Tothers_Response_singleR2_ref3_dimplot.pdf")
DimPlot(pbmc, reduction = "tsne", group.by = "singleR2",split.by = 'Response')
dev.off()


print("start write pdf..." )

pdf("umap_Tothers_singleR2_ref3_dimplot_fine.pdf")
DimPlot(pbmc, reduction = "umap", group.by = "singleR.main2")
dev.off()
pdf("tsne_Tothers_singleR2_ref3_dimplot_fine.pdf")
DimPlot(pbmc, reduction = "tsne", group.by = "singleR.main2")
dev.off()

pdf("umap_Response_Tothers_singleR2_ref3_dimplot_fine.pdf")
DimPlot(pbmc, reduction = "umap", group.by = "singleR.main2",split.by = 'Response')
dev.off()
pdf("tsne_Response_Tothers_singleR2_ref3_dimplot_fine.pdf")
DimPlot(pbmc, reduction = "tsne", group.by = "singleR.main2",split.by = 'Response')
dev.off()





