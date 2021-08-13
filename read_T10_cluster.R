library(Seurat)
library(cowplot)
library(Matrix)
library(dplyr)
library(ggplot2)
library(reshape2)
setwd("singlecell/P-MATRIX")


print("start...\n")

pbmc<- readRDS("IntegrateSample_T10.rds")

sample_info = as.data.frame(colnames(pbmc))
colnames(sample_info) = c('ID')
rownames(sample_info) = sample_info$ID
sample_info$sample = pbmc@meta.data$orig.ident
samples=read.csv("../meta.csv",header=T,sep=",")
sample_info = dplyr::left_join(sample_info,samples)
rownames(sample_info) = sample_info$ID
pbmc = AddMetaData(object = pbmc, metadata = sample_info)

print("pbmc assays:")
names(pbmc@assays)
print("pbmc active assays:")
pbmc@active.assay

pbmc.int=pbmc

print("pbmc runing...")
DefaultAssay(pbmc) <- "RNA"
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
Atlas.integrated <- CellCycleScoring(Atlas.integrated, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Atlas.integrated <- ScaleData(Atlas.integrated, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(Atlas.integrated))
Atlas.integrated <- FindVariableFeatures(object = Atlas.integrated, mean.function = ExpMean, dispersion.function = LogVMR)

pbmc[['percent.mito']] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 1e4)
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 2000,verbose = FALSE)

pbmc <- ScaleData(pbmc, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito"))
pbmc <- RunPCA(pbmc, verbose = FALSE,npcs = 100)
pbmc <- FindNeighbors(object = pbmc, dims = 1:50)
pbmc <- FindClusters(object = pbmc, resolution = 1)
pbmc <- RunTSNE(object = pbmc, dims = 1:50)
pbmc  <- RunUMAP(pbmc  , dims = 1:50)
saveRDS(pbmc, file="pbmc_RNA_res1_tsne_umap_T10.rds")

print("pbmc.int runing...")
DefaultAssay(pbmc.int) <- "integrated"

print("pbmc.int active assays:")
pbmc.int@active.assay
names(pbmc.int@meta.data)
#pbmc.int[['percent.mito']] <- PercentageFeatureSet(pbmc.int, pattern = "^MT-")
pbmc.int <- NormalizeData(object = pbmc.int, normalization.method = "LogNormalize", scale.factor = 1e4)
pbmc.int <- FindVariableFeatures(object = pbmc.int, selection.method = "vst", nfeatures = 2000,verbose = FALSE)
pbmc.int <- ScaleData(pbmc.int, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito"))
pbmc.int <- RunPCA(pbmc.int, verbose = FALSE,npcs = 100)
pbmc.int <- FindNeighbors(object = pbmc.int, dims = 1:50)
pbmc.int <- FindClusters(object = pbmc.int, resolution = 1)
pbmc.int <- RunTSNE(object = pbmc.int, dims = 1:50)
pbmc.int  <- RunUMAP(pbmc.int  , dims = 1:50)
saveRDS(pbmc.int, file="pbmc_int_res1_tsne_umap_T10.rds")

pbmc<- readRDS("pbmc_RNA_res1_tsne_umap_T10.rds")
DefaultAssay(pbmc) <- "RNA"

print("write figures...")
plot_grid(ncol = 3,
          DimPlot(pbmc, reduction = "pca", group.by = "orig.ident")+NoAxes()+ggtitle("PCA raw_data"),
          DimPlot(pbmc, reduction = "tsne", group.by = "orig.ident")+NoAxes()+ggtitle("tSNE raw_data"),
          DimPlot(pbmc, reduction = "umap", group.by = "orig.ident")+NoAxes()+ggtitle("UMAP raw_data"),
          
          DimPlot(pbmc.int, reduction = "pca", group.by = "orig.ident")+NoAxes()+ggtitle("PCA integrated"),
          DimPlot(pbmc.int, reduction = "tsne", group.by = "orig.ident")+NoAxes()+ggtitle("tSNE integrated"),
          DimPlot(pbmc.int, reduction = "umap", group.by = "orig.ident")+NoAxes()+ggtitle("UMAP integrated")
)
