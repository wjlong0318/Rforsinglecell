library(Seurat)
library(cowplot)
library(Matrix)
library(dplyr)
library(ggplot2)
library(reshape2)
setwd("singlecell/P-MATRIX")

#setwd("H:\\shanghai\\singlecell\\P-MATRIX")


print("start...\n")

nCoV.integrated<- readRDS("IntegrateSample_T10.rds")

sample_info = as.data.frame(colnames(nCoV.integrated))
colnames(sample_info) = c('ID')
rownames(sample_info) = sample_info$ID
sample_info$sample = nCoV.integrated@meta.data$orig.ident
samples=read.csv("../meta.csv",header=T,sep=",")
sample_info = dplyr::left_join(sample_info,samples)
rownames(sample_info) = sample_info$ID
nCoV.integrated = AddMetaData(object = nCoV.integrated, metadata = sample_info)

DefaultAssay(nCoV.integrated) <- "RNA"
nCoV.integrated[['percent.mito']] <- PercentageFeatureSet(nCoV.integrated, pattern = "^MT-")
nCoV.integrated <- NormalizeData(object = nCoV.integrated, normalization.method = "LogNormalize", scale.factor = 1e4)
nCoV.integrated <- FindVariableFeatures(object = nCoV.integrated, selection.method = "vst", nfeatures = 2000,verbose = FALSE)
nCoV.integrated <- ScaleData(nCoV.integrated, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito"))
#VlnPlot(object = nCoV.integrated, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
#FeatureScatter(object = nCoV.integrated, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
nCoV.integrated <- ScaleData(nCoV.integrated, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito"))
nCoV.integrated <- RunPCA(nCoV.integrated, verbose = FALSE,npcs = 100)
nCoV.integrated <- ProjectDim(object = nCoV.integrated)
#ElbowPlot(object = nCoV.integrated,ndims = 100)

nCoV.integrated <- FindNeighbors(object = nCoV.integrated, dims = 1:50)
nCoV.integrated <- FindClusters(object = nCoV.integrated, resolution = 1.2)
nCoV.integrated <- RunTSNE(object = nCoV.integrated, dims = 1:50)

nCoV.integrated  <- RunUMAP(nCoV.integrated  , dims = 1:50)
saveRDS(nCoV.integrated, file="IntegrateSample_tsne_umap_T10.rds")
# pdf("umap_dimplot.pdf")
# DimPlot(nCoV.integrated , reduction = "umap", label = TRUE)
# dev.off()
# 
# pdf("tsne_dimplot.pdf")
# DimPlot(nCoV.integrated , reduction = "tsne", label = TRUE)
# dev.off()

pdf("tchain_tsne_dimplot.pdf")
DimPlot(nCoV.integrated , reduction = "tsne", group.by = "t_chain")
dev.off()

pdf("bchain_tsne_dimplot.pdf")
DimPlot(nCoV.integrated , reduction = "tsne", group.by = "b_chain")
dev.off()

# pdf("tchain_dimplot.pdf")
# DimPlot(nCoV.integrated ,group.by = "t_chain")
# dev.off()
# 
# pdf("bchain_dimplot.pdf")
# DimPlot(nCoV.integrated ,group.by = "b_chain")
# dev.off()
#T?????????marker???????????????:
pdf("tmarker_fplot.pdf")
t_cell_markers <- c("CD3D","CD3E")
FeaturePlot(nCoV.integrated , features = t_cell_markers)
dev.off()



# DimPlot(object = nCoV.integrated, reduction = 'tsne',label = TRUE)
# 
# markers = c("CD4",
#             "CD8A",
#             "FCGR3A",
#             "CD14",
#             "MS4A1",
#             "CD3D",
#             "KLRF1",
#             "KLRD1",
#             "FCER1A",
#             "PF4",
#             "CD68",
#             "IGHG4"
# )
# 

# pdf("DotPlot.pdf")
# pp = DotPlot(nCoV.integrated, features = rev(markers),cols = c('white','#F8766D'),dot.scale =5) + RotatedAxis()
# pp = pp + 
#   theme(axis.text.x = element_text(size = 8),
#         axis.text.y = element_text(size = 8)) + 
#   labs(x='',y='') + 
#   guides(color = guide_colorbar(title = 'Scale expression'),
#          size = guide_legend(title = 'Percent expressed')) + 
#   theme(axis.line = element_line(size = 0.6))
# pp
# dev.off()

