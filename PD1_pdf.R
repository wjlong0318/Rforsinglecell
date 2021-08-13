library(Seurat)
library(cowplot)
library(Matrix)
library(dplyr)
library(ggplot2)
library(reshape2)
setwd("singlecell/P-MATRIX")

print("start...\n")

prepdata=function(data,label,group){
  print("load function predata...\n")
  # data=M164
  # label="M164"
  # group="PR"
  Mdata<-CreateSeuratObject(counts = data, project =label,min.cells = 3, min.features = 200)
  Mdata[["percent.mt"]]<-PercentageFeatureSet(Mdata,pattern = "^MT")
  #  p164<-VlnPlot(M164, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
  Mdata.filter <- subset(Mdata, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10 & nCount_RNA >1000)
  Mdata.norm <- NormalizeData(Mdata.filter)
  Mdata.genes <- FindVariableFeatures(Mdata.norm, selection.method = "vst", nfeatures = 2000)
  Mdata.genes
}

print("read files...\n")

M164<-Read10X("P-164")
M165<-Read10X("P-165")
M174<-Read10X("P-174")
M175<-Read10X("P-175")
M176<-Read10X("P-176")
M177<-Read10X("P-177")
M178<-Read10X("P-178")
M179<-Read10X("P-179")
M180<-Read10X("P-180")
M181<-Read10X("P-181")
M183<-Read10X("P-183")
M184<-Read10X("P-184")
M185<-Read10X("P-185")
M186<-Read10X("P-186")
M187<-Read10X("P-187")
M189<-Read10X("P-189")
M191<-Read10X("P-191")
M192<-Read10X("P-192")
M193<-Read10X("P-193")
M197<-Read10X("P-197")
M198<-Read10X("P-198")
M199<-Read10X("P-199")
M202<-Read10X("P-202")
M164.gene=prepdata(M164,"M164","PR")
M165.gene=prepdata(M165,"M165","PR")
M174.gene=prepdata(M174,"M174","PD")
M175.gene=prepdata(M175,"M175","PD")
M176.gene=prepdata(M176,"M176","PR")
M177.gene=prepdata(M177,"M177","SD")
M178.gene=prepdata(M178,"M178","PD")
M179.gene=prepdata(M179,"M179","SD")
M180.gene=prepdata(M180,"M180","PD")
M181.gene=prepdata(M181,"M181","PR")
M183.gene=prepdata(M183,"M183","SD")
M184.gene=prepdata(M184,"M184","SD")
M185.gene=prepdata(M185,"M185","SD")
M186.gene=prepdata(M186,"M186","SD")
M187.gene=prepdata(M187,"M187","SD")
M189.gene=prepdata(M189,"M189","PR")
M191.gene=prepdata(M191,"M191","PR")
M192.gene=prepdata(M192,"M192","PD")
M193.gene=prepdata(M193,"M193","PD")
M197.gene=prepdata(M197,"M197","PD")
M198.gene=prepdata(M198,"M198","PR")
M199.gene=prepdata(M199,"M199","PR")
M202.gene=prepdata(M202,"M202","PD")


sample.anchors <- FindIntegrationAnchors(object.list = list(
  M164.gene,
  M165.gene,
  M174.gene,
  M175.gene,
  M176.gene,
  M177.gene,
  M178.gene,
  M179.gene,
  M180.gene,
  M181.gene,
  M183.gene,
  M184.gene,
  M185.gene,
  M186.gene,
  M187.gene,
  M189.gene,
  M191.gene,
  M192.gene,
  M193.gene,
  M197.gene,
  M198.gene,
  M199.gene,
  M202.gene), dims = 1:50)
sample.combined <- IntegrateData(anchorset = sample.anchors, dims = 1:50)

saveRDS(sample.combined, file="sample.combined164_202.rds")
# nCoV.integrated<- readRDS("sample.combined164_202.rds")
nCoV.integrated=sample.combined
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
DimPlot(object = nCoV.integrated, reduction = 'tsne',label = TRUE)

markers = c("CD4",
            "CD8A",
            "FCGR3A",
            "CD14",
            "MS4A1",
            "CD3D",
            "KLRF1",
            "KLRD1",
            "FCER1A",
            "PF4",
            "CD68",
            "IGHG4"
)
pp_temp = FeaturePlot(object = nCoV.integrated, features = markers,cols = c("lightgrey","#ff0000"),combine = FALSE)
plots <- lapply(X = pp_temp, FUN = function(p) p + 
                  theme(axis.title = element_text(size = 8),
                        axis.text = element_text(size = 8),
                        plot.title = element_text(family = 'sans',face='italic',size=8),
                        legend.text = element_text(size = 8),
                        legend.key.height = unit(1.8,"line"),
                        legend.key.width = unit(1.2,"line") ))
pdf("combineplot.pdf")
CombinePlots(plots = plots,ncol = 4,legend = 'right')
dev.off()

pdf("DimCluster1.pdf")
DimPlot(object = nCoV.integrated, reduction = 'tsne',label = TRUE)
dev.off()

pdf("DotPlot.pdf")
pp = DotPlot(nCoV.integrated, features = rev(markers),cols = c('white','#F8766D'),dot.scale =5) + RotatedAxis()
pp = pp + 
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8)) + 
  labs(x='',y='') + 
  guides(color = guide_colorbar(title = 'Scale expression'),
         size = guide_legend(title = 'Percent expressed')) + 
  theme(axis.line = element_line(size = 0.6))
dev.off()


pdf("DimCluster_Response.pdf")
DimPlot(object = nCoV.integrated, reduction = 'tsne',label = FALSE, group.by = 'Response')
dev.off()
pdf("DimCluster_Treatment.pdf")
DimPlot(object = nCoV.integrated, reduction = 'tsne',label = FALSE, group.by = 'Treatment')
dev.off()
pdf("DimCluster_Gender.pdf")
DimPlot(object = nCoV.integrated, reduction = 'tsne',label = FALSE, group.by = 'Gender')
dev.off()
pdf("DimSplit_Response.pdf")
DimPlot(object = nCoV.integrated, reduction = 'tsne',label = FALSE, split.by = 'Response')
dev.off()
pdf("DimSplit_Treatment.pdf")
DimPlot(object = nCoV.integrated, reduction = 'tsne',label = FALSE, split.by = 'Treatment')
dev.off()
pdf("DimSplit_Gender.pdf")
DimPlot(object = nCoV.integrated, reduction = 'tsne',label = FALSE, split.by = 'Gender')
dev.off()

