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
  Mdata$Group=group
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

M164.gene=prepdata(M164,"M164","PR")
M165.gene=prepdata(M165,"M165","PR")
M174.gene=prepdata(M174,"M174","PD")
M175.gene=prepdata(M175,"M175","PD")
M176.gene=prepdata(M176,"M176","PR")
M177.gene=prepdata(M177,"M177","SD")



sample.anchors <- FindIntegrationAnchors(object.list = list(
  M164.gene,
  M165.gene,
  M174.gene,
  M175.gene,
  M176.gene,
  M177.gene,
), dims = 1:30)
sample.combined <- IntegrateData(anchorset = sample.anchors, dims = 1:30)

saveRDS(sample.combined, file="sample.combined164_177.rds")
