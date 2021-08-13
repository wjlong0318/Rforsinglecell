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


M193<-Read10X("P-193")
M197<-Read10X("P-197")
M198<-Read10X("P-198")
M199<-Read10X("P-199")
M202<-Read10X("P-202")

M193.gene=prepdata(M193,"M193","PD")
M197.gene=prepdata(M197,"M197","PD")
M198.gene=prepdata(M198,"M198","PR")
M199.gene=prepdata(M199,"M199","PR")
M202.gene=prepdata(M202,"M202","PD")


sample.anchors <- FindIntegrationAnchors(object.list = list(
  M193.gene,
  M197.gene,
  M198.gene,
  M199.gene,
  M202.gene), dims = 1:50)
sample.combined <- IntegrateData(anchorset = sample.anchors, dims = 1:50)

saveRDS(sample.combined, file="sample.combined193_202.rds")
