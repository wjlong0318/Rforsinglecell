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


M178<-Read10X("P-178")
M179<-Read10X("P-179")
M180<-Read10X("P-180")
M181<-Read10X("P-181")
M183<-Read10X("P-183")
M184<-Read10X("P-184")


M178.gene=prepdata(M178,"M178","PD")
M179.gene=prepdata(M179,"M179","SD")
M180.gene=prepdata(M180,"M180","PD")
M181.gene=prepdata(M181,"M181","PR")
M183.gene=prepdata(M183,"M183","SD")
M184.gene=prepdata(M184,"M184","SD")



sample.anchors <- FindIntegrationAnchors(object.list = list(
  
  M178.gene,
  M179.gene,
  M180.gene,
  M181.gene,
  M183.gene,
  M184.gene
  ), dims = 1:30)
sample.combined <- IntegrateData(anchorset = sample.anchors, dims = 1:30)

saveRDS(sample.combined, file="sample.combined178_184.rds")
