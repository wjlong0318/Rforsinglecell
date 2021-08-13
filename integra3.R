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


M185<-Read10X("P-185")
M186<-Read10X("P-186")
M187<-Read10X("P-187")
M189<-Read10X("P-189")
M191<-Read10X("P-191")
M192<-Read10X("P-192")


M185.gene=prepdata(M185,"M185","SD")
M186.gene=prepdata(M186,"M186","SD")
M187.gene=prepdata(M187,"M187","SD")
M189.gene=prepdata(M189,"M189","PR")
M191.gene=prepdata(M191,"M191","PR")
M192.gene=prepdata(M192,"M192","PD")



sample.anchors <- FindIntegrationAnchors(object.list = list(
  
 
  M185.gene,
  M186.gene,
  M187.gene,
  M189.gene,
  M191.gene,
  M192.gene
), dims = 1:30)
sample.combined <- IntegrateData(anchorset = sample.anchors, dims = 1:30)

saveRDS(sample.combined, file="sample.combined185_192.rds")
