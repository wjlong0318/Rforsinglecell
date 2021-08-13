library(Seurat)
library(cowplot)
library(Matrix)
library(dplyr)
library(ggplot2)
library(reshape2)
library(hdf5r)
setwd("H:\\shanghai\\singlecell\\matrix.h5")
M164<-Read10X_h5("164\\filtered_feature_bc_matrix.h5")
M165<-Read10X_h5("165\\filtered_feature_bc_matrix.h5")
M174<-Read10X_h5("174\\filtered_feature_bc_matrix.h5")
M175<-Read10X_h5("175\\filtered_feature_bc_matrix.h5")
M164.gene=prepdata(M164,"M164","PR")
M165.gene=prepdata(M165,"M165","PR")
M174.gene=prepdata(M174,"M174","PD")
M175.gene=prepdata(M175,"M175","PD")


# saveRDS(M175.gene, file="M175.gene.Rds")
# M175.gene.rds<- readRDS("M175.gene.Rds")

sample.anchors <- FindIntegrationAnchors(object.list = list(M164.gene,M165.gene), dims = 1:50)
sample.combined <- IntegrateData(anchorset = sample.anchors, dims = 1:50)

prepdata=function(data,label,group){
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

M164<-Read10X_h5("164\\filtered_feature_bc_matrix.h5")
Mdata<-CreateSeuratObject(counts = M164, project ="M164",min.cells = 3, min.features = 200)
Mdata[["percent.mt"]]<-PercentageFeatureSet(Mdata,pattern = "^MT")
p164<-VlnPlot(Mdata, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)



