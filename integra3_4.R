library(Seurat)
library(cowplot)
library(Matrix)
library(dplyr)
library(ggplot2)
library(reshape2)
setwd("singlecell/P-MATRIX")

print("start...\n")


in1<- readRDS("sample.combined164_177.rds")
in2<- readRDS("sample.combined178_184.rds")
sample.anchors <- FindIntegrationAnchors(object.list = list(
  in1,in2
), dims = 1:30)
sample.combined <- IntegrateData(anchorset = sample.anchors, dims = 1:30)

saveRDS(sample.combined, file="sample.combined164_184.rds")
