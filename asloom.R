library(dplyr)
library(Seurat)
print(Sys.time())
print("load library complete...")

setwd("singlecell/P-MATRIX")


print("load RDS complete...")
pbmc<- readRDS("sample.combined_T10_others_singleR_ref3_main_fine_celltype.rds")
#write.csv(pbmc@meta.data,file="sample.combined_T10_others_singleR_ref3_main_fine_celltype.csv")
scdata1=subset(pbmc,Treatment!=5)
scdata2=subset(scdata1,Response!="Na")
print(Sys.time())
pbmc.loom <- as.loom(scdata2, filename = "Tall_celltype.loom", verbose = FALSE)
pbmc.loom

pbmc.loom$close_all()