library(Seurat)
library(celldex)
library(SingleR)
library(dplyr)
library(ggplot2)
library(reshape2)
library(tibble)
print("load library complete...")



setwd("singlecell/P-MATRIX")

# scdata<- readRDS("sample.combined_T10_others_singleR_ref3_main_fine.rds")
# 
# scdata<- RenameIdents(object = scdata,
#                       "0"="CD14+ monocytes (M1)",
#                       "1"="naive CD4+ T cells (T1)",
#                       "2"="naive CD4+ T cells (T1)",
#                       "3"="C56-CD16+ NK cells (NK2)",
#                       "4"="CD16+ monocytes (M2)",
#                       "5"="C56-CD16+ NK cells (NK2)",
#                       "6"="effector memory CD8+T cells (T5, CD8 Tm)",
#                       "7"="naive B cells (B1)",
#                       "8"="CD14+ monocytes (M1)",
#                       "9"="CD14+ monocytes (M1)",
#                       "10"="C56-CD16+ NK cells (NK2)",
#                       "11"="CD14+ monocytes (M1)",
#                       "12"="CD14+ monocytes (M1)",
#                       "13"="CD14+ monocytes (M1)",
#                       "14"="immature B cells (B2)",
#                       "15"="cytotoxic CD8+ lymphocytes (CD8+ CTL) (T6)",
#                       "16"="CD4+ T cells (T3)",
#                       "17"="cytotoxic CD8+ lymphocytes (CD8+ CTL) (T6)",
#                       "18"="C56-CD16+ NK cells (NK2)",
#                       "19"="naive CD8+ T cells (T4)",
#                       "20"="regulatory CD4+ T cells (T2, Treg)",
#                       "21"="Undefined T cells",
#                       "22"="CD14+ monocytes (M1)",
#                       "23"="γδ T cells (T9) ",
#                       "24"="CD1C+ dendritic cells (M3)",
#                       "25"="Proliferating T cells (T7, Tprol) ",
#                       "26"="CD14+ monocytes (M1)",
#                       "27"="naive CD4+ T cells (T1)",
#                       "28"="CD56+CD16- NK cells (NK1)",
#                       "29"="cytotoxic CD8+ lymphocytes (CD8+ CTL) (T6)",
#                       "30"="C56-CD16+ NK cells (NK2)",
#                       "31"="CD16+ monocytes (M2)",
#                       "32"="Platelets",
#                       "33"="cytotoxic CD8+ lymphocytes (CD8+ CTL) (T6)",
#                       "34"="memory B cells (B3)",
#                       "35"="Undefined T cells",
#                       "36"="CD14+ monocytes (M1)",
#                       "37"="CD14+ monocytes (M1)",
#                       "38"="cytotoxic CD8+ lymphocytes (CD8+ CTL) (T6)",
#                       "39"="cytotoxic CD8+ lymphocytes (CD8+ CTL) (T6)",
#                       "40"="Undefined T cells",
#                       "41"="cytotoxic CD8+ lymphocytes (CD8+ CTL) (T6)",
#                       "42"="mucosal-associated invariant T (MAIT) cells (T8) ",
#                       "43"="cytotoxic CD8+ lymphocytes (CD8+ CTL) (T6)",
#                       "44"="plasmacytoid dendritic cells (CLEC4C+ pDC) (M4)",
#                       "45"="cytotoxic CD8+ lymphocytes (CD8+ CTL) (T6)",
#                       "46"="CD16+ monocytes (M2)",
#                       "47"="cytotoxic CD8+ lymphocytes (CD8+ CTL) (T6)",
#                       "48"="Platelets",
#                       "49"="cytotoxic CD8+ lymphocytes (CD8+ CTL) (T6)",
#                       "50"="CD14+ monocytes (M1)",
#                       "51"="CD14+ monocytes (M1)",
#                       "52"="naive B cells (B1)",
#                       "53"="cytotoxic CD8+ lymphocytes (CD8+ CTL) (T6)",
#                       "54"="C56-CD16+ NK cells (NK2)",
#                       "55"="HSC",
#                       "56"="CLEC9A+ dendritic cells (M5)",
#                       "57"="C56-CD16+ NK cells (NK2)",
#                       "58"="naive B cells (B1)"
# )
# 
# scdata@meta.data$celltype = Idents(scdata)
# saveRDS(scdata, file="sample.combined_T10_others_singleR_ref3_main_fine_celltype.rds")

scdata0<- readRDS("sample.combined_T10_others_singleR_ref3_main_fine_celltype.rds")
scdata1=subset(scdata0,Treatment!=5)
scdata=subset(scdata1,Response!="Na")

# pbmc=subset(scdata,singleR.main2=="Monocyte")
# 
# 
# pdf("monocytes_Tall_celltype_dimplot_tsne_split_Response.pdf")
# DimPlot(pbmc , reduction = "tsne", group.by = "celltype", split.by = 'Response',pt.size=1,ncol=1,label = T)
# dev.off()
# 
# pdf("monocytes_Tall_celltype_dimplot_tsne.pdf")
# DimPlot(pbmc , reduction = "tsne", group.by = "celltype",pt.size=1,ncol=1,label = T)
# dev.off()


#monocytes
# markers = c(
#   "LYZ",
#   "CD14",
#   "CD16",
#   "FCGR3A",
#   "CD1C",
#   "LILRA4",
#   "CLEC9A",
#   "CLEC4C"
#   )

# print("write DotPlot_markers_T10 complete...")
# 
# pdf("Monocyte_DotPlot_markers_T10_Tothers_cluster1.pdf")
# pp = DotPlot(pbmc, features = markers,cols = c('white','#F8766D'),dot.scale =5,group.by ="seurat_clusters" ) + RotatedAxis()
# pp = pp +
#   theme(axis.text.x = element_text(size = 8),
#         axis.text.y = element_text(size = 8)) +
#   labs(x='',y='') +
#   guides(color = guide_colorbar(title = 'Scale expression'),
#          size = guide_legend(title = 'Percent expressed')) +
#   theme(axis.line = element_line(size = 0.6))
# pp
# dev.off()
# 
# pdf("Monocyte_FeaturePlot_markers_T10_Tothers1.pdf")
# FeaturePlot(pbmc , features = markers,cols = c('gray','#F8766D'),reduction = "tsne")
# dev.off()


# ###Tcell
# clusters_files=c(1, 2, 6, 15, 16, 17, 19, 20, 21,
#                  23, 25, 27, 29, 33, 35, 38, 39,
#                  40, 41, 42, 43, 45, 47, 49, 53
# )
# 
# pbmc=subset(scdata,seurat_clusters %in% clusters_files)
# 
# 
# 
# pdf("Tcell_Tall_celltype_dimplot_tsne_split_Response.pdf")
# DimPlot(pbmc , reduction = "tsne", group.by = "celltype", split.by = 'Response',pt.size=1,ncol=1,label = T)
# dev.off()
# 
# pdf("Tcell_Tall_celltype_dimplot_tsne.pdf")
# DimPlot(pbmc , reduction = "tsne", group.by = "celltype",pt.size=1,ncol=1,label = T)
# dev.off()
# markers=c(
#  "CD3E",  "CD4", "CD8A",  "CD8B", "CCR7",
#   "LEF1", "TCF7",  "AQP3", "CD69",  "CCR6",
#  "CXCR6",  "CCL5", "PRDM1",  "FOXP3", "GZMK",
#   "GZMB", "GNLY",  "PRF1", "TYMS",  "MKI67",
#   "SLC4A10", "TRGV9"
# 
# )

# print("write DotPlot_markers_T10 complete...")
# 
# pdf("Tcell_DotPlot_markers_T10_Tothers_cluster1.pdf")
# pp = DotPlot(pbmc, features = markers,cols = c('white','#F8766D'),dot.scale =5,group.by ="seurat_clusters" ) + RotatedAxis()
# pp = pp +
#   theme(axis.text.x = element_text(size = 8),
#         axis.text.y = element_text(size = 8)) +
#   labs(x='',y='') +
#   guides(color = guide_colorbar(title = 'Scale expression'),
#          size = guide_legend(title = 'Percent expressed')) +
#   theme(axis.line = element_line(size = 0.6))
# pp
# dev.off()
# 
# pdf("Tcell_FeaturePlot_markers_T10_Tothers1.pdf")
# FeaturePlot(pbmc , features = markers,cols = c('gray','#F8766D'),reduction = "tsne")
# dev.off()
# 
# 
# # 
pbmc=subset(scdata,singleR.main2=="NK cells")

pdf("NKcell_Tall_celltype_dimplot_tsne_split_Response.pdf")
DimPlot(pbmc , reduction = "tsne", group.by = "celltype", split.by = 'Response',pt.size=0.1,ncol=1,label = T)
dev.off()

pdf("NKcell_Tall_celltype_dimplot_tsne.pdf")
DimPlot(pbmc , reduction = "tsne", group.by = "celltype",pt.size=0.1,ncol=1,label = T)
dev.off()
# # 
# #NK cells
# markers = c(
#   "KLRF1",
#   "KLRD1",
#   "KLRC1",
#   "NCAM1",#"CD56"
#   "FCGR3A" #"CD16"
# 
#   
#   )
# 
# print("write DotPlot_markers_T10 complete...")
# 
# pdf("NKcells_DotPlot_markers_T10_Tothers_cluster.pdf")
# pp = DotPlot(pbmc, features = markers,cols = c('white','#F8766D'),dot.scale =5,group.by ="seurat_clusters" ) + RotatedAxis()
# pp = pp +
#   theme(axis.text.x = element_text(size = 8),
#         axis.text.y = element_text(size = 8)) +
#   labs(x='',y='') +
#   guides(color = guide_colorbar(title = 'Scale expression'),
#          size = guide_legend(title = 'Percent expressed')) +
#   theme(axis.line = element_line(size = 0.6))
# pp
# dev.off()
# 
# pdf("NKcells_FeaturePlot_markers_T10_Tothers.pdf")
# FeaturePlot(pbmc , features = markers,cols = c('gray','#F8766D'),reduction = "tsne")
# dev.off()


pbmc=subset(scdata,singleR.main2=="B-cells")
pdf("Bcell_Tall_celltype_dimplot_tsne_split_Response.pdf")
DimPlot(pbmc , reduction = "tsne", group.by = "celltype", split.by = 'Response',pt.size=0.1,ncol=1,label = T)
dev.off()

pdf("Bcell_Tall_celltype_dimplot_tsne.pdf")
DimPlot(pbmc , reduction = "tsne", group.by = "celltype",pt.size=0.1,ncol=1,label = T)
dev.off()

#B-cells
# markers = c(
#   "CD19",
#   "MS4A1",
#   "IGHD",
#   "IGHM",
#   "IL4R",
#   "TCL1A",
#   "CD27",
#   "CD38",
#   "IGHG",
#   "XBP1",
#   "MZB1"
#   
# )
# 
# print("write DotPlot_markers_T10 complete...")
# 
# pdf("Bcells_DotPlot_markers_T10_Tothers_cluster.pdf")
# pp = DotPlot(pbmc, features = markers,cols = c('white','#F8766D'),dot.scale =5,group.by ="seurat_clusters" ) + RotatedAxis()
# pp = pp +
#   theme(axis.text.x = element_text(size = 8),
#         axis.text.y = element_text(size = 8)) +
#   labs(x='',y='') +
#   guides(color = guide_colorbar(title = 'Scale expression'),
#          size = guide_legend(title = 'Percent expressed')) +
#   theme(axis.line = element_line(size = 0.6))
# pp
# dev.off()
# 
# pdf("Bcells_FeaturePlot_markers_T10_Tothers.pdf")
# FeaturePlot(pbmc , features = markers,cols = c('gray','#F8766D'),reduction = "tsne")
# dev.off()



