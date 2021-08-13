library(Seurat)
library(cowplot)
library(hdf5r)

setwd("H:\\shanghai\\singlecell\\P-MATRIX")
#加载Cellranger output file
balbc_pbmc <- Read10X("P-175")
s_balbc_pbmc <- CreateSeuratObject(counts = balbc_pbmc, min.cells = 3, min.features = 200, project = "cellranger")
##counts: Either a matrix-like object with unnormalized data with cells as columns and features as rows or an Assay-derived object
##project：Project name for the Seurat object
##min.cells: Include features detected in at least this many cells. Will subset the counts matrix as well. To reintroduce excluded features, create a new object with a lower cutoff.
##min.features: Include cells where at least this many features are detected.

#提取线粒体基因
s_balbc_pbmc$percent.mito <- PercentageFeatureSet(s_balbc_pbmc, pattern = "^mt-")

#增加T和B细胞的克隆信息
add_clonotype <- function(tcr_prefix, seurat_obj, type="t"){
  #读入contig_annotations.csv
  tcr <- read.csv(paste(tcr_prefix,"filtered_contig_annotations.csv", sep="/"))
  # Remove the -1 at the end of each barcode
  # tcr$barcode <- gsub("-1", "", tcr$barcode)
  
  # Subsets so only the first line of each barcode is kept,as each entry for given barcode will have same clonotype.
  tcr <- tcr[!duplicated(tcr$barcode), ]
  names(tcr)[names(tcr) == "raw_clonotype_id"] <- "clonotype_id"
  # Clonotype-centric info.
  clono <- read.csv(paste(tcr_prefix,"clonotypes.csv", sep="/"))
  # Slap the AA sequences onto our original table by clonotype_id.
  tcr <- merge(tcr, clono)
  #Reorder so barcodes are first column and set them as rownames.
  rownames(tcr) <- tcr[,2]
  tcr[,2] <- NULL
  colnames(tcr) <- paste(type, colnames(tcr), sep="_")
  # Add to the Seurat object's metadata.
  clono_seurat <- AddMetaData(object=seurat_obj, metadata=tcr)
  return(clono_seurat)
}
s_balbc_pbmc <- add_clonotype("P-175-T", s_balbc_pbmc, "t")
s_balbc_pbmc <- add_clonotype("P-175-B", s_balbc_pbmc, "b")
head(s_balbc_pbmc[[]])
table(!is.na(s_balbc_pbmc$t_clonotype_id),!is.na(s_balbc_pbmc$b_clonotype_id))
s_balbc_pbmc <- subset(s_balbc_pbmc, cells = colnames(s_balbc_pbmc)[!(!is.na(s_balbc_pbmc$t_clonotype_id) & !is.na(s_balbc_pbmc$b_clonotype_id))])
table(!is.na(s_balbc_pbmc$t_clonotype_id),!is.na(s_balbc_pbmc$b_clonotype_id))
#进行常规workflow
s_balbc_pbmc <- subset(s_balbc_pbmc, percent.mito <= 10)
s_balbc_pbmc <- subset(s_balbc_pbmc, nCount_RNA >= 500 & nCount_RNA <= 40000)
s_balbc_pbmc <- NormalizeData(s_balbc_pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
s_balbc_pbmc <- FindVariableFeatures(s_balbc_pbmc, selection.method = "vst", nfeatures = 2000)

all.genes <- rownames(s_balbc_pbmc)
s_balbc_pbmc <- ScaleData(s_balbc_pbmc, features = all.genes)
s_balbc_pbmc <- RunPCA(s_balbc_pbmc, features = VariableFeatures(object = s_balbc_pbmc))
use.pcs = 1:30
s_balbc_pbmc <- FindNeighbors(s_balbc_pbmc, dims = use.pcs)
s_balbc_pbmc <- FindClusters(s_balbc_pbmc, resolution = c(0.5))
s_balbc_pbmc <- RunUMAP(s_balbc_pbmc, dims = use.pcs)
p1 <- DimPlot(s_balbc_pbmc, reduction = "umap", label = TRUE)
p2 <-DimPlot(s_balbc_pbmc,group.by = "t_chain")

#T细胞的marker的表达情况：
t_cell_markers <- c("CD3D","CD4")
p3 <- FeaturePlot(s_balbc_pbmc, features = t_cell_markers)

#在UMAP 图上表示指定蛋白序列，eg:IGH:CARWGGYGYDGGYFDYW;IGK:CGQSYSYPYTF，然后：
p4 <- DimPlot(s_balbc_pbmc, cells.highlight = Cells(subset(s_balbc_pbmc, subset = b_cdr3s_aa == "IGH:CAREAGSSGKAGWFDPW;IGK:CQQYNSYPVTF")))

library(patchwork)
p1+p2+p3+p4
library(sankeywheel)
sankeywheel(
  from = s_balbc_pbmc@meta.data$t_chain, to = s_balbc_pbmc@meta.data$t_v_gene,
  weight = s_balbc_pbmc@meta.data$nCount_RNA,#type = "sankey",
  title = "sankeywheel",
  subtitle = "chain & t_v_gene",
  seriesName = "", width = "100%", height = "600px")
#计算不同链的差异基因
library(DT)
DT::datatable(FindMarkers(s_balbc_pbmc,group.by = "b_chain",ident.1 = 'IGK'))
mk<-FindMarkers(s_balbc_pbmc,group.by = "b_chain",ident.1 = 'IGK')
mk_max10 <- rownames(mk[order(mk$p_val,decreasing = T)[1:10],])
DotPlot(s_balbc_pbmc,features=mk_max10,group.by = group) + RotatedAxis()
#有问题，无法分组
RidgePlot(s_balbc_pbmc,features=mk_max10,group.by = "b_chain")
RidgePlot(s_balbc_pbmc,features=t_cell_markers)
DotPlot(s_balbc_pbmc,features=t_cell_markers) + RotatedAxis()
