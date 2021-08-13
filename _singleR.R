rm(list=ls())
options(stringsAsFactors = F)
library(Seurat)
library(ggplot2)
load(file = 'first_sce.Rdata')
# Specify genes  
genes_to_check = c("PTPRC","EPCAM","CD3G","CD3E", "CD79A", "BLNK","MS4A1", "CD68", "CSF1R", 
                   "MARCO", "CD207", "PMEL", "ALB", "C1QB", "CLDN5", "FCGR3B", "COL1A1")
# All on Dotplot 
p <- DotPlot(sce, features = genes_to_check) + coord_flip()
p

dat=p$data 
cd45=dat[dat$features.plot=='PTPRC',]
fivenum(cd45$avg.exp.scaled)
imm=cd45[cd45$avg.exp.scaled > -0.5,]$id
imm
sce@meta.data$immune_annotation <-ifelse(sce@meta.data$RNA_snn_res.0.5  %in% imm ,'immune','non-immune')
# MAke a table 
table(sce@meta.data$immune_annotation)

p <- TSNEPlot(object = sce, group.by = 'immune_annotation')
p 

sce@meta.data$immune_annotation <-ifelse(sce@meta.data$RNA_snn_res.0.5  %in% imm ,'immune','non-immune')
# MAke a table 
table(sce@meta.data$immune_annotation)
phe=sce@meta.data
save(phe,file = 'phe-of-immune-or-not.Rdata')

rm(list=ls())
options(stringsAsFactors = F)
library(Seurat)
library(ggplot2)
load(file = 'first_sce.Rdata')
sce <- FindClusters(sce, resolution = 0.5)
table(sce@meta.data$RNA_snn_res.0.5)  
load(file = 'phe-of-immune-or-not.Rdata')
table(phe$immune_annotation)
cells.use <- row.names(sce@meta.data)[which(phe$immune_annotation=='immune')]
length(cells.use)
sce <-subset(sce, cells=cells.use)  
sce

sce
sce <- NormalizeData(sce, normalization.method =  "LogNormalize", 
                     scale.factor = 10000)
GetAssay(sce,assay = "RNA")
sce <- FindVariableFeatures(sce, 
                            selection.method = "vst", nfeatures = 2000) 
sce <- ScaleData(sce) 
sce <- RunPCA(object = sce, pc.genes = VariableFeatures(sce)) 
res.used <- 0.7
sce <- FindClusters(object = sce, verbose = T, resolution = res.used)
set.seed(123)
sce <- RunTSNE(object = sce, dims = 1:15, do.fast = TRUE)
DimPlot(sce,reduction = "tsne",label=T)
DimPlot(sce,reduction = "tsne",label=T, group.by = "patient_id")
table(sce@meta.data$seurat_clusters) 

sce_for_SingleR <- GetAssayData(sce, slot="data")
sce_for_SingleR
library(SingleR)
hpca.se <- HumanPrimaryCellAtlasData()
hpca.se
clusters=sce@meta.data$seurat_clusters
pred.hesc <- SingleR(test = sce_for_SingleR, ref = hpca.se, labels = hpca.se$label.main,
                     method = "cluster", clusters = clusters, 
                     assay.type.test = "logcounts", assay.type.ref = "logcounts")
table(pred.hesc$labels)
celltype = data.frame(ClusterID=rownames(pred.hesc), celltype=pred.hesc$labels, stringsAsFactors = F) 
sce@meta.data$singleR=celltype[match(clusters,celltype$ClusterID),'celltype']
DimPlot(sce, reduction = "tsne", group.by = "singleR")
phe=sce@meta.data
table(phe$singleR)
save(phe,file = 'phe-of-subtypes-Immune-by-singleR.Rdata')

as.data.frame(sort(table(phe$singleR)))

free_annotation <- c("T-cells","MF-Monocytes", "MF-Monocytes", "B-cells-PB", "MF-Monocytes", "T-cells", "T-cells", "Neutrophils", "Dendritic", "Mast-cells", "MF-Monocytes", "T-cells", "B-cells-M", "Unknown", "T-cells", "pDCs", "B-cells-M", "MF-Monocytes")
