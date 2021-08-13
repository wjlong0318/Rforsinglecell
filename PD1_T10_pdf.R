library(Seurat)
library(cowplot)
library(Matrix)
library(dplyr)
library(ggplot2)
library(reshape2)
#setwd("singlecell/P-MATRIX")

setwd("H:\\shanghai\\singlecell\\P-MATRIX")


print("start...\n")

add_clonotype <- function(tcr_prefix, seurat_obj, type="t"){
  print("load function add_clonotype...\n")
  tcr <- read.csv(paste(tcr_prefix,"filtered_contig_annotations.csv", sep="/"))
   tcr <- tcr[!duplicated(tcr$barcode), ]
  names(tcr)[names(tcr) == "raw_clonotype_id"] <- "clonotype_id"
   clono <- read.csv(paste(tcr_prefix,"clonotypes.csv", sep="/"))
  tcr <- merge(tcr, clono)
  rownames(tcr) <- tcr[,2]
  tcr[,2] <- NULL
  colnames(tcr) <- paste(type, colnames(tcr), sep="_")
  clono_seurat <- AddMetaData(object=seurat_obj, metadata=tcr)
  return(clono_seurat)
}

prepdata=function(data,label){
  print("load function predata...\n")
  Mdata<-CreateSeuratObject(counts = data, project =label,min.cells = 3, min.features = 200)
  tname=paste(label,"-T",sep="")
  bname=paste(label,"-B",sep="")
  Mdata<- add_clonotype(tname, Mdata, "t")
  Mdata<- add_clonotype(bname, Mdata, "b")

  s_balbc_pbmc <- subset(Mdata, cells = colnames(Mdata)[!(!is.na(Mdata$t_clonotype_id) & !is.na(Mdata$b_clonotype_id))])
  
  Mdata[["percent.mt"]]<-PercentageFeatureSet(Mdata,pattern = "^MT")
  Mdata.filter <- subset(Mdata, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt < 10 & nCount_RNA >1000)
  Mdata.norm <- NormalizeData(Mdata.filter)
  Mdata.genes <- FindVariableFeatures(Mdata.norm, selection.method = "vst", nfeatures = 2000)
  Mdata.genes
}

print("read files...\n")
v=c("P-165",
    #"P-174",
    # "P-175",
    # "P-177",
    # "P-179",
    # "P-180",
    # "P-183",
    # "P-185",
    # "P-186",
    # "P-191",
    # "P-197",
    # "P-198",
    # "P-199",
    # "P-202"
    "P-164"
)
all=list()
i=0
for(n in v ){
  i=i+1
  file_10x=Read10X(n)
   M_gene=prepdata(file_10x,n)
   all[[i]]=M_gene
}


sample.anchors <- FindIntegrationAnchors(object.list = all, dims = 1:30)
sample.combined <- IntegrateData(anchorset = sample.anchors, dims = 1:30)

saveRDS(sample.combined, file="IntegrateSample_T10.rds")
# nCoV.integrated<- readRDS("sample.combined164_202.rds")
nCoV.integrated=sample.combined
sample_info = as.data.frame(colnames(nCoV.integrated))
colnames(sample_info) = c('ID')
rownames(sample_info) = sample_info$ID
sample_info$sample = nCoV.integrated@meta.data$orig.ident
samples=read.csv("../meta.csv",header=T,sep=",")
sample_info = dplyr::left_join(sample_info,samples)
rownames(sample_info) = sample_info$ID
nCoV.integrated = AddMetaData(object = nCoV.integrated, metadata = sample_info)

DefaultAssay(nCoV.integrated) <- "RNA"
nCoV.integrated[['percent.mito']] <- PercentageFeatureSet(nCoV.integrated, pattern = "^MT-")
nCoV.integrated <- NormalizeData(object = nCoV.integrated, normalization.method = "LogNormalize", scale.factor = 1e4)
nCoV.integrated <- FindVariableFeatures(object = nCoV.integrated, selection.method = "vst", nfeatures = 2000,verbose = FALSE)
nCoV.integrated <- ScaleData(nCoV.integrated, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito"))
#VlnPlot(object = nCoV.integrated, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
#FeatureScatter(object = nCoV.integrated, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
nCoV.integrated <- ScaleData(nCoV.integrated, verbose = FALSE, vars.to.regress = c("nCount_RNA", "percent.mito"))
nCoV.integrated <- RunPCA(nCoV.integrated, verbose = FALSE,npcs = 100)
nCoV.integrated <- ProjectDim(object = nCoV.integrated)
#ElbowPlot(object = nCoV.integrated,ndims = 100)

nCoV.integrated <- FindNeighbors(object = nCoV.integrated, dims = 1:50)
nCoV.integrated <- FindClusters(object = nCoV.integrated, resolution = 1.2)
nCoV.integrated <- RunTSNE(object = nCoV.integrated, dims = 1:50)
DimPlot(object = nCoV.integrated, reduction = 'tsne',label = TRUE)

markers = c("CD4",
            "CD8A",
            "FCGR3A",
            "CD14",
            "MS4A1",
            "CD3D",
            "KLRF1",
            "KLRD1",
            "FCER1A",
            "PF4",
            "CD68",
            "IGHG4"
)
pp_temp = FeaturePlot(object = nCoV.integrated, features = markers,cols = c("lightgrey","#ff0000"),combine = FALSE)
plots <- lapply(X = pp_temp, FUN = function(p) p + 
                  theme(axis.title = element_text(size = 8),
                        axis.text = element_text(size = 8),
                        plot.title = element_text(family = 'sans',face='italic',size=8),
                        legend.text = element_text(size = 8),
                        legend.key.height = unit(1.8,"line"),
                        legend.key.width = unit(1.2,"line") ))
pdf("combineplot.pdf")
CombinePlots(plots = plots,ncol = 4,legend = 'right')
dev.off()

pdf("DimCluster1.pdf")
DimPlot(object = nCoV.integrated, reduction = 'tsne',label = TRUE)
dev.off()

pdf("DotPlot.pdf")
pp = DotPlot(nCoV.integrated, features = rev(markers),cols = c('white','#F8766D'),dot.scale =5) + RotatedAxis()
pp = pp + 
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8)) + 
  labs(x='',y='') + 
  guides(color = guide_colorbar(title = 'Scale expression'),
         size = guide_legend(title = 'Percent expressed')) + 
  theme(axis.line = element_line(size = 0.6))
dev.off()


pdf("DimCluster_Response.pdf")
DimPlot(object = nCoV.integrated, reduction = 'tsne',label = FALSE, group.by = 'Response')
dev.off()
pdf("DimCluster_Treatment.pdf")
DimPlot(object = nCoV.integrated, reduction = 'tsne',label = FALSE, group.by = 'Treatment')
dev.off()
pdf("DimCluster_Gender.pdf")
DimPlot(object = nCoV.integrated, reduction = 'tsne',label = FALSE, group.by = 'Gender')
dev.off()
pdf("DimSplit_Response.pdf")
DimPlot(object = nCoV.integrated, reduction = 'tsne',label = FALSE, split.by = 'Response')
dev.off()
pdf("DimSplit_Treatment.pdf")
DimPlot(object = nCoV.integrated, reduction = 'tsne',label = FALSE, split.by = 'Treatment')
dev.off()
pdf("DimSplit_Gender.pdf")
DimPlot(object = nCoV.integrated, reduction = 'tsne',label = FALSE, split.by = 'Gender')
dev.off()

