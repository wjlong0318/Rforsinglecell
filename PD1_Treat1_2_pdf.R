library(Seurat)
library(cowplot)
library(Matrix)
library(dplyr)
library(ggplot2)
library(reshape2)
setwd("singlecell/P-MATRIX")

#setwd("H:\\shanghai\\singlecell\\P-MATRIX")


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
v=c("P-164",
   "P-165",
   "P-174",
   "P-175",
   "P-176",
   "P-177",
   "P-180",
   "P-181",
   "P-184",
   "P-186",
   "P-187",
   "P-191",
   "P-193",
   "P-197",
   "P-198",
   "P-199",
   "P-202"
    
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

saveRDS(sample.combined, file="IntegrateSample_Treat1_2.rds")

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
nCoV.integrated <- RunPCA(nCoV.integrated, verbose = FALSE,npcs = 100)
nCoV.integrated <- ProjectDim(object = nCoV.integrated)

nCoV.integrated <- FindNeighbors(object = nCoV.integrated, dims = 1:50)
nCoV.integrated <- FindClusters(object = nCoV.integrated, resolution = 1.2)
nCoV.integrated <- RunTSNE(object = nCoV.integrated, dims = 1:50)
nCoV.integrated  <- RunUMAP(nCoV.integrated  , dims = 1:50)

saveRDS(nCoV.integrated, file="IntegrateSample_tsne_umap_Treat1_2.rds")

print("SaveData IntegrateSample_tsne_umap_Treat1_2.rds complete...\n")
pdf("umap_Treat1_2_dimplot.pdf")
DimPlot(nCoV.integrated , reduction = "umap", label = TRUE)
dev.off()

pdf("tsne_Treat1_2_dimplot.pdf")
DimPlot(nCoV.integrated , reduction = "tsne", label = TRUE)
dev.off()

pdf("tchain_tsne_Treat1_2_dimplot.pdf")
DimPlot(nCoV.integrated , reduction = "tsne", group.by = "t_chain")
dev.off()

pdf("bchain_tsne_Treat1_2_dimplot.pdf")
DimPlot(nCoV.integrated , reduction = "tsne", group.by = "b_chain")
dev.off()

pdf("tchain_umap_Treat1_2_dimplot.pdf")
DimPlot(nCoV.integrated , reduction = "umap", group.by = "t_chain")
dev.off()

pdf("bchain_umap_Treat1_2_dimplot.pdf")
DimPlot(nCoV.integrated , reduction = "umap", group.by = "b_chain")
dev.off()

pdf("umap_Response_Treat1_2_dimplot.pdf")
DimPlot(nCoV.integrated , reduction = "umap", label = TRUE,split.by = 'Response')
dev.off()

pdf("tsne_Response_Treat1_2_dimplot.pdf")
DimPlot(nCoV.integrated , reduction = "tsne", label = TRUE,split.by = 'Response')
dev.off()

pdf("tchain_tsne_Response_Treat1_2_dimplot.pdf")
DimPlot(nCoV.integrated , reduction = "tsne", group.by = "t_chain",split.by = 'Response')
dev.off()

pdf("bchain_tsne_Response_Treat1_2_dimplot.pdf")
DimPlot(nCoV.integrated , reduction = "tsne", group.by = "b_chain",split.by = 'Response')
dev.off()

pdf("tchain_umap_Response_Treat1_2_dimplot.pdf")
DimPlot(nCoV.integrated , reduction = "umap", group.by = "t_chain",split.by = 'Response')
dev.off()

pdf("bchain_umap_Response_Treat1_2_dimplot.pdf")
DimPlot(nCoV.integrated , reduction = "umap", group.by = "b_chain",split.by = 'Response')
dev.off()

print("start singleR...")
library(celldex)
library(SingleR)
pbmc=readRDS("IntegrateSample_tsne_umap_Treat1_2.rds")

pbmc_for_SingleR <- GetAssayData(pbmc, slot="data")
print("read pbmc_for_SingleR...")

hpca.se <- celldex::HumanPrimaryCellAtlasData()
Blue.se=celldex::BlueprintEncodeData() 
Immune.se=celldex::DatabaseImmuneCellExpressionData()
print("read hpca.se,blue.se,Immune.se...")


clusters=pbmc@meta.data$seurat_clusters
print("read clusters...")
head(clusters)

pred.hesc <- SingleR(test = pbmc_for_SingleR, ref = list(BP=Blue.se, HPCA=hpca.se,Im=Immune.se), labels =list(Blue.se$label.fine, hpca.se$label.fine,Immune.se$label.fine),
                     method = "cluster", clusters = clusters, 
                     assay.type.test = "logcounts", assay.type.ref = "logcounts")
print("[result]pred.hesc...")
pred.hesc$labels

celltype = data.frame(ClusterID=rownames(pred.hesc), celltype=pred.hesc$labels, stringsAsFactors = F) 
print("[result]celltype...")
celltype
pbmc@meta.data$singleR=celltype[match(clusters,celltype$ClusterID),'celltype']
#write.csv(pbmc@meta.data,file="pbmc_Treat1_2_meta.data.ref3_fine.csv")


pred.main <- SingleR(test = pbmc_for_SingleR, ref = list(BP=Blue.se, HPCA=hpca.se,Im=Immune.se), labels =list(Blue.se$label.main, hpca.se$label.main,Immune.se$label.main),
                     method = "cluster", clusters = clusters, 
                     assay.type.test = "logcounts", assay.type.ref = "logcounts")
print("[result]pred.main...")
pred.main$labels

celltype.main = data.frame(ClusterID=rownames(pred.main), celltype.main=pred.main$labels, stringsAsFactors = F) 
print("[result]celltype.main...")
celltype.main
pbmc@meta.data$singleR.main=celltype.main[match(clusters,celltype.main$ClusterID),'celltype.main']
print("[result]write metadata fine and main...")

tem=pbmc@meta.data$singleR
tem[which(tem=="Monocytes")] <-"Monocyte"
pbmc@meta.data$singleR2=tem

tem2=pbmc@meta.data$singleR.main
tem2[which(tem2=="Monocytes")] <-"Monocyte"
pbmc@meta.data$singleR.main2=tem2

write.csv(pbmc@meta.data,file="pbmc_Treat1_2_meta.data.ref3_main_fine.csv")
saveRDS(pbmc, file="sample.combined_Treat1_2_singleR_ref3_main_fine.rds")

print("start write pdf..." )
pdf("umap_Treat1_2_singleR2_ref3_dimplot.pdf")
DimPlot(pbmc, reduction = "umap", group.by = "singleR2")
dev.off()
pdf("tsne_Treat1_2_singleR2_ref3_dimplot.pdf")
DimPlot(pbmc, reduction = "tsne", group.by = "singleR2")
dev.off()

pdf("umap_Treat1_2_Response_singleR2_ref3_dimplot.pdf")
DimPlot(pbmc, reduction = "umap", group.by = "singleR2",split.by = 'Response')
dev.off()
pdf("tsne_Treat1_2_Response_singleR2_ref3_dimplot.pdf")
DimPlot(pbmc, reduction = "tsne", group.by = "singleR2",split.by = 'Response')
dev.off()


print("start write pdf..." )

pdf("umap_Treat1_2_singleR2_ref3_dimplot_fine.pdf")
DimPlot(pbmc, reduction = "umap", group.by = "singleR.main2")
dev.off()
pdf("tsne_Treat1_2_singleR2_ref3_dimplot_fine.pdf")
DimPlot(pbmc, reduction = "tsne", group.by = "singleR.main2")
dev.off()

pdf("umap_Response_Treat1_2_singleR2_ref3_dimplot_fine.pdf")
DimPlot(pbmc, reduction = "umap", group.by = "singleR.main2",split.by = 'Response')
dev.off()
pdf("tsne_Response_Treat1_2_singleR2_ref3_dimplot_fine.pdf")
DimPlot(pbmc, reduction = "tsne", group.by = "singleR.main2",split.by = 'Response')
dev.off()





