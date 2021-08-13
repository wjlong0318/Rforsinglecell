library(Seurat)
library(celldex)
library(SingleR)
library(dplyr)
library(ggplot2)
library(reshape2)
library(tibble)
print("load library complete...")



setwd("singlecell/P-MATRIX")

pbmc<- readRDS("sample.combined_T10_others_singleR_ref3_main_fine.rds")

# clusters_files=c(0, 4, 8, 9, 11, 12, 13,
#                  22, 24, 26, 31, 32, 36,
#                  37, 44, 46, 48, 50, 51, 56)
clusters_files=c(0, 3, 8, 48, 56)

scdata=subset(pbmc,seurat_clusters %in% clusters_files)

markers_genes=NULL
for(x in clusters_files){
  filename=paste("T10_Tothers_cluster_",x,"_diffgene.csv",sep="")
  print(filename)
  markers_tem=read.csv(filename,header=T,sep=",")
  if(length(markers_genes)==0){
    markers_genes=markers_tem
  }else{
    markers_genes=rbind(markers_genes,markers_tem)    
  }
  
}
dim(markers_genes)
names(markers_genes)
top5 <- markers_genes %>% group_by(cluster_id) %>% top_n(5, avg_log2FC)



# create a scale.data slot for the selected genes
scdata <- ScaleData(scdata, features = as.character(unique(top5$gene)), assay = "RNA")
pdf('heatmap_monocytes_group_cluster.pdf', width = 14, height = 10)
DoHeatmap(scdata, features = as.character(unique(top5$gene)), group.by = "seurat_clusters",
          size = 3, angle = -50, hjust=0.8,assay = "RNA")+
  scale_fill_gradientn(colors = c("dodgerblue", "white", "firebrick1"))
dev.off()

pdf('heatmap_monocytes.pdf', width = 14, height = 10)
DoHeatmap(scdata, features = as.character(unique(top5$gene)),
          size = 3, angle = -50, hjust=0.8,assay = "RNA")+
  scale_fill_gradientn(colors = c("dodgerblue", "white", "firebrick1"))
dev.off()