library(Seurat)
library(celldex)
library(SingleR)
library(dplyr)
library(ggplot2)
library(reshape2)
library(tibble)
print("load library complete...")



setwd("singlecell/P-MATRIX")

scdata0<- readRDS("sample.combined_T10_others_singleR_ref3_main_fine_celltype.rds")
scdata1=subset(scdata0,Treatment!=5)
scdata=subset(scdata1,Response!="Na")
pbmc=scdata
# pbmc<- readRDS("sample.combined_T10_others_singleR_ref3_main_fine.rds")
# 
# print("load data complete...")
# tem=pbmc@meta.data$singleR
# tem[which(tem=="Monocytes")] <-"Monocyte"
# pbmc@meta.data$singleR2=tem

#tem2=pbmc@meta.data$singleR.main
#tem2[which(tem2=="Monocytes")] <-"Monocyte"
#pbmc@meta.data$singleR.main2=tem2
# 
# print("start write pdf..." )
# pdf("umap_singleR2_ref3_dimplot.pdf")
# DimPlot(pbmc, reduction = "umap", group.by = "singleR2")
# dev.off()
# pdf("tsne_singleR2_ref3_dimplot.pdf")
# DimPlot(pbmc, reduction = "tsne", group.by = "singleR2")
# dev.off()
# 
# pdf("umap_Response_singleR2_ref3_dimplot.pdf")
# DimPlot(pbmc, reduction = "umap", group.by = "singleR2",split.by = 'Response')
# dev.off()
# pdf("tsne_Response_singleR2_ref3_dimplot.pdf")
# DimPlot(pbmc, reduction = "tsne", group.by = "singleR2",split.by = 'Response')
# dev.off()


# pdf("tchain_Response_umap_dimplot.pdf")
# DimPlot(pbmc , reduction = "umap", group.by = "t_chain",split.by = 'Response')
# dev.off()
# 
# pdf("bchain_Response_umap_dimplot.pdf")
# DimPlot(pbmc , reduction = "umap", group.by = "b_chain",split.by = 'Response')
# dev.off()
# 
# pdf("tchain_Response_tsne_dimplot.pdf")
# DimPlot(pbmc , reduction = "tsne", group.by = "t_chain",split.by = 'Response')
# dev.off()
# 
# pdf("bchain_Response_tsne_dimplot.pdf")
# DimPlot(pbmc , reduction = "tsne", group.by = "b_chain",split.by = 'Response')
# dev.off()
# # 
# pdf("cluster_sample_tsne_size_dimplot.pdf")
# DimPlot(pbmc , reduction = "tsne", group.by = "singleR2", split.by = 'sample',pt.size=5,ncol=4)
# dev.off()
# 
# pdf("cluster_sample_umap_size_dimplot.pdf")
# DimPlot(pbmc , reduction = "umap", group.by = "singleR2", split.by = 'sample',pt.size=5,ncol=4)
# dev.off()

# markers = c("MS4A1",
#             "CD3D",
#             "KLRF1",
#             "KLRD1",
#             "CD14",
#             "CD8A",
#             "PF4",
#             "CD68"
# )
# #rev(markers)
# pdf("heatmap_markers_T10.pdf")
# DoHeatmap(pbmc,features = rev(markers),label = F,group.by ="singleR")
# dev.off()
# 
# 
#  markers_genes <- FindAllMarkers(alldata,logfc.threshold = 0.2, test.use = "wilcox", assay = "RNA",
#                                  min.pct = 0.1, min.diff.pct = 0.2, only.pos = TRUE, max.cells.per.ident = 50)
# top25 <- markers_genes %>% group_by(cluster) %>% top_n(-25, p_val_adj)
# mypar(2, 5, mar = c(4, 6, 3, 1))
# for (i in unique(top25$cluster)) {
#   barplot(sort(setNames(top25$avg_logFC, top25$gene)[top25$cluster == i], F), horiz = T, 
#           las = 1, main = paste0(i, " vs. rest"), border = "white", yaxs = "i")
#   abline(v = c(0, 0.25), lty = c(1, 2))
# }
# top5 <- markers_genes %>% group_by(cluster) %>% top_n(-5, p_val_adj)
# 
# # create a scale.data slot for the selected genes
# alldata <- ScaleData(alldata, features = as.character(unique(top5$gene)), assay = "RNA")
# DoHeatmap(alldata, features = as.character(unique(top5$gene)), group.by = sel.clust, assay = "RNA")

 annotations <- read.csv("annotation.csv")
# # markers_all <- FindAllMarkers(object = pbmc, 
# #                               only.pos = TRUE,
# #                               logfc.threshold = 0.25)
# # write.csv(markers_all,file="markers_all.csv")
# # print("write markers_all  complete...")
# # 
# get_conserved <- function(seurat_integrated,cluster,annotations){
#   FindMarkers(seurat_integrated,
#                        ident.1 = cluster,
#                       min.pct = 0.25,
#                        only.pos = TRUE) %>%
#     tibble::rownames_to_column(var = "gene") %>%
#     left_join(y = unique(annotations[, c("gene_name", "description")]),
#               by = c("gene" = "gene_name")) %>%
#     cbind(cluster_id = cluster, .)
# }
# for(n in c(0:58)){
# diffgenes=get_conserved(pbmc,n,annotations)
# #n=1
# filename=paste("T10_Tothers_cluster_",n,"_diffgene.csv",sep="")
# write.csv(diffgenes,file=filename)
# print("write complete...")
# print(filename)
# }
# genes_C41=get_conserved(pbmc,41,annotations)
# write.csv(genes_C41,file="gene_C41.csv")
# print("write gene_C41 complete...")

# markers = c("CD3D",
#             "CCR7",
#             "CD8A",
#             "MS4A1",
#             "KLRF1",
#             "KLRD1",
#             "PF4",
#             "CD14",
#             "CD68"
#             
# )
 markers=c("TNF",
           "HLA-E",
           "HLA-F")
print("write DotPlot_markers_T10 complete...")

# pdf("DotPlot_markers_T10_Tothers_cluster.pdf")
# pp = DotPlot(pbmc, features = rev(markers),cols = c('white','#F8766D'),dot.scale =5,group.by ="seurat_clusters" ) + RotatedAxis()
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
pdf("DotPlot_Tall_TNF_celltype.pdf")
pp = DotPlot(pbmc, features = markers,cols = c('white','#F8766D'),dot.scale =5,group.by ="celltype") + RotatedAxis()
pp = pp +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8)) +
  labs(x='',y='') +
  guides(color = guide_colorbar(title = 'Scale expression'),
         size = guide_legend(title = 'Percent expressed')) +
  theme(axis.line = element_line(size = 0.6))
pp
dev.off()
pdf("DotPlot_Tall_TNF_sample.pdf")
pp = DotPlot(pbmc, features = markers,cols = c('white','#F8766D'),dot.scale =5,group.by ="sample") + RotatedAxis()
pp = pp +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8)) +
  labs(x='',y='') +
  guides(color = guide_colorbar(title = 'Scale expression'),
         size = guide_legend(title = 'Percent expressed')) +
  theme(axis.line = element_line(size = 0.6))
pp
dev.off()

pdf("FeaturePlot_Tall_TNF.pdf")
FeaturePlot(pbmc , features = markers,cols = c('gray','#F8766D'),reduction = "tsne")
dev.off()
pdf("FeaturePlot_Tall_TNF_celltype.pdf")
FeaturePlot(pbmc , features = markers,cols = c('gray','#F8766D'),reduction = "tsne",split.by ="celltype")
dev.off()
pdf("FeaturePlot_Tall_TNF_sample.pdf")
FeaturePlot(pbmc , features = markers,cols = c('gray','#F8766D'),reduction = "tsne",split.by ="sample")
dev.off()
pdf("FeaturePlot_Tall_TNF_response.pdf")
FeaturePlot(pbmc , features = markers,cols = c('gray','#F8766D'),reduction = "tsne",split.by ="response")
dev.off()


##############################################################################
# par(mar = c(4, 8, 2, 1))
# pdf("boxplot_most_express_T10.pdf")
# C <- pbmc@assays$RNA@counts
# C <- Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
# most_expressed <- order(apply(C, 1, median), decreasing = T)[20:1]
# 
# boxplot(as.matrix(t(C[most_expressed, ])), cex = 0.1, las = 1, xlab = "% total count per cell", 
#         col = (scales::hue_pal())(20)[20:1], horizontal = TRUE)
# 
# dev.off()
##########################################################################
# cell_selection=subset(pbmc,b_is_cell=="true")
# 
# #cell_selection <- subset(alldata, cells = colnames(alldata)[alldata@meta.data[, sel.clust] ==  2])
# cell_selection <- SetIdent(cell_selection, value = "Response")

# Compute differentiall expression
# DGE_cell_selection <- FindAllMarkers(cell_selection, logfc.threshold = 0.2, test.use = "wilcox", 
#                                      min.pct = 0.1, min.diff.pct = 0.2, only.pos = TRUE, max.cells.per.ident = 50, 
#                                      assay = "RNA")
# write.csv(DGE_cell_selection,file="DGE_b_cell.csv")
