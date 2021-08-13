library('scCATCH')
# install.packages(pkgs = 'devtools')
# devtools::install_github('ZJUFanLab/scCATCH')
#load("mouse_kidney_203_Seurat.RData")

setwd("singlecell/P-MATRIX")


print("load library complete...")
pbmc<- readRDS("IntegrateSample_tsne_umap_T10.rds")
print("load data complete...")
clu_markers <- findmarkergenes(object = pbmc,
                               species = 'Human',
                               cluster = 'All',
                               match_CellMatch = FALSE,
                               cancer = NULL,
                               tissue = NULL,
                               cell_min_pct = 0.25,
                               logfc = 0.25,
                               pvalue = 0.05)
print("start scCATCH..." )
clu_markers$clu_markers
clu_ann <- scCATCH(object = clu_markers$clu_markers,
                   species = 'Human',
                   cancer = NULL,
                   tissue = c('Blood','Peripheral blood','Bone marrow'))

new.cluster.ids <- clu_ann$cell_type
names(new.cluster.ids) <- clu_ann$cluster
pbmc <- RenameIdents(pbmc, new.cluster.ids)#进行重新命名

print("start write pdf..." )
pdf("umap_scCATH_dimplot.pdf")
DimPlot(pbmc, reduction = "umap", label = TRUE, pt.size = 0.5) + NoLegend()
dev.off()