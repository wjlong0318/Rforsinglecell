# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("monocle")
library(ggplot2)
library("DDRTree")
library("pheatmap")
library("monocle")
# setwd("F:\\signalcell\\P-MATRIX")
# scdata<-readRDS("sample.combined.Rds") 
print(Sys.time())
print("load library complete...")
print(Sys.time())
setwd("singlecell/P-MATRIX")

pbmc<- readRDS("sample.combined_T10_singleR_ref3.rds")
tem=pbmc@meta.data$singleR
tem[which(tem=="Monocytes")] <-"Monocyte"
pbmc@meta.data$singleR2=tem
print("load RDS complete...")


#scdata=subset(pbmc,b_is_cell=="true")
scdata=subset(pbmc,singleR2=="Monocyte")
# dim(scdata@assays$RNA@counts)
print("pbmc@assays$RNA@counts...")
str(pbmc@assays$RNA@counts)
print("scdata@assays$RNA@counts...")
str(scdata@assays$RNA@counts)
class(scdata@assays$RNA@counts)
# dim(scdata@meta.data)
# tem.mat=as.matrix(scdata@assays$RNA@counts)
# print("create tem.mat complete...")
# data<- as(tem.mat, 'sparseMatrix')
data=scdata@assays$RNA@counts

# print("row...")
# 
# row_t=apply(data,1,function(x) sum(x>0) )
# fivenum(row_t)
# pdf("T10_Tcell_row.pdf")
# boxplot(row_t)
# dev.off()
# 
# print("col...")
# col_t=apply(data,2,function(x) sum(x>0) )
# fivenum(col_t)
# pdf("T10_Tcell_col.pdf")
# hist(col_t)
# dev.off()


print(Sys.time())
print("create data complete...")

pd <- new('AnnotatedDataFrame', data = scdata@meta.data)
print("create pd complete...")

fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data)) 
fd <- new('AnnotatedDataFrame', data = fData)    
print("create fd complete...")

mycds<- newCellDataSet(data,
                        phenoData = pd,
                        featureData = fd, 
                        lowerDetectionLimit = 0.5, 
                        expressionFamily = negbinomial.size())
print("create CellDataset complete...")

mycds
saveRDS(mycds, file="T10_Monocyte_mycds.rds")
###filter cell and genes#########################################################
print("filter genes...")
mycds <- detectGenes(mycds, min_expr = 0.1)
expressed_genes <- row.names(subset(mycds@featureData@data,
                                    num_cells_expressed >= 100))
mycds <- mycds[expressed_genes,]
mycds
print("filter cells...")

cell_anno <- mycds@phenoData@data
valid_cells <- row.names(cell_anno[cell_anno$num_genes_expressed>2000,] )
mycds <- mycds[,valid_cells]
mycds 

saveRDS(mycds, file="T10_Monocyte_mycds_filters.rds")
##########################################################
print("estimateSizeFactors...")
 mycds <- estimateSizeFactors(mycds)
 mycds <- estimateDispersions(mycds, cores=4, relative_expr = TRUE)
 
 saveRDS(mycds, file="T10_Monocyte_mycds_filters_estimate.rds")
print("reduceDimension start...")

disp_table <- dispersionTable(mycds)
unsup_clustering_genes <- subset(disp_table, 
                                 mean_expression >= 0.1)
mycds  <- setOrderingFilter(mycds, unsup_clustering_genes$gene_id)
pdf("T10_Monocyte_plot_ordering_genes1.pdf")
plot_ordering_genes(mycds) 
dev.off()

mycds <- reduceDimension(mycds, max_components = 2, num_dim = 10,
                       reduction_method = 'tSNE', verbose = T)


mycds <- clusterCells(mycds, num_clusters = 10) 
saveRDS(mycds, file="T10_Monocyte_singleR_ref3_monocle_tSNE1.rds")

print("gene detect...")
start=Sys.time()
diff_test_res <- differentialGeneTest(mycds,
                                      fullModelFormulaStr = "~Response")
end=Sys.time()
end-start

sig_genes <- subset(diff_test_res, qval < 0.1)

write.csv(sig_genes,file="T10_Monocyte_sig_gene_monocle21.csv")

ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
mycds <- setOrderingFilter(mycds, ordering_genes)
pdf("T10_Monocyte_plot_sig_genes1.pdf")
plot_ordering_genes(mycds)
dev.off()

 mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
 mycds <- orderCells(mycds)
 
 saveRDS(mycds, file="T10_Monocyte_singleR_ref3_monocle_tSNE_DDR1.rds")
 print("save mycds complete...")
 
 #State轨迹分布�?
 pdf("T10_Monocyte_trajectory_state1.pdf")
 plot_cell_trajectory(mycds, color_by = "State")
 dev.off()
 
 pdf("T10_Monocyte_trajectory_state_split_Response1.pdf")
 plot_cell_trajectory(mycds, color_by = "State") + facet_wrap(~Response, nrow = 1)
 dev.off()

 pdf("T10_Monocyte_trajectory_singleR21.pdf")
 plot_cell_trajectory(mycds, color_by = "singleR2")
 dev.off()
 
 pdf("T10_Monocyte_trajectory_singleR2_split_Response1.pdf")
 plot_cell_trajectory(mycds, color_by = "singleR2") + facet_wrap(~Response, nrow = 1)
 dev.off()

 pdf("T10_Monocyte_trajectory_b_chain.pdf")
 plot_cell_trajectory(mycds, color_by = "b_chain")
 dev.off()
 
 pdf("T10_Monocyte_trajectory_b_chain_split_Response1.pdf")
 plot_cell_trajectory(mycds, color_by = "b_chain") + facet_wrap(~Response, nrow = 1)
 dev.off()
 
 pdf("T10_Monocyte_trajectory_Pseudotime1.pdf")
 plot3 <- plot_cell_trajectory(mycds, color_by = "Pseudotime")
 dev.off()
 
 pdf("T10_Monocyte_trajectory_Pseudotime_split_Response1.pdf.pdf")
 plot3 <- plot_cell_trajectory(mycds, color_by = "Pseudotime") + facet_wrap(~Response, nrow = 1)
 dev.off()
 
 
 pdf("T10_Monocyte_gene_Pseudotime1.pdf")
 s.genes <- c("ITGB1","CCR7","KLRB1","GNLY")

 p1 <- plot_genes_jitter(mycds[s.genes,], grouping = "Response", color_by = "State")

 p2 <- plot_genes_violin(mycds[s.genes,], grouping = "Response", color_by = "State")

 p3 <- plot_genes_in_pseudotime(mycds[s.genes,], color_by = "State")
 
 p1|p2|p3
 dev.off()
 