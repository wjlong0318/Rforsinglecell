# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("monocle")
library(ggplot2)
library("DDRTree")
library("pheatmap")
library("monocle")
library(dplyr)
library(Seurat)
# setwd("F:\\signalcell\\P-MATRIX")
# scdata<-readRDS("sample.combined.Rds") 
print(Sys.time())
print("load library complete...")
print(Sys.time())
setwd("singlecell/P-MATRIX")


print("load RDS complete...")
pbmc<- readRDS("sample.combined_T10_others_singleR_ref3_main_fine_celltype.rds")
#write.csv(pbmc@meta.data,file="sample.combined_T10_others_singleR_ref3_main_fine_celltype.csv")
scdata1=subset(pbmc,Treatment!=5)
scdata2=subset(scdata1,Response!="Na")
dir.create("time")
names=unique(scdata2@meta.data$celltype)
#unique(meta.data$celltype)
#1:#9#(10:15):16:21
for(n in names[c(1,4)]){
  print("start...")
  print(n)
 # scdata=subset(scdata2,celltype==n)
#  clusters_files=unique(scdata@meta.data$seurat_clusters)
  ctype=n
  
  # 
  # #monotype 
  # #ctype="monocytes"
  # # clusters_files=c(0, 4, 8, 9, 11, 12, 13,
  # #                  22, 24, 26, 31, 32, 36,
  # #                  37, 44, 46, 50, 51, 56)
  # # #Bcell cluster
  # # clusters_files=c(7, 14, 34, 52, 58)
  # 
  # #scdata=subset(scdata2,singleR.main2=="Bcell")
  # #scdata=subset(scdata2,b_is_cell=="true")
  # #scdata=subset(scdata2,seurat_clusters %in% clusters_files)
  # 
  # # print("pbmc@assays$RNA@counts...")
  # # str(pbmc@assays$RNA@counts)
  # # print("scdata@assays$RNA@counts...")
  # # str(scdata@assays$RNA@counts)
  # # class(scdata@assays$RNA@counts)
  # 
  # data=scdata@assays$RNA@counts
  # 
  # # print("row...")
  # #
  # # row_t=apply(data,1,function(x) sum(x>0) )
  # # fivenum(row_t)
  # # pdf("Tall_Bcell_row.pdf")
  # # boxplot(row_t)
  # # dev.off()
  # #
  # # print("col...")
  # # col_t=apply(data,2,function(x) sum(x>0) )
  # # fivenum(col_t)
  # # pdf("Tall_Bcell_col.pdf")
  # # hist(col_t)
  # # dev.off()
  # 
  # 
  # print(Sys.time())
  # print("create data complete...")
  # 
  # pd <- new('AnnotatedDataFrame', data = scdata@meta.data)
  # print("create pd complete...")
  # 
  # fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
  # fd <- new('AnnotatedDataFrame', data = fData)
  # print("create fd complete...")
  # 
  # mycds<- newCellDataSet(data,
  #                        phenoData = pd,
  #                        featureData = fd,
  #                        lowerDetectionLimit = 0.5,
  #                        expressionFamily = negbinomial.size())
  # print("create CellDataset complete...")
  # 
  # mycds
  # ###filter cell and genes#########################################################
  # print("filter genes...")
  # mycds <- detectGenes(mycds, min_expr = 0.1)
  # expressed_genes <- row.names(subset(mycds@featureData@data,
  #                                     num_cells_expressed >= 100))
  # mycds <- mycds[expressed_genes,]
  # mycds
  # print("filter cells...")
  # 
  # cell_anno <- mycds@phenoData@data
  # valid_cells <- row.names(cell_anno[cell_anno$num_genes_expressed>1000,] )
  # #valid_cells <- row.names(cell_anno[cell_anno$num_genes_expressed>4000,] )
  # mycds <- mycds[,valid_cells]
  # mycds
  # 
  # saveRDS(mycds, file=paste("time/Tall_",ctype,"_mycds_filters_100_1000.rds",sep=""))
  # ##########################################################
  # 
  # #mycds<- readRDS("Tall_Bcell_mycds_filters.rds")
  # print("load Tall_Bcell_mycds_filters.rds complete...")
  # 
  # print("estimateSizeFactors...")
  # mycds <- estimateSizeFactors(mycds)
  # mycds <- estimateDispersions(mycds, cores=4, relative_expr = TRUE)
  # print("reduceDimension start...")
  # 
  # disp_table <- dispersionTable(mycds)
  # unsup_genes <- subset(disp_table, 
  #                       mean_expression >= 0.1)
  # mycds  <- setOrderingFilter(mycds, unsup_genes$gene_id)
  # 
  # 
  # # clusters_files=c(0, 3, 8, 48, 56)
  # 
  # markers_genes=NULL
  # for(x in clusters_files){
  #   filename=paste("T10_Tothers_cluster_",x,"_diffgene.csv",sep="")
  #   print(filename)
  #   markers_tem=read.csv(filename,header=T,sep=",")
  #   if(length(markers_genes)==0){
  #     markers_genes=markers_tem
  #   }else{
  #     markers_genes=rbind(markers_genes,markers_tem)    
  #   }
  #   
  # }
  # 
  # top50 <- markers_genes %>% group_by(cluster_id) %>% top_n(100, avg_log2FC)
  # ordering_genes=unique(top50$gene)
  # print("unique genes...")
  # mycds <- setOrderingFilter(mycds, ordering_genes)
  # mycds
  # pdf(paste("time/Tall_",ctype,"_plot_ordering_genes2.pdf",sep=""))
  # p=plot_ordering_genes(mycds) 
  # print(p)
  # dev.off()
  # # print("gene detect...")
  # # start=Sys.time()
  # # diff_test_res <- differentialGeneTest(mycds,
  # #                                       fullModelFormulaStr = "~seurat_clusters")
  # # end=Sys.time()
  # # end-start
  # # 
  # # sig_genes <- subset(diff_test_res, qval < 0.1)
  # # 
  # # write.csv(sig_genes,file="Tall_Bcells_sig_gene_monocle.csv")
  # # 
  # # ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
  # # mycds <- setOrderingFilter(mycds, ordering_genes)
  # 
  # 
  # print("reduceDimension...")
  # 
  # mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
  # print("orderCells...")
  # mycds <- orderCells(mycds)
  # 
  # saveRDS(mycds, file=paste("time/Tall_",ctype,"_singleR_ref3_monocle_tSNE_DDR1.rds",sep=""))
  # print("save mycds complete...")
  # # 
  # setwd("H:\\shanghai\\singlecell\\P-MATRIX\\")
  filename=paste("time/Tall_",ctype,"_singleR_ref3_monocle_tSNE_DDR1.rds",sep="")
  mycds<- readRDS(filename)
  print("load mycds complete...")
  
  write.csv(phenoData(mycds)@data,file=paste("time/Tall_",ctype,"_fdata.csv",sep=""))
  
  pdf(paste("time/Tall_",ctype,"_trajectory_state.pdf",sep=""))
  p=plot_cell_trajectory(mycds, color_by = "State")
  print(p)
  dev.off()
  
  pdf(paste("time/Tall_",ctype,"_trajectory_state_split_Response.pdf",sep=""))
  p=plot_cell_trajectory(mycds, color_by = "State") + facet_wrap(~Response, nrow = 1)
  print(p)
  dev.off()
  
  pdf(paste("time/Tall_",ctype,"_trajectory_celltype.pdf",sep=""))
  p=plot_cell_trajectory(mycds, color_by = "celltype")
  print(p)
  dev.off()
  
  pdf(paste("time/Tall_",ctype,"_trajectory_singleR2_split_celltype.pdf",sep=""))
  p=plot_cell_trajectory(mycds, color_by = "celltype") + facet_wrap(~Response, nrow = 1)
  print(p)
  dev.off()
  
  pdf(paste("time/Tall_",ctype,"_trajectory_Pseudotime.pdf",sep=""))
  p=plot_cell_trajectory(mycds, color_by = "Pseudotime")
  print(p)
  dev.off()
  
  pdf(paste("time/Tall_",ctype,"_trajectory_Pseudotime_split_Response.pdf.pdf",sep=""))
  p=plot_cell_trajectory(mycds, color_by = "Pseudotime") + facet_wrap(~Response, nrow = 1)
  print(p)
  dev.off()
  
}

# pdf(paste("Tall_","Tall_Bcell_gene_Pseudotime.pdf",sep=""))
# s.genes <- c("CD19",
#              "MS4A1",
#              "IGHD",
#              "IGHM",
#              "IL4R",
#              "TCL1A",
#              "CD27",
#              "CD38"
# )
# 
# p1 <- plot_genes_jitter(mycds[s.genes,], grouping = "Response", color_by = "State")
# 
# p2 <- plot_genes_violin(mycds[s.genes,], grouping = "Response", color_by = "State")
# 
# p3 <- plot_genes_in_pseudotime(mycds[s.genes,], color_by = "State")
# 
# p1|p2|p3
# dev.off()
