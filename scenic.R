# if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
# BiocManager::version()
# # If your bioconductor version is previous to 4.0, see the section bellow
# 
# ## Required
# BiocManager::install(c("AUCell", "RcisTarget"))
# BiocManager::install(c("GENIE3")) # Optional. Can be replaced by GRNBoost
# 
# ## Optional (but highly recommended):
# # To score the network on cells (i.e. run AUCell):
# BiocManager::install(c("zoo", "mixtools", "rbokeh"))
# # For various visualizations and perform t-SNEs:
# BiocManager::install(c("DT", "NMF", "ComplexHeatmap", "R2HTML", "Rtsne"))
# # To support paralell execution (not available in Windows):
# BiocManager::install(c("doMC", "doRNG"))
# # To export/visualize in http://scope.aertslab.org
# if (!requireNamespace("devtools", quietly = TRUE)) 
# install.packages("devtools")
# devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)
# 
# if (!requireNamespace("devtools", quietly = TRUE)) 
# install.packages("devtools")
# library(devtools)
# devtools::install_github("aertslab/SCENIC")
# packageVersion("SCENIC")
# 
# dbFiles <- c("https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-500bp-upstream-7species.mc9nr.feather",
#              "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-7species.mc9nr.feather")
# setwd("F:\\signalcell")
# dir.create("cisTarget_databases"); setwd("cisTarget_databases") # if needed
# for(featherURL in dbFiles)
# {
#   download.file(featherURL, destfile=basename(featherURL)) # saved in current dir
# }
#rm(list = ls()) 
library(Seurat) 
library("AUCell")
library("RcisTarget")
library("GENIE3") 
library(SCENIC)
#devtools::install_github('satijalab/seurat-data')
#SeuratData::InstallData("pbmc3k")
# library(SeuratData)
# AvailableData()
# # # InstallData("pbmc3k") #  (89.4 MB) 
# 
#   attachNamespace('pbmc3k.SeuratData')
#   data("pbmc3k")
# scdata<-readRDS("sample.combined.Rds") 
print(Sys.time())
print("load library complete...")
print(Sys.time())
setwd("singlecell/P-MATRIX")
all_pbmc<- readRDS("sample.combined_T10_others_singleR_ref3_main_fine.rds")
pbmc=subset(all_pbmc,b_is_cell=="true")
exprMat  <-  as.matrix(pbmc@assays$RNA@counts)
dim(exprMat)
exprMat[1:4,1:4] 
cellInfo <-  data.frame(pbmc@meta.data)
# colnames(cellInfo)=c('CellType', 'nGene' ,'nUMI')
# head(cellInfo)
# table(cellInfo$CellType)

### Initialize settings

# 保证cisTarget_databases 文件夹下面有下载好2个1G的文件
scenicOptions <- initializeScenic(org="hgnc", 
                                  dbDir="cisTarget_databases", nCores=4)
scenicOptions@inputDatasetInfo$cellInfo=cellInfo
#dir.create("Scenic_bcell")
#setwd("Scenic_bcell") 
saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

### Co-expression network
#genesKept <- geneFiltering(exprMat, scenicOptions)
genesKept <-geneFiltering(exprMat, scenicOptions,
                          minCountsPerGene = 3 * 0.01 *ncol(exprMat), 
                          minSamples =ncol(exprMat) * 0.01)
exprMat_filtered <- exprMat[genesKept, ]
exprMat_filtered[1:4,1:4]
dim(exprMat_filtered)
print("runCorrelation...")
runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered_log <- log2(exprMat_filtered+1) 
print("runGenie3...")
runGenie3(exprMat_filtered_log, scenicOptions)
saveRDS(scenicOptions, file="runGenie3_scenicOptions.Rds") 
### Build and score the GRN
exprMat_log <- log2(exprMat+1)
scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # Toy run settings
print("runSCENIC_1_coexNetwork2modules...")
scenicOptions <- runSCENIC_1_coexNetwork2modules(scenicOptions)
print("runSCENIC_2_createRegulons...")
scenicOptions <- runSCENIC_2_createRegulons(scenicOptions,
                                            coexMethod=c("top5perTarget"))
print("runSCENIC_3_scoreCells...")# Toy run settings
scenicOptions <- runSCENIC_3_scoreCells(scenicOptions, exprMat_log ) 
print("runSCENIC_4_aucell_binarize...")
scenicOptions <- runSCENIC_4_aucell_binarize(scenicOptions)
saveRDS(scenicOptions, file="runSCENIC_runGenie3_scenicOptions.Rds") 
print("tsneAUC...")
tsneAUC(scenicOptions, aucType="AUC") # choose settings` 
saveRDS(scenicOptions, file="tsneAUC_runSCENIC_runGenie3_scenicOptions.Rds") 
