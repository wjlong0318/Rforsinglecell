library(Seurat)
library(SeuratData)
library(tidyverse)
library(CellChat)
library(NMF)
library(ggalluvial)
library(patchwork)
library(ggplot2)
library(svglite)
options(stringsAsFactors = FALSE)
print(Sys.time())
print("load library complete...")

setwd("singlecell/P-MATRIX")

dir.create("CellChat")
print("load RDS complete...")
pbmc<- readRDS("sample.combined_T10_others_singleR_ref3_main_fine_celltype.rds")
#write.csv(pbmc@meta.data,file="sample.combined_T10_others_singleR_ref3_main_fine_celltype.csv")
scdata1=subset(pbmc,Treatment!=5)
scdata2=subset(scdata1,Response!="PD")

cellchat <- createCellChat(object=scdata2,group.by = "celltype")
#cellchat <- createCellChat(pbmc3k.final@assays$RNA@data, meta = pbmc3k.final@meta.data, group.by = "cell_type")
print("seurat to cellchat down... ")
cellchat
summary(cellchat)
str(cellchat)
levels(cellchat@idents)
#cellchat <- setIdent(cellchat, ident.use = "celltype")
groupSize <- as.numeric(table(cellchat@idents))  
print("groupSize...")
groupSize
# [1] 711 480 472 344 27

CellChatDB <- CellChatDB.human
#导入小鼠是CellChatDB <- CellChatDB.mouse
print("cellchatdb...")
str(CellChatDB) #查看数据库信息
#包含interaction、complex、cofactor和geneInfo这4个dataframe
print("colnames(CellChatDB$interaction) ...")
colnames(CellChatDB$interaction)
print("CellChatDB$interaction[1:4,1:4] ...")
CellChatDB$interaction[1:4,1:4]
print("head(CellChatDB$cofactor) ...")
head(CellChatDB$cofactor)
print("head(CellChatDB$complex) ...")
head(CellChatDB$complex)
print("head(CellChatDB$geneInfo) ...")
head(CellChatDB$geneInfo)
print("showDatabaseCategory(CellChatDB) ...")
showDatabaseCategory(CellChatDB)

unique(CellChatDB$interaction$annotation)#查看可以选择的侧面，也就是上图左中的三种
#选择"Secreted Signaling"进行后续细胞互作分析
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB
# set the used database in the object
cellchat@DB <- CellChatDB.use # set the used database in the object

## 在矩阵的所有的基因中提取signaling gene，结果保存在data.signaling(13714个基因，过滤完只有270个）
cellchat <- subsetData(cellchat)
future::plan("multiprocess", workers = 4)
#相当于Seurat的FindMarkers，找每个细胞群中高表达的配体受体
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat) #Identify over-expressed ligand-receptor interactions (pairs) within the used CellChatDB
#上一步运行的结果储存在cellchat@LR$LRsig
cellchat <- projectData(cellchat, PPI.human) 

#根据表达值推测细胞互作的概率（cellphonedb是用平均表达值代表互作强度）。
cellchat <- computeCommunProb(cellchat, raw.use = FALSE, population.size = TRUE) #如果不想用上一步PPI矫正的结果，raw.use = TRUE即可。
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
write.csv(df.net, "PD_net_lr.csv")

cellchat <- computeCommunProbPathway(cellchat)
saveRDS(cellchat, file=paste("cellchat/PD_",ctype,"_cellchat.rds",sep=""))
df.netp <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.netp, "PD_net_pathway.csv")
#找到配体受体关系后，projectData将配体受体对的表达值投射到PPI上，来对@data.signaling中的表达值进行校正。结果保存在@data.project

cellchat <- aggregateNet(cellchat)
#计算每种细胞各有多少个
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
pdf(paste("cellchat/PD_",ctype,"Number_interactions.pdf",sep=""))
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Number of interactions")
dev.off()
pdf(paste("cellchat/PD_",ctype,"Interaction_weights.pdf",sep=""))
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, 
                 label.edge= F, title.name = "Interaction weights/strength")
dev.off()