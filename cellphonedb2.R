library(dplyr)
library(Seurat)

print(Sys.time())
print("load library complete...")
setwd("singlecell/P-MATRIX/cellphonedb")

# as_matrix <- function(mat){
#   
#   tmp <- matrix(data=0L, nrow = mat@Dim[1], ncol = mat@Dim[2])
#   
#   row_pos <- mat@i+1
#   col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])+1
#   val <- mat@x
#   
#   for (i in seq_along(val)){
#     tmp[row_pos[i],col_pos[i]] <- val[i]
#   }
#   
#   row.names(tmp) <- mat@Dimnames[[1]]
#   colnames(tmp) <- mat@Dimnames[[2]]
#   return(tmp)
# }
# 
# write_count=function(pbmc,label){
#   write.table(as_matrix(pbmc@assays$RNA@data), paste("cellphonedb/",label,"_cellphonedb_count.txt",sep=""), sep='\t', quote=F)
#   meta_data <- cbind(rownames(pbmc@meta.data), pbmc@meta.data[,'celltype', drop=F])  
#   meta_data <- as.matrix(meta_data)
#   meta_data[is.na(meta_data)] = "Unkown" #  细胞类型中不能有NA
#   write.table(meta_data, paste("cellphonedb/",label,"_cellphonedb_meta.txt",sep=""), sep='\t', quote=F, row.names=F)
# }
#meta.data$sample
# scdata<- readRDS("sample.combined_T10_others_singleR_ref3_main_fine_celltype.rds")
# scdata1=subset(scdata,Treatment!=5)
# sample=unique(scdata1@meta.data$sample)
# for(s in sample){
#   pbmc=subset(scdata1,sample==s)
#   print(Sys.time())
#   print(s)
#   print("load RDS complete...") 
#   write_count(pbmc,s)
# }
library(hash)
h = hash(
  "Undefined T cells"="Tx",
  "CD56+CD16- NK cells (NK1)"="NK1",
  "naive CD4+ T cells (T1)"="T1",
  "C56-CD16+ NK cells (NK2)"="NK2",
  "immature B cells (B2)"="B2",
  "cytotoxic CD8+ lymphocytes (CD8+ CTL) (T6)"="T6",
  "Platelets"="P1",
  "mucosal-associated invariant T (MAIT) cells (T8) "="T8",
  "CD14+ monocytes (M1)"="M1",
  "effector memory CD8+T cells (T5, CD8 Tm)"="T5",
  "Î³Î´ T cells (T9) "="T9",
  "naive B cells (B1)"="B1",
  "CD16+ monocytes (M2)"="M2",
  "regulatory CD4+ T cells (T2, Treg)"="T2",
  "CD1C+ dendritic cells (M3)"="M3",
  "Proliferating T cells (T7, Tprol) "="T7",
  "naive CD8+ T cells (T4)"="T4",
  "CD4+ T cells (T3)"="T3",
  "HSC"="HSC",
  "plasmacytoid dendritic cells (CLEC4C+ pDC) (M4)"="M4",
  "memory B cells (B3)"="B3",
  "CLEC9A+ dendritic cells (M5)"="M5")

file2line=function(label,f){
  # label="P-181"
  #f="significant_means.txt"
  # f="means.txt"
  filename=paste(label,"\\",f,sep="")
  print(Sys.time())
  print(filename)
  data=read.table(filename,header=F,sep="\t")
  name=data[1,-c(1:12)]
  all_name=NULL
  for(n in name){
    ab=strsplit(n,"|",fixed=T)
    all_name=c(all_name,ab[[1]][1],ab[[1]][2])
  }
  ctype=unique(all_name)
  
  all_line=NULL
  for(x in c(13:dim(data)[2])){
    #x=13
    ab=strsplit(data[1,x],"|",fixed=T)
    # s1=ab[[1]][1]
    # s2=ab[[1]][2]
    a=ab[[1]][1]
    b=ab[[1]][2]
    s1=ifelse(is.null(h[[a]]),a,h[[a]])
    s2=ifelse(is.null(h[[b]]),b,h[[b]])
    
    data1=data[-1,]
    data1[which(data1[,5]==""),5]=data1[which(data1[,5]==""),3]
    data1[which(data1[,6]==""),6]=data1[which(data1[,6]==""),4]
    rownum=which(data1[,x]!=""&data1[,5]!=""&data1[,6]!="")
    pp=data1[rownum,c(5,6,x)]
    line=cbind(rep(s1,dim(pp)[1]),rep(s2,dim(pp)[1]),pp)
    colnames(line)=c("sec1","sec2","p1","p2","mean")
    if(length(all_line)==0){
      all_line=line
    }else{
      all_line=rbind(all_line,line)
    }
  }
  res=all_line[all_line$mean!=0,]
  res1=cbind(res,rep(label,dim(res)[1]))
  res1
}

allh=function(h,keys){
  res=NULL
  for(x in keys){
    res=c(res,h[[x]])
  }
  res
}


###########################################################################################
#setwd("H:\\shanghai\\singlecell\\P-MATRIX\\cellphonedb")
sample=c("P-164", "P-165", "P-174", "P-175", "P-176", "P-177", "P-178", "P-179", "P-180",
         "P-181", "P-182", "P-184", "P-185", "P-186", "P-187", "P-189", "P-191", "P-192",
         "P-193", "P-197", "P-198", "P-199", "P-201", "P-202", "P-203")
all=NULL
for(s in sample){
  print(s)
  ppi=file2line(s,"pvalues.txt")
  write.csv(ppi,file=paste(s,"_interaction_pvalues.csv",sep=""))
  if(is.null(all)){
    all=ppi
  }else{
    all=rbind(all,ppi)
  }
}
# all=NULL
# for(s in sample){
#   print(s)
#   ppi=read.csv(paste(s,"_interaction.csv",sep=""),sep=",")
#   if(is.null(all)){
#     all=ppi
#   }else{
#     all=rbind(all,ppi)
#     )
#   }
# }
write.csv(all,file=paste("all","_interaction_pvalues.csv",sep=""))
print("down...")




