#R 数据可视化 —— circlize 简单介绍
library(circlize)
library(ggplot2)
setwd("H:\\shanghai\\singlecell\\P-MATRIX\\cellphonedb")
# filename="PD\\significant_means.txt"
# PD=read.table(filename,header=F,sep="\t")

library(stringr)
str_extract("CD56+CD16- NK cells (NK1)", "\\(.*\\d\\)")
type="Platelets"
type="CD56+CD16- NK cells (NK1)"
type="plasmacytoid dendritic cells (CLEC4C+ pDC) (M4)"
ifelse(is.na(str_extract(type, "\\(.*\\d\\)")),type,
       str_extract(type, "\\(\\w*\\d\\)"))
h[[type]]
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
file2count=function(label){
  #label="PD"
  filename=paste(label,"\\means.txt",sep="")
 # filename=paste(label,"\\significant_means.txt",sep="")
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
    a=ab[[1]][1]
    b=ab[[1]][2]
    # s1=ifelse(is.na(str_extract(a, "\\(.*\\d\\)")),a,
    #           str_extract(a, "\\(\\w*\\d\\)"))
    # s2=ifelse(is.na(str_extract(b, "\\(.*\\d\\)")),b,
    #           str_extract(b, "\\(\\w*\\d\\)"))
    s1=ifelse(is.null(h[[a]]),"T9",h[[a]])
    s2=ifelse(is.null(h[[b]]),"T9",h[[b]])
    data1=data[-1,]
    rownum=which(data1[,x]!="0.0"&data1[,5]!=""&data1[,6]!="")
  #  pp=data1[rownum,c(5,6,x)]
    line=cbind(s1,s2,length(rownum))
    colnames(line)=c("sec1","sec2","count")
    if(length(all_line)==0){
      all_line=line
    }else{
      all_line=rbind(all_line,line)
    }
  }
  all_line
}
PD=file2count("PD")
#head(PD)
SD=file2count("SD")
PR=file2count("PR")
line_merge=rbind(cbind(PD,label=rep("PD",dim(PD)[1])),
               cbind(PR,label=rep("PR",dim(PR)[1])),
               cbind(SD,label=rep("SD",dim(SD)[1]))
)
#write.csv(line_merge,file="line_merge.csv")
dim(line_merge)

line_merge[1:20,1:4]
line_merge=data.frame(line_merge)
line_merge$label_color=ifelse(line_merge$label=="PD","#93F49FFF",ifelse(line_merge$label=="PD","#00723AFF","#BD151AFF"))
library(ggplot2)
line_merge1=line_merge[line_merge$sec1=="B1)",]
#table(line_merge1$label)
ggplot(data = line_merge,aes(x = sec2, y = as.numeric(count))) + 
  geom_bar(stat = 'identity',aes(fill=line_merge$label),position = position_dodge(0.9))+
  facet_wrap(.~sec1,scales="free_y",ncol=4)+
  theme_bw() +
  scale_y_continuous(name = "Interaction pairs (count)") +
  scale_x_discrete(name = "Cell type") +
  
  theme(plot.title = element_text(size = 11, face =  "bold"),
        text = element_text(size = 8, color="black",face="bold",),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 8, color="black",face="bold",angle = 45, hjust = 0.5, vjust = 0.5),
        axis.text.y=element_text(size = 10, color="black",face="bold"))

PD=data.frame(PD)
vec2max=function(v){
  #v=PD
  name=unique(v[,1])
  line=NULL
  for(n in name){
    #n="NK2"
    if(is.null(line)){
      line=as.numeric(v[v[,1]==n,3])
      names(line)=v[v[,1]==n,2]
    }else{
      line=cbind(line,v[v[,1]==n,3])
    }
    
  }
  colnames(line)=name
  line
}
max2vec=function(v){
  #v=dcols
  class(v)
  name=col.names(v)
  line=NULL
  for(n in name){
    #n="NK2"
    if(is.null(line)){
      line=as.numeric(v[v[,1]==n,3])
      names(line)=v[v[,1]==n,2]
    }else{
      line=cbind(line,v[v[,1]==n,3])
    }
    
  }
  colnames(line)=name
  line
}
line_max=apply(vec2max(line_merge),2,as.numeric)
colnames(line_max)=unique(line_merge[,1])

PD_max=apply(vec2max(PD),2,as.numeric)
rownames(PD_max)=unique(PD[,1])
colnames(PD_max)=paste(colnames(PD_max),"_PD",sep="")
PR_max=apply(vec2max(PR),2,as.numeric)
rownames(PR_max)=unique(PR[,1])
colnames(PR_max)=paste(colnames(PR_max),"_PR",sep="")
SD_max=apply(vec2max(SD),2,as.numeric)
rownames(SD_max)=unique(SD[,1])
colnames(SD_max)=paste(colnames(SD_max),"_SD",sep="")
all_merge=cbind(PD_max,PR_max,SD_max)
all_merge_ordered=all_merge[,order(colnames(all_merge),decreasing=F)]
library(pheatmap)
p=pheatmap(all_merge_ordered,cluster_rows=T,cluster_cols=T)
dd=NULL
for(x in c(0:21)){
  #x=1
  id1=x*3+1
  id2=x*3+2
  id3=x*3+3
  cellname1=colnames(all_merge_ordered[,id1:id2])
  dcols1 = dist(t(all_merge_ordered[,id1:id2]), method = "euclidean")
    cellname2=colnames(all_merge_ordered[,c(id1,id3)])
  dcols2 = dist(t(all_merge_ordered[,c(id1,id3)]), method = "euclidean")
  
  cellname3=colnames(all_merge_ordered[,id2:id3])
  dcols3 = dist(t(all_merge_ordered[,id2:id3]), method = "euclidean")
  min_d=min(dcols1[1:1],dcols2[1:1],dcols3[1:1])
  d1f=abs(dcols1[1:1]-min_d)
  d2f=abs(dcols2[1:1]-min_d)
  d3f=abs(dcols3[1:1]-min_d)
  d1=c(cellname1,dcols1[1:1],d1f)
  d2=c(cellname2,dcols2[1:1],d2f)
  d3=c(cellname3,dcols3[1:1],d3f)
  if(is.null(dd)){dd=rbind(d1,d2,d3)
  }else{
    dd=rbind(dd,d1,d2,d3)  
  }
}
write(dd,file="dist_euclidean.csv")

line_merge=data.frame(line_merge)
p=ggplot(line_merge,aes(sec1,sec2,colour=as.numeric(count)))+
  geom_point(size=2,shape=16)+
 # scale_color_gradient(low = "cyan",high = "red")
  scale_color_gradient(low="#5098EF",high="#F21952")+
# facet_wrap(.~label)+ 
  theme_bw()+
  # scale_y_discrete(name = "Interaction pairs") +
  # scale_x_discrete(name = "Cell type") +
  scale_size_continuous(range=c(0,3))+
  theme(plot.title = element_text(size = 11, face =  "bold"),
        text = element_text(size = 8, color="black",face="bold",),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 8, color="black",face="bold",angle = 45, hjust = 0.5, vjust = 0.5),
        axis.text.y=element_text(size = 8, color="black",face="bold"))


###############################################################################3
gene2id=function(col){
  all_genes=unique(c(data1[,5],data1[,6]))
  num=NULL
  for(x in col){
    num=c(num,which(all_genes==x))  
  }
  num
}
file2line=function(label,f){
 # label="P-181"
  #f="significant_means.txt"
 # f="means.txt"
  filename=paste(label,"\\",f,sep="")
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
    s1=ifelse(is.null(h[[a]]),"T9",h[[a]])
    s2=ifelse(is.null(h[[b]]),"T9",h[[b]])
    
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
  # p1num=gene2id(all_line[,3])
  # p2num=gene2id(all_line[,4])
  # all_line=cbind(all_line,p1num=p1num,p2num=p2num)
  all_line
}
PD=file2line("PD","significant_means.txt")
SD=file2line("SD","significant_means.txt")
PR=file2line("PR","significant_means.txt")
# PD=file2line("PD","means.txt")
# SD=file2line("SD","means.txt")
# PR=file2line("PR","means.txt")
line_merge=rbind(cbind(PD,label=rep("PD",dim(PD)[1])),
                 cbind(PR,label=rep("PR",dim(PR)[1])),
                 cbind(SD,label=rep("SD",dim(SD)[1]))
)

line_merge[1:4,1:8]
line_merge[line_merge$sec1=='P1'&line_merge$sec2=="NK1",]
all_genes=unique(c(data1[,5],data1[,6]))
#########################################################################
n <- length(all_genes)
allh=function(h,keys){
  res=NULL
  for(x in keys){
    res=c(res,h[[x]])
  }
  res
}
df=data.frame(
  sec=rep(allh(h,ctype),each=n),x=rep(c(1:n),length(ctype))
)

# df <- data.frame(
#   sectors = sample(letters[1:8], n, replace = TRUE),
#   x = rnorm(n), y = runif(n)
# )
library(circlize)
col22=rand_color(22)

circos.par("track.height" = 0.05)
circos.initialize(df$sec, x = df$x)
circos.track(df$sec, y = df$x,
             panel.fun = function(x, y) {
               circos.text(CELL_META$xcenter, 
                           CELL_META$cell.ylim[2] + mm_y(7), 
                           CELL_META$sector.index)
             #  circos.axis(labels.cex = 0.2)
             })
#bgcol <- rep(c("#fb8072", "#80b1d3"), 4)

circos.trackPlotRegion(ylim = c(0, 1), track.height = 0.05, 
                       bg.col = col22, bg.border = NA) 
#circos.trackHist(df$sectors, df$x, bin.size = 0.2, bg.col = bgcol, col = NA)
#
all_line1=line_merge[line_merge$sec1=='P1'&line_merge$label=="PD",]
all=dim(all_line1)[1]
#all_line[1:4,1:6]
#x=2

for(x in c(1:all)){
circos.link(all_line1[x,1], all_line1[x,6], all_line1[x,2], all_line1[x,7], h = 0.4,
            col = col22[which(allh(h,ctype)==all_line1[x,2])], border = NA)
 }
circos.clear()
###newdata--dotplot##############################################################
line_merge$p1p2=paste(line_merge$p1,line_merge$p2,sep="-")

name=unique(new_data$s1s2)
for(x in name){
#pdf(paste(x,"_dotplot_sigmean.pdf",sep=""),paper="a4")
#dot=line_merge[line_merge$sec1==x,]
  unique(new_data$sample)
#SD  P-177  P-179   P-184  P-185  P-186  P-187
#PR  P-164  P-165   P-176   P-181   P-189   P-191   P-198  P-199
#PD  P-174   P-175   P-178   P-180   P-192   P-193   P-197   P-202
library(ggplot2)
setwd("H:\\shanghai\\singlecell\\P-MATRIX\\cellphonedb")
new_data=read.csv("allinteraction_mean_pvalue_newdata2.csv")
#new_data$sec2=ifelse(new_data$sec2=="γδ T cells (T9) ","T9",new_data$sec2)
#unique(new_data$s1s2)
#new_data$s1s2=paste(new_data$sec1,new_data$sec2,sep="-")
#write.csv(new_data,file="allinteraction_mean_pvalue_newdata2.csv")
#sec1name=unique(new_data$sec1)
all=NULL
for(s1 in sec1name){
  #s1="T3"
  print(s1)
  sec1p=sec.test(new_data,"sec1",s1)
  sec2p=sec.test(new_data,"sec2",s1)
  if(is.null(all)){
    all=rbind(sec1p,sec2p)
  }else{
    all=rbind(all,sec1p,sec2p)    
  }
}

write.csv(all,file="sec1sec2_X2_test.csv")

sec.test=function(new_data,sec1,s1){
  res=NULL
 # new_data$sec1
  dot1=new_data[new_data[[sec1]]==s1&new_data$Response!="Na",]
  ppname1=unique(dot1$pp)
  for(p1 in ppname1){
    #p1="CD48???CD244"
    #sec1="sec2"
    #s1="T3"
    dot1p1=dot1[dot1$pp==p1&dot1$pvalue<0.05,]
    if(dim(dot1p1)[1]!=0){
    #dim(dot1p1)
    max_len1=max(table(dot1p1[dot1p1$Response=="PD",]$s1s2))
    max_len2=max(table(dot1p1[dot1p1$Response=="PR",]$s1s2))
    max_len3=max(table(dot1p1[dot1p1$Response=="SD",]$s1s2))
    num1=dim(dot1p1[dot1p1$Response=="PD",])[1]
    num2=dim(dot1p1[dot1p1$Response=="PR",])[1]
    num3=dim(dot1p1[dot1p1$Response=="SD",])[1]
    PDPR=chisq.test(matrix(c(num1,num2,176-num1,176-num2),nrow=2,ncol=2))$p.value
    PDSD=chisq.test(matrix(c(num1,num3,176-num1,132-num3),nrow=2,ncol=2))$p.value
    PRSD=chisq.test(matrix(c(num2,num3,176-num2,132-num3),nrow=2,ncol=2))$p.value
    out=c(sec1,s1,p1,PDPR,PDSD,PRSD,max_len1,max_len2,max_len3)
   # print(paste(out,collapse="|"))
    if(is.null(res)){
      res=out
    }else{
      res=rbind(res,
                out
      )
    }
    }
  }
  res
}

p=secpp.test(new_data)
write.csv(p,file="sspp_fishertest.csv")
secpp.test=function(new_data0){
  res=NULL
  # new_data$sec1
  new_data0=new_data
  print("remove Na...")
  new_data=new_data0[new_data0$Response!="Na",]
  print("generate sspp...")
  new_data$sspp=paste(new_data$sec1,new_data$sec2,new_data$p1,new_data$p2,sep="-")
  print("unique sspp...")
  ppname1=unique(new_data$sspp)
  print("for sspp...")
  for(p1 in ppname1){
    #p1="CD48???CD244"
    #sec1="sec2"
    #s1="T3"
    dot1p1=new_data[new_data$sspp==p1&new_data$pvalue<0.05,]
    if(dim(dot1p1)[1]!=0){
      num1=dim(dot1p1[dot1p1$Response=="PD",])[1]
      num2=dim(dot1p1[dot1p1$Response=="PR",])[1]
      num3=dim(dot1p1[dot1p1$Response=="SD",])[1]
      PDPR=fisher.test(matrix(c(8-num1,num1,8-num2,num2),nrow=2,ncol=2))$p.value
      PDSD=fisher.test(matrix(c(8-num1,num1,6-num3,num3),nrow=2,ncol=2))$p.value
      PRSD=fisher.test(matrix(c(8-num2,num2,6-num3,num3),nrow=2,ncol=2))$p.value
      out=c(p1,PDPR,PDSD,PRSD)
      print(paste(out,collapse="|"))
      if(is.null(res)){
        res=out
      }else{
        res=rbind(res,
                  out
        )
      }
    }
  }
  res
}



dot=new_data[new_data$sec1=="NK2"&new_data$pp=="CD74-COPA"&new_data$Response!="Na",]
# dot=new_data[new_data$pp=="CD74-APP",]
# unique(dot$sec1)
# dot$sec1=ifelse(grepl("T9",dot$sec1),"T9",dot$sec1)
# dot$sec2=ifelse(grepl("T9",dot$sec2),"T9",dot$sec2)
dot$RS=paste(dot$Response,dot$sample,sep="-")
dot$logp=ifelse(dot$pvalue==0,-log(0.0001),-log(dot$pvalue))
#dot$s1s2=ifelse(dot$s1s2=="T9-γδ T cells (T9) ","T9-T9",dot$s1s2)
p=ggplot(dot,aes(RS,s1s2))+
  geom_point(aes(size=as.numeric(logp),color=Response))+
#  scale_color_gradient(low="gray",high="red")+
 #facet_grid(.~sample)+ 
  ggtitle("CD74-COPA") +
  theme_bw()+
  scale_y_discrete(name = "sec1") +
  scale_x_discrete(name = "sec2",position = "top") +
  scale_size_continuous(range=c(0,3))+
  theme(plot.title = element_text(size = 11, face =  "bold"),
        text = element_text(size = 8, color="black",face="bold",),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 8, color="black",face="bold",angle = 90),# hjust = 0.5, vjust = 0.5),
        axis.text.y=element_text(size = 8, color="black",face="bold"))
print(p)
dev.off()
}
###########################################################################################
setwd("H:\\shanghai\\singlecell\\P-MATRIX\\cellphonedb")
sample=c("P-164", "P-165", "P-174", "P-175", "P-176", "P-177", "P-178", "P-179", "P-180",
"P-181", "P-182", "P-184", "P-185", "P-186", "P-187", "P-189", "P-191", "P-192",
"P-193", "P-197", "P-198", "P-199", "P-201", "P-202", "P-203")

#SD  P-177  P-179   P-184  P-185  P-186  P-187
#PR  P-164  P-165   P-176   P-181   P-189   P-191   P-198  P-199
#PD  P-174   P-175   P-178   P-180   P-192   P-193   P-197   P-202
for(s in sample){
  print(s)
  ppi=file2line(s,"significant_means.txt")
  write.csv(ppi,file=paste(s,"interaction_significant.csv",sep=""))
}
all=NULL
for(s in sample){
  print(s)
 # ppi=read.csv(paste(s,"interaction_significant.csv",sep=""),sep=",")
  ppi=read.csv(paste(s,"_interaction.csv",sep=""),sep=",")
  
  if(is.null(all)){
    all=cbind(ppi,rep(s,dim(ppi)[1]))
  }else{
    all=rbind(all,
              cbind(ppi,rep(s,dim(ppi)[1]))
    )
  }
}
#write.csv(all,file=paste("all","interaction_significant.csv",sep=""))
#data=read.csv("allinteraction_significant.csv",sep=",")
data=read.csv("allinteraction_mean_pvalue.csv",sep=",")
data=read.csv("PD_PR_SDinteraction_mean_pvalue.csv",sep=",")
# p=read.csv("all_interaction_pvalues.csv",sep=",")
# data$sspp=paste(data$sec1,data$sec2,data$p1,data$p2,sep="_")
# p$sspp=paste(p$sec1,p$sec2,p$p1,p$p2,sep="_")
# data_p = dplyr::left_join(data,p,by="sspp")
dim(data)
#colnames(data)=c("id","X","sec1","sec2","p1","p2","mean", "sample")
colnames(data)=c("id","sec1","sec2","p1","p2","mean", "sample","pvalue")
data[1:4,1:8]
samples=read.csv("../../meta.csv",header=T,sep=",")
dim(samples)
samples[1:4,1:6]
new_data = dplyr::left_join(data,samples)
#new_data=data
new_data[1:4,1:12]
new_data$pp=paste(new_data$p1,new_data$p2,sep="-")
new_data[1:4,1:13]
new_data$s1s2=paste(new_data$sec1,new_data$sec2,sep="-")
new_data[1:4,1:10]
new_data$logp=-log10(new_data$pvalue)
new_data=data.frame(new_data)
write.csv(new_data,file="allinteraction_mean_pvalue_newdata.csv")
new_data[new_data$s1s2=="T5-M3",]

cellcell=unique(new_data$s1s2)
P_value=NULL
res=NULL
for(n in cellcell){
  d1=new_data[new_data$s1s2==n&new_data$Response!="Na",]
  pairs=unique(d1$pp)
  for(m in pairs) {
    d2=d1[d1$pp==m,]
    
    #if(length(which((d2$mean>0)))>12){
      #   p =wilcox.test(d1[[m]]  ~ Response,data = d,paired = FALSE)
    tryCatch({ 
     fit<-kruskal.test(d2$mean ~ Response,data = d2)
      out=c(n,m,fit$p.value)
      P_value=rbind(P_value,out)
    },error=function(e){
      cat("ERROR :",conditionMessage(e),"\n")
     
    },finally = {
      out=c(n,m,paste(table(d2$Response),collapse="|"))
     print(out)
      P_value=rbind(P_value,out)
    }
    )
   # }
  }
}
P_value=data.frame(P_value)
colnames(P_value)=c("cell_pairs","gene_pairs","Pvalue")
#P_value$fdr=p.adjust(P_value[,3], "BH")

#write.csv(P_value,file="all_interaction_significant_pvalue.min.csv")  
unique(new_data$sec1)
dd=new_data[new_data$sec1=="T7",]
#write.csv(unique(dd$pp),file="ftr.csv")
dd=new_data[new_data$pp=="CD74-COPA"& new_data$sec1=="T7",]
dim(dd)
library(ggpubr)
my_comparisons=list(c("PD","SD"),c("PD","PR"),c("PR","SD"))
dd1=dd[dd$Response!="Na",]
p<-ggplot(dd1,aes(x=Response,y = as.numeric(mean),fill=Response)) #DURATION HBA1C_Ratio
p+# stat_boxplot(geom = "errorbar",width=0.3)+
  geom_boxplot(outlier.shape = NA, alpha=0.7)+ # alpha=0.7
#  coord_cartesian(ylim =c(0.5,1.3))+
  # scale_y_continuous(name = "Carnitine.C10.OH (Log10, �M)")+
  scale_y_continuous(name = "mean") +
  scale_x_discrete(name = "Response") +
  ggtitle("HLA-E-KLRK1") +
  facet_wrap(.~s1s2,scales="free_y")+
  # geom_line()
  geom_point(position = position_jitterdodge(),size=0.75)+
  theme_bw() +
  theme(plot.title = element_text(size = 14, face =  "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 9, color="black",face="bold"),
        axis.text.y=element_text(size = 9, color="black",face="bold"))+ 
  stat_compare_means(comparisons=my_comparisons,size = 2,color="black" )#+
  stat_compare_means(size = 4,color="black",label="p" )
  
###post_fishertest#############################################################
sspp=c("B1-M3-COPA-P2RY6",
     #  "B1-M5-TNFRSF13C-TNFSF13B",
       "B1-NK1-LGALS9-SORL1",
       "B1-T7-LGALS9-SLC1A5",
     #  "B1-T9-CD47-SIRPG",
     #  "B2-M5-CD40-TNFSF13B",
      # "B2-M5-TNFRSF13C-TNFSF13B",
       #"B2-T9-CD47-SIRPG",
      "B2-T9-LTA-TNFRSF1A",
       "M1-M3-TGFB1-complex:TGFbeta receptor1",
       "M1-T3-ICAM1-complex:aMb2 complex",
       "M1-T9-ICAM1-complex:aLb2 complex",
       "M1-T9-ICAM1-ITGAL",
       "M2-M3-TGFB1-complex:TGFbeta receptor1",
       "M2-T3-ICAM1-complex:aMb2 complex",
       "M2-T8-IL15-IL15RA",
       "M3-T9-ICAM1-complex:aLb2 complex",
       "M3-T9-ICAM1-ITGAL",
       #"M5-NK2-HLA-E-KLRK1",
      # "M5-T2-MIF-TNFRSF14",
     #  "M5-T5-ICAM3-complex:aLb2 complex",
    #   "M5-T6-MIF-TNFRSF14",
     "M5-T8-BTLA-TNFRSF14",
      "P1-M3-TGFB1-complex:TGFbeta receptor1",
       "T2-B2-CD70-TNFRSF17",
    #   "T3-T7-LAMP1-FAM3C",
      # "T3-Tx-LAMP1-FAM3C",
       "T4-M3-COPA-P2RY6",
       "T4-T9-LTA-TNFRSF14",
       "T6-P1-CD74-APP",
   #    "T6-T5-ICAM2-complex:aLb2 complex",
       "T7-M1-CD74-APP",
       "T7-M1-CD74-COPA",
       "T7-M2-CD74-APP",
       "T7-M2-CD74-COPA",
       "T7-M3-CD74-COPA",
       "T7-M3-CD74-MIF",
       "T7-M4-HLA-DPB1-TNFSF13B",
       "T7-NK1-CD74-MIF",
       "T7-NK2-CD74-MIF",
     #  "T7-NK2-FAM3C-CLEC2D",
       "T7-T7-CD74-COPA",
       "T7-Tx-CD74-APP",
       "T7-Tx-CD74-MIF",
       "T8-M2-TNF-RIPK1",
       "T8-M3-COPA-P2RY6",
       "T8-M3-TNF-RIPK1",
       "T8-NK2-TNF-FAS",
       "T8-T1-TNF-FAS",
       "T8-T3-CD48-CD244",
       "T9-T4-TNFSF12-TNFRSF25"#,

       #"T9-T7-CD48-CD244",
      # "Tx-B2-MIF-TNFRSF14"
       
)
 
  sspp=c(
         "B1-M5-TNFRSF13C-TNFSF13B",
         "B1-T9-CD47-SIRPG",
         "B2-M5-CD40-TNFSF13B",
         "B2-M5-TNFRSF13C-TNFSF13B",
         "B2-T9-CD47-SIRPG",
        # "B2-T9-LTA-TNFRSF1A",
         "M5-NK2-HLA-E-KLRK1",
         "M5-T2-MIF-TNFRSF14",
          "M5-T5-ICAM3-complex:aLb2 complex",
         "M5-T6-MIF-TNFRSF14",
        
           

          "T3-T7-LAMP1-FAM3C",
          "T3-Tx-LAMP1-FAM3C",

            "T6-T5-ICAM2-complex:aLb2 complex",
          "T7-NK2-FAM3C-CLEC2D",
        # "T9-T7-CD48-CD244",
          "Tx-B2-MIF-TNFRSF14"

  ) 
library(ggplot2)
setwd("H:\\shanghai\\singlecell\\P-MATRIX\\cellphonedb")
#new_data=read.csv("allinteraction_mean_pvalue_newdata2.csv")
new_data0=new_data
new_data=new_data0[new_data0$Response!="Na",]
new_data$sspp=paste(new_data$sec1,new_data$sec2,new_data$p1,new_data$p2,sep="-")
ppname1=unique(new_data$sspp)
all_dot=NULL
for(pp in sspp){
    dot=new_data[new_data$sspp==pp,]
    all_dot=rbind(all_dot,dot)
}
all_dot$logp=ifelse(all_dot$pvalue==0,-log(0.0001),-log(all_dot$pvalue))
all_dot$RS=paste(all_dot$Response,all_dot$sample,sep="-")
#"#D6604D","#4393C3","#5AAE61"
all_dot$color_label=ifelse(all_dot$pvalue>0.05,"gray",
                           ifelse(all_dot$Response=="PD","#D6604D",
                                  ifelse(all_dot$Response=="PR","#4393C3",
                                         ifelse(all_dot$Response=="SD","#5AAE61","gray"))))
p=ggplot(all_dot,aes(RS,sspp))+
  geom_point(aes(size=as.numeric(logp)),color=all_dot$color_label)+
  #  scale_color_gradient(low="gray",high="red")+
  #facet_grid(.~sample)+ 
#  ggtitle("CD74-COPA") +
  theme_bw()+
  scale_y_discrete(name = "sec1") +
  scale_x_discrete(name = "sec2",position = "top") +
  scale_size_continuous(range=c(0,3))+
  theme(plot.title = element_text(size = 11, face =  "bold"),
        text = element_text(size = 8, color="black",face="bold",),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 8, color="black",face="bold",angle = 90),# hjust = 0.5, vjust = 0.5),
        axis.text.y=element_text(size = 8, color="black",face="bold"))
print(p)
all_dot[,c(10,15,21)]

