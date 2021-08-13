
vdj_cgene=function(n,tb_label){
  result=NULL
  # n="P-164"
  # tb_label="B"
  name=paste(n,"-",tb_label,"/filtered_contig_annotations.csv",sep="")
  d=read.csv(name,sep=",")
 # names(d)
  new_line=NULL
 if(tb_label=="B"){
  IGH=data.frame(prop.table(table(d[d$chain=="IGH",]$c_gene)))
  IGL=data.frame(prop.table(table(d[d$chain!="IGH",]$c_gene)))
  new_line=rbind(cbind(IGH,sample=rep(n,dim(IGH)[1]),cell_type=rep(tb_label,dim(IGH)[1]),IG_label=rep("IGH",dim(IGH)[1])),
            cbind(IGL,sample=rep(n,dim(IGL)[1]),cell_type=rep(tb_label,dim(IGL)[1]),IG_label=rep("IGL",dim(IGL)[1]))
           )
 }else{
   TR=data.frame(prop.table(table(d$c_gene)))
   new_line=cbind(TR,sample=rep(n,dim(TR)[1]),cell_type=rep(tb_label,dim(TR)[1]),IG_label=rep("TR",dim(TR)[1]))
            
 }
  if(length(result)==0){
    result=new_line
  }else{
    result=rbind(result,new_line)
  }
  result
}




setwd("H:\\shanghai\\singlecell\\P-MATRIX")
meta=read.csv("..\\meta.csv",sep=",")


v=c("P-164",
    "P-165",
    "P-174",
    "P-175",
    "P-176",
    "P-177",
    "P-178",
    "P-179",
    "P-180",
    "P-181",
    "P-183",
    "P-184",
    "P-185",
    "P-186",
    "P-187",
    "P-189",
    "P-191",
    "P-192",
    "P-193",
    "P-197",
    "P-198",
    "P-199",
    "P-202"
)
res=NULL
for(n in v){
  r1=vdj_cgene(n,"T")
  r2=vdj_cgene(n,"B")
  
  if(length(res)==0){
    res=rbind(r1,r2)
  }else{
    res=rbind(res,r1,r2)
  }
}
res_meta=dplyr::left_join(res,meta)
write.csv(res_meta,file="vdj_c_gene.csv")

setwd("H:\\shanghai\\singlecell\\P-MATRIX")
meta1=read.csv("vdj_c_gene.csv",sep=",")
library(ggplot2)
library(ggpubr)
meta=meta1[meta1$IG_label=="IGL",]
x=factor(meta$Response)
y=meta$Freq
summary(meta)
p<-ggplot(meta,aes(x=factor(Var1),y = y,fill=x)) #DURATION HBA1C_Ratio
p+# stat_boxplot(geom = "errorbar",width=0.3)+
  geom_boxplot(outlier.shape = NA, alpha=0.7)+ # alpha=0.7
 # coord_cartesian(ylim =c(0.5,1.3))+
  # scale_y_continuous(name = "Carnitine.C10.OH (Log10, µM)")+
  #scale_y_continuous(name = "HbA1c (baseline %)") +
 # scale_x_discrete(name = "Duration") +
  ggtitle("") +
 facet_grid(.~meta$Gender)+
  # geom_line()
  geom_point(position = position_jitterdodge(0.1),size=1)+
  theme_bw() +
  theme(plot.title = element_text(size = 14, face =  "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 9, color="black",face="bold"),
        axis.text.y=element_text(size = 9, color="black",face="bold"))+ 
  stat_compare_means(method = "t.test",size = 4,color="black",label.y=0.8 )
  #stat_compare_means(aes(label =round(..p..,3)),size = 4,color="black" )


