# setwd("singlecell/P-MATRIX")
# pbmc<- readRDS("sample.combined_T10_singleR_hesc_dice.rds")
# write.csv(pbmc@meta.data,file="pbmc_meta.data.csv")
setwd("H:\\shanghai\\singlecell\\P-MATRIX\\all")
#meta.data=read.csv("pbmc_meta.data.ref3_fine.csv",sep=",")
#meta.data=read.csv("pbmc_Tothers_meta.data.ref3_main_fine.csv",sep=",")
# 
# 
names(meta.data)
table(meta.data$orig.ident)
unique(meta.data$sample)
meta.data=read.csv("sample.combined_T10_others_singleR_ref3_main_fine_celltype.csv",sep=",")
#meta.data=meta.data0[meta.data0$singleR.main2=="Monocyte",]
table(meta.data$orig.ident)
summary(meta.data$percent.mito)

table(meta.data$Treatment!=5)
table(meta.data$Response!="Na")
# tem=meta.data$singleR
# tem[which(tem=="Monocytes")] <-"Monocyte"
# meta.data$singleR2=tem

# table(meta.data[meta.data$sample=="P-199",]$seurat_clusters)
# 620/12161
# dim(meta.data)
# names(meta.data)
# summary(meta.data)
library(ggplot2)
library(ggpubr)
# pdf("bar_response.pdf")
table(meta.data$singleR.main2)
sample_names=names(table(meta.data$sample))
new_data=NULL
for(n in sample_names){
 # n="P-165"
 # population=data.frame(prop.table(table(meta.data[meta.data$sample==n,]$singleR.main2)))
 # population=data.frame(prop.table(table(meta.data[meta.data$sample==n,]$singleR2)))
 #population=data.frame(prop.table(table(meta.data[meta.data$sample==n,]$seurat_clusters)))
  population=data.frame(prop.table(table(meta.data[meta.data$sample==n,]$celltype)))
  col=cbind(population,sample=rep(n,dim(population)[1]))
  if(length(new_data)==0){
    new_data=col
  }else{
    new_data=rbind(new_data,col)
  }
}


samples=read.csv("../../meta.csv",header=T,sep=",")
new_data = dplyr::left_join(new_data,samples)
table(new_data$Var1)
names(new_data)
my_comparisons=list(c("PD","SD"),c("PD","PR"),c("PR","SD"))
new_data1=new_data[new_data$Response!='Na'
                #   &new_data$Var1 %in% c(0,3,8,48,56)
                   &new_data$Treatment!=5
                #& new_data$sample!="P-175"&
                ,]
#factor(Response,levels=c("PD","SD","PR"))
p<-ggplot(new_data1,aes(x=Response,y = round(as.numeric(Freq)*100,3),fill=Response)) #DURATION HBA1C_Ratio
p+ #stat_boxplot(geom = "errorbar",width=0.3)+
  geom_boxplot(outlier.shape = NA)+ # alpha=0.7
 # coord_cartesian(ylim =c(0.5,1.3))+
  # scale_y_continuous(name = "Carnitine.C10.OH (Log10, µM)")+
  scale_y_continuous(name = "Population (%)") +
# scale_x_discrete(name = "T_cells") +
  ggtitle("treatment1234") +
 facet_wrap(.~new_data1$Var1,scales="free_y")+
  # geom_line()
  geom_point(position = position_jitterdodge(0.01),size=0.5)+
  theme_bw() +
  theme(plot.title = element_text(size = 11, face =  "bold"),
        text = element_text(size = 11),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 8, color="black",face="bold",angle = 45, hjust = 0.5, vjust = 0.5),
        axis.text.y=element_text(size = 10, color="black",face="bold"))+
 stat_compare_means(comparisons=my_comparisons,size = 2,color="black")
  # stat_compare_means(size = 1,color="black",label="p")


p<-ggplot(meta.data,aes(x=sample,y =percent.mito)) #DURATION HBA1C_Ratio
p+ #stat_boxplot(geom = "errorbar",width=0.3)+
  geom_boxplot(outlier.shape = NA)+ # alpha=0.7
  # coord_cartesian(ylim =c(0.5,1.3))+
  # scale_y_continuous(name = "Carnitine.C10.OH (Log10, µM)")+
  scale_y_continuous(name = "Population (%)") +
  # scale_x_discrete(name = "T_cells") +
  #  ggtitle("T_cells") +
  # facet_grid(.~new_data$Gender)+
  # geom_line()
  geom_point(position = position_jitterdodge(0.01),size=1)+
  theme_bw() +
  theme(plot.title = element_text(size = 8, face =  "bold"),
        text = element_text(size = 8),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 5, color="black",face="bold"),
        axis.text.y=element_text(size = 5, color="black",face="bold"))+ 
  stat_compare_means(comparisons=my_comparisons,size = 4,color="black",method = "t.test")
# stat_compare_means(size = 3,color="black",label="p")





ggplot(new_data,aes(Response,fill=var1))+
  geom_bar(position="fill")+
  ggtitle("")+
  labs(x="Response",y="percent")+
  # coord_flip()+
  theme_bw()+
  theme(legend.position="right",
        plot.title = element_text(size = 14, face =  "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 9, color="black",face="bold"),
        axis.text.y=element_text(size = 9, color="black",face="bold"))
# dev.off()
# pdf("bar_sample.pdf")
new_bar$sample_orderd=factor(new_bar$sample,levels=c(
  "P-174",
  "P-175",
  "P-180",
  "P-197",
  "P-202",
  "P-165",
  "P-191",
  "P-198",
  "P-199",
  "P-177",
  "P-179",
  "P-183",
  "P-185",
  "P-186"
))
ggplot(new_bar,aes(sample_orderd,fill=celltype))+
  geom_bar(position="fill")+
  ggtitle("")+
  labs(x="Sample",y="percent")+
  # coord_flip()+
  theme_bw()+
  theme(legend.position="right",
        plot.title = element_text(size = 14, face =  "bold"),
        text = element_text(size = 12),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 9, color="black",face="bold"),
        axis.text.y=element_text(size = 9, color="black",face="bold"))
# dev.off()