library(ggplot2)
library(ggpubr)
library(meta)
setwd("H:\\shanghai\\singlecell\\P-MATRIX\\time")

meta.data=read.csv("sample.combined_T10_others_singleR_ref3_main_fine_celltype.csv",sep=",")
unique(meta.data$celltype)
# "Undefined T cells"                                
# [2] "CD56+CD16- NK cells (NK1)"                        
# [3] "naive CD4+ T cells (T1)"                          
# [4] "C56-CD16+ NK cells (NK2)"                         
# [5] "immature B cells (B2)"                            
# [6] "cytotoxic CD8+ lymphocytes (CD8+ CTL) (T6)"       
# [7] "Platelets"                                        
# [8] "mucosal-associated invariant T (MAIT) cells (T8) "
# ###[9] "CD14+ monocytes (M1)"####19.6G                             
# [10] "effector memory CD8+T cells (T5, CD8 Tm)"         
# [11] "Î³Î´ T cells (T9) "                               
# [12] "naive B cells (B1)"                               
# [13] "CD16+ monocytes (M2)"                             
# [14] "regulatory CD4+ T cells (T2, Treg)"               
#  [15] "CD1C+ dendritic cells (M3)"                       
# [16] "Proliferating T cells (T7, Tprol) "               
# [17] "naive CD8+ T cells (T4)"                          
# [18] "CD4+ T cells (T3)"                                
# [19] "HSC"                                              
# [20] "plasmacytoid dendritic cells (CLEC4C+ pDC) (M4)"  
# [21] "memory B cells (B3)"                              
# [22] "CLEC9A+ dendritic cells (M5)"                 
ctype="CD4+ T cells (T3)"
# fdata=read.csv(paste("Tall_",ctype,"_fdata.csv",sep=""),sep=",")
myfiles <- list.files(pattern = "*_fdata.csv")
for(f in myfiles){
  fdata=read.csv(f,sep=",")
  xGroup=fdata$Response
  ydata=fdata$Pseudotime
  my_comparisons=list(c("PD","SD"),c("PD","PR"),c("PR","SD"))
  p<-ggplot(fdata,aes(x=xGroup,y = ydata,fill=xGroup))+
    geom_boxplot(outlier.shape = NA)+ # alpha=0.7
     scale_y_continuous(name = "Pseudotime") +
      geom_point(position = position_jitterdodge(0.1),size=0.5)+
    ggtitle(f)+
    theme_bw() +
    theme(plot.title = element_text(size = 11, face =  "bold"),
          text = element_text(size = 11),
          axis.title = element_text(face="bold"),
          axis.text.x=element_text(size = 8, color="black",face="bold",angle = 45, hjust = 0.5, vjust = 0.5),
          axis.text.y=element_text(size = 10, color="black",face="bold"))+
    stat_compare_means(comparisons=my_comparisons,size = 4,color="black")
  print(p)
}


sample_names=names(table(fdata$sample))
new_data=NULL
for(n in sample_names){
   population=data.frame(prop.table(table(fdata[fdata$sample==n,]$State)))
   count_data=data.frame(count=table(fdata[fdata$sample==n,]$State))
   all_count=length(fdata[fdata$sample==n,]$State)
  col=cbind(population,count_data,allcount=rep(all_count,dim(population)[1]),sample=rep(n,dim(population)[1]))
  if(length(new_data)==0){
    new_data=col
  }else{
    new_data=rbind(new_data,col)
  }
}
names(new_data)

samples=read.csv("../../meta.csv",header=T,sep=",")
new_data = dplyr::left_join(new_data,samples)
table(new_data$Var1)
names(new_data)
for(x in unique(new_data$Var1)){
  for(y in unique(new_data$Response)){
    new_data1=new_data[new_data$Var1==x&new_data$Response==y,]
    meta1 <- metaprop (count.Freq ,allcount, data=new_data1,sm="PFT")
    
  }
}
str(meta1)



my_comparisons=list(c("PD","SD"),c("PD","PR"),c("PR","SD"))

p<-ggplot(new_data,aes(x=Response,y = round(as.numeric(Freq)*100,3),fill=Response)) #DURATION HBA1C_Ratio
p+ #stat_boxplot(geom = "errorbar",width=0.3)+
  geom_boxplot(outlier.shape = NA)+ # alpha=0.7
  # coord_cartesian(ylim =c(0.5,1.3))+
  # scale_y_continuous(name = "Carnitine.C10.OH (Log10, ÂµM)")+
  scale_y_continuous(name = "Population (%)") +
  # scale_x_discrete(name = "T_cells") +
  ggtitle("treatment1234") +
  facet_wrap(.~new_data$Var1,scales="free_y")+
  # geom_line()
  geom_point(position = position_jitterdodge(0.01),size=0.5)+
  theme_bw() +
  theme(plot.title = element_text(size = 11, face =  "bold"),
        text = element_text(size = 11),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 8, color="black",face="bold",angle = 45, hjust = 0.5, vjust = 0.5),
        axis.text.y=element_text(size = 10, color="black",face="bold"))+
  stat_compare_means(comparisons=my_comparisons,size = 2,color="black")


xGroup=fdata$Response
ydata=fdata$Pseudotime
my_comparisons=list(c("PD","SD"),c("PD","PR"),c("PR","SD"))
p<-ggplot(fdata,aes(x=xGroup,y = ydata,fill=xGroup)) #DURATION HBA1C_Ratio
p+ #stat_boxplot(geom = "errorbar",width=0.3)+
  geom_boxplot(outlier.shape = NA)+ # alpha=0.7
  # coord_cartesian(ylim =c(0.5,1.3))+
  # scale_y_continuous(name = "Carnitine.C10.OH (Log10, ÂµM)")+
  scale_y_continuous(name = "Pseudotime") +
  # scale_x_discrete(name = "T_cells") +
 # ggtitle("treatment1234") +
 # facet_wrap(.~new_data1$Var1,scales="free_y")+
  # geom_line()
  geom_point(position = position_jitterdodge(0.1),size=0.5)+
  theme_bw() +
  theme(plot.title = element_text(size = 11, face =  "bold"),
        text = element_text(size = 11),
        axis.title = element_text(face="bold"),
        axis.text.x=element_text(size = 8, color="black",face="bold",angle = 45, hjust = 0.5, vjust = 0.5),
        axis.text.y=element_text(size = 10, color="black",face="bold"))+
#stat_compare_means(size = 2,color="black")
 stat_compare_means(comparisons=my_comparisons,size = 4,color="black")
