install.packages("ggseqlogo")
#从GitHub中安装
devtools::install.github("omarwagih/ggseqlogo")

#1RidgePlot
#install.packages("entropy")
##加载entropy这个包
library(entropy)
##我们举一个简单的例子来计算信息熵，例如我有一个位点，在这个位点有400个A，200个G，10个C，计算该位点的信息熵
##首先我们生成一个“400个A，200个G，10个C”的长字符变量
pos_A <- rep("A",times=400)
pos_G <- rep("G",times=200)
pos_C <- rep("C",times=10)
pos <- c(pos_A,pos_G,pos_C)
##对数据进行统计
pos_stat<- as.data.frame(table(pos))
##展示统计结果，A:400,C:10,G:200
pos_stat
'''
 pos Freq
1   A  400
2   C   10
3   G  200
'''
##计算信息熵，使用第二列统计数目进行计算
pos_entropy <- entropy(pos_stat[,2],unit = "log2") ##[1] 1.023923
##根据信息熵计算bits值
pos_bits <- log2(4)-pos_entropy

