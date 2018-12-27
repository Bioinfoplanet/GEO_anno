#3-1 reshape
#输入文件：前两步的Rdata
#输出文件：韦恩图：总和分
rm(list = ls())
load('hgnc_family_type.Rdata')
load('ABM_gene_probe.Rdata')
if(!require(tidyverse))install.packages('tidyverse')
library(tidyverse)
count_family <- count(hgnc_family,gene_family,sort = T)
family_top20 <- as.character((head(count_family,20))$gene_family)

family_s <- filter(hgnc_family,gene_family %in% family_top20)
Aff_family <- merge(Aff,family_s,
                    by.x='ensembl_id',
                    by.y ='ensembl_gene_id')[,-1]
Bio_family <- merge(Bio,family_s,
                    by.x='ensembl_id',
                    by.y ='ensembl_gene_id')[,-1]
Mine_family <- merge(Mine,family_s,
                     by.x='ensembl_id',
                     by.y ='ensembl_gene_id')[,-1]

#画图
#包装韦恩图函数
venn <- function(x,y,z,name){
  if(!require(VennDiagram))install.packages('VennDiagram')
  library (VennDiagram)
  venn.diagram(x= list(Aff = x,Bio = y,Mine = z),
               filename=paste(name,".png"),
               lwd=1,#圈线粗度
               lty=1, #圈线类型
               col=c('#0099CC','#FF6666','#FFCC99'), #圈线颜色
               fill=c('#0099CC','#FF6666','#FFCC99'), #填充颜色
               cat.col=c('#0099CC','#FF6666','#FFCC99'),#A和B的颜色
               cat.cex = 1.5,# A和B的大小
               rotation.degree = 0,#旋转角度
               main = name,#主标题内容
               main.cex = 1.5,#主标题大小
               cex=1.5,#里面交集字的大小
               alpha = 0.5,#透明度
               reverse=TRUE)
}
uni <- function(x){
  (unite(x,"x1",c(colnames(x)[1],colnames(x)[2]),sep = " "))[,1]
}
venn(uni(Aff_family),uni(Bio_family),uni(Mine_family),"top20_all")
save(family_top20,Bio_family,Aff_family,Mine_family,file = 'ABM_top20family.Rdata')


# 3-2 分开画图
load('ABM_top20family.Rdata')
ABM_family_list <- list(Aff_family=Aff_family,Bio_family=Bio_family,Mine_family=Mine_family)
#需要改掉特殊字符
family_top20[6] <- 'Small nucleolar RNAs C or D box'
family_top20[10] <- 'Small nucleolar RNAs H or ACA box'
output_list <- list()
for (j in 1:3){
  output <- list()
  for (i in 1 :length(family_top20)){
    output[[i]] <- (filter(ABM_family_list[[j]],
                           gene_family==family_top20[i])
    )$probe_id
  }
  names(output) <- paste(1:20,family_top20,sep = '.')
  output_list[[j]] = output
}
names(output_list) <- names(ABM_family_list)


for (i in 1:20){
  venn(output_list[[1]][[i]],
       output_list[[2]][[i]],
       output_list[[3]][[i]],
       name = names(output_list[[1]][i]))
}
