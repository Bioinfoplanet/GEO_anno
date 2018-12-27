# step4 type作图
#输入文件：前两步的Rdata
#输出文件：韦恩图
#4-1 idchange and 4 venn
rm(list = ls())
load('hgnc_family_type.Rdata')
load('ABM_gene_probe.Rdata')
library(tidyverse)
count_type <- count(hgnc_type,locus_group,sort = T)
type_top4 <- as.character((count_type)$locus_group)

Aff_type <- merge(Aff,hgnc_type,
                  by.x='ensembl_id',
                  by.y ='ensembl_gene_id')[,-1]
Bio_type <- merge(Bio,hgnc_type,
                  by.x='ensembl_id',
                  by.y ='ensembl_gene_id')[,-1]
Mine_type <- merge(Mine,hgnc_type,
                   by.x='ensembl_id',
                   by.y ='ensembl_gene_id')[,-1]


#画图
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
venn(uni(Aff_type),uni(Bio_type),uni(Mine_type),"top4type_all")
save(type_top4,Bio_type,Aff_type,Mine_type,file = 'ABM_top4type.Rdata')


# 4-2 分开画图
load('ABM_top4type.Rdata')
ABM_type_list <- list(Aff_type=Aff_type,Bio_type=Bio_type,Mine_type=Mine_type)

output_list <- list()
for (j in 1:3){
  output <- list()
  for (i in 1 :length(type_top4)){
    output[[i]] <- (filter(ABM_type_list[[j]],
                           locus_group==type_top4[i])
    )$probe_id
  }
  names(output) <- paste(1:length(type_top4),type_top4,sep = ".")
  output_list[[j]] = output
}
names(output_list) <- names(ABM_type_list)
            
for (i in c(1:length(type_top4))){
  venn(output_list[[1]][[i]],
       output_list[[2]][[i]],
       output_list[[3]][[i]],
       name = names(output_list[[1]][i]))
}
