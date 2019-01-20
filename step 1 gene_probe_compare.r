#输入文件：官网注释文件HG-U133_Plus_2.na36.annot.csv，自由注释文件GPL570_probe2ensemb.csv'
#输出文件：三个注释的Rdata，格式两列分别是id和ensemblid，两个韦恩图，分别是probe和gene
#http://www.affymetrix.com/support/technical/byproduct.affx?product=hg-u133-plus,
#下载地址http://www.affymetrix.com/Auth/analysis/downloads/na36/ivt/HG-U133_Plus_2.na36.annot.csv.zip

#step 1 gene_probe_compare
if(!require("hgu133plus2.db"))BiocManager::install("hgu133plus2.db")
if(!require("clusterProfiler"))BiocManager::install("clusterProfiler")
if(!require("org.Hs.eg.db"))BiocManager::install("org.Hs.eg.db")
if(!require(tidyverse))install.packages('tidyverse')
library(hgu133plus2.db)
library(tidyverse)
empty.omit <- function(x){
  x[x=="---"] <- NA
  na.omit(x)
}
zz1 = read.csv('HG-U133_Plus_2.na36.annot.csv',
               comment.char = "#")[,c(1,19)]
zz1 <- empty.omit(zz1)
zz2<- apply(zz1,
            1,
            function(x){
              paste(x[1],
                    str_split(x[2],'///',simplify=T),
                    sep = "///")
            })
zz3 <- tibble(unlist(zz2))

colnames(zz3) <- "lala" 
zz4 <- separate(zz3,lala,c("probe","gene"),sep = "///")
library(clusterProfiler) 
library(org.Hs.eg.db)
zz5<- bitr(zz4$gene, fromType = "ENTREZID",
           toType = c("ENSEMBL"),
           OrgDb = org.Hs.eg.db)
zz6 <- merge(zz4,zz5,by.x = "gene",by.y='ENTREZID')[,-1]
#官网注释
Aff <- zz6
rm(list =paste0("zz",1:6))
#biocductor包
Bio = toTable(hgu133plus2ENSEMBL)
#新流程注释
Mine =  read.csv('GPL570_probe2ensemb.csv')[,c(6,12)]
#统一列名
colnames(Aff) <- colnames(Bio)
colnames(Mine)<- colnames(Bio)
dumd <- function(x){
  colname <- vector("character")
  count <- vector("integer")
  for(i in 1:ncol(x)){
    colname[i] = colnames(x)[[i]]
    count[i]=nrow(x[!duplicated(x[,i]),])
  }
  df <- tibble(colname,count) %>%
    arrange(desc(count))
  df
}
multidumd <- function(...){
  input <- list(...)
  output <- list()
  for (i in 1:length(input)){
    output[[i]] <- dumd(noquote((input)[[i]]))
  }
  output
}
#三个注释的基因与探针汇总
sum_all <- multidumd(Aff,Bio,Mine)
#画韦恩图
#包装韦恩图函数
venn <- function(x,y,z,name){
  if(!require(VennDiagram))install.packages('VennDiagram')
  library (VennDiagram)
  venn.diagram(x= list(Aff = x,Bio = y,Mine = z),
               filename=paste(name,".tiff"),
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
# #一张
venn(Aff$probe_id,Bio$probe_id,Mine$probe_id,"mapped_probe")
# #两张
venn(Aff$ensembl_id,Bio$ensembl_id,Mine$ensembl_id,"mapped_gene")
# 保存基因与探针数据
save(Aff,Bio,Mine,file = 'mapped_gene_probe.Rdata')


