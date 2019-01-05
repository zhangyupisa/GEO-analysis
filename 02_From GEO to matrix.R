##SOP for GEO analysis##

#02_From GEO to matrix

library(GEOquery)
library(xlsx)
library(openxlsx)

##----------------------------------------------------------------------------------------
#Get expression matrix
gset = getGEO(filename="./GSEXXXX_series_matrix.txt.gz",AnnotGPL=T)
pdata = pData( gset ) ## “GEOquery”包中的pData函数用来取出样本信息

exprSet_raw<-read.xlsx("GSExxxx_series_matrix.xlsx",sheet=1,rowNames = TRUE,colNames = TRUE)
#提前对xlsx进行了处理，删去了表达矩阵之前的内容，直接读入

##----------------------------------------------------------------------------------------
#Get group list data and new expression matrix
group_list = as.character( pdata[, 35] )#在样例中，pdata中第35列记录了样本类型
#表达矩阵只提供了样本编号，而没有样本类型，所以需要group_list数据对应

n_expr = exprSet_raw[ , grep( 'non-tumor',group_list )]#使用通配符匹配
g_expr = exprSet_raw[ , grep( 'glioblastoma',group_list )]
o_expr = exprSet_raw[ , grep( 'oligodendroglioma',group_list )]
a_expr = exprSet_raw[ , grep( 'astrocytoma',group_list )]

exprSet = log(cbind( n_expr, g_expr ))
exprSet_boxplot=log(cbind(n_expr,o_expr,a_expr,g_expr))
#直接获取了对数值

group_list = c(rep('normal',ncol( n_expr)), rep('gbm',ncol(g_expr)))
group_list_boxplot = c(rep('NC',ncol( n_expr)), rep('OD',ncol( o_expr)),
                       rep('A',ncol(a_expr)),rep('GBM',ncol( g_expr)))
#同时新建了后续做boxplot要用到的表达矩阵与group_list

##----------------------------------------------------------------------------------------
#Select lncRNA from ref and change gene id to gene name, as well as select the highest expression rows

ref<-read.xlsx2(file="lncRNA.xlsx",sheetIndex=1)
dim(ref)
#读入lncRNA reference矩阵，共有1950个lncRNA
colnames(ref)
#我们需要的是probeID与Gene.Symbol这两列

lncref<-ref[,c(1,4)]#取出这两列
colnames(lncref)<-c("ID","gene")




#筛选exprSet，去除其中没有在ref构建成的lncref内的ID
#ids的第一列是芯片所用的ID，第二列是基因名
expr = exprSet[rownames(exprSet) %in% lncref[,1],]
expr_boxplot = exprSet_boxplot[rownames(exprSet_boxplot) %in% lncref[,1],]

dim(expr)#54613-->1950
dim(expr_boxplot)




MAX = by(expr,lncref[,2], 
         function(x) rownames(x)[ which.max( rowMeans(x) ) ] )
MAX = as.character(MAX)
#相同基因表达数据取最大值
expr = expr[ rownames(expr) %in% MAX,]
#筛选出exprSet矩阵里在MAX中的ID
rownames( expr ) = lncref[ match( rownames( expr ), lncref[ , 1 ] ), 2 ]
#将筛选出的ID改为基因名，获得以基因名记录的表达矩阵

MAX_boxplot = by(expr_boxplot,lncref[,2], 
                 function(x) rownames(x)[ which.max( rowMeans(x) ) ] )
MAX_boxplot = as.character(MAX_boxplot)
expr_boxplot = expr_boxplot[ rownames(expr_boxplot) %in% MAX_boxplot,]
rownames( expr_boxplot ) = lncref[ match( rownames( expr_boxplot ), lncref[ , 1 ] ), 2 ]

#现在我们从最一开始的gset得到了3个基因表达矩阵，分别是原始未过滤数据exprSet_raw，对样本过滤后的、用对数记录表达量的exprSet,
#以及用lncRNA过滤后的、选出了最大表达量的expr。需要保存的是expr和后续boxplot要用的expr_boxplot。另外两个group_list也要保存。
save(expr,group_list,group_list_boxplot, file = 'expr_grouplist.Rdata' )