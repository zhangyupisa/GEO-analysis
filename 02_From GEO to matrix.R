##SOP for GEO analysis##

#02_From GEO to matrix

library(GEOquery)
library(xlsx)
library(openxlsx)

##----------------------------------------------------------------------------------------
#Get expression matrix
gset = getGEO(filename="./GSEXXXX_series_matrix.txt.gz",AnnotGPL=T)
pdata = pData( gset ) ## pData in  “GEOquery” package is used to extract sample information

exprSet_raw<-read.xlsx("GSExxxx_series_matrix.xlsx",sheet=1,rowNames = TRUE,colNames = TRUE)
#If you dealt with this file beforehand and deleted the information ahead of the expression matrix, you can read it in like that directly

#----------------------------------------------------------------------------------------
#Get group list data and new expression matrix
group_list = as.character( pdata[, 35] )#Column No.35 records sample information in the example
#The expression provides only the number of the samples but not their disease type, so group_list is needed to match them

n_expr = exprSet_raw[ , grep( 'non-tumor',group_list )]#wildcard capture
g_expr = exprSet_raw[ , grep( 'glioblastoma',group_list )]
o_expr = exprSet_raw[ , grep( 'oligodendroglioma',group_list )]
a_expr = exprSet_raw[ , grep( 'astrocytoma',group_list )]

exprSet = log(cbind( n_expr, g_expr ))
exprSet_boxplot=log(cbind(n_expr,o_expr,a_expr,g_expr))


group_list = c(rep('normal',ncol( n_expr)), rep('gbm',ncol(g_expr)))
group_list_boxplot = c(rep('NC',ncol( n_expr)), rep('OD',ncol( o_expr)),
                       rep('A',ncol(a_expr)),rep('GBM',ncol( g_expr)))


##----------------------------------------------------------------------------------------
#Select lncRNA from ref and change gene id to gene name, as well as select the highest expression rows

ref<-read.xlsx2(file="lncRNA.xlsx",sheetIndex=1)
dim(ref)
#Read lncRNA reference matrix, 1950 lncRNAs in all
colnames(ref)
#We only need columns "probeID" and "Gene.Symbol"

lncref<-ref[,c(1,4)]#Extract these two columns
colnames(lncref)<-c("ID","gene")



#Screen exprSet and delete the IDs that are not in lncRNA reference matrix
expr = exprSet[rownames(exprSet) %in% lncref[,1],]
expr_boxplot = exprSet_boxplot[rownames(exprSet_boxplot) %in% lncref[,1],]

dim(expr)#54613-->1950
dim(expr_boxplot)




MAX = by(expr,lncref[,2], 
         function(x) rownames(x)[ which.max( rowMeans(x) ) ] )
MAX = as.character(MAX)
#Get the maximum expression number for the same gene
expr = expr[ rownames(expr) %in% MAX,]
#Screen the IDs in exprSet that are in MAX
rownames( expr ) = lncref[ match( rownames( expr ), lncref[ , 1 ] ), 2 ]
#Change the IDs that are selected to gene name, and we get a gene expression matrix written in gene names

MAX_boxplot = by(expr_boxplot,lncref[,2], 
                 function(x) rownames(x)[ which.max( rowMeans(x) ) ] )
MAX_boxplot = as.character(MAX_boxplot)
expr_boxplot = expr_boxplot[ rownames(expr_boxplot) %in% MAX_boxplot,]
rownames( expr_boxplot ) = lncref[ match( rownames( expr_boxplot ), lncref[ , 1 ] ), 2 ]

save(expr,group_list,group_list_boxplot, file = 'expr_grouplist.Rdata' )#For next time use
