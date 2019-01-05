##SOP for GEO analysis##

#03_Differential Expression Analysis

##-------------------
#分组矩阵

library(limma)

load('expr_grouplist.Rdata')

#仅对expr进行差异分析
design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
length(group_list)
rownames(design)=colnames(expr)
head(design)
#design是分组矩阵，记录每个样本的分组情况

##-------------------
#差异比较矩阵

contrast.matrix<-makeContrasts(paste0(unique(group_list),collapse = "-"),levels = design)
contrast.matrix #这个矩阵声明，我们要把gbm组跟normal组进行差异分析比较

##-------------------
#差异分析

#step1
fit <- lmFit(expr,design)
#step2
fit2 <- contrasts.fit(fit, contrast.matrix) 
fit2 <- eBayes(fit2)  
#step3
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput) 
head(nrDEG)