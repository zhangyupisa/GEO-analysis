##SOP for GEO analysis##

#03_Differential Expression Analysis

##-------------------
#Grouping matrix

library(limma)

load('expr_grouplist.Rdata')

design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
length(group_list)
rownames(design)=colnames(expr)
head(design)

##------------------
#Differential comparation matrix

contrast.matrix<-makeContrasts(paste0(unique(group_list),collapse = "-"),levels = design)
contrast.matrix #The matrix compares gbm with normal. When comparing group A with group B, A=1, B=-1. In "group_list",A=1, B=0.


##-------------------
#Differential analysis

#step1
fit <- lmFit(expr,design)
#step2
fit2 <- contrasts.fit(fit, contrast.matrix) 
fit2 <- eBayes(fit2)  
#step3
tempOutput = topTable(fit2, coef=1, n=Inf)
nrDEG = na.omit(tempOutput) 
head(nrDEG)
