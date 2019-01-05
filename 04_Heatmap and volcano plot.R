##SOP for GEO analysis##

#04_Heatmap and volcano plot
#Heatmap
library("gplots")
#先构建绘制热图需要的矩阵，我们只对上调的84个基因与下调的78个基因绘图，并根据样本类型分类。
#该矩阵行为UP与DOWN的基因名称，列为样本编号，每个值为FC
heatexpr = expr[nrDEG$change != 'NOT' ,]
heatexpr<-as.matrix(heatexpr)

plot_color=c(rep("green",23),rep("blue",77))
heatmap.2(heatexpr, col=greenred(75), scale="row",
          ColSideColors = plot_color,margins=c(3,2),
          main="Non-tumor vs Glioblastoma",dendrogram = "col",Colv = T,
          key=T, symkey=FALSE, density.info="none", trace="none",
          labRow = NA,labCol = NA)
#在保存文件时若使用Rmarkdown则会出现一些小问题，不如直接右键单击图片，copy之后复制到其他文件中保存


#Volcano
library(ggplot2)
logFC_cutoff = 0.56
#原文说82个下调，81个上调。但是若按原文设置logFC_cutoff=1仅有15个下调，35个上调。
nrDEG$change = as.factor( ifelse( nrDEG$P.Value <0.05 & abs(nrDEG$logFC) > logFC_cutoff,
                                  ifelse( nrDEG$logFC > logFC_cutoff , 'UP', 'DOWN' ), 'NOT' ) )
sort(table(nrDEG$change))
head(nrDEG)

volcano = ggplot(data = nrDEG, aes( x = logFC, y = -log10(P.Value), color = change)) +
  geom_point( alpha = 0.4, size = 1.75) +xlim(-4,4)+ylim(0,22)+
  labs(title="Volcanoplot")+
  xlab( 'Log2(Fold Change)' ) + ylab( '-Log10(P Value)' ) +
  scale_colour_manual( values = c('red','grey','blue') )+
  coord_fixed(ratio=0.5)
#添加阈值线
volcano<-volcano+geom_hline(yintercept=1.3,linetype=3)+geom_vline(xintercept=c(-0.56,0.56),linetype=3)
ggsave(volcano,filename = "Volcano.png",width = 5,height = 7)

