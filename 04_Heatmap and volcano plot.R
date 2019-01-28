##SOP for GEO analysis##

#04_Heatmap and volcano plot
#Heatmap
library("gplots")
#Create the matrix needed for drawing heatmap, row=gene name, column=sample number

heatexpr = expr[nrDEG$change != 'NOT' ,]
heatexpr<-as.matrix(heatexpr)

plot_color=c(rep("green",23),rep("blue",77))
heatmap.2(heatexpr, col=greenred(75), scale="row",
          ColSideColors = plot_color,margins=c(3,2),
          main="Non-tumor vs Glioblastoma",dendrogram = "col",Colv = T,
          key=T, symkey=FALSE, density.info="none", trace="none",
          labRow = NA,labCol = NA)
#May encounter problems if saving the image directly in R markdown. Right-click and copy the image.


#Volcano
library(ggplot2)
logFC_cutoff = 0.56

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
#Add threshold line
volcano<-volcano+geom_hline(yintercept=1.3,linetype=3)+geom_vline(xintercept=c(-0.56,0.56),linetype=3)
ggsave(volcano,filename = "Volcano.png",width = 5,height = 7)

