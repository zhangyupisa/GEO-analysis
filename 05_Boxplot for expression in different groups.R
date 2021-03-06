##SOP for GEO analysis##

#05_Boxplot for expression in different groups

library(ggplot2)
library(magrittr) 
library(ggpubr)
library(ggbeeswarm)

load('expr_grouplist.Rdata')

TMP = expr_boxplot[ rownames( expr_boxplot ) == 'PVT1', ]
data = as.data.frame(TMP)
data=t(data)
group_list_boxplot=as.data.frame(group_list_boxplot)
data<-cbind(data,group_list_boxplot)
colnames(data)<-c('PVT1','group')

p1<-ggplot(data,aes(x = factor(data$group, 
                               level = c("NC", "OD", "A", "GBM")),
                    y=data$PVT1,
                    fill=data$group,
                    shape=data$group,colour=data$group))+
  labs(x='',y='PVT1 expression')+#change axis names
  theme_classic()+#change bg color to white/transparency
  theme(axis.title.y=element_text(size=18),axis.text.x=element_text(size=18,
                                                                    color="black",
                                                                    face="bold"))+
  geom_boxplot(fill="white",color="black",lwd=1,outlier.color = NULL, 
               outlier.fill = NULL)+#create boxplot and change the width of line
  geom_quasirandom(size=1.5,width=0.3,varwidth=T)+#create beeswarm
  scale_colour_manual(values=c("NC"<-"black","OD"<-"black","A"<-"black","GBM"<-"black"))+
  scale_shape_manual(values=c("NC"=23,"OD"=24,"A"=25,"GBM"=21))+
  scale_fill_manual(values=c("NC"<-"black","OD"<-"black","A"<-"black","GBM"<-"black"))+
  theme(legend.position="none")+#delete graphic symbol
  annotate("text",x=2.5,y=7.43,label="One-way P<0.001",size=5)+
  annotate("text",x=2.5,y=7.78,label="GSE4290",size=7,face="bold")+
  geom_segment(aes(x=1,y=7.3,xend=4,yend=7.3))

ggsave(p1,filename = "p1.png",width = 7,height = 5)
