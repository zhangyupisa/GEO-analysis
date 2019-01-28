##SOP for GEO analysis##

#06_Survival analysis

##--------------------------------------------------
#Get data needed from GEO
library(GEOquery)
library(openxlsx)
sset = getGEO(filename="./GSE43378_series_matrix.txt.gz",AnnotGPL=T)
#"s" denotes survive 
sexpr<-read.xlsx("GSE43378_series_matrix.xlsx",sheet=1,rowNames = TRUE,colNames = TRUE)
dim(sexpr)#50 samples

sdata<-pData(sset)
#survival time is recorded in line No.46
colnames(sdata)[46]<-c("survival time")
colnames(sdata)[2]<-c('ID')

##--------------------------------------------------
#Creatr survival object
library(survival)
library(dplyr)
library(survminer)

#According to file lncRNA.xlsx, PVT1 and HAR1A' IDs are 1558290_a_at,1557098_s_at.

row2<-as.data.frame(as.numeric(sdata[,46]))
colnames(row2)<-c("time")

row3<-t(as.data.frame(sexpr[ '1558290_a_at', ]))
colnames(row3)<-c('PVT1')

row4<-t(as.data.frame(sexpr[ '1557098_s_at', ]))
colnames(row4)<-c('HAR1A')

surv<-cbind(row2,row3,row4)

surv$PVT1_status<-as.factor(ifelse(surv$PVT1<mean(surv$PVT1),'Low PVT1 expression',
                                   'High PVT1 expression'))
surv$HAR1A_status<-as.factor(ifelse(surv$HAR1A<mean(surv$HAR1A),'Low HAR1A expression',
                                    'High HAR1A expression'))
status<-as.factor(ifelse(sdata$characteristics_ch1.6=="outcome: ALIVE",0,1))
surv$status<-as.numeric(status)
#1 denotes survive, 2 denotes dead
surv

##--------------------------------------------------
#Draw Kaplan-Meier survival curve
s_PVT1 <- Surv(time = surv$time, event = surv$status)
km1<-survfit(s_PVT1~PVT1_status,data = surv)
s1<-ggsurvplot(km1,legend =c(0.78,0.91), 
               xlab="Time(day)",ylab="Percent survival",
               title="GSE43378",
               pval=TRUE,pval.size=3.6,pval.coord=c(2500,0.8),
               legend.title = '',
               legend.labs = c("High PVT1 expression", "Low PVT1 expression"),
               palette=c("red","blue"))+
  guides(colour=guide_legend(nrow=2))
s1<-s1$plot + theme(plot.title = element_text(hjust = 0.5, face = "bold"))

ggsave(s1,filename = "PVT1.png",width = 7,height = 5)
