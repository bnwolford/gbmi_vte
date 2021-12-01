
###revised from flagship github repo

library(data.table)
library(reshape)
library(ggplot2)
library(magrittr)
library(plyr)
library(dplyr)
library(optparse)
library(RColorBrewer)
library(stringr)

#data <- fread("VTE_bothsex_ss.txt") #this file is part of this sample size form
#data = data[,c(4,5,6,7)]
#data<-data %>% filter(row_number() <= n()-1)
#colnames(data) = c("pop","biobankNames", "case","control")
#nfe      UKBB 31169 379656

#set pheno
#pheno<-"VTE"
#pheno<-"Stroke"
pheno<-"HF"
library(googlesheets4)
#https://www.tidyverse.org/blog/2020/05/googlesheets4-0-2-0/
data<-read_sheet("https://docs.google.com/spreadsheets/d/1zn-vD7Xtg6GYFoRSjNKqWObcRajVwaRgcoNpZJ9HjJ0/edit#gid=797843764",sheet=2)
data<-data %>% filter(Endpoint==pheno) %>% filter(Sex=="Bothsex") %>% filter(row_number() <= n()-1)
data = data[,c(4,5,6,7)]
colnames(data) = c("pop","biobankNames", "case","control")
data$pop<-unlist(data$pop)
data$biobankNames<-unlist(data$biobankNames) #weird bug being read in as list

library(scales)
show_col(hue_pal()(7))

group.colors <- c(AFR = hue_pal()(7)[1],
                  MID = hue_pal()(7)[2],
                  AMR = hue_pal()(7)[3],
                  EAS = hue_pal()(7)[4],
                  FIN = hue_pal()(7)[5],
                  NFE = hue_pal()(7)[6],
                  SAS = hue_pal()(7)[7])

group.colors2 <- c(AFR = '#941494',
                  MID = '#33CC33',
                  AMR = '#ED1E24',
                  EAS = '#108C44',
                  NFE = '#6AA5CD',
                  FIN = '#000080',
                  SAS = '#FF9912')





bbknames_st=c("BBJ"  ,     "BioMe" ,    "BioVU"  ,   "CCPM"  ,    "CKB"    ,   "DECODE"  ,  "ESTBB"  ,
              "FinnGen" ,  "GNH"   ,    "GS"   ,     "HUNT"  ,    "Lifelines" ,"MGB"    ,   "MGI",
              "QSKIN"   ,  "TWB","UCLA"  ,    "UKBB" )
sampletype=c("Hospital/Health center-based","Hospital/Health center-based","Hospital/Health center-based",
             "Mixed","Population-based","Population-based","Population-based","Mixed","Population-based",
             "Population-based","Population-based","Population-based","Hospital/Health center-based","Hospital/Health center-based",
             "Mixed","Population-based", "Hospital/Health center-based","Population-based")
data_sampletype=data.frame(biobankNames=bbknames_st,sampletype=sampletype)


data$pop[which(data$pop == "eas")] = "EAS"
data$pop[which(data$pop == "nfe")] = "NFE"
data$pop[which(data$pop == "sas")] = "SAS"
data$pop[which(data$pop == "afr")] = "AFR"
data$pop[which(data$pop == "amr")] = "AMR"
data$pop[which(data$pop == "fin")] = "FIN"
data$pop[which(data$pop == "ami")] = "MID"
data = data[which(data$case >= 50),]
totalCase = sum(data$case)
#data$perc=100*data$case/(data$case + data$control)
#pdffile1= paste0(pheno,"_case_by_biobank.pdf")
data$totaln = (data$case + data$control)

bbklist=unique(data$biobankNames)
bbklist = bbklist[order(bbklist, decreasing=T)]

prevL = NULL
#totalN = NULL
for(j in bbklist){
  prevL = c(prevL, 100*sum(data$case[which(data$biobankNames == j)])/(sum(data$control[which(data$biobankNames == j)]) + sum(data$case[which(data$biobankNames == j)])))
  #totalN = c(totalN, (sum(data$control[which(data$biobankNames == j)]) + sum(data$case[which(data$biobankNames == j)])))
}
prevdata = data.frame(prev = prevL, bbk=bbklist)
#prevdata = prevdata[order(prevdata$prev),]
#prevdata = prevdata[order(prevdata$prev, decreasing=T),] #no order by prev
#prevdata$bbk <- factor( prevdata$bbk, levels=prevdata$bbk)

prevdata$bbk <- factor(prevdata$bbk, levels = prevdata$bbk[order(prevdata$prev)])

poplist=unique(data$pop)
prevL=NULL
nL=NULL
for (i in poplist){
  prevL = c(prevL, 100*sum(data$case[which(data$pop == i)])/(sum(data$control[which(data$pop == i)]) + sum(data$case[which(data$pop == i)])))
  nL=c(nL,sum(data$control[which(data$pop == i)]) + sum(data$case[which(data$pop == i)]))
}
prevpop=data.frame(prev=prevL,pop=poplist,n=nL)
prevpop$pop<-factor(prevpop$pop,levels=prevpop$pop[order(prevpop$prev)])
prevpop$se<-sqrt((prevpop$prev/100*(1-prevpop$prev/100))/prevpop$n)*100
pdf(file=paste0(pheno,"_ancestry_prev.pdf"),height=3,width=2)
print(ggplot(prevpop,aes(x=prev,y=pop,fill=pop)) +  geom_bar(stat = "identity")  + theme_classic()  +
  scale_fill_manual(name="Ancestry", values=group.colors2) + labs(x="VTE Prevalence",y="Ancestry") +
  geom_errorbarh(aes(xmin=prev-se, xmax=prev+se), height=.2,
                position=position_dodge(.9)) + theme(legend.position="none"))
dev.off()

caseN=NULL
for(j in bbklist){
  caseN = c(caseN, sum(data$case[which(data$biobankNames == j)]))
}
casess = data.frame(bbk=bbklist, caseN = caseN)
data$biobankNames<- factor(data$biobankNames, levels=prevdata$bbk[order(prevdata$prev)])
data = merge(data, prevdata, by.x="biobankNames", by.y="bbk",all.x=T)
data = data[order(data$prev),]

library(grid)

library(gridExtra)
data2=data.frame(biobankNames = unique(data$biobankNames))
g.mid<-ggplot(data2,aes(x=1,y=biobankNames))+geom_text(size=8, aes(label=biobankNames))+
  geom_segment(aes(x=0.94,xend=0.95,yend=biobankNames))+
  geom_segment(aes(x=1.05,xend=1.06,yend=biobankNames))+
  ggtitle("")+
  ylab(NULL)+
  scale_x_continuous(expand=c(0,0),limits=c(0.94,1.065))+
  theme(axis.title=element_blank(),
        panel.grid=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.background=element_blank(),
        axis.text.x=element_text(color=NA),
        axis.ticks.x=element_line(color=NA),
        plot.margin = unit(c(1,-1,1,-1), "mm"))


prevdata = merge(prevdata, data_sampletype, by.x="bbk", by.y="biobankNames", all.x=T)

prevdata = prevdata[order(prevdata$prev),]
#prevdata$sampletype = factor(prevdata$sampletype, levels = c("Population-based","Mixed","Hospital/Health center-based"))
prevdata$sampletype = factor(prevdata$sampletype, levels = c("Hospital/Health center-based","Mixed","Population-based"))


#my_colors <- RColorBrewer::brewer.pal(4, "Greys")[2:4]
my_colors <- RColorBrewer::brewer.pal(4, "Greys")[4:2]
g1 <-  ggplot(data = prevdata, aes(y = prev, x = bbk, fill=sampletype))+xlab(NULL) +
  geom_bar(stat = "identity") + ggtitle(paste0(pheno, ": Prevalence (%)")) +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=20),
        legend.text = element_text(size=22),legend.title = element_text(size=24),
        plot.margin = unit(c(1,-1,1,0), "mm"), plot.title = element_text(size = 24, face = "bold",hjust = 0.5)) +
  scale_y_reverse() + coord_flip() +   guides(colour = guide_legend(nrow = 3)) + expand_limits(y=10) +
  scale_x_discrete(limits = prevdata$bbk,
                   labels = function(x) str_wrap(x, width = 10)) + scale_fill_manual(values=my_colors,name="Sampling Strategy")

#

g2 <- ggplot(data = data, aes(x = biobankNames , y = case, fill = pop)) +
  geom_bar(stat = "identity") + ggtitle(paste0(pheno, ": Number of cases")) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=20),
        legend.text = element_text(size=22),legend.title = element_text(size=24),
        plot.margin = unit(c(1,0,1,-1), "mm"), text = element_text(size=20), plot.title = element_text(size = 24, face = "bold",hjust = 0.5))  +
  coord_flip()  + scale_x_discrete(limits = prevdata$bbk,
 labels = function(x) str_wrap(x, width = 10)) + scale_fill_manual(name="Ancestry", values=group.colors2)


library(cowplot)
my_legend1 <- get_legend(g1)
library(ggpubr)
pdf(file="legends_left.pdf",height=6,width=6)
print(as_ggplot(my_legend1))
dev.off()

my_legend2 <- get_legend(g2)
pdf(file="legends_right.pdf",height=6,width=6)
print(as_ggplot(my_legend2))
dev.off()

pdf(file="both_legends.pdf",height=6,width=10)
x1<-as_ggplot(my_legend1)
x2<-as_ggplot(my_legend2)
grid.arrange(x1,x2,nrow=1)
dev.off()


g1=g1+theme(panel.background = element_blank())+theme(legend.position="none")
g2=g2+theme(panel.background = element_blank())+theme(legend.position="none")
gg1 <- ggplot_gtable(ggplot_build(g1))
gg2 <- ggplot_gtable(ggplot_build(g2))
gg.mid <- ggplot_gtable(ggplot_build(g.mid))

p3 = grid.arrange(gg1,gg.mid,gg2,ncol=3,widths=c(4/9,1/9,4/9))
pdffile1= paste0(pheno,"_prevalence_Ncase_by_biobank.pdf")
ggsave(p3,filename=pdffile1, dpi=300,width = 14, height = 8)
dev.off()
