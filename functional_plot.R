library(ggplot2)
library(ggpubr)
library(dplyr)
library(readxl)
library(stringr)
library(RColorBrewer)

df<-read_excel("Raw laser Injury Data.xlsx",skip=2)
df2<-df %>% pivot_longer(everything())
df2$name<-gsub("2nd Control","2nd_Control",df2$name)
df2$name<-gsub("2nd Injected","2nd_Injected",df2$name)
df2$name<-gsub("First Injected","1st_Injected",df2$name)
df2$name<-gsub("First Control","1st_Control",df2$name)
df3<-cbind(df2,str_split_fixed(df2$name," ", 2))
names(df3)<-c("name","value","status","gene")
df3$inj<-as.factor(df3$status)
levels(df3$inj)<-c("Control","Injected","Control","Injected","Control","Injected","Control")

df3[df3$gene=="control",]$gene<-"no sgRNA"
df3[df3$gene=="Cas9",]$gene<-"no sgRNA"

#reorder
df3$gene<-factor(df3$gene,levels=c("no sgRNA","F7","RASIP1","TC2N","STAB2","TSPAN15","PLCG2"))           

#paired
pdf(file="paired.pdf",height=4,width=8)
ggplot(df3,aes(x=gene,y=value,fill=inj)) +  geom_point(position=position_jitterdodge(),alpha=0.3,aes(color=inj)) +
 geom_boxplot(outlier.shape=NA) +
  labs(x="Targeted locus",y="Time to occlusion (seconds)") + 
  theme_bw() + scale_fill_manual(values=c("dark blue","cyan4"),name="Experimental\ncondition") + 
  scale_color_manual(values=c("dark blue","cyan4"),name="Experimental\ncondition") +
  theme(axis.text.x = element_text(angle = 45,vjust=0.75,face="italic"))
dev.off()