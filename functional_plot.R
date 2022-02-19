library(ggplot2)
library(ggpubr)
library(dplyr)
library(readxl)
library(stringr)

df<-read_excel("Raw laser Injury Data.xlsx",skip=2)
df2<-df %>% pivot_longer(everything())
df2$name<-gsub("2nd Control","2nd_Control",df2$name)
df2$name<-gsub("2nd Injected","2nd_Injected",df2$name)
df2$name<-gsub("First Injected","1st_Injected",df2$name)
df2$name<-gsub("First Control","1st_Control",df2$name)
df3<-cbind(df2,str_split_fixed(df2$name," ", 2))
names(df3)<-c("name","value","status","gene")


pdf(file="",height=4,width=4)

dev.off()