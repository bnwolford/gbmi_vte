library(data.table)
library(googlesheets4)
library(googledrive)
library(dplyr)
library(reshape2)
library(gridExtra)
library(forcats)
library(RColorBrewer)
library(cowplot)

setwd("~/2021_analysis/vte/")

#df<-fread("Supplementary_Table3.txt")

#id<-"842021993"
#saved .xlsx as google sheet
df<-read_sheet("https://docs.google.com/spreadsheets/d/1ycGtr2s0LpYF87f6w6oGmKPP59d3IpeRpzX6D2sHshE/edit?pli=1#gid=842021993",
               sheet=3)
#gs4_auth()

df2<-df %>% select(DEPICT, POPS, "Closest Gene", eQTL, "Credible Set", Clinvar, PWAS, Integrative_priortization_gene) 
melt_df<-reshape2::melt(df2)
melt_df$value<-as.factor(melt_df$value)

sum<-df%>% select(Sum,Integrative_priortization_gene) %>%
  left_join(melt_df,by="Integrative_priortization_gene") %>% 
  filter(Integrative_priortization_gene!="NA")


pdf(file="integrative_heatmap.pdf",height=5,width=10,useDingbats=TRUE)
  x2<-ggplot(sum,aes(x=reorder(Integrative_priortization_gene,Sum,decreasing=TRUE),y=variable,fill=value)) + geom_tile() + theme_bw() +
  theme(axis.text.x = element_text(angle = 45,hjust=1),legend.position="none") + 
  scale_fill_manual(values=c("white","grey")) + 
  labs(x="Prioritized Gene",y="Evidence")
  x1<-ggplot(sum,aes(x=reorder(Integrative_priortization_gene,Sum,decreasing=TRUE),y=as.factor(0),fill=as.factor(Sum))) + 
    geom_tile() + theme_bw() + 
    theme(legend.position="top",axis.title.x=element_blank(), 
               axis.text.y=element_blank(),
                                     axis.ticks.y = element_blank(), 
          axis.title.y=element_blank(),axis.text.x=element_blank()) + 
    scale_fill_manual(values=RColorBrewer::brewer.pal(7, "YlGnBu"),name="Number of supporting\nlines of evidence")
  #grid.arrange(x1,x2,heights=c(0.3, 0.7),ncol=1) 
  plot_grid(x1, x2, align = "v", ncol=1,rel_heights=c(0.3,1))
  dev.off()
  
  
  
  ##### can't get this to work, something wrong with vector for font
  
  gold_std<-c("ADAMTS13","PTPN11","VWF","F10","F2","F3","F5","F7","F8","F9","FGA","FGG","FGB",
              "GP5","GP6","GP9","PLG","PROC","PROS1",
              "PROZ","SERPINC1","SERPINE1","SERPINE2","THBD","VKORC1")
 
  
sum<-sum%>% mutate(gene_reorder=fct_reorder(Integrative_priortization_gene,desc(Sum)))
setface<-ifelse(levels(sum$gene_reorder) %in% gold_std,"bold.italic","italic")
                

sum$variable<-recode_factor(sum$variable, "Credible Set"="Nonsynonymous variant\nin Credible Set",
                            "Clinvar"="ClinVar", "PWAS"="PWAS & Colocalization")  
  pdf(file="integrative_heatmap_bold.pdf",height=5,width=10,useDingbats=TRUE)
  x2<-ggplot(sum,aes(gene_reorder,y=variable,fill=value)) + geom_tile() + theme_bw() +
    theme(axis.text.x = element_text(face=setface,angle=45,hjust=1),
          legend.position="none") + 
    scale_fill_manual(values=c("white","grey")) + 
    labs(x="Prioritized Gene",y="Evidence")
  x1<-ggplot(sum,aes(x=gene_reorder,y=as.factor(0),fill=as.factor(Sum))) + 
    geom_tile() + theme_bw() + 
    theme(legend.position="top",axis.title.x=element_blank(), 
          axis.text.y=element_blank(),
          axis.ticks.y = element_blank(), 
          axis.title.y=element_blank(),axis.text.x=element_blank()) + 
    scale_fill_manual(values=RColorBrewer::brewer.pal(7, "YlGnBu"),
                      guide = guide_legend(reverse = TRUE),
                      name="Number of supporting\nlines of evidence")
  #grid.arrange(x1,x2,heights=c(0.3, 0.7),ncol=1) 
  plot_grid(x1, x2, align = "v", ncol=1,rel_heights=c(0.3,1))
  dev.off()
  
  ### TODO: add rectangles
  
  
  