###format tables
library(data.table)
library(dplyr)

setwd("/net/hunt/disk2/bwolford/GBMI")

##replication
df<-fread(file="GBMI_lookup.txt")
names(df)[13]<-"pval"
df<-df[!(pval==".")]
df<-df[,c(5,6,13)] #only cols of interst aka pvalue
table(as.numeric(df$pval) < 0.05/nrow(df)) #how many replicate

##known and novel
main<-fread("/net/hunt/disk2/wukh/GBMI_bio/VTE/GBMI_VTE_IndexVariants_KnownTrait.txt",na.strings=".",fill=TRUE)

main$beta_pos=abs(main$beta)
main<-main %>% mutate(freq_pos=case_when(beta < 0~(1-freq), beta>0~freq),
                      effect_allele=case_when(beta<0~ref,beta>0~alt))

main$OR<-exp(main$beta_pos)
main$LB<-exp(main$beta_pos-(1.96*main$se))
main$UB<-exp(main$beta_pos+(1.96*main$se))
main$CI<-paste0("[",formatC(main$LB,digits=2),",",formatC(main$UB,digits=2),"]")
#main<-main[,c(1,2,3,4,5,14,16,20,21,22,23,24,25,26)]
main$OR<-formatC(main$OR,digits=2)
main$p<-formatC(main$p,digits=3)
main$VTE<-ifelse(main$VTE==0,"Potentially Novel","Known")

df2<-left_join(main,df,by=c("pos"="pos_hg38"))
df2$chr_hg38<-NULL
df2$p<-ifelse(as.numeric(df2$p)==0,2.22E-308,df2$p)

## to do: how to fix the numerical overflow and print "." in cells

#subset columns and write
df3<-df2[,c(1,2,3,4,5,10,11,14,6,15,7)]
names(df3)<-c("Chromosome","Position (hg38)","Reference Allele", "Alternate Allele","Effect Allele","Lead SNP rsID","Odds Ratio","95% CI",
              "P-value","Replication p-value","Previously identified in GWAS")
write.table(df3,file="Formatted_Table1.txt",sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)

library(ggplot2)
library(ggrepel)
anno<-fread("VTE_tophits.onetophitperlocus.txt_hg38_meta_annovar_500kb_05252021_heterPval_ancestry_noami_withaf.txt")
anno<-anno[,c(1,2,11)]
tag<-fread("/net/hunt/disk2/wukh/GBMI_bio/VTE/RaceSpecific/data/ALL_tagSNP.txt")
names(tag)[1]<-"CHR"
tag<-left_join(tag,df2,by=c("POS"="pos"))
df4<-left_join(tag,anno,by=c("POS"="pos"))
df4$gene <- gsub('NONE;', '', df4$gene)
df4$gene <- gsub(';NONE', '', df4$gene)
df4$gene <- gsub('NONE', '', df4$gene)

pdf(file="smile_plot.pdf",height=5,width=8)
ggplot(df4,aes(x=freq_pos,y=beta_pos,label=gene,color=VTE,size=-log10(inv_var_meta_p))) + scale_size_continuous(name="-log10(p-value)",range=c(3,6)) +
  theme_bw() +  labs(x="Meta-analysis Effect Allele Frequency",y="Observed GBMI\nMeta-analysis effect size") + geom_text_repel(show.legend=FALSE,max.overlaps=16) + geom_point(alpha=0.7) + coord_cartesian(xlim=c(0,1)) +
  theme(legend.position="bottom") + scale_color_manual(values=c("blue","red"),name="")
dev.off()

df4$rep<-as.numeric(df4$pval)<0.05/38
df5<-df4[df4$rep==TRUE,]
pdf(file="smile_plot_rep.pdf",height=5,width=8)
ggplot(df5,aes(x=freq_pos,y=beta_pos,label=gene,color=VTE,size=-log10(inv_var_meta_p))) + scale_size_continuous(name="-log10(p-value)",range=c(3,6)) +
  theme_bw() +  labs(x="Meta-analysis Effect Allele Frequency",y="Observed GBMI\nMeta-analysis effect size") + geom_text_repel(show.legend=FALSE,max.overlaps=16) + geom_point(alpha=0.7) + coord_cartesian(xlim=c(0,1)) +
  theme(legend.position="bottom") + scale_color_manual(values=c("blue","red"),name="")
dev.off()

##### Kristin code
library(readr)
library(readxl)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(data.table)
library(ggrepel)

# smile plot function
plot_af_effect <- function(input_data, pt_size, alpha_level, text_size){
  plot <- ggplot(input_data, aes(x = risk_allele_freq, y = risk_allele_beta, color = Category)) +
    geom_point(size=pt_size, alpha=alpha_level) +
    geom_smooth(se = FALSE, size=pt_size) +
    scale_color_manual(values = c('PotentiallyNovel' = 'orangered3', 'PreviouslyReported' = 'steelblue3'), labels=c('Potentially Novel', 'Previously Reported')) +
    labs(color='Category') +
    xlab('Risk Allele Frequency') + ylab('Risk Allele Effect Size') +
    theme_bw() +
    theme(text=element_text(size=text_size),panel.grid = element_blank()) #+
  geom_label_repel(data = input_data,
                   aes(x=risk_allele_freq, y=risk_allele_beta, color=Category, label=Gene.refGene),
                   box.padding=0.5,
                   show.legend=FALSE)
  return(plot)
}

# add risk allele columns function
add_risk_allele <- function(input_data, beta_col, AF_col, ALT_col, REF_col){
  df <- input_data %>% mutate(risk_allele = ifelse(get(beta_col) > 0, get(ALT_col), get(REF_col)),
                              risk_allele_freq = ifelse(get(beta_col) > 0, get(AF_col), 1 - get(AF_col)),
                              risk_allele_beta = ifelse(get(beta_col) > 0, get(beta_col), -get(beta_col)))
  return(df)
}

# plot

top_hits<-fread("VTE_tophits.onetophitperlocus.txt_hg38_meta_annovar_500kb_05252021_heterPval_ancestry_noami_withaf.txt")
top_hits = read_excel("~/Bothsex_inv_var_meta.tophits.onetophitperlocus.txt_hg38_meta_annovar_500kb_05252021_newKnowAsthmaList.xlsx",
                      sheet = "Asthma")
top_hits <- add_risk_allele(top_hits, 'all_inv_var_meta_beta', 'all_meta_AF', 'ALT', 'REF')

pt = 3
alpha = 0.5
text = 16
plot <- plot_af_effect(top_hits, pt, alpha, text)

plot_output_name = 'smileplot.jpeg'
ggsave(plot, height = 8, width = 12, dpi = 300, filename=plot_output_name)



# geom_text_repel(
#   aes(point.size = -log10(inv_var_meta_p)), # data point size
#   size = 5, # font size in the text labels
#   point.padding = 0, # additional padding around each point
#   min.segment.length = 0, # draw all line segments
#   max.time = 1, max.iter = 1e5, # stop after 1 second, or after 100,000 iterations
#   box.padding = 0.3 # additional padding around each text label
# )
