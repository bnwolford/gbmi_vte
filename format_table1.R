###format tables
library(data.table)
library(dplyr)

#setwd("/net/hunt/disk2/bwolford/GBMI")
setwd("~/2021_analysis/vte")

##replication from INVENT
df<-fread(file="GBMI_lookup.txt")
names(df)[13]<-"pval"
df$proxy<-factor(df$proxy,levels=c("no","yes"))
df<-df[!(pval==".")] #remove the variant that needed a proxy
#n<-ncol(df)
df<-df[,c(5:19)] #remove hg19 coords
df$replicate<-ifelse(as.numeric(df$pval)<5e-8,1,0) #TO DO FIX the p-value overflow for scatter plot
bonferroni<-0.05/nrow(df)
df$replicate_bonferroni<-ifelse(as.numeric(df$pval)<bonferroni,1,0)
table(as.numeric(df$pval) < 0.05/nrow(df)) #how many replicate
df<-df %>% rename_with(~paste0("replication_",.))#renamecolumns as replication_X

## this has the old known and novel designation (22 novel)
main<-fread("GBMI_VTE_IndexVariants_KnownTrait.txt",na.strings=".",fill=TRUE)
main$beta_pos=abs(main$beta)
main<-main %>% mutate(freq_pos=case_when(beta < 0~(1-freq), beta>0~freq),
                      effect_allele=case_when(beta<0~ref,beta>0~alt))

main$OR<-exp(main$beta_pos)
main$LB<-exp(main$beta_pos-(1.96*main$se))
main$UB<-exp(main$beta_pos+(1.96*main$se))
main$CI<-paste0("[",formatC(main$LB,digits=3),",",formatC(main$UB,digits=3),"]")
#main<-main[,c(1,2,3,4,5,14,16,20,21,22,23,24,25,26)]
main$OR<-formatC(main$OR,digits=3)
main$p<-formatC(main$p,digits=3)
#main$VTE<-ifelse(main$VTE==0,"Potentially Novel","Known")
main$Category<-NULL

#join
df2<-left_join(main,df,by=c("pos"="replication_pos_hg38"))

#supp table 6 from GBMI flagship (best)
st<-fread("VTE_variants_flagship_ST6.csv")
st$PMID<-as.character(st$PMID)
df3<-left_join(df2,st,by=c("pos"="POS(hg38)"))

#all the evidence about loci
info<-fread("integrative_prioritization_2021_09_24.csv")
#2 clinvar entries for rs6025 COMBINE
info[info$rsid=="rs6025",]$CLINVAR[1]<-paste0(info[info$rsid=="rs6025",]$CLINVAR[1],info[info$rsid=="rs6025",]$CLINVAR[2])
rowtodrop<-which(info$rsid=="rs6025")[2]
names(info)[which(names(info)=="gene")]<-"gene_pwas"
info<-info[-rowtodrop,]
info$LocusIndex<-seq(1:nrow(info)) #fix mis numbering of locus index
df4<-left_join(df3,info,by=c("pos"="POS"))

#this file has number of lines of evidence manually curated
info2<-fread("integrative_prioritization_2021_09_24_toplot.csv")
df5<-left_join(df4,info2,by=c("pos"="POS")) %>% dplyr::select(-c("LocusIndex.y","rsid","CHR","REF","ALT.y"))

#1 to 25 is original file from discovery, 26 to 41 is replication, 42 to 67 is the flagship table, 
#68 to 121 is the itnegrative prioritization table 
##### subset to relevant columns for Supplementary Table 3

supp<-df5 %>% dplyr::select(68,45:48,43,113,67,53:66,29:41,18,78:112,114:116)
                            
#selectively rename
names(supp)[1]<-"Locus Index"
names(supp)[13]<-"StdError"
supp$direction_ALT<-paste0("'", supp$direction_ALT,"'")
supp$PMID<-paste0("'",supp$PMID,"'")
supp$`replication_Direction(INVENT_EA,INVENT_AA,MVP_EA,MVP_AA,MVP_HIS)`<-paste0("'", supp$`replication_Direction(INVENT_EA,INVENT_AA,MVP_EA,MVP_AA,MVP_HIS)`,"'")

supp$OR<-exp(supp$`beta (ALT)`)
supp$UB<-exp(supp$`beta (ALT)`+1.96*supp$StdError)
supp$LB<-exp(supp$`beta (ALT)`-1.96*supp$StdError)

supp$replication_OR<-exp(as.numeric(supp$replication_Effect))
supp$replication_UB<-exp(as.numeric(supp$replication_Effect)+1.96*as.numeric(supp$replication_StdErr))
supp$replication_LB<-exp(as.numeric(supp$replication_Effect)-1.96*as.numeric(supp$replication_StdErr))


## to do: how to fix the numerical overflow and print "." in cells
#write
write.table(supp,file="Supplementary_Table3.txt",sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)


######################### Smile plot
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

##### Kristin Tsuo code
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
