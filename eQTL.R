library(data.table)
library(tidyr)
library(dplyr)
library(gcookbook)
library(ggpubr)
library(gtable)

### script to merge eQTL with other annotation data, merge PWAS, and make heatmap

#df<-fread("/net/hunt/disk2/wukh/GBMI_bio/result/VTE/VTE_locus_cs.txt")
#df<-fread("/net/hunt/disk2/wukh/GBMI_bio/result/VTE/VTE_locus_cs_052621.txt")
df<-fread("/net/hunt/disk2/wukh/GBMI_bio/VTE/VTE_locus_cs_OldLocusFile_091021.txt")
df$CHR<-as.character(df$CHR)
df$POS<-as.character(df$POS)

tissues<-c("Whole_Blood","Heart_Atrial_Appendage","Heart_Left_Ventricle","Artery_Aorta","Artery_Coronary","Artery_Tibial","Cells_EBV-transformed_lymphocytes")

for (i in 1:length(tissues)){
    fn<-paste(sep="","/net/hunt/disk2/bwolford/GBMI/GTEx_Analysis_v8_eQTL/",tissues[i],".v8.signif_variant_gene_pairs.txt.gz")
    eqtl<-fread(fn)
    fn2<-paste(sep="","/net/hunt/disk2/bwolford/GBMI/GTEx_Analysis_v8_eQTL/",tissues[i],".v8.egenes.txt.gz")
    gene_info<-fread(fn2,select=c("gene_id","gene_name"))
    eqtl2<-eqtl%>% left_join(gene_info,by=c("gene_id"="gene_id"))
    eqtl3<-eqtl2 %>% separate(variant_id,into=c("CHR","POS","A1","A2")) %>% mutate(CHR=gsub("chr","",CHR)) %>% right_join(df,by=c("CHR"="CHR","POS"="POS"))
    eqtl4<-eqtl3 %>% group_by(SNP) %>% summarise(genes=paste(gene_name, collapse=", ")) %>% separate(SNP,into=c("CHR","POS"))
    col_name<-paste(sep=".","eQTL_eGENE",tissues[i])
    df<-df %>% left_join(eqtl4,by=c("CHR"="CHR","POS"="POS")) %>% rename("{col_name}":=genes)
}

write.table(df,"VTE_locus_cs_eqtls.txt",quote=FALSE,sep="\t",col.names=TRUE,row.names=FALSE)
print("eqtl")
#clinvar
vcf<-fread("/net/hunt/bwolford/annot/clinvar_20210517.vcf.gz",skip="##",select=c("#CHROM","POS","INFO"))
vcf$POS<-as.character(vcf$POS)
df2<-left_join(df,vcf,by=c("CHR"="#CHROM","POS"="POS")) %>% rename("CLINVAR":=INFO)
#handle repeat Clinvar
write.table(df2,"VTE_locus_cs_eqtls_clinvar.txt",quote=FALSE,sep="\t",col.names=TRUE,row.names=FALSE)
print("clinvar")


#note rs6025 has 2 entries
df2<-fread(file="VTE_locus_cs_eqtls_clinvar.txt")

#PWAS
pwas<-fread("PWAS-VTE-results.csv")
pwas<-pwas[pwas$evidence=="BOTH",] #want both pqtl and coloc evidence
pwas$rsid<-sapply(strsplit(as.character(pwas$ld_snp),"_",fixed=TRUE), '[', 1)
write.table(pwas,"PWAS-VTE-results-rsID.tab",quote=FALSE,sep="\t",col.names=TRUE,row.names=FALSE)
print("pwas")


###intersected with dbsnp to get positions
pwas<-fread("PWAS-VTE-results-rsID.tab")
rsid<-fread("rsid_lookup.txt")
rsid$CHR<-as.numeric(sapply(strsplit(sapply(strsplit(rsid$V1,".",fixed=TRUE),'[',1),"_",fixed=TRUE),'[',2))
rsid2<-rsid[,c("CHR","V2","V3")]
names(rsid2)<-c("CHR","POS","RSID") #hg38
pwas2<-left_join(pwas,rsid2,by=c("rsid"="RSID"))
pwas2<-pwas2[pwas2$CHR %in% seq(1,22),]
pwas2<-pwas2[order(pwas2$CHR)]

closest<-function(xv,sv){
    xv[which(abs(xv-sv)==min(abs(xv-sv)))]
}
for (i in 1:nrow(df2)){
    res<-pwas2[df2[i]$CHR==pwas2$CHR,c("Gene_gene","rsid","CHR","POS")]
    if (nrow(res)>0){
        res<-res[res$POS==closest(res$POS,df2[i]$POS),] #pick closest
        df2[i,"gene"]<-res$Gene_gene
        df2[i,"rsid_pwas"]<-res$rsid
        df2[i,"chr_pwas"]<-res$CHR
        df2[i,"pos_pwas"]<-res$POS
    }
}

write.table(df2,"VTE_locus_cs_eqtls_clinvar_pwas.txt",quote=FALSE,sep="\t",col.names=TRUE,row.names=FALSE)

### go through by hand and mark binary yes/no

############# make figure with heat map

dat<-fread("/net/hunt/disk2/bwolford/GBMI/integrative_prioritization_2021_09_24_toplot.csv")
dat$label<-paste0(dat$CHR,":",dat$POS,"-",dat$`Prioritized gene`)
dat<-dat[!is.na(dat$`Prioritized gene`),]
dat2<-dat[,c(7,8,9,10,11,12,13,14,19)]
dat3<-dat2 %>% gather("method","value",-c(label,Sum)) #wide to long
dat3$method<-recode(dat3$method,"Credible Set"="Nonsynonymous variant\nin Credible Set","PWAS"="PWAS & Colocalization")

dat3$label<-reorder(dat3$label,-dat3$Sum)
dat4<-dat3[dat3$method!="Sum",]
dat4[dat4$method=="POPS",]$method<-"PoPS" #reformat name

#variants included too
pdf(file="heatmap.pdf",height=4,width=8)
ggplot(dat4,aes(x=label,y=method,fill=as.factor(value))) + geom_tile() + theme_bw() +
    scale_fill_manual(values=c("white","dark grey")) + theme(axis.text.x = element_text(angle = 45, hjust=1),legend.position="none") +
    labs(y="",x="Prioritized Gene")
dev.off()

#just gene names
dat4$gene<-sapply(strsplit(as.character(dat4$label),"-",fixed=TRUE), '[',2)
dat4$gene<-reorder(dat4$gene,-dat3$Sum)
pdf(file="heatmap_gene.pdf",height=4,width=8)
ggplot(dat4,aes(x=gene,y=method,fill=as.factor(value))) + geom_tile() + theme_bw() +
    scale_fill_manual(values=c("white","dark grey")) + theme(axis.text.x = element_text(angle = 45, hjust=1,face="italic"),legend.position="none") +
    labs(y="",x="Prioritized Gene")
dev.off()

dat5<-unique(dat4[c("Sum","label")])
dat5$Sum<-as.numeric(dat5$Sum)
pdf(file="heatmap_sum.pdf",height=1,width=8)
ggplot(dat5,aes(x=label,y=1,fill=as.factor(Sum))) + geom_tile() + theme_classic() +
 theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), line=element_blank(),
       axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
       legend.position="none") +
    scale_fill_viridis_d() +
    labs(y="",x="")
dev.off()

pdf(file="legend.pdf",height=3,width=3)
p<-get_legend(ggplot(dat5,aes(x=label,y=1,fill=as.factor(Sum))) + geom_tile() + theme_classic() +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), line=element_blank(),
          axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank()) +
    scale_fill_viridis_d() + labs(y="",x=""))
as_ggplot(p)
dev.off()


#which genes are gold standard or known by GWAS
#italicize genes

known_genes<-c("F5","PLEK","PROS1","F11","F2","FGG","PROC","VWF","ABO","F7;F10","PROCR","SLC44A2","STAB2")

### combine into 1 plot
pdf(file="heatmap_all.pdf",height=6,width=8)
g1<-ggplotGrob(ggplot(dat5,aes(x=label,y=1,fill=as.factor(Sum))) + geom_tile() + theme_classic() +
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(), line=element_blank(),
          axis.title.y=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(),
          legend.position="top") +
    scale_fill_viridis_d(name="Sum") + guides(fill=guide_legend(reverse = TRUE)) +
    labs(y="",x=""))
g2<-ggplotGrob(ggplot(dat4,aes(x=gene,y=method,fill=as.factor(value))) + geom_tile() + theme_bw() +
    scale_fill_manual(values=c("white","dark grey")) +
    theme(legend.position="none",
          axis.text.x = element_text(angle = 45, hjust=1,
                                face=ifelse(levels(dat4$gene) %in% known_genes,"bold.italic","italic"))) +
     labs(y="",x="Prioritized Gene"))
g<-rbind(g1,g2,size="last")
g$widths<-unit.pmax(g1$widths,g2$widths)
panels <- g$layout$t[grep("panel", g$layout$name)]
g$heights[panels] <- unit(c(1,5), "null")
#grid.newpage()
grid.draw(g)
dev.off()
