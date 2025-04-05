library(tidyverse)
library(data.table)

#if (!requireNamespace("BiocManager", quietly = TRUE))
 # install.packages("BiocManager")

#BiocManager::install("GenomicRanges")

library(GenomicRanges)

####### oops actually we want to make sure these variants within 50Kb are significant first, so let's do intersectbed on the comand line with the summ statss

df2<-fread("VTE_Bothsex_leave_BioVU_inv_var_meta_GBMI_052021_signif_50kb.tab.gz",header=FALSE)
names(df2)<-c("tagSNP","gene_gwas","chr","posS","posE","pval")
df2<-df2 %>% filter(pval!=".") %>% select(-chr)


#read in variants
df<-fread("GBMI_VTE_IndexVariants_KnownTrait.txt",na.strings=".",fill=TRUE) %>% select(-PreviouslyReportedTraits)
df$var<-paste0(df$chr,":",df$pos)

dat<-left_join(df,df2,by=c("var"="tagSNP"))

dat<-dat %>% mutate(posS=ifelse(is.na(posS),pos,posS)) %>% mutate(posE=ifelse(is.na(posE),pos+1,posE)) %>% mutate(chr=ifelse(chr==23,"X",chr))

#granges object with 50KB upstream and downstream of the start position
obj<-makeGRangesFromDataFrame(dat,start.field="posS",end.field="posE",keep.extra.columns=TRUE)

#clinvar
vcf<-fread("/net/hunt/bwolford/annot/clinvar_20210517.vcf.gz",skip="##",select=c("#CHROM","POS","INFO"))
vcf<-vcf %>% rename(chr="#CHROM") %>% mutate(posE=POS+1)
vcf_obj<-makeGRangesFromDataFrame(vcf,start.field="POS",end.field="posE",keep.extra.columns=TRUE)

# Find overlaps
overlaps <- findOverlaps(obj, vcf_obj)

# Extract overlapping regions
intersected_regions <- pintersect(obj[queryHits(overlaps)], vcf_obj[subjectHits(overlaps)])

#pull INFO column out of metadata
metadata <- mcols(vcf_obj)[subjectHits(overlaps),]

#put clivnar metadata in the intersected file
elementMetadata(intersected_regions)[["INFO"]] <- metadata

df3<-data.frame(intersected_regions)

df4<-df3 %>% separate_wider_delim(INFO,delim=";GENEINFO=",names=c("junk","geneINFO"),too_few="align_start") %>% 
  separate_wider_delim(geneINFO,delim=":",names=c("clinvar_gene","junk2"),too_few="align_start",too_many="merge") %>% 
  select(-c("junk","hit","junk2")) %>% mutate(snp=paste0(as.character(seqnames),":",as.character(pos))) %>%
  group_by(snp) %>% select(snp,pos,gene,clinvar_gene) %>% unique() 

df5<- df4 %>% summarise(clinvar_genes=paste(clinvar_gene, collapse=", ")) %>%
  left_join(df4,by="snp") %>% select(snp,clinvar_genes,pos,gene) %>% unique()


### previous results
dat2<-fread("/net/csgspare2/spare2/hunt/bwolford/2021_analysis/vte/integrative_prioritization_2021_09_24.csv")
dat3<-dat2 %>% select(LocusIndex,CHR,POS,CLINVAR) %>% separate_wider_delim(CLINVAR,delim=";GENEINFO=",names=c("junk","geneINFO"),too_few="align_start") %>% 
  separate_wider_delim(geneINFO,delim=":",names=c("clinvar_gene_v1","junk2"),too_few="align_start",too_many="merge") %>% 
  select(CHR,POS,clinvar_gene_v1)

#intersect with new clinvar results 
clinvar_v2<-dat3%>% left_join(df5,by=c("POS"="pos")) %>% select(CHR,POS,clinvar_gene_v1,clinvar_genes,gene)

#write file 
write.csv(clinvar_v2,"clinvar_v2.csv",row.names=FALSE)



##### eQTL 

tissues<-c("Whole_Blood","Heart_Atrial_Appendage","Heart_Left_Ventricle","Artery_Aorta","Artery_Coronary","Artery_Tibial","Cells_EBV-transformed_lymphocytes")
tissue_list<-list()
for (i in 1:length(tissues)){
  fn<-paste(sep="","/net/csgspare2/spare2/hunt/bwolford/2021_analysis/vte/GTEx_Analysis_v8_eQTL/",tissues[i],".v8.signif_variant_gene_pairs.txt.gz")
  eqtl<-fread(fn)
  fn2<-paste(sep="","/net/csgspare2/spare2/hunt/bwolford/2021_analysis/vte/GTEx_Analysis_v8_eQTL/",tissues[i],".v8.egenes.txt.gz")
  gene_info<-fread(fn2,select=c("gene_id","gene_name"))
  eqtl2<-eqtl%>% left_join(gene_info,by=c("gene_id"="gene_id"))
  eqtl3<-eqtl2 %>% separate(variant_id,into=c("CHR","POS","A1","A2")) %>% mutate(CHR=gsub("chr","",CHR))
  eqtl3$tissue<-tissues[i]
  eqtl3<-eqtl3 %>% mutate(POS=as.numeric(POS)) %>% mutate(posE=POS+1)
  eqtl_obj<-makeGRangesFromDataFrame(eqtl3,start.field="POS",end.field="posE",keep.extra.columns=TRUE)
  
  # Find overlaps
  overlaps <- findOverlaps(obj, eqtl_obj)
  
  # Extract overlapping regions
  intersected_regions <- pintersect(obj[queryHits(overlaps)], eqtl_obj[subjectHits(overlaps)])
  
  #pull INFO column out of metadata
  metadata <- mcols(eqtl_obj)[subjectHits(overlaps),]
  
  metadata_2<-elementMetadata(intersected_regions)
  
  #put all metadata in the intersected file
  elementMetadata(intersected_regions)<-cbind(metadata,metadata_2)
  
  df2<-data.frame(intersected_regions) %>% mutate(snp=paste0(seqnames,":",pos))
  orig<-df%>%select(gene,chr,pos) %>%  mutate(snp=paste0(chr,":",pos))
  df3<- df2 %>% group_by(snp) %>% select(snp,pos,gene_name,tissue) %>% unique() 
  
  df4 <- df3 %>% summarise(eqtl_genes=paste(gene_name, collapse=", ")) %>% unique() 
  col_name<-paste(sep=".","eQTL_eGENE",tissues[i])
  eqtl_out<-orig %>% left_join(df4,by="snp") %>% unique() %>% select(snp,gene,eqtl_genes)  %>% rename("{col_name}":=eqtl_genes)
  tissue_list[[i]]<-eqtl_out
}

my_eqtl<-do.call(cbind,tissue_list) %>% distinct() %>% select(snp,gene,contains("eQTL_eGENE")) %>%
  unite(all_genes, starts_with("eQTL"), remove = F, sep = ";") %>% select(snp,gene,all_genes)

#write file 
write.csv(my_eqtl,"eqtl_v2.csv",row.names=FALSE)

############################
####### old code #######
############################
#check the assumption that eQTL/clinvar results are the same whether we use the lead SNP only or 50 kb +/-

#read in variants
df<-fread("GBMI_VTE_IndexVariants_KnownTrait.txt",na.strings=".",fill=TRUE) %>% select(-PreviouslyReportedTraits)

#look upstream and downstream 50KB
df$end50<-df$pos+50000
df$start50<-df$pos-50000


df <- df %>% mutate(chr=ifelse(chr==23,"X",chr))

#granges object with 50KB upstream and downstream of the start position
obj<-makeGRangesFromDataFrame(df,start.field="start50",end.field="end50",keep.extra.columns=TRUE)

#clinvar
vcf<-fread("/net/hunt/bwolford/annot/clinvar_20210517.vcf.gz",skip="##",select=c("#CHROM","POS","INFO"))
vcf<-vcf %>% rename(chr="#CHROM") %>% mutate(posE=POS+1)
vcf_obj<-makeGRangesFromDataFrame(vcf,start.field="POS",end.field="posE",keep.extra.columns=TRUE)

# Find overlaps
overlaps <- findOverlaps(obj, vcf_obj)

# Extract overlapping regions
intersected_regions <- pintersect(obj[queryHits(overlaps)], vcf_obj[subjectHits(overlaps)])

#pull INFO column out of metadata
metadata <- mcols(vcf_obj)[subjectHits(overlaps),]

#put clivnar metadata in the intersected file
elementMetadata(intersected_regions)[["INFO"]] <- metadata

df2<-data.frame(intersected_regions)

df3<-df2 %>% separate_wider_delim(INFO,delim=";GENEINFO=",names=c("junk","geneINFO"),too_few="align_start") %>% 
  separate_wider_delim(geneINFO,delim=":",names=c("clinvar_gene","junk2"),too_few="align_start",too_many="merge") %>% 
  select(-c("junk","hit","junk2")) %>% mutate(snp=paste0(as.character(seqnames),":",as.character(pos))) %>%
  group_by(snp) %>% select(snp,pos,gene,clinvar_gene) %>% unique() 

df4<- df3 %>% summarise(clinvar_genes=paste(clinvar_gene, collapse=", ")) %>%
  left_join(df3,by="snp") %>% select(snp,clinvar_genes,pos,gene) %>% unique()


### previous results
dat<-fread("/net/csgspare2/spare2/hunt/bwolford/2021_analysis/vte/integrative_prioritization_2021_09_24.csv")
dat2<-dat %>% select(LocusIndex,CHR,POS,CLINVAR) %>% separate_wider_delim(CLINVAR,delim=";GENEINFO=",names=c("junk","geneINFO"),too_few="align_start") %>% 
  separate_wider_delim(geneINFO,delim=":",names=c("clinvar_gene_v1","junk2"),too_few="align_start",too_many="merge") %>% 
  select(CHR,POS,clinvar_gene_v1)

#intersect with new clinvar results 
clinvar_v2<-dat2 %>% left_join(df4,by=c("POS"="pos")) %>% select(CHR,POS,clinvar_gene_v1,clinvar_genes,gene)

#write file 
write.csv(clinvar_v2,"clinvar_v2.csv",row.names=FALSE)


######### 


#eQTL
#redownload from https://www.gtexportal.org/home/downloads/adult-gtex/qtl
tissues<-c("Whole_Blood","Heart_Atrial_Appendage","Heart_Left_Ventricle","Artery_Aorta","Artery_Coronary","Artery_Tibial","Cells_EBV-transformed_lymphocytes")
tissue_list<-list()
for (i in 1:length(tissues)){
  fn<-paste(sep="","/net/csgspare2/spare2/hunt/bwolford/2021_analysis/vte/GTEx_Analysis_v8_eQTL/",tissues[i],".v8.signif_variant_gene_pairs.txt.gz")
  eqtl<-fread(fn)
  fn2<-paste(sep="","/net/csgspare2/spare2/hunt/bwolford/2021_analysis/vte/GTEx_Analysis_v8_eQTL/",tissues[i],".v8.egenes.txt.gz")
  gene_info<-fread(fn2,select=c("gene_id","gene_name"))
  eqtl2<-eqtl%>% left_join(gene_info,by=c("gene_id"="gene_id"))
  eqtl3<-eqtl2 %>% separate(variant_id,into=c("CHR","POS","A1","A2")) %>% mutate(CHR=gsub("chr","",CHR))
  eqtl3$tissue<-tissues[i]
  eqtl3<-eqtl3 %>% mutate(POS=as.numeric(POS)) %>% mutate(posE=POS+1)
  eqtl_obj<-makeGRangesFromDataFrame(eqtl3,start.field="POS",end.field="posE",keep.extra.columns=TRUE)
  
  # Find overlaps
  overlaps <- findOverlaps(obj, eqtl_obj)
  
  # Extract overlapping regions
  intersected_regions <- pintersect(obj[queryHits(overlaps)], eqtl_obj[subjectHits(overlaps)])
  
  #pull INFO column out of metadata
  metadata <- mcols(eqtl_obj)[subjectHits(overlaps),]
  
  metadata_2<-elementMetadata(intersected_regions)
  
  #put all metadata in the intersected file
  elementMetadata(intersected_regions)<-cbind(metadata,metadata_2)
  
  df2<-data.frame(intersected_regions) %>% mutate(snp=paste0(seqnames,":",pos))
  orig<-df%>%select(gene,chr,pos) %>%  mutate(snp=paste0(chr,":",pos))
  df3<- df2 %>% group_by(snp) %>% select(snp,pos,gene_name,tissue) %>% unique() 
  
  df4 <- df3 %>% summarise(eqtl_genes=paste(gene_name, collapse=", ")) %>% unique() 
  col_name<-paste(sep=".","eQTL_eGENE",tissues[i])
  eqtl_out<-orig %>% left_join(df4,by="snp") %>% unique() %>% select(snp,gene,eqtl_genes)  %>% rename("{col_name}":=eqtl_genes)
  tissue_list[[i]]<-eqtl_out
}


