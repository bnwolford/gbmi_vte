rm(list=ls())
library(data.table)
library(forestplot)
setwd("/Users/wei/Documents/research/projects/GBMI/Flagship/052021/meta/summaryloci/ancestrySpecific/AFinAncestry/byancestry")


phenoData = fread("/Users/wei/Documents/research/projects/GBMI/Flagship/052021/meta/createInput/pheno_sex.list.Bothsex", header=F, data.table=F)
nloci=0
for (j in 1:nrow(phenoData)){
  pheno=phenoData[j,1]
  if(pheno!="Appendectomy" & pheno != "UtC" ){
    
  sex=phenoData[j,2]
  data = fread(paste0(pheno,"_",sex,"_all_notin_eur.intersect.bed_onetophitperlocus.eurAF.noneurAF.byancestry.txt"), data.table=F, sep="\t", header=T)
  data = data[which(data$fc > 10 | is.na(data$fc)), ]
  nloci = nloci+nrow(data[which(data$Category == "PotentiallyNovel"),])
  
  }
} 




rm(list=ls())
library(data.table)
library(forestplot)
setwd("/Users/wei/Documents/research/projects/GBMI/Flagship/052021/meta/summaryloci/ancestrySpecific/AFinAncestry/byancestry")


phenoData = fread("/Users/wei/Documents/research/projects/GBMI/Flagship/052021/meta/createInput/pheno_sex.list.Bothsex", header=F, data.table=F)

for (j in 1:nrow(phenoData)){


 
#pheno="ThC"  
#sex = "Bothsex"

pheno=phenoData[j,1]
sex=phenoData[j,2]
#if(pheno!="Appendectomy" & pheno!="Stroke" & pheno!="ThC" & pheno != "UtC" & pheno != "POAG"){
  if(pheno!="Appendectomy" & pheno != "UtC" ){
    
#pheno="Asthma"
#sex="Bothsex"
#AAA_Bothsex_all_notin_eur.intersect.bed_onetophitperlocus.eurAF.noneurAF.txt
data = fread(paste0(pheno,"_",sex,"_all_notin_eur.intersect.bed_onetophitperlocus.eurAF.noneurAF.byancestry.txt"), data.table=F, sep="\t", header=T)
data = data[which(data$fc > 10 | is.na(data$fc)), ]

if(nrow(data) > 0){
for(i in 1:nrow(data)){  
  print(i)
  print(pheno)
locusIndex = data$LocusIndex[i]
chr=data$CHR[i]
pos=data$POS[i]
gene=data$Gene.refGene[i]
gene = gsub(";","_",gene)
#pdf("COPD_2_55839380.pdf", width=10)
#data1 = data[which(data$chr == 2 & data$pos == 55839380 & data$trait == "COPD"),]

#pdf(paste0(pheno,"_",sex,"_locusIndex_", locusIndex,".pdf"), width=10)
#data1 = data[which(data$chr == 7 & data$pos == 146660507 & data$trait == "ThC"),]
data1 = data[i,]


#setwd("/Users/wei/Documents/research/projects/GBMI/Flagship/052021/meta/summaryloci")
#data = fread("COPD_Bothsex_inv_var_meta.gz_inv_var_meta.tophits.txt_hg38_meta.txt",  data.table=F, sep="\t")
ss=fread(paste0("/Users/wei/Documents/research/projects/GBMI/Flagship/052021/samplesize/byPhenobySexByPopByBBK/",pheno,"_",sex,".txt"),  data.table=F, sep="\t", header=F)
ss = ss[which(ss[,4] != 3), ]
ss = ss[order(ss[,2]),]

poplist= paste0(ss[,3], "_", ss[,2])
#poplist = poplist[which (poplist != "CCPM_nfe")]
#pdf("COPD_2_55839380_More.pdf", width=12)
#data1 = data[which(data[,1] == 2 & data$POS == 55839380),]

#pdf("ThC_7_146660507.pdf", width=10)
#data1 = data[which(data$chr == 7 & data$pos == 146660507 & data$trait == "ThC"),]

#poplist=c("afr", "amr", "eas", "fin", "nfe", "sas")
meanvec = data1[,paste0(poplist,"_beta")]
poplist=poplist[which(!is.na(meanvec))]
meanvec = c(t(as.vector(meanvec[which(!is.na(meanvec))])), data1$all_inv_var_meta_beta)
sevec = data1[,c(paste0(poplist,"_sebeta"), "all_inv_var_meta_sebeta")]
sevec = t(as.vector(sevec))
pvec = t(as.vector(data1[,c(paste0(poplist,"_pval"), "all_inv_var_meta_p")]))
pvec = formatC(pvec, format = "e", digits = 2)


lowervec = meanvec - 1.96*sevec
highervec = meanvec + 1.96*sevec
freqvec = t(as.vector(data1[,paste0(poplist,"_AF_Allele2")]))
freqvec = c(freqvec, data1$all_meta_AF)
freqvec = round(freqvec, digits=4)
cochrane_from_rmeta = data.frame(mean = c(NA,as.vector(meanvec)), lower = c(NA, as.vector(lowervec)), upper = c(NA, as.vector(highervec)))
#cochrane_from_rmeta = exp(cochrane_from_rmeta)




# if(sum(highervec > 1.5) > 0){
# 
# 
# 
#   orvec = meanvec
# 
# #orvec = exp(meanvec)
# #lowervec = exp(lowervec)
# #highervec = exp(highervec)
# xlogval=FALSE
# orvec = round(orvec, digits=4)
# 
# 
# 
# heightval = round(7*(length(freqvec)/10), digits=0)
# if(heightval < 7){heightval = 7}
# 
# 
# pdf(paste0(pheno,"_",sex,"_locusIndex_", locusIndex,".pdf"), width=10, height=heightval)
# 
# tabletext <- cbind(
#   c( "Ancestry",poplist, "Summary"),
#   c("Freq",freqvec),
#   c("BETA",orvec),
#   c("P-value", pvec))
# forestplot(tabletext, 
#            #hrzl_lines = gpar(col = "#444444"),
#            cochrane_from_rmeta,new_page = TRUE,
#            is.summary = c(TRUE,rep(FALSE,nrow(tabletext)-2),TRUE),
#            xlog=xlogval,
#            xlab="Effect size",
#            clip=c(-2,2),
#            col = fpColors(box = "royalblue",
#                           line = "darkblue",
#                           summary = "royalblue"), 
#            txt_gp = fpTxtGp(ticks=gpar(cex=2), summary = gpar(cex=2), label= gpar(cex=2), xlab= gpar(cex=2))
# 
# )
# 
# dev.off()
# }else{
  
  #orvec = exp(meanvec)
orvec=meanvec
  #lowervec = exp(lowervec)
  #highervec = exp(highervec)
  #xlogval=TRUE
xlogval=FALSE
  orvec = round(orvec, digits=4)
  
  
  
  heightval = round(7*(length(freqvec)/10), digits=0)
  if(heightval < 7){heightval = 7}
  
  
  #pdf(paste0(pheno,"_",sex,"_chr_", chr,"_pos_", pos,".pdf"), width=10, height=heightval)
  pdf(paste0(pheno,"_",sex,"_gene_", gene,".pdf"), width=11, height=heightval)
  
  tabletext <- cbind(
    c( "Biobank_Ancestry",poplist, "Summary"),
    c("Freq",freqvec),
    c("Beta",orvec),
    c("P-value", pvec))
  forestplot(tabletext, 
             #hrzl_lines = gpar(col = "#444444"),
             cochrane_from_rmeta,new_page = FALSE,
             is.summary = c(TRUE,rep(FALSE,nrow(tabletext)-2),TRUE),
             xlog=xlogval,
             xlab="Beta:log(Odds Ratio)",
            # clip=c(0.01,5),
            clip=c(-2,2),
             col = fpColors(box = "royalblue",
                            line = "darkblue",
                            summary = "royalblue"), 
            boxsize=0.1,
             txt_gp = fpTxtGp(ticks=gpar(cex=2), summary = gpar(cex=2), label= gpar(cex=2), xlab= gpar(cex=2))
             
  )
  
  dev.off()
  
  
  
  
#}
}
}
}
}