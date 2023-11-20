#!/usr/bin/env Rscript

#from https://github.com/globalbiobankmeta/PLOTS/blob/master/plot_scripts/ManhattanPlot.r

options(stringsAsFactors=F,bitmapType='cairo')
Sys.setlocale("LC_CTYPE", "C.UTF-8")

req_packages <- c("optparse", "data.table",  "RColorBrewer", "plotrix")
#req_packages <- c("RColorBrewer")
for (pack in req_packages) {
  if(!require(pack, character.only = TRUE)) {
    #install.packages(pack, repos = "https://cloud.r-project.org")
    install.packages(pack, repos='http://cran.us.r-project.org')
  }
}


library(data.table)
library(optparse)
library(plotrix)

# Load extra R functions
source("/net/snowwhite/home/bwolford/2021_analysis/vte/gbmi_vte/helperFunctions.r")

option_list <- list(
  make_option(c("-i", "--input"), type="character", default="",
              help="full path to input file"),
  make_option(c("-p", "--prefix"), type="character", default="",
              help="prefix for output files"),
  make_option(c("-t", "--title"), type="character", default=character(0),
              help="Plot Title"),
  make_option("--pheno", type="character", default="",
              help="phenotype"),
  make_option("--SNP", type="character", default="snpid",
              help="colnames for SNP ID"),
  make_option("--CHR", type="character", default="chrom",
              help="colnames for chromosome"),
  make_option("--ALLELE1", type="character", default="Allele1",
              help="colnames for allele1"),
  make_option("--ALLELE2", type="character", default="Allele2",
              help="colnames for allele2"),
  make_option("--POS", type="character", default="pos",
              help="colnames for genome position"),
  make_option("--AC", type="character", default="ac",
              help="colnames for allele count"),
  make_option("--MAC", type="character", default="mac",
              help="colnames for minor allele count"),
  make_option("--N", type="character", default="N",
              help="samples size"),
  make_option("--PVAL", type="character", default="",
              help="p values"),
  make_option("--minMAC", type="numeric", default=2,
              help="minimal minor allele count threshold [default=2]"),
  make_option("--sigThreshold", type="numeric", default=5E-8,
              help="significance threshold for p values [default=5E-8]"),
  make_option("--knownGWASList", type="character", default="",
              help="known GWAS list provided by the user"),
  make_option("--includeGWASCatinPlot", type="logical", default=FALSE,
              help="is GWAS catalog also used for plots"),
  make_option("--knownRegionFlank", type="numeric", default=500000,
              help="known region flank [default=500000]"),
  make_option("--knownGWASRegion", type="character", default="",
              help="full path to known GWAS loci Region file"),	
  make_option("--useknownGWASRegion", type="logical", default=FALSE,
              help="whether use known GWAS Region for QQ plot by Region or not"),
  make_option("--isqqplot", type="logical", default=FALSE,
              help="whether generate qq plots"),
  make_option("--iscallambda", type="logical", default=FALSE,
              help="whether calculate lambda"),
  make_option("--issummarytable", type="logical", default=FALSE,
              help="whether generate summary table"),
  make_option("--isannovar", type="logical", default=FALSE,
              help="whether annotate table using annovar"),
  make_option("--ismanhattanplot", type="logical", default=FALSE,
              help="whether generate manhattan plot"),
  make_option("--top.size", type="numeric", default=0.125,
              help="top size = proportion of total length y axis [default=0.125]"),
  make_option("--break.top", type="numeric", default=15,
              help="set axis break at -log10(P) [default=15]"),
  make_option("--qq.width", type="numeric", default=900,
              help="Width QQ plot in pixel [default=900]"),
  make_option("--qq.height", type="numeric", default=900,
              help="Height QQ plot in pixel [default=900]"),
  make_option("--mh.width", type="numeric", default=1600,
              help="Width Manhattan plot in pixel [default=1600]"),
  make_option("--mh.height", type="numeric", default=900,
              help="Height Manhattan plot in pixel [default=900]"),
  make_option("--pointsize", type="numeric", default=16,
              help="Point size of plots [default=16]")
  
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options

print(opt)

sigThreshold <- opt$sigThreshold
yLine <- -log10(sigThreshold)
chrcol <- opt$CHR
poscol <- opt$POS
A1col <- opt$ALLELE1
A2col <- opt$ALLELE2
isqqplot <- opt$isqqplot
iscallambda <- opt$iscallambda
ismanhattanplot <- opt$ismanhattanplot
issummarytable <- opt$issummarytable
isannovar <- opt$isannovar
pheno <- opt$pheno
pvalcol <- opt$PVAL

# Settings for plot
plotLine <- T
colLine <- "red"
ycol <- "log10P"
file_imp <- opt$input
prefix <- opt$prefix
maintitle <- gsub("_"," ",opt$title)
top.size <- opt$top.size
break.top <- opt$break.top
qq.width <- opt$qq.width
qq.height <- opt$qq.height
mh.width <- opt$mh.width
mh.height <- opt$mh.height
pointsize <- opt$pointsize
freqthres <- opt$minMAF
phewas.extract <- opt$phewas.extract
knownGWASRegionFile <- opt$knownGWASRegion
useknownGWASRegion <- opt$useknownGWASRegion
knownGWASList <- opt$knownGWASList
includeGWASCatinPlot <- opt$includeGWASCatinPlot
flank <- as.numeric(opt$knownRegionFlank)


#outdir <- gsub("(.+)/[^/]+$","\\1",prefix)
filename <- gsub(".+/([^/]+$)","\\1",prefix)

#if(!file.exists(outdir)) dir.create(outdir)
#file_qq_byMAF <- paste0(outdir,"/",filename,".byMAF_QQ_Plot.png")
#file_qq_byRegion <- paste0(outdir,"/",filename,".byRegion_QQ_Plot.png")
#file_mh <- paste0(outdir,"/",filename,".Manhattan.png")
#file_lambda <- paste0(outdir,"/",filename,".Lambda.txt")
#file_summarytable <- paste0(outdir,"/",filename,".regions.txt")
#file_tophits <- paste0(outdir,"/",filename,".tophits.txt")
#file_tophits_anno <- paste0(outdir,"/",filename,".tophits.annovar.txt")
#file_maf <- paste0(outdir,"/MAF/",filename,".MAF.txt")
#file_pgc <- paste0(outdir,"/LOG10PGC/",filename,".LOG10PGC.txt")
file_mh <- paste0(filename,".Manhattan.png")
file_tophits <- paste0(filename,".tophits.txt")
file_summarytable <- paste0(filename,".regions.txt")


if( grepl(".gz$",file_imp) | grepl(".bgz$",file_imp) ) {
  resIn_1 = fread(cmd=paste0("gunzip -c ", file_imp), header=T, select=c(chrcol, poscol, A1col, A2col,pvalcol))
} else {
  resIn_1 <- fread(file_imp, header=T, select=c(chrcol, poscol, A1col, A2col,pvalcol))
}


#resIn_1 <- fread(file_imp, header=T)
names(resIn_1)



print("names(resIn)")


#resIn = resIn_1[, c(SNPcol,chrcol, poscol, A1col, A2col,pvalcol), with=FALSE]
#setnames(resIn, c(SNPcol, chrcol, poscol, A1col, A2col,pvalcol), c("SNP","CHR", "BP", "ALLELE1", "ALLELE2","PVAL"))
resIn = resIn_1[, c(chrcol, poscol, A1col, A2col,pvalcol), with=FALSE]
setnames(resIn, c(chrcol, poscol, A1col, A2col,pvalcol), c("CHR", "BP", "ALLELE1", "ALLELE2","PVAL"))
print("OKKKKKKKKKKKK2")
resIn$ALLELE1 = toupper(resIn$ALLELE1)
resIn$ALLELE2 = toupper(resIn$ALLELE2)
print("OKKKKKKKKKKKK1")
resIn$SNP = paste0(resIn$CHR,"_", resIn$BP, "_", resIn$ALLELE1,"/",resIn$ALLELE2)
resIn$CHR <- gsub("chr","",resIn$CHR)
resIn$CHR <- as.numeric(gsub("X","23",resIn$CHR))
print("OKKKKKKKKKKKK")
minMAC = opt$minMAC

# needed second position to create data.table interval
resIn$BP2 <- resIn$BP

Nmarkers = nrow(resIn)
cat("nrow(resIn) ", nrow(resIn), "\n")
resIn$PVAL = as.numeric(resIn$PVAL)
resIn$PVAL[which(resIn$PVAL < 1*10^-300)] = 1*10^-300
print(min(resIn$PVAL))
print(max(resIn$PVAL))
resIn$log10P = -log10(resIn$PVAL)
print(length(resIn$log10P))


# Wei: remove the poorly imputed variants
####################################
# extract candidate hits and candidate regions before / for lambda correction
#####################################

#print("OK4")
if(useknownGWASRegion == FALSE){
  if(knownGWASList!=""){
    GWAScat1 = data.frame(read.table(knownGWASList, sep="\t", header=F, colClasses = c("character")))
  }else{
    GWAScat1 = data.frame(CHROMPOS=c("1:-999999"))
  }
  #GWAScat = cbind(GWAScat1, as.character(pheno))
  GWAScat = cbind(GWAScat1, "pheno")
  GWAScat = data.frame(GWAScat)
  colnames(GWAScat) = c("CHROMPOS","GWAScat_Trait")					
  GWAScat$CHROMPOS <- gsub("X","23",GWAScat$CHROMPOS)
  splitGWAS = matrix(as.numeric(do.call('rbind', strsplit(as.character(GWAScat$CHROMPOS),':',fixed=TRUE))), ncol=2)
  GWAScat$CHROM = splitGWAS[,1]
  GWAScat$POS = splitGWAS[,2]
  GWAScat = GWAScat[complete.cases(GWAScat),]
  
  knownGWAS <- data.table(
    CHROM=GWAScat$CHROM,
    BEGIN=GWAScat$POS-flank,
    END=GWAScat$POS+flank,key=c("CHROM", "BEGIN", "END"))
  
  
  print("OK7")
}else{
  knownGWASRegion = data.frame(read.table(knownGWASRegionFile, sep="\t", header=T, colClasses = c("numeric", "numeric", "numeric")))
  knownGWAS <- data.table(
    CHROM=knownGWASRegion$CHROM,
    BEGIN=knownGWASRegion$BEGIN,
    END=knownGWASRegion$END,key=c("CHROM", "BEGIN", "END"))
}


withinKnownRegion_v1 <- !is.na(foverlaps(resIn, knownGWAS,
                                         by.x=c("CHR", "BP", "BP2"),
                                         by.y=c("CHROM", "BEGIN", "END"),
                                         type="any", which=TRUE, mult="first"))

print("withinKnownRegion")
print(dim(resIn[withinKnownRegion_v1,]))

sum(!withinKnownRegion_v1)


if(isqqplot == TRUE){
  plotdata1b<-qqplotdata(as.numeric(resIn$log10P),Lambda=1)
  #print("OK10")
  cat(dim(plotdata1b[[1]]), "dim(plotdata1b)\n")
  
  plotdata2b<-qqplotdata(resIn$log10P[which(!withinKnownRegion_v1)],Lambda=1)
  cat(dim(plotdata2b[[1]]), "dim(plotdata2b)\n")
  legendtext <- paste(c("All","outside of +- 1Mb of Known GWAS loci"),"; N = ",c(
    format(length(resIn$log10P),big.mark=","),
    format(length(resIn$log10P[which(!withinKnownRegion_v1)]),big.mark=",")),sep="")
  
  #SNP,CHR,BP,ALLELE1,ALLELE2,MAF
  print(legendtext)
  png(filename = file_qq_byRegion, width = qq.width, height = qq.height, pointsize = pointsize)
  plotqq(plotdata1b[[1]],plotdata2b[[1]],legendtext=legendtext,maintitle=maintitle,break.top=break.top,top.size=top.size)
  #plotqq(plotdata1b[[1]],plotdata2b[[1]],legendtext=legendtext,maintitle=maintitle,break.top=break.top,top.size=top.size, plotdata1v2b[[1]], plotdata2v2b[[1]])
  
  dev.off()
}


if(ismanhattanplot == TRUE){
  issummarytable = TRUE
}

if(issummarytable == TRUE){
  tophits <- which(resIn$log10P > -log10(sigThreshold))
  if(length(tophits)>0){
    tophits <- resIn[tophits,]
    tophits$numCHR <- as.numeric(gsub("X","23",tophits$CHR))
    x <- as.numeric(tophits$BP)
    y <- tophits$numCHR
    start = c(1, which(diff(y) != 0 | diff(x) <= (-1)*flank | diff(x) >= flank) + 1)
    end = c(start - 1, length(x))	
    
    highlightRegions <- data.frame(
      'Nearest_Gene'=NA,
      'CHROM'=tophits$CHR[start],
      'START'=tophits$BP[start] - flank ,
      'END'=tophits$BP[end] + flank,
      'COLOR'="green3",
      'STATUS'="potentially_novel",
      'TOP_SNP'=NA,
      'minPVALUE'=NA,
      'Reported_GWAS_Catalog'=NA,
      'Reported_disGen'=NA
    )
    print("highlightRegions")
    
    # Extract top SNPs from each hit region
    hits <- data.table(
      CHROM=highlightRegions$CHROM,
      BEGIN=highlightRegions$START,
      END=highlightRegions$END, key = c("CHROM", "BEGIN", "END"))
    
    print("hits")
    
    regionHITS = NULL
    
    
    for(a in 1:dim(highlightRegions)[1]){
      withinKnownRegion <- !is.na(foverlaps(resIn, hits[a,], 
                                            by.x=c("CHR", "BP", "BP2"),
                                            by.y=c("CHROM", "BEGIN", "END"),
                                            type="any", which=TRUE, mult="first"))
      
      regionhits <- resIn[which(withinKnownRegion & resIn$log10P > -log10(1E-4)),]
      
      regionhits$CHROMPOS = paste0(regionhits$CHR, ":", regionhits$BP)
      knowncate = rep("potential_novel", nrow(regionhits))
      
      knownHit = which(withinKnownRegion & resIn$log10P > -log10(1E-4) & withinKnownRegion_v1)
      
      if(length(knownHit) > 0){
        highlightRegions$'STATUS'[a] <- "known"
        highlightRegions$'COLOR'[a] <- "blue"
        knowncate = rep("known", nrow(regionhits))
      }
      
      regionhits$knowncate = knowncate
      regionHITS = rbind(regionHITS, regionhits)
    }
    
    write.table(highlightRegions, file_summarytable, col.names=T, row.names=F, quote=F)
    
    print("head(regionHITS)")
    print(head(regionHITS))
    #SNP CHR BP ALLELE1 ALLELE2 MAF N PVAL MAC BP2 log10P
    resInkeepAnno = regionHITS[,c("CHR", "BP", "BP", "ALLELE1", "ALLELE2","PVAL", "log10P", "knowncate")]
    write.table(resInkeepAnno, file_tophits, col.names=F, row.names=F, quote=F, sep="\t")
    if(isannovar){	
      cmd = paste0("perl /net/dumbo/home/zhowei/tools/annovar/table_annovar.pl ",file_tophits," /net/dumbo/home/zhowei/tools/annovar/humandb/ -buildver hg19 -out ", file_tophits_anno," -remove -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a -operation gx,r,f,f,f -nastring NA  -polish -xref /net/dumbo/home/zhowei/tools/annovar/example/gene_xref.txt -otherinfo")
      system(cmd)
    }
  }else{
    resInkeepAnno = c("CHR", "BP", "BP", "ALLELE1", "ALLELE2","PVAL", "log10P", "knowncate")
    write.table(t(resInkeepAnno), file_tophits, col.names=F, row.names=F, quote=F, sep="\t")
    highlightRegions <- data.frame(
      'Nearest_Gene'=NA,
      'CHROM'=NA,
      'START'=NA,
      'END'=NA,
      'COLOR'="NA",
      'STATUS'="NA",
      'TOP_SNP'=NA,
      'minPVALUE'=NA,
      'Reported_GWAS_Catalog'=NA,
      'Reported_disGen'=NA
    )	
    write.table(highlightRegions, file_summarytable, col.names=T, row.names=F, quote=F)
    rm(resInkeepAnno)
    rm(highlightRegions)
    
    
  }
}



if(ismanhattanplot == TRUE){
  # keep all SNPs within the candidate regions
  if(exists("hits")){
    withinKnownRegion <- !is.na(foverlaps(resIn, hits,
                                          by.x=c("CHR", "BP", "BP2"),
                                          by.y=c("CHROM", "BEGIN", "END"),
                                          type="any", which=TRUE, mult="first"))
    keep <- which(resIn[[ycol]] > -log10(1E-4) | withinKnownRegion)
  } else {
    keep <- which(resIn[[ycol]] > -log10(1E-4))
  }
  
  
  # Thinning 
  prethin0 <- which(resIn[[ycol]] < -log10(0.1))
  prethin0 <- sample(prethin0,length(prethin0)/10)
  prethin1 <- which(resIn[[ycol]] < -log10(0.05) & resIn[[ycol]] >= -log10(0.1))
  prethin2 <- which(resIn[[ycol]] < -log10(1E-4) & resIn[[ycol]] >= -log10(0.05))
  
  
  
  chrcol="CHR"
  poscol="BP"
  thinned0 <- prethin0[which(!duplicated(data.frame(resIn[[chrcol]][prethin0],round(resIn[[poscol]][prethin0]/100000),round(resIn[[ycol]][prethin0],1))))]
  thinned1 <- prethin1[which(!duplicated(data.frame(resIn[[chrcol]][prethin1],round(resIn[[poscol]][prethin1]/100000),round(resIn[[ycol]][prethin1],1))))]
  thinned2 <- prethin2[which(!duplicated(data.frame(resIn[[chrcol]][prethin2],round(resIn[[poscol]][prethin2]/10000),round(resIn[[ycol]][prethin2],2))))]
  thinned <- sort(unique(c(thinned0,thinned1,thinned2,keep)))
  
  #thinned = 1:nrow(resIn)
  
  
  CHR <- resIn[[chrcol]][thinned]
  POS <- resIn[[poscol]][thinned]
  log10P <- resIn[[ycol]][thinned]
  
  chrs <- c(1:22,"X",23,"Y",24,"XY",25,"MT",26)
  chrs <- chrs[which(chrs %in% CHR)]
  chrNr <- as.numeric(as.character(factor(CHR,levels=chrs,labels=1:length(chrs))))
  chrColors <- c("grey40","grey60")[(1:length(chrs)%%2)+1]
  
  names(chrColors) <- chrs
  
  plotdata <- data.frame(
    CHR,
    POS,
    log10P,
    'plotPos'=NA,
    'chrNr'= chrNr,
    pch=20,
    highlightColor=NA,
    pcol=chrColors[CHR],
    'POS2'=POS,
    check.names=F)
  
  endPos <- 0
  plotPos <- numeric(0)
  chrLab <- numeric(0)
  chrGAP <- 1E7
  for(chr in chrs){
    chrTemp <- which(CHR == chr)
    chrPOS <- POS[chrTemp]-min(POS[chrTemp],na.rm=T)+endPos+1
    chrLab <- c(chrLab,mean(chrPOS,na.rm=T))
    endPos <- max(chrPOS,na.rm=T)+chrGAP
    plotPos <- c(plotPos,chrPOS)
  }
  
  chrName <- c(1:22,"X","Y","XY","MT","X","Y","XY","MT")
  names(chrName) <- as.character(c(1:26,"X","Y","XY","MT"))
  chrs <- chrName[chrs]
  
  plotdata$plotPos <- plotPos
  
  print("OK14b")
  
  if(exists("highlightRegions")){
    for(a in 1:dim(highlightRegions)[1]){ 
      # Set colors for plot
      setDT(plotdata,key,c("CHR","POS","POS2"))
      overlap <- which(!is.na(foverlaps(plotdata, hits[a,], 
                                        by.x=c("CHR", "POS", "POS2"),
                                        by.y=c("CHROM", "BEGIN", "END"),
                                        type="any", which=TRUE, mult="first")))
      plotdata$highlightColor[overlap] <- highlightRegions$COLOR[a]
      plotdata$pcol[overlap] <- NA
    }
  }
  
  
  print("OK15")
  
  
  #print(highlightRegions)
  png(filename = file_mh, width = mh.width, height = mh.height, pointsize = pointsize)
  par(mar=c(5.1,5.1,4.1,1.1),las=1)
  x = plotdata$plotPos
  y = plotdata$log10P
  maxY <- max(y,na.rm=T)
  #if(setYmax != F) maxY <- setYmax
  
  if(maxY > break.top/0.75){
    # Manhattan plot with two different y axis scales
    
    # set axis labels of both scales
    lab1 <- pretty(c(0,break.top),n=ceiling(12 * (1-top.size)))
    lab1 <- c(lab1[lab1 < break.top],break.top)
    lab2 <- pretty(c(break.top,maxY),n=max(3,floor(12 * top.size)))
    lab2 <- lab2[lab2 > max(lab1)]
    
    # resulting range of top scale in bottom scale units
    top.range = break.top/(1 - top.size) - break.top
    #top.data = max(y,na.rm=T)-break.top
    top.data = max(lab2)-break.top
    # function to rescale the top part
    rescale = function(y) { break.top+(y-break.top)/(top.data/top.range)}
    
    # plot bottom part / rescaled
    #print("plotdata")
    #print( plotdata[1:1000,])
    
    plot(x[y<break.top],y[y<break.top],ylim=c(0,break.top+top.range),axes=FALSE,
         pch=plotdata$pch[y<break.top], cex=1.5,cex.lab=2.0,cex.axis=2.0, xaxt="n",
         col=plotdata$pcol[y<break.top], ylab=expression(-log[10]*italic(P)), xlab="",bty="n",
         main=paste(maintitle," / ",format(Nmarkers,big.mark=",",scientific=F)," Variants",sep=""))
    #col=plotdata$pcol[y<break.top], ylab=expression(-log[10]*italic(P[GC])), xlab="",bty="n",
    # plot top part
    points(x[y>break.top],rescale(y[y>break.top]),pch=plotdata$pch[y>break.top],
           col=plotdata$pcol[y>break.top],cex=1.5)
    
    # plot highlighted regions
    for(hcol in unique(plotdata$highlightColor)){
      topDot <- plotdata[which(plotdata$highlightColor == hcol & y>break.top),]
      if(length(topDot)>0) {
        points(topDot$plotPos,rescale(topDot$log10P),pch=20,col=hcol, cex=1.5)
      }
      bottomDot <- plotdata[which(plotdata$highlightColor == hcol & y<=break.top),]
      if(length(bottomDot)>0) {
        points(bottomDot$plotPos,bottomDot$log10P,pch=20,col=hcol, cex=1.5)
      }
    }
    
    # add axes and axis labels
    axis(1,at=chrLab[seq(1,length(chrLab),by=2)],
         labels=chrs[1:length(chrLab)][seq(1,length(chrLab),by=2)],
         las=1,tick=F,cex.axis=2.0,cex.lab=2.0,line=2)
    axis(1,at=chrLab[seq(2,length(chrLab),by=2)],
         labels=chrs[1:length(chrLab)][seq(2,length(chrLab),by=2)],
         las=1,tick=F,cex.axis=2.0,cex.lab=2.0,line=0)
    axis(side=2,at=lab1,cex.axis=2.0,cex.lab=2.0)
    axis(side=2,at=rescale(lab2),labels=lab2,cex.axis=2.0,cex.lab=2.0)
    
    # plot axis breaks and indicate line of axis break
    box()
    axis.break(axis=2,breakpos=break.top,style="zigzag",brw=0.02)
    axis.break(axis=4,breakpos=break.top,style="zigzag",brw=0.02)
    abline(h=break.top,lwd=2.0,lty=2,col="grey")
    if(plotLine) {
      for(rl in 1:length(yLine)){
        if(yLine[rl] <= break.top) {
          abline(h=yLine[rl],lwd=2.0,col=colLine[rl],lty=2)
        } else {
          abline(h=rescale(yLine[rl]),lwd=2.0,col=colLine[rl],lty=2)
        }
      }
    }
  } else {
    # Standard Manhattan plot if no association signal above break
    plot(x,y,xaxt="n",pch=plotdata$pch,cex=1.5,cex.lab=2.0,cex.axis=2.0,xaxt="n",
         col=plotdata$pcol,ylab=expression(-log[10]*italic(P)),xlab="",bty="o",
         main=paste(maintitle," / ",format(Nmarkers,big.mark=",",scientific=F)," Variants",sep=""),ylim=c(0,maxY/0.9))
    #col=plotdata$pcol,ylab=expression(-log[10]*italic(P[GC])),xlab="",bty="o",
    #main=paste(maintitle," / ",format(Nmarkers,big.mark=",",scientific=F)," Variants",sep=""),ylim=c(0,break.top/0.75))
    axis(1,at=chrLab[seq(1,length(chrLab),by=2)],
         labels=chrs[1:length(chrLab)][seq(1,length(chrLab),by=2)],
         las=1,tick=F,cex.axis=2.0,cex.lab=2.0,line=2)
    axis(1,at=chrLab[seq(2,length(chrLab),by=2)],
         labels=chrs[1:length(chrLab)][seq(2,length(chrLab),by=2)],
         las=1,tick=F,cex.axis=2.0,cex.lab=2.0,line=0)
    # plot highlighted regions
    for(hcol in unique(plotdata$highlightColor)){
      extraDot <- plotdata[which(plotdata$highlightColor == hcol),]
      points(extraDot$plotPos,extraDot$log10P,pch=20,col=hcol, cex=1.5)
    }
    if(plotLine) abline(h=yLine,lwd=2.0,col=colLine,lty=2)
  }
  
  addLegend <- which(c("blue","green3") %in% unique(plotdata$highlightColor))
  
  #if(length(addLegend) > 0){	
  #	legend("topleft",c("Known Loci","Potentially Novel Loci")[addLegend],col=c("blue","green3")[addLegend],pch=15,bty="n")
  #}
  dev.off()
}
