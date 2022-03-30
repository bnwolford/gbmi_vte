rm(list=ls())
library(data.table)
library(ggthemes)
library(ggplot2)
library(MASS)
library(IsoplotR)
#fit is a weigted linear regression model
#lowH, highH, lowV, highV for error bars
#xMaxVal, yMaxVal, xMinVal, yMinVal for axis limits

ggplotRegression <- function (fit, dat, invVar, lowH, highH, lowV, highV, xlabtext, ylabtext, pdffile, xMaxVal, yMaxVal, xMinVal, yMinVal) {
  require(ggplot2)
  #pdf(pdffile)
  #invVarF= as.factor(1/invVar)
  invVarF= as.factor(invVar)
  yorkfit=york(dat)
  intercept1=yorkfit$a[1]
  slope1=yorkfit$b[1]
  p=ggplot(fit$model, aes_string(x = names(fit$model)[2], y = names(fit$model)[1])) +
    #geom_errorbar(aes(ymin=lowV, ymax=highV, colour = invVarF), width=.005) +
    #geom_errorbarh(aes( xmin=lowH,xmax=highH, colour = invVarF), height=.005) +
    #geom_point(aes(size=4, colour = invVarF) ) +
    geom_errorbar(aes(ymin=lowV, ymax=highV), width=.005) +
    geom_errorbarh(aes( xmin=lowH,xmax=highH), height=.005) +
    #geom_smooth(method = "lm", mapping = aes(weight = invVar), fullrange=TRUE,  se=TRUE,
    #            color = "red", show.legend = FALSE)
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5,size=24), axis.text=element_text(size=20),axis.title=element_text(size=20),
          plot.margin=grid::unit(c(5,5,5,5), "mm")) +
    #geom_point(aes(size=((1/invVar)/15)) ) +
    geom_abline(intercept = intercept1, slope = slope1,  color = "red", show.legend = T, size=1.2) +
    geom_point() +
    #scale_colour_grey(start=0,end=0.6) +
    geom_abline(intercept = 0, slope=1, col="grey3", show.legend = T, linetype = "dashed") +
    xlab(xlabtext) + ylab(ylabtext)+
    #ylim(0, 0.25) + xlim(0, 0.25) +
    #geom_rangeframe() +
    scale_x_continuous(expand=c(0,0), limits=c(xMinVal-0.1,xMaxVal+0.1)) +
    scale_y_continuous(expand=c(0,0), limits=c(yMinVal-0.1,yMaxVal+0.1)) +
    coord_cartesian(xlim=c(xMinVal,xMaxVal), ylim=c(yMinVal,yMaxVal)) +
    labs(title = paste(
      "Intercept =",signif(intercept1,2 ),
      " Slope =",signif(slope1, 2)

    )
    )
  #ggsave(p,filename=pdffile, dpi=300,width = 5, height = 5) # ID will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).
  ggsave(p,filename=pdffile,dpi=300, width=6.5, height=6.5)
  return(p)
}

#setwd("/net/hunt/disk2/bwolford/GBMI")
setwd("~/2021_analysis/vte")
##replication from INVENT
inv<-fread(file="GBMI_lookup.txt")
inv$IDtemp<-paste0(inv$chr_hg38,":",inv$pos_hg38)
## all biobank meta-analysis
data<-fread("GBMI_VTE_IndexVariants_KnownTrait.txt",na.strings=".",fill=TRUE)
data$IDtemp<-paste0(data$chr,":",data$pos)
 ###sensitivity analysis meta without VU
#d2<-
data2 = merge(data, inv, by="IDtemp")

#check alleles match, they do except for the proxy variant 
which((toupper(data2$Allele1) != data2$ref | toupper(data2$Allele2) != data2$alt) & (toupper(data2$Allele2) != data2$ref | toupper(data2$Allele1) != data2$alt))

#drop proxy and the unreplicated variant
data2<-data2[data2$proxy!="yes" & data2$Effect!="."]

#keep certain columns
data3<-data2[,c("IDtemp","ref","alt","beta","se","p","Allele1","Allele2","Effect","StdErr","P-value")]
#Allele1 is effect allele from INVENT
#alt is effect allele from GBMI
data3$new_effect<-ifelse(toupper(data3$Allele2)==data3$alt,as.numeric(data3$Effect),(-1*as.numeric(data3$Effect)))
data3$Effect<-as.numeric(data3$Effect)
data3$StdErr<-as.numeric(data3$StdErr)
data3$GBMI_LB<-data3$beta-(1.96*data3$se)
data3$GBMI_UB<-data3$beta+(1.96*data3$se)
data3$GBMI_invse<-1/data3$se

data3$INV_LB<-data3$Effect-(1.96*data3$StdErr)
data3$INV_LB<-data3$Effect+(1.96*data3$StdErr)
data3$INV_invse<-1/data3$StdErr
  #betacol=paste0(b,"_beta")
  #secol=paste0(b,"_sebeta")
  #bse=secol
  #diffindex=which(data3[,betacol] < 0)
  #data3[diffindex,"beta"] = (-1)* data3[diffindex,"beta"]

  #data3[diffindex,betacol] = (-1)* data3[diffindex,betacol]

  #binvse = paste0(b,"_INVsebeta")
  #blowV = paste0(b, "_lowV")
  #data3[[blowV]] =  data3[,betacol] - data3[[secol]]*1.96
  #bhighV = paste0(b, "_highV")
  #data3[[bhighV]] =  data3[,betacol] + data3[[secol]]*1.96
  #data3[[binvse]] = 1/data3[[secol]]


  #lbose="standard_error"
  #lboinvse="invse_GWAS"
  #lbolowH = paste0("GWAS","_lowH")
  #data3[[lbolowH]] =  data3$beta - data3[[lbose]]*1.96
  #lbohighH = paste0("lbo_",b, "_highH")
  #lbohighH = paste0("GWAS", "_highH")
  #data3[[lbohighH]] =  data3$beta + data3[[lbose]]*1.96
  #data3[[lboinvse]] = 1/data3[[lbose]]

  #formula3 = as.formula(paste0("Effect_GWAS ~ ", b, "_beta"))
  formula3 = as.formula(paste0(b,"new_effect ~ beta"))
  #fit3 <- lm(formula3, weights=data3[[binvse]], data = data3)
  fit3 <- lm(formula3, data = data3)
  minVal=min(0,min(data3[[blowV]])-0.01,min(data3[[lbolowH]])-0.01)
  maxVal=max(max(data3[[bhighV]])+0.01,max(data3[[lbohighH]])+0.01 )
  cat("maxVal: ", maxVal, "\n")
  cat("minVal: ", minVal, "\n")
  #ggplotRegression(fit3, data3[[bse]],data3[[blowH]], data3[[bhighH]], data3[[lbolowV]],data3[[lbohighV]],
  #                 ylabtext=paste0("previous GWAS"),
  #                xlabtext=paste0(d),
  #                 pdffile=paste0("effect_size_",b,"_vs_previous_GWAS.png"),
  #                 xMaxVal = maxVal, yMaxVal = maxVal,
  #                 xMinVal= minVal, yMinVal= minVal)
  #dat=data3[,c("beta","standard_error",betacol,secol)]

 ggplotRegression(fit3, data3, data3[[binvse]], data3[[lbolowH]],data3[[lbohighH]], data3[[blowV]], data3[[bhighV]],
                 xlabtext=paste0("Effect sizes reported in TAGC multiancestry"),
                  ylabtext=paste0("Effect sizes in all-biobank meta-analyses"),
                 pdffile=paste0("effect_size_",b,"_vs_TAGC_multiancestry_harmonized.png"),
                  xMaxVal = maxVal, yMaxVal = maxVal,
                  xMinVal= minVal, yMinVal= minVal)
