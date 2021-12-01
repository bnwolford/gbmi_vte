library(plyr)
library(dplyr)
library(ggplot2)
library(see)
library(data.table)
#library(ggridges)
library(gridExtra)
library(RNOmni)
library(optparse)
library(ResourceSelection)
library(Hmisc)
library(corrplot)
#library(GenomicRanges)
library(ltm)
library(pROC)
library(RColorBrewer)
library(DescTools)
library(rcompanion)
library(LaCroixColoR)
library(cvAUC)
library(grid)
library(brant)
library(MASS)
library(nnet)
library(lmtest)

print(Sys.time())
print(sessionInfo())

#pwd("/net/hunt/disk2/bwolford/UKBB_HUNT_proxy/UKBB/GPS/CAD")
cad<-fread("UKBB_CoronaryArteryDisease_PRS_LDpred_rho0.001_allchr.scores_pheno.txt")
t2d<-fread("/net/hunt/disk2/bwolford/UKBB_HUNT_proxy/UKBB/GPS/T2D/UKBB_Type2Diabetes_PRS_LDpred_rho0.01_allchr.scores_pheno.txt")

batch<-6 #use array
fh<-34
birthYear<-8
pcs<-c(9,10,11,12)
sex<-7
pheno<-21 #ischemic heart disease
out<-"UKBB.CAD"
grs<-3
dat<-cad


batch<-6 #use array
fh<-69
birthYear<-8
pcs<-c(9,10,11,12)
sex<-7
pheno<-26
out<-"UKBB.T2D"
grs<-3
#dat<-cad
dat<-t2d


cut_age<-65
cutpt<-0.9
quantiles<-c(20) ##  To do: test with more
censor<-FALSE
young<-FALSEout<-out
strat_col<-fh
pheno_col<-pheno
grs_col<-grs
dig<-3
print(out)

dat$partAge<-apply(X=data.frame(dat$age0,dat$age1,dat$age2), MARGIN=1, FUN=max,na.rm=TRUE)
part_age<-which(names(dat)=="partAge")

subset<-dat[!is.na(dat[[strat_col]])] #remove if NA for stratum
subset[[strat_col]]<-as.factor(subset[[strat_col]])
print(paste("Data dimensions after removing samples with NA stratum:",dim(subset)[1],dim(subset)[2]))

subset<-subset[!is.na(subset[[pheno_col]])] #remove if NA for pheno
print(paste("Data dimensions after removing samples with NA phenotype:", dim(subset)[1],dim(subset)[2]))


subset$invNormGRS<-RankNorm(subset[[grs_col]])
igrs_col<-which(names(subset)=="invNormGRS")
formula<-as.formula(paste(colnames(subset)[strat_col], "~",
                          paste(colnames(subset)[igrs_col], collapse = "+"),
                          sep = ""))
obj<-glm(formula,data=subset,family="binomial")
print(format(coef(summary(obj))[8],digits=3)) #pvalue
OR<-exp(coef(summary(obj))[2]) #OR
CI<-c(exp(coef(summary(obj))[2]-(1.96*coef(summary(obj))[4])),exp(coef(summary(obj))[2]+(1.96*coef(summary(obj))[4])))
print(format(OR,digits=3))
print(format(CI,digits=3))

qsub<-subset

qsub$age<-2021-qsub[[birthYear]]
  qsub$age_std<-data.frame(scale(qsub$age))
  qsub$part_age_sq<-qsub[[part_age]]^2
part_age_sq<-which(names(qsub)=="part_age_sq") #participation age squared


 qsub$part_age_std<-data.frame(scale(qsub[[part_age]])) #standardize so comparable to binary and inverse normalized GRS
  qsub$part_age_sq_std<-data.frame(scale(qsub[[part_age_sq]]))
  qsub$birthYear_std<-data.frame(scale(qsub[[birthYear]]))
  part_age_sq_std<-which(names(qsub)=="part_age_sq_std")
  part_age_std<-which(names(qsub)=="part_age_std")
  birthYear_std<-which(names(qsub)=="birthYear_std")
age_std<-which(names(qsub)=="age_std")

birthYear_std<-age_std

  covar<-c(part_age_std,birthYear_std,pcs,sex,batch,part_age_sq_std)
  ## model for family history + GRS with interaction term
  qsub$invNormGRS<-RankNorm(qsub[[grs_col]])
  igrs_col<-which(names(qsub)=="invNormGRS") #new GRS col
  formula<-as.formula(paste(colnames(qsub)[pheno_col], "~",
                            paste(colnames(qsub)[c(covar)], collapse = "+"), "+" ,
                            paste(colnames(qsub)[c(strat_col,igrs_col)],collapse="*"),
                            sep = ""))
  best<-glm(formula=formula,data=qsub,family="binomial")
  best_df<-data.frame(summary(best)$coefficients)
  best_df$OR<-exp(best_df$Estimate)
  best_df$LB<-exp(best_df$Estimate-(1.96*best_df$Std..Error))
  best_df$UB<-exp(best_df$Estimate+(1.96*best_df$Std..Error))
  options(scipen=99)
format(best_df,digits=3)


 base_mod<-glm(formula=as.formula(paste(colnames(qsub)[pheno_col], "~",
                            paste(colnames(qsub)[c(covar)], collapse = "+"),
                            sep = "")),data=qsub,family="binomial")
  grs_mod<-glm(as.formula(paste(colnames(qsub)[pheno_col], "~",
                   paste(colnames(qsub)[c(covar,igrs_col)], collapse = "+"),
                   sep = "")),data=qsub,family="binomial")
  fh_mod<-glm(as.formula(paste(colnames(qsub)[pheno_col], "~",
                       paste(colnames(qsub)[c(covar,strat_col)], collapse = "+"),
                       sep = "")),data=qsub,family="binomial")
  add_mod<-glm(as.formula(paste(colnames(qsub)[pheno_col], "~",
                                paste(colnames(qsub)[c(covar,strat_col,igrs_col)], collapse = "+"),
                                sep = "")),data=qsub,family="binomial")
  #grs vs base
  formatC(anova(grs_mod,base_mod,test="LRT")$`Pr(>Chi)`[2],digits=3)
  PseudoR2(grs_mod,which="all")[["Nagelkerke"]]-PseudoR2(base_mod,which="all")[["Nagelkerke"]]
  #fh vs base
  formatC(anova(fh_mod,base_mod,test="LRT")$`Pr(>Chi)`[2],digits=3)
  PseudoR2(fh_mod,which="all")[["Nagelkerke"]]-PseudoR2(base_mod,which="all")[["Nagelkerke"]]
  #additive vs grs
  formatC(anova(add_mod,grs_mod,test="LRT")$`Pr(>Chi)`[2],digits=3)
  PseudoR2(add_mod,which="all")[["Nagelkerke"]]-PseudoR2(grs_mod,which="all")[["Nagelkerke"]]
  #additive vs fh
  formatC(anova(add_mod,fh_mod,test="LRT")$`Pr(>Chi)`[2],digits=3)
  PseudoR2(add_mod,which="all")[["Nagelkerke"]]-PseudoR2(fh_mod,which="all")[["Nagelkerke"]]
  #best (Full) vs additive
  formatC(anova(best,add_mod,test="LRT")$`Pr(>Chi)`[2],digits=3)
  PseudoR2(best,which="all")[["Nagelkerke"]]-PseudoR2(add_mod,which="all")[["Nagelkerke"]]
