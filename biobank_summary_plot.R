library(googlesheets4)
library(googledrive)
library(readxl)
library(rcartocolor)
library(dplyr)
library(ggplot2)
library(stringi)
library(tidyverse)

#I think this doesn't work because it's a .xlsx
#file<-"https://docs.google.com/spreadsheets/d/1DNKd1KzI8WOIfG2klXWskbCSyX6h5gTu/edit?usp=sharing&ouid=115040676596308698055&rtpof=true&sd=true"

#download file to use
file<-drive_find(pattern="Intervene_flagship_endpoint_collection.xlsx") #takes a bit of time..would be better to know path
drive_download(file$name,overwrite=TRUE)

#read in names
names<-read_xlsx(file$name,sheet=1)
names2<-names[c("Endpoint","Endpoint family","PhenoID")]
my_colors = carto_pal(length(unique(names2$`Endpoint family`)), "Safe")
names2$Group<-names2$`Endpoint family`
names2$`Endpoint family`<-as.factor(names2$`Endpoint family`)
levels(names2$`Endpoint family`)<-my_colors #set colorblind frienldy colors

#list sheets
sheets<-excel_sheets(file$name)
indices<-grep("Metrics",sheets) #select sheets with metrics

drop_rows_all_na <- function(x, pct=0.90) x[!rowSums(is.na(x)) >= ncol(x)*pct,]

### function to read and parse each metric sheet into consistent formatting
my_read_func<-function(sheet){
  df<-read_xlsx(file$name,sheet=sheet,na=c("NA","na"))
  df<-drop_rows_all_na(df) #handle finngen formatting with merged columns
  df<-df[c(1:39),c(1:10)] #cut extra columns that some sheets have, manually checked that base 10 are the same (moved extra to the end)
  names(df)<-c("endpoint","def","cases","controls","age_baseline","follow_dist","age_onset_dist","age_corr","sex_corr","female")

  print(sheet)

  df$cases<-as.numeric(df$cases)
  df$controls<-as.numeric(df$controls)

  #format the IQR columns
  df$age_baseline_median<-trimws(sapply(stri_split_regex(as.character(df$age_baseline),pattern="\\s|\n|\\(",n=2),"[",1))
  df$age_baseline_iqr<- trimws(sapply(stri_split_regex(as.character(df$age_baseline),pattern="\\s|\n|\\(",n=2),"[",2))

  df$age_onset_dist_median<-trimws(sapply(stri_split_regex(as.character(df$age_onset_dist),pattern="\\s|\n|\\(",n=2),"[",1))
  df$age_onset_dist_iqr<-trimws(sapply(stri_split_regex(as.character(df$age_onset_dist),pattern="\\s|\n|\\(",n=2),"[",2))

  df$follow_dist_median<-trimws(sapply(stri_split_regex(as.character(df$follow_dist),pattern="\\s|\n|\\(",n=2),"[",1))
  df$follow_dist_iqr<-trimws(sapply(stri_split_regex(as.character(df$follow_dist),pattern="\\s|\n|\\(",n=2),"[",2))

  df$sex_corr<-trimws(sapply(stri_split_regex(as.character(df$sex_corr),pattern="\\s|\n|\\(",n=2),"[",1))
  df$age_corr<-trimws(sapply(stri_split_regex(as.character(df$age_corr),pattern="\\s|\n|\\(",n=2),"[",1)) #pearson?

  df$female_est<-as.numeric(trimws(sapply(stri_split_regex(as.character(df$female),pattern="\\s|\n|\\(",n=2),"[",1)))
  df$female_ci<-trimws(sapply(stri_split_regex(as.character(df$female),pattern="\\s|\n|\\d\\(",n=2),"[",2))

  df$biobank<-unlist(strsplit(sheets[sheet],"-"))[[2]] #label biobank

  df<-df %>% select("endpoint","cases","controls","age_baseline_median","age_baseline_iqr","age_onset_dist_median","age_onset_dist_iqr","follow_dist_median","follow_dist_iqr","age_corr","sex_corr","female_est","female_ci","biobank")
  return(df)
}
#loop over sheets
list_of_biobanks<-lapply(indices,my_read_func)

dat<-bind_rows(list_of_biobanks) %>% left_join(names2,by=c("endpoint"="Endpoint"))
dat<-dat[!is.na(dat$endpoint),] #remove when endpoint is NA
dat<-dat[!is.na(dat$Group),]
dat<-dat[dat$endpoint!="BMI (obesity)",]#remove BMI since it is quant
dat$endpoint<-as.factor(dat$endpoint)
dat$Group<-as.factor(dat$Group)

######### number of cases
options(scipen=10)
pdf(file="cases.pdf",height=10,width=12)
ggplot(dat,aes(y=endpoint,x=cases)) + geom_bar(aes(fill = `Endpoint family`),stat="identity") + facet_wrap(~biobank,nrow=2) +
  scale_fill_manual(values = levels(dat$`Endpoint family`), labels = levels(dat$Group)) + theme_bw() +
  theme(legend.position="bottom") + scale_y_discrete(limits=rev(levels(dat$endpoint)))
dev.off()

## median follow up
dat$follow_dist_low<-as.numeric(gsub("[^0-9.-]", "",sapply(stri_split_regex(as.character(dat$follow_dist_iqr),pattern="-",n=2),"[",1)))
dat$follow_dist_high<-as.numeric(gsub("[^0-9.-]", "",sapply(stri_split_regex(as.character(dat$follow_dist_iqr),pattern="-",n=2),"[",2)))
pdf(file="followup.pdf",height=10,width=12)
ggplot(dat[!is.na(dat$follow_dist_median),],aes(y=endpoint,x=as.numeric(follow_dist_median))) + geom_point(aes(color = `Endpoint family`),stat="identity") + facet_wrap(~biobank,nrow=2) +
  scale_color_manual(values = levels(dat$`Endpoint family`), labels = levels(dat$Group)) +
  theme_bw() + labs(x="Median and interquartile range of follow-up from baseline (yrs)") + geom_linerange(aes(xmin=follow_dist_low,xmax=follow_dist_high,color=`Endpoint family`)) +
  theme(legend.position="bottom")  + scale_y_discrete(limits=rev(levels(dat$endpoint)))
dev.off()

##median age of onset
dat$age_onset_dist_low<-as.numeric(gsub("[^0-9.-]", "",sapply(stri_split_regex(as.character(dat$age_onset_dist_iqr),pattern="-",n=2),"[",1)))
dat$age_onset_dist_high<-as.numeric(gsub("[^0-9.-]", "",sapply(stri_split_regex(as.character(dat$age_onset_dist_iqr),pattern="-",n=2),"[",2)))
pdf(file="onset.pdf",height=10,width=12)
ggplot(dat[!is.na(dat$age_onset_dist_median),],aes(y=endpoint,x=as.numeric(age_onset_dist_median))) + geom_point(aes(color = `Endpoint family`),stat="identity") + facet_wrap(~biobank,nrow=2) +
  scale_color_manual(values = levels(dat$`Endpoint family`), labels = levels(dat$Group)) +   theme_bw() + labs(x="Median age of onset and interquartile range") +
  geom_linerange(aes(xmin=age_onset_dist_low,xmax=age_onset_dist_high,color=`Endpoint family`)) +
theme(legend.position="bottom") + scale_y_discrete(limits=rev(levels(dat$endpoint)))
dev.off()

##### % female
dat$female_est_perc<-as.numeric(ifelse(dat$female_est<=1,dat$female_est*100,dat$female_est))
pdf(file="female.pdf",height=10,width=12)
ggplot(dat[!is.na(dat$female_est_perc),],aes(y=endpoint,x=female_est_perc)) + geom_point(aes(color = `Endpoint family`),stat="identity") + facet_wrap(~biobank,nrow=2) +
  scale_color_manual(values = levels(dat$`Endpoint family`), labels = levels(dat$Group)) + theme_bw() + labs(x="% Female") + geom_vline(xintercept=50,linetype="dashed",color="red") +
  theme(legend.position="bottom") +  scale_y_discrete(limits=rev(levels(dat$endpoint)))
dev.off()

###correlations
dat_long<-melt(dat[,c("endpoint","age_corr","sex_corr","biobank","Endpoint family","Group")],id.vars=c("endpoint","biobank","Endpoint family","Group"))
dat_long$value<-as.numeric(dat_long$value)
pdf(file="corr.pdf",height=10,width=12)
ggplot(dat_long[!is.na(dat_long$value),],aes(y=endpoint,x=value,shape=variable)) + geom_point(aes(color = `Endpoint family`),stat="identity") + facet_wrap(~biobank,nrow=2) + geom_vline(linetype="dashed",color="red",xintercept=0) +
  scale_color_manual(values = levels(dat_long$`Endpoint family`), labels = levels(dat_long$Group)) + theme_bw()  + theme(legend.position="bottom") +  scale_y_discrete(limits=rev(levels(dat_long$endpoint)))
dev.off()


##### total cases across biobanks
tot<-read_xlsx(file$name,sheet=sheets[sheets=="Case-Totals"],na=c("NA","na"))
tot<-tot[tot$Endpoints!="BMI (obesity)",]
tot<-tot %>% left_join(names2,by=c("Endpoints"="Endpoint"))
tot$Endpoints<-as.factor(tot$Endpoints)
tot$Group<-as.factor(tot$Group)
pdf(file="all_biobanks_cases.pdf",height=6,width=10)
ggplot(tot,aes(y=Endpoints,x=`N cases`)) + geom_bar(aes(fill = `Endpoint family`),stat="identity") +
  scale_fill_manual(values = levels(tot$`Endpoint family`), labels = levels(tot$Group)) + theme_bw() +
  theme(legend.position="bottom") + scale_y_discrete(limits=rev(levels(tot$Endpoints)))
dev.off()


#would be good to plot the disease prevalence across biobank one vs the others. Like scatter plots.
#and also the correlation with age and sex
#to understand the consistency across endpoints
