###format tables
library(data.table)
library(dplyr)

##replication
df<-fread(file="GBMI_lookup.txt")
names(df)[13]<-"pval"
df<-df[!(pval==".")]
df$Effect<-as.numeric(df$Effect)
df$StdErr<-as.numeric(df$StdErr)

##GBMI
gb<-fread("/net/hunt/disk2/wukh/GBMI_bio/VTE/VTE_locus_cs_OldLocusFile_072021.txt")

df2<-left_join(gb,df,by=c("POS"="pos_hg38"))
df2$chr_hg38<-NULL
df2$p<-ifelse(as.numeric(df2$p)==0,2.22E-308,df2$p)

df3<-df2[,c("CHR","POS","REF","ALT","all_inv_var_meta_beta","all_inv_var_meta_sebeta","Allele1","Allele2","Effect","StdErr","pval","all_inv_var_meta_p")]
df3$rep_UB<-df3$Effect+1.96*df3$StdErr
df3$rep_LB<-df3$Effect-1.96*df3$StdErr
df3$gbmi_UB<-df3$all_inv_var_meta_beta + 1.96*df3$all_inv_var_meta_sebeta
df3$gbmi_LB<-df3$all_inv_var_meta_beta - 1.96*df3$all_inv_var_meta_sebeta

cor<-cor.test(df3$Effect,df3$all_inv_var_meta_beta) #0.846

pdf(file="replication_scatter.pdf",height=3,width=3)
ggplot(df3,aes(x=Effect,y=all_inv_var_meta_beta)) + geom_point() + theme_bw() +
  labs(x="INVENT+MVP Replication",y="GBMI") +
  geom_errorbarh(aes(xmin=rep_LB,xmax=rep_UB)) + geom_errorbar(aes(ymin=gbmi_LB,ymax=gbmi_UB))
dev.off()

##known and novel
main<-fread("/net/hunt/disk2/wukh/GBMI_bio/VTE/GBMI_VTE_IndexVariants_KnownTrait.txt",na.strings=".",fill=TRUE)

main$beta_pos=abs(main$beta)
main<-main %>% mutate(freq_pos=case_when(beta < 0~(1-freq), beta>0~freq),
                      effect_allele=case_when(beta<0~ref,beta>0~alt))

main$OR<-exp(main$beta_pos)
main$LB<-exp(main$beta_pos-(1.96*main$se))
main$UB<-exp(main$beta_pos+(1.96*main$se))
main$CI<-paste0("[",formatC(main$LB,digits=2),",",formatC(main$UB,digits=2),"]")
main<-main[,c(1,2,3,4,5,14,16,20,21,22,23,24,25,26)]
main$OR<-formatC(main$OR,digits=2)
main$p<-formatC(main$p,digits=3)
main$VTE<-ifelse(main$VTE==0,"Potentially Novel","Known")
sub<-main[,c("chr","pos","VTE")]

#merge
df4<-left_join(df3,sub,by=c("POS"="pos"))

#direction of effect
table(df4$Effect/df4$all_inv_var_meta_beta < 0)

pdf(file="replication_scatter_known_novel.pdf",height=4,width=4)
ggplot(df4,aes(x=Effect,y=all_inv_var_meta_beta,color=VTE)) +geom_abline(intercept=0,slope=1,linetype="dashed",color="grey",alpha=0.5) +
geom_point() + theme_bw() +
  labs(x="INVENT+MVP Replication Effect Size",y="GBMI Effect Size") + scale_color_manual(values=c("blue","red")) +
  geom_errorbarh(aes(xmin=rep_LB,xmax=rep_UB)) + geom_errorbar(aes(ymin=gbmi_LB,ymax=gbmi_UB)) + geom_vline(xintercept=0,color="black") + geom_hline(yintercept=0,color="black") +
  theme(legend.position="bottom") + geom_text(x = -0.35, y = 1, label = paste("R^2==",format(cor$estimate,digits=2)),parse=T,color="black")
dev.off()


#whiskers
#geom_errorbarh(aes(xmin=rep_LB,xmax=rep_UB),height=0.05) + geom_errorbar(aes(ymin=gbmi_LB,ymax=gbmi_UB),width=0.05) +
