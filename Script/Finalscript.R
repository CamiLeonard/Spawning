library(hms)
library(lubridate)
library(tidyverse)
library(car)
library(lme4)
library(glmmTMB)
library(rstatix)
library(ggpubr)
library(DataCombine)
library(multcomp)
library(plyr)
library(readxl)
library(writexl)
library(corrplot)
library(FactoMineR)
library(reshape2)
library(FSA)
library(Rmisc)
library(vegan)
library(tibble)

sem <- function(x) sd(x)/sqrt(length(x))

#timing
meta <- read.csv("Data/Metadata.csv",sep=";")
meta[meta == ''] <- NA
meta$FEC<-as.numeric(meta$FEC)
meta$TASS<-as_hms(meta$TASS)
meta$timing<-hour(meta$TASS)*60 + minute(meta$TASS)
meta$DAFM<-as.numeric(meta$DAFM)
meta<-data.frame(meta)
hist(meta$DAFM)
shapiro.test(meta$DAFM)
mod1<-glmer(DAFM~Species*Type+(1|Genotype),meta, family=nbinom2)
Anova(mod1)

meta %>%
  group_by(Species) %>%
  summarise(Days=mean(DAFM))

model<-meta %>% filter(Species=="A. cytherea")
model$Year<-as.factor(model$Year)
mod1<-glmer(DAFM~Type*Year+(1|Genotype),model, family=nbinom2)
Anova(mod1)
aggregate(DAFM~Type, model, mean)
kruskal.test(model$DAFM,model$Type)
shapiro.test(model$timing)
kruskal.test(model$timing,model$Type)

model<-meta %>% filter(Species=="A. pulchra")
shapiro.test(model$DAFM)
kruskal.test(model$DAFM,model$Type)
aggregate(DAFM~Type, model, mean)

model<-meta %>% filter(Species=="A. hyacinthus")
shapiro.test(model$DAFM)
kruskal.test(model$DAFM,model$Type)

model<-meta %>% filter(Species=="A. nasuta")
shapiro.test(model$DAFM)
kruskal.test(model$DAFM,model$Type)

model<-meta %>% filter(Species=="A. retusa")
shapiro.test(model$DAFM)
kruskal.test(model$DAFM,model$Type)

##Fig 2
total<-ddply(meta,c("Species","Type"),summarise,count=length(Type))
total$name<-paste(total$Species,total$Type,sep="_")
month<-ddply(meta,c("Species","Type","Month"),summarise,count=length(month))
daygraph<-ddply(meta,c("Species","Type","DAFM"),summarise,count=length(DAFM))
daygraph$name<-paste(daygraph$Species,daygraph$Type,sep="_")
daygraph$total<-daygraph$name
daygraph<-FindReplace(daygraph,"total",total,"name","count",exact=TRUE)
daygraph$percent<-daygraph$count*100/as.numeric(daygraph$total)
daygraph$DAFM<-as.character(daygraph$DAFM)
daygraph$DAFM<-factor(daygraph$DAFM,c("6","7","8","9","10","11","12","13"))
daygraph$Species<-factor(daygraph$Species,c("A. cytherea","A. pulchra", "A. hyacinthus" ,"A. nasuta","A. retusa")) 

dafm<-ggplot(daygraph, aes(x = DAFM, y = Type, fill = percent)) + 
  geom_tile()+
  scale_fill_gradientn(colours=c("#FFFFCC","#A1DAB4","#41B6C4", "#253494"), na.value = "grey98",limits = c(0, 100))+
  xlab("Days AFM")+
  guides(fill=guide_legend(title="% of colonies"))+  
  facet_grid(Species~., scales = "free", space = "free_y", switch = "y") +
  scale_x_discrete(position = "top")+
  scale_y_discrete(limits=rev)+
  theme_minimal()+
  theme(strip.text = element_text(face = "italic"),
        axis.title.y = element_blank(),
        panel.spacing = unit(0, "lines"),
        legend.position = c(0.85, 0.2))
dafm

meta$timing<-as.numeric(meta$timing)
time<-ddply(meta[meta$timing>0,],c("Species","Type"),summarise,time=mean(timing),error=sem(timing))
time<-left_join(total,time, by=c("Species","Type"))
time$time[5]<-240
time$time[6]<-240
time$time[7]<-180
time$time[8]<-180
time$Species<-factor(time$Species,c("A. cytherea","A. pulchra", "A. hyacinthus" ,"A. nasuta","A. retusa")) 
timing2<-ggplot(time, aes(x = time, y = Type)) + 
  geom_pointrange(aes(xmin=time-error, xmax=time+error))+
  xlab("Minutes after sunset")+ 
  facet_grid(Species~., scales = "free", space = "free_y", switch = "y") +
  scale_x_continuous(position = "top",limits=c(180,260))+
  scale_y_discrete(limits=rev)+
  theme_minimal()+
  theme(strip.text = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        panel.spacing = unit(0, "lines"))
timing2

graph<-ggarrange(dafm,timing2,widths = c(1,0.5))
ggsave("Output/timing.jpg",graph,width = 5.5, height = 4.5)
ggsave("Output/timing.svg",graph,width = 5.5, height = 4.5)

##split spawning
kruskal.test(meta$multiple_spawning~meta$Type)
meta$split<-0
meta$split[meta$multiple_spawning=="Y"]<-1
mod1<-glmer(split~Species*Type+(1|Genotype),meta,family=binomial)
Anova(mod1)
mod1<-glmer(split~Species+(1|Genotype),meta,family=binomial)
summary(glht(mod1, mcp(Species="Tukey")))
pairwise.wilcox.test(meta$split,meta$Species)

model<-meta %>% filter(Species=="A. cytherea")
model$Year<-as.factor(model$Year)
mod1<-glmer(DAFM~Type*Year+(1|Genotype),model, family=nbinom2)
Anova(mod1)
cyt<-ggplot(model,aes(x=Year,y=DAFM,color=Type))+
  geom_jitter(height=0.1,width = 0.1)+
  stat_compare_means(method="kruskal.test",label.y=12.5)+
  theme_bw()+
  scale_color_manual(values=c("#999722","#41B6C4"))
ggsave("Output/cyttiming.jpg",cyt,width = 8, height = 3)

##Fecundity
model<-meta %>% filter(FEC>=0)
mean(model$FEC)
sem(model$FEC)
hist(model$FEC)
shapiro.test(model$FEC)
mod1<-glmer(FEC~Species*Type+(1|Genotype),data=model,family=nbinom2)
plot(mod1)
Anova(mod1)
model$title<-"A. Fecundity"
figfec<-ggplot(model,aes(y=FEC,x=Species,fill=Type))+
  geom_boxplot()+
  theme_bw()+
  guides(x= guide_axis(angle=90)) +
  theme(axis.title.x = element_blank(),legend.position = "none",
        axis.text.x = element_text(lineheight=0.75,face="italic"))+  
  ylab(expression(paste("Fecundity (eggs/cm"^{2},")")))+
  scale_fill_manual(values=c("#FFFFCC","#41B6C4"))+
  facet_wrap(~title)
figfec

##Eggs per bundle
model<-meta %>% filter(EPB>=0)
shapiro.test(model$EPB)
mod1<-glmer(EPB~Species*Type+(1|Genotype),data=model,family=nbinom2)
plot(mod1)
Anova(mod1)
mean(model$EPB)
sem(model$EPB)
pairwise.wilcox.test(model$EPB,model$Species)
model<-meta %>% filter(Species=="A. cytherea" & EPB>=0)
kruskal.test(model$EPB,model$Type)
model<-meta %>% filter(Species=="A. hyacinthus" & EPB>=0)
kruskal.test(model$EPB,model$Type)
model<-meta %>% filter(Species=="A. nasuta" & EPB>=0)
kruskal.test(model$EPB,model$Type)
model<-meta %>% filter(Species=="A. pulchra" & EPB>=0)
kruskal.test(model$EPB,model$Type)
model<-meta %>% filter(Species=="A. retusa" & EPB>=0)
kruskal.test(model$EPB,model$Type)

model<-meta %>% filter(EPB>=0)
model$title<-"B. Eggs per bundle"
model$text[model$Species=="A. retusa"]<-"*"
figepb<-ggplot(model,aes(y=EPB,x=Species,fill=Type))+
  geom_boxplot()+
  theme_bw()+
  geom_signif(comparisons=list(c("A. retusa", "A. hyacinthus")),y_position = 15,annotation = "*")+
  geom_signif(comparisons=list(c("A. retusa", "A. nasuta")),y_position = 13,annotation = "*")+
  geom_signif(comparisons=list(c("A. pulchra","A. hyacinthus")),y_position = 11,annotation = "*")+
  guides(x= guide_axis(angle=90)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(lineheight=0.75,face="italic"),legend.position = "none")+
  ylab("Eggs per bundle")+
  scale_fill_manual(values=c("#FFFFCC","#41B6C4"))+
  geom_text(aes(label = text,y=12.5 ),size=4)+
  facet_wrap(~title)+
  ylim(4,16)
figepb

## Egg size

model<-meta %>% filter(ES>=0)
hist(model$ES,breaks=30)
shapiro.test(model$ES)
mod1<-lmer(ES~Species*Type+(1|Genotype),data=model)
Anova(mod1)
mean(model$ES)
sem(model$ES)
pairwise.wilcox.test(model$ES,model$Species)
ddply(model,c("Species"),summarize,Es=mean(ES),sem=sem(ES))
model$title<-"C. Eggs size"
figes<-ggplot(model,aes(y=ES,x=Species,fill=Type))+
  geom_boxplot()+
  theme_bw()+
  guides(x= guide_axis(angle=90)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(lineheight=0.75,face="italic"))+
  ylab("Eggs size (mm)")+
  geom_signif(comparisons=list(c("A. retusa", "A. hyacinthus")),y_position = 0.72,annotation = "*")+
  geom_signif(comparisons=list(c("A. pulchra", "A. hyacinthus")),y_position = 0.7,annotation = "*")+
  scale_fill_manual(values=c("#FFFFCC","#41B6C4"))+
  labs(fill = "Population")+
  facet_wrap(~title)+
  ylim(0.52,0.75)
figes

fig2<-ggarrange(figfec,figepb,figes,ncol=3,widths=c(1,1,1.4))
ggsave("Output/Fig2.pdf",fig2,width = 9, height = 3)
ggsave("Output/Fig2.jpg",fig2,width = 9, height = 3)

##Fertilization 
fert<- read_excel("Data/fullfertrate.xlsx")
fert2<-ddply(fert,c("Year","Species","Tube.Number","Population","Type","Temp","Male","Female","nParents"),
             summarise,NotFert = sum(NotFert), Fert.Eggs = sum(Fert.Eggs),Deformed = sum(Deformed),Total=sum(Total))
fert2$UF = fert2$NotFert/fert2$Total*100
fert2$NF = fert2$Fert.Eggs/fert2$Total*100
fert2$De = fert2$Deformed/fert2$Total*100
fert2$TF= 100-fert2$UF
fert2$nParents<-as.factor(fert2$nParents)
model<-fert2 %>% filter(Species!="A. pulchra")
shapiro.test(fert2$NF)
mod1<-glmer(NF~Temp*Species*Type+(1|nParents),data=model,family=nbinom2)
Anova(mod1)
plot(mod1)
anova<- as.data.frame(Anova(mod1)) 
anova1<-cbind("Explanatory"=rownames(anova),anova)
anova1<-anova1 %>% mutate(Explanatory = str_replace_all(Explanatory,":", " * "))
anova1$`Pr(>Chisq)`<-ifelse(anova1$`Pr(>Chisq)`<0.001,"<0.001",round(anova1$`Pr(>Chisq)`,digits=3))
anova1$Chisq<-round(anova1$Chisq,digits=2)
anova1 %>% write_xlsx("Fertanova2.xlsx")

model2<-model%>% filter(Species=="A. retusa" & Temp=="27°C")
kruskal.test(model2$NF,model2$Type)
model2<-model%>% filter(Species=="A. retusa" & Temp=="31°C")
kruskal.test(model2$NF,model2$Type)
model2<-model%>% filter(Temp=="27°C")
pairwise.wilcox.test(model2$NF,model2$Species)
model2<-model%>% filter(Temp=="31°C")
pairwise.wilcox.test(model2$NF,model2$Species)
ddply(model,c("Species"),summarize,Fert=mean(NF),sem=sem(NF))
model$text[model$Species=="A. retusa"]<-"*"

full<-ggplot(model,aes(x = Species, y = NF, fill=Type)) +
  geom_boxplot(position = position_dodge(preserve = "single"))+
  scale_fill_manual(values=c("#FFFFCC","#41B6C4")) +
  facet_grid(~Temp) +
  labs(y = "Normal fertilization %",fill = "Population") +
  theme_bw()+
  geom_text(aes(label = text,y=75 ),size=4)+
  guides(x= guide_axis(angle=90)) +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(lineheight=0.75,face="italic"))+
  ylim(0,110)
ggsave("Output/Fig4.pdf",full,width = 6, height = 3)
ggsave("Output/Fig4.svg",full,width = 6, height = 3)

##multivariate
meta2<-meta %>% filter(size>0)
hist(meta2$size)
shapiro.test(meta2$size)
mod1<-glmer(size~Type*Species+(1|Genotype),meta2,family=nbinom2)
Anova(mod1)
ggplot(meta2,aes(x=Species,y=size,fill=Type))+geom_boxplot()
ggplot(meta2,aes(x=Type,y=size,fill=Type))+geom_boxplot()+stat_compare_means(method="kruskal.test")
mod1<-glmer(EPB~size*Species+(1|Genotype),meta,family=nbinom2)
Anova(mod1)
mod1<-lmer(ES~size*Species+(1|Genotype),meta)
Anova(mod1)
mod1<-lmer(ES~EPB*Species+(1|Genotype),meta)
Anova(mod1)
mod1<-glmer(FEC~size*Species+(1|Genotype),meta,family=nbinom2)
Anova(mod1)
ggplot(meta,aes(x=size,y=FEC))+geom_point()+geom_smooth(method="lm")
cor.test(meta$FEC,meta$size)
meta2<-meta %>% filter(size<1000)
cor.test(meta2$FEC,meta2$size)
mod1<-glmer(FEC~size*Species+(1|Genotype),meta2,family=nbinom2)
Anova(mod1)
ggplot(meta2,aes(x=size,y=FEC))+geom_point()+geom_smooth(method="lm")

fert3<-ddply(fert2,c("Species","Type","Temp"),summarise,fert=mean(NF))
fert4<-spread(fert3[fert3$Species!="A. pulchra",],Temp,fert)
colnames(fert4)[3]<-"Fert27°C"
colnames(fert4)[4]<-"Fert31°C"
meta2<-meta[,c(3,4,5,7,13,17,19,21,22,24)]
meta2<-left_join(meta2,fert4, by=c("Species","Type"))
meta2$size[meta2$size==1000]<-NA
meta2$size<-as.numeric(meta2$size)
M=cor(meta2[,5:12],use="pairwise.complete.obs")
res1 <- cor.mtest(meta2[,5:12], conf.level = 0.95)
pvalues<-as.data.frame(res1$p)
corrplot(M, method = 'color',tl.col = 'black',p.mat = res1$p,sig.level = 0.05,insig = "blank",type='lower')
##export jpeg 470x410

spec<-c("#0A2F5B","#307BB5","#75C6BD","#EA9579","#7A0822")
pca1<-PCA(meta2[,5:12])
pca.vars <- pca1$var$coord %>% data.frame
pca.vars$vars <- rownames(pca.vars)
pca.vars.m <- melt(pca.vars, id.vars = "vars")
meta2$pc1 <- pca1$ind$coord[, 1]
meta2$pc2 <- pca1$ind$coord[, 2]
Pca<-ggplot(data = meta2, aes(x = pc1, y = pc2,color=Species,shape=Type)) +
  geom_point()+
  geom_segment(data = pca.vars, inherit.aes =F, aes(x = 0, xend = Dim.1*4, y = 0, 
                                                    yend = Dim.2*4),
               arrow = arrow(length = unit(0.025, "npc"), type = "open"), 
               lwd = 0.5)+
  stat_ellipse( aes(fill = Species), 
                show.legend = FALSE, 
                level = 0.95)+
  geom_text(data = pca.vars, inherit.aes =F,
            aes(x = Dim.1*4.5, y =  Dim.2*4.5,label=vars), 
            check_overlap = F, size = 3) +
  labs(x="PCA1 (21.34% variance)",y="PCA2 (21.08% variance)")+
  geom_hline(yintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
  geom_vline(xintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
  theme_minimal()+theme(legend.position="right")+
  scale_color_manual(values=spec,guide =
                       guide_legend(label.theme = element_text(angle = 0,size=10, face = "italic")))
Pca
ggsave("Output/PCAsp.jpg",Pca,width=7,height=4)
ggsave("Output/PCAsp.svg",Pca,width=7,height=4)

##selective breeding
all<-read_xlsx("Data/allindividualcrosses.xlsx")
all$Temp<-as.factor(all$Temp)
all$Male<-as.factor(all$Male)
all$Female<-as.factor(all$Female)
all$cross<-paste(all$Male,all$Female,sep="x")
hist(all$NF)
mod1<-glmer(NF~Temp*cross+(1|Year),data=all,family=nbinom2)
plot(mod1)
Anova(mod1)
aggregate(NF~Temp,all, mean)
aggregate(NF~Temp,all, sem)
all31<-all %>% filter(Temp=="31°C")
mod1<-glmer(NF~Female+Male+cross+(1|Year),data=all31,family=nbinom2)
Anova(mod1)
all27<-all %>% filter(Temp=="27°C")
mod1<-glmer(NF~Female+Male+cross+(1|Year),data=all27,family=nbinom2)
Anova(mod1)
all31mean<-ddply(all31,c("cross"),summarise,mean=mean(NF),sem=sem(NF))
dunn<-dunnTest(all31$NF, all31$cross,method="bonferroni")
dunn2<-data.frame(dunn$res)
dunn3<-dunn2 %>% filter(P.adj<0.05)
dunn3$P.unadj<-ifelse(dunn3$P.unadj<0.001,"<0.001",round(dunn3$P.unadj,digits=3))
dunn3$P.adj<-ifelse(dunn3$P.adj<0.001,"<0.001",round(dunn3$P.adj,digits=3))
dunn3$Z<-round(dunn3$Z,digits=2)
dunn3%>%write_xlsx("Output/bestcrosses.xlsx")
all31mean<-ddply(all31,c("Female"),summarise,mean=mean(NF),sem=sem(NF))
dunn<-dunnTest(all31$NF, all31$Female,method="bonferroni")
dunn2<-data.frame(dunn$res)
dunn3<-dunn2 %>% filter(P.adj<0.05)
dunn3$Comparison
all31mean<-ddply(all31,c("Male"),summarise,mean=mean(NF),sem=sem(NF))
dunn<-dunnTest(all31$NF, all31$Male,method="bonferroni")
dunn2<-data.frame(dunn$res)
dunn3<-dunn2 %>% filter(P.adj<0.05)
dunn3$Comparison

sum=summarySE(all, measurevar="NF", groupvars=c("Temp", "Male", "Female"))
ind<-ggplot(sum, aes(x = Female, y = Male, fill = NF)) + 
  geom_tile()+
  facet_wrap(~ Temp)+
  scale_fill_gradientn(colours=c("#FFFFCC","#A1DAB4","#41B6C4", "#253494"), na.value = "grey98",limits = c(0, 100))+
  ylab("Sperm from colony#")+
  xlab("Eggs from colony#")+
  guides(fill=guide_legend(title="Successful\nfertilization\nrate (%)"))+  theme_bw()
ind
ggsave("Output/Fig5.pdf",ind,width = 5.5, height = 2.5)
ggsave("Output/Fig5.svg",ind,width = 5.5, height = 2.5)
ggsave("Output/Fig5.jpg",ind,width = 5.5, height = 2.5)

##link with symbionts
symbionts<-read_xlsx("Data/postmedsummary2.xlsx")
symbionts<-symbionts %>% filter(genotype!=10)
symbionts$year<-"2021"
symbionts$year[symbionts$time=="T12"]<-"2022"
symbionts$genotype[symbionts$genotype==4]<-1
symbionts2<-gather(symbionts,"strain","percent",c(6:25))

symb<-ddply(symbionts2,c("year","genotype","strain"),summarize,percent=mean(percent))
symsum<-ddply(symb,c("year","genotype"),summarize,percent=sum(percent))
symbiont<-spread(symb,strain, percent)
symbiont21<-symbiont %>% filter(year=="2021")
symbiont21$CladeA<-rowSums(symbiont21[,c(3,4,20)])
symbiont21$CladeC<-rowSums(symbiont21[,c(5:16,21)])
symbiont21$CladeD<-rowSums(symbiont21[,c(17:19,22)])
colnames(symbiont21)[2]<-"Genotype"

##link between diversity and individual crosses
bcdist<-vegdist(symbiont[1:7,3:22], method="bray")
dist<-as.matrix(bcdist)
dist2<-as.data.frame(dist)
colnames(dist2)[4]<-5
colnames(dist2)[5]<-7
colnames(dist2)[6]<-8
colnames(dist2)[7]<-9

rownames(dist2)[7]<-9
rownames(dist2)[6]<-8
rownames(dist2)[5]<-7
rownames(dist2)[4]<-5

dist2 <- tibble::rownames_to_column(dist2, "gen")
dist3<-gather(dist2,gen2,dissimilarity,2:8)
dist3$cross<-paste(dist3$gen,dist3$gen2,sep="x")
allmean<-ddply(all,c("cross","Temp"),summarise,mean=mean(NF),sem=sem(NF))
allmean2<-left_join(allmean,dist3,by="cross")
allmean3<-allmean2 %>% filter(dissimilarity>=0)
Correlation<-ddply(allmean3,c("Temp"),summarise,
                   cor=cor(mean,dissimilarity),
                   pval=cor.test(mean,dissimilarity)$p.value)
allmean31<-allmean3 %>% filter(Temp=="31°C")
cor.test(allmean31$mean,allmean31$dissimilarity)
allmean27<-allmean3 %>% filter(Temp=="27°C")
cor.test(allmean27$mean,allmean27$dissimilarity)

fathersym<-all %>% filter(Male!=11)
fathersym$Genotype<-as.factor(fathersym$Male)
fathersym<-left_join(fathersym,symbiont21,by="Genotype")
fathersym$group<-"Sperm_27C"
fathersym$group[fathersym$Temp=="31°C"]<-"Sperm_31C"
mothersym<-all %>% filter(Female!=11)
mothersym$Genotype<-as.factor(mothersym$Female)
mothersym<-left_join(mothersym,symbiont21,by="Genotype")
mothersym$group<-"Eggs_27C"
mothersym$group[mothersym$Temp=="31°C"]<-"Eggs_31C"

fullfert<-rbind(fathersym[,c(5,8,30:33)],mothersym[,c(5,8,30:33)])
cormean<-ddply(fullfert,c("group"),summarise,
               CladeA=cor(NF,CladeA),
               CladeC=cor(NF,CladeC),
               CladeD=cor(NF,CladeD))

corsig<-ddply(fullfert,c("group"),summarise,
              pvalA=cor.test(NF,CladeA)$p.value,pvalC=cor.test(NF,CladeC)$p.value,
              pvalD=cor.test(NF,CladeD)$p.value)

corsig$CladeC[corsig$pvalC<0.05]<-"*"
corsig$CladeD[corsig$pvalD<0.05]<-"*"
cormean2<-gather(cormean,"Clade","cor",c(2:4))
corsig2<-gather(corsig,"Clade","pval",c(2:4))
corsig2$sig[corsig2$pval<0.05]<-"*"
cormean2$sig<-corsig2$sig

cor<-ggplot(cormean2, aes(x = Clade, y = group, fill = cor)) + 
  geom_tile()+
  scale_fill_gradientn(colours=c("#68001F","#E27C63","#FFFFFF","#67AACF", "#053061"), 
                       na.value = "grey98",limits = c(-1, 1))+
  theme_bw()+
  geom_text(aes(label=sig))+
  theme(axis.title=element_blank() )
cor
ggsave("Output/indcrossescor.jpg",cor,width=4,height=3)

eggs<-ddply(all,c("Temp","Female"),summarise,mean=mean(NF))
eggs1<-spread(eggs,Temp,mean)
colnames(eggs1)[1]<-"Genotype"
colnames(eggs1)[2]<-"Eggs_27C"
colnames(eggs1)[3]<-"Eggs_31C"
sperm<-ddply(all,c("Temp","Male"),summarise,mean=mean(NF))
sperm1<-spread(sperm,Temp,mean)
colnames(sperm1)[1]<-"Genotype"
colnames(sperm1)[2]<-"Sperm_27C"
colnames(sperm1)[3]<-"Sperm_31C"

summary<-left_join(eggs1,sperm1,by="Genotype")
summary<-left_join(summary,symbiont21[,c(2,23,24,25)],by="Genotype")
summary$Genotype<-paste("CYT",summary$Genotype,sep="")



##use all metadata
cyt<-meta %>% filter(Species=="A. cytherea" & Type=="Nursery")
cytmeta<-cyt[,c(4,5,13,17,19,21,22,24)]
cytmeta<-left_join(cytmeta,summary,by="Genotype")
M=cor(cytmeta[,3:15],use="pairwise.complete.obs")
res1 <- cor.mtest(cytmeta[,3:15], conf.level = 0.95)
pvalues<-as.data.frame(res1$p)
corrplot(M, method = 'color',tl.col = 'black',p.mat = res1$p,sig.level = 0.05,insig = "blank",type='lower')
##save svg 500x500
cytmeta2<-cytmeta %>% filter(Genotype!="CYT11" & Genotype!="CYT5")
pairwise.wilcox.test(cytmeta2$ES,cytmeta2$Genotype)
mod1<-glmer(EPB~size*CladeA*CladeC*CladeD+(1|Genotype),cytmeta,family=nbinom2)
Anova(mod1)
mod1<-lmer(ES~size*CladeA*CladeC*CladeD+(1|Genotype),cytmeta)
Anova(mod1)
mod1<-glmer(FEC~size*CladeA*CladeC*CladeD+(1|Genotype),cytmeta,family=nbinom2)
Anova(mod1)
ggplot(cytmeta,aes(x=size,y=ES))+geom_point()+geom_smooth(method="lm")
cor.test(cytmeta$FEC,cytmeta$size)

##PCA
geno<-c("#0A2F5B","#307BB5","#BFDBEB","#47AB5A","#D4D42C","#EA9579","#BE3137","#860B25")
pca1<-PCA(cytmeta[,3:15])
pca.vars <- pca1$var$coord %>% data.frame
pca.vars$vars <- rownames(pca.vars)
pca.vars.m <- melt(pca.vars, id.vars = "vars")
cytmeta$pc1 <- pca1$ind$coord[, 1]
cytmeta$pc2 <- pca1$ind$coord[, 2]

Pca<-ggplot(data = cytmeta, aes(x = pc1, y = pc2,color=Genotype)) +
  geom_point()+
  geom_segment(data = pca.vars, inherit.aes =F, aes(x = 0, xend = Dim.1*4.5, y = 0, 
                                                    yend = Dim.2*4.5),
               arrow = arrow(length = unit(0.025, "npc"), type = "open"), 
               lwd = 0.5)+
  stat_ellipse( aes(fill = Genotype), 
                show.legend = FALSE, 
                level = 0.95)+
  geom_text(data = pca.vars, inherit.aes =F,
            aes(x = Dim.1*5.5, y =  Dim.2*5,label=vars), 
            check_overlap = F, size = 3) +
  labs(x="PCA1 (25.94% variance)",y="PCA2 (19.57% variance)")+
  geom_hline(yintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
  geom_vline(xintercept = 0, lty = 2, color = "grey", alpha = 0.9) +
  theme_minimal()+theme(legend.position="right")+
  scale_color_manual(values=geno)
Pca
ggsave("Output/PCA.jpg",Pca,width=6,height=4)

##Larval seeding
rec<-read_xlsx("Data/RECRUITS.xlsx")
rec$den2<-as.numeric(rec$'12 months')*100/rec$hemisphere_surface_cm2
rec$net<-"Seeding"
rec$net[rec$TRMT=="CTRL"]<-"Control"
rec$title<-"B. Recruit density after 1 year"
plot1<-ggplot(rec,aes(x=net,y=den2,fill=net))+geom_boxplot()+
  ylab(expression(paste("Recruits/100 cm"^{2})))+ theme_bw()+
  theme(
    legend.position="none",axis.title.x = element_blank())+
  xlab("")+
  scale_fill_manual(values=c("#FFFFCC","#FFFFCC"))+
  facet_grid(cols=vars(title))+
  geom_signif(comparisons = list(c("Control", "Seeding")),annotation = "*")+
  ylim(0,0.33)
plot1

mean(rec$den2[rec$TRMT=="NET"])
sem(rec$den2[rec$TRMT=="NET"])
hist(rec$den2)
kruskal.test(rec$den2,rec$net)

tiles<-read_xlsx("Data/FILETS.xlsx")
tiles$traitement<-gsub('FILET','NET',tiles$traitement)

##if we assume plugs are 2cm in diameter: density
tiles$rec<-tiles$recrues*100/(tiles$`nbr plots`*3.14)
tiles$title<-"A. Recruit density after 10 days"
tiles$net<-"Seeding"
tiles$net[tiles$traitement=="CONTROL"]<-"Control"

##whats the mean initial recruit density on each bommie
density<-ddply(tiles,c("patate", "net"),summarise,dens=mean(rec))
mean(density$dens[density$dens!=0])

plot<-ggplot(tiles,aes(x=net,y=rec,fill=net))+geom_boxplot()+
  ylab(expression(paste("Recruits/100 cm"^{2})))+theme_bw()+
  theme(
    legend.position="none",axis.title.x = element_blank())+
  xlab("")+
  scale_fill_manual(values=c("#FFFFCC","#FFFFCC"))+
  facet_grid(cols=vars(title))+
  geom_signif(comparisons = list(c("Control", "Seeding")),annotation = "*")+
  ylim(0,55)
plot

mean(tiles$rec[tiles$net=="Net"])
aggregate(rec~net,tiles,mean)
aggregate(rec~net,tiles,sem)
kruskal.test(tiles$rec~tiles$net)

size<-read_xlsx("Data/SIZE.xlsx")
aggregate(SIZE~AGE,size,mean)
aggregate(SIZE~AGE,size,sem)

fig<-ggarrange(plot,plot1)
fig
ggsave("Output/Fig6.pdf",fig,width = 5, height = 2.5)
ggsave("Output/Fig6.svg",fig,width = 5, height = 2.5)
ggsave("Output/Fig6.jpg",fig,width = 5, height = 2.5)

