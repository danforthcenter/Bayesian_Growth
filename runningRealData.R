library(brms)
library(ggplot2)
library(stringr)
library(patchwork)
library(cowplot)
library(plyr)
library(tidyverse)
library(bayesplot)
library(viridis)
library(viridisLite)

# get functions
source('sourcingFunctions.R')

# Read in Data

area_convert <- 13.2*3.7/46856


# here SynComA is going to be relabelled a and control is going to be relabelled b 
#sv_shapes <- read.table("b51632508ce74bdd9ceffcd073b871e6-shapes_outliersRM.tsv",sep="\t",header=T,stringsAsFactors = F)

sv_shapes <- read.table("data.tsv",sep="\t",header=T,stringsAsFactors = F)
sv_shapes <- sv_shapes[sv_shapes$Drought == "ABB" & sv_shapes$Microbes %in% c("SynCom B","SynCom A"),] #add "SynCom B"
sv_shapes$area_c <- sv_shapes$area * area_convert

df <- aggregate(data=sv_shapes[sv_shapes$DAP > 4,], area_c ~ Microbes+Barcodes+DAP, FUN=function(i)mean(i))
df <- df[(df$Barcodes %in% df$Barcodes[df$DAP == 25]) & (df$Barcodes %in% df$Barcodes[df$DAP == 5]),]
df$Microbes <- gsub(" ","",df$Microbes)

df<-df%>%rename(treatment=Microbes, y=area_c, sample=Barcodes)%>%
  mutate(treatment=ifelse(treatment=="SynComA", "a", "b"))
head(df)

## Version with separate brm() calls and df's


df_5_9<-df%>%filter(DAP<=9)
df_5_11<-df%>%filter(DAP<=11)
df_5_13<-df%>%filter(DAP<=13)
df_5_15<-df%>%filter(DAP<=15)
df_5_17<-df%>%filter(DAP<=17)
df_5_19<-df%>%filter(DAP<=19)
df_5_21<-df%>%filter(DAP<=21)
df_5_23<-df%>%filter(DAP<=23)
df_5_25<-df%>%filter(DAP<=25)

prior1 <-
  prior(lognormal(log(125), .25),nlpar = "phi1",coef="treatmentb") + #MicrobesControl
  prior(lognormal(log(125), .25),nlpar = "phi1",coef="treatmenta") + #MicrobesSynComA
  prior(lognormal(log(15), .25), nlpar = "phi2",coef="treatmentb") +
  prior(lognormal(log(15), .25), nlpar = "phi2",coef="treatmenta") +
  prior(lognormal(log(5), .25), nlpar = "phi3",coef="treatmentb") +
  prior(lognormal(log(5), .25), nlpar = "phi3",coef="treatmenta") +
  prior(student_t(3,0,5), dpar="sigma") +
  prior(gamma(2,0.1), class="nu")

fit_5_to_9 <- brm(bf(y ~ phi1/(1+exp((phi2-DAP)/phi3)),
                      sigma~s(DAP,by=treatment, k=4),
                      phi1 + phi2 + phi3 ~ 0+treatment,
                      autocor = ~arma(~DAP|sample:treatment,1,1),nl = TRUE),
                   family = student, prior = prior1, data = df_5_9, iter = 1000,
                   cores = 2, chains = 2, backend = "cmdstanr", #threads = threading(4),
                   control = list(adapt_delta = 0.999,max_treedepth = 20),
                   inits = function(){list(b_phi1=rgamma(3,1),b_phi2=rgamma(3,1),b_phi3=rgamma(3,1))})

fit_5_to_11 <- brm(bf(y ~ phi1/(1+exp((phi2-DAP)/phi3)),
                      sigma~s(DAP,by=treatment, k=6),
                      phi1 + phi2 + phi3 ~ 0+treatment,
                      autocor = ~arma(~DAP|sample:treatment,1,1),nl = TRUE),
                   family = student, prior = prior1, data = df_5_11, iter = 1000,
                   cores = 2, chains = 2, backend = "cmdstanr", #threads = threading(4),
                   control = list(adapt_delta = 0.999,max_treedepth = 20),
                   inits = function(){list(b_phi1=rgamma(3,1),b_phi2=rgamma(3,1),b_phi3=rgamma(3,1))})

fit_5_to_13 <- brm(bf(y ~ phi1/(1+exp((phi2-DAP)/phi3)),
                      sigma~s(DAP,by=treatment, k=8),
                      phi1 + phi2 + phi3 ~ 0+treatment,
                      autocor = ~arma(~DAP|sample:treatment,1,1),nl = TRUE),
                   family = student, prior = prior1, data = df_5_13, iter = 1000,
                   cores = 2, chains = 2, backend = "cmdstanr", #threads = threading(4),
                   control = list(adapt_delta = 0.999,max_treedepth = 20),
                   inits = function(){list(b_phi1=rgamma(3,1),b_phi2=rgamma(3,1),b_phi3=rgamma(3,1))})

fit_5_to_15 <- brm(bf(y ~ phi1/(1+exp((phi2-DAP)/phi3)),
                      sigma~s(DAP,by=treatment),
                      phi1 + phi2 + phi3 ~ 0+treatment,
                      autocor = ~arma(~DAP|sample:treatment,1,1),nl = TRUE),
                   family = student, prior = prior1, data = df_5_15, iter = 1000,
                   cores = 2, chains = 2, backend = "cmdstanr", #threads = threading(4),
                   control = list(adapt_delta = 0.999,max_treedepth = 20),
                   inits = function(){list(b_phi1=rgamma(3,1),b_phi2=rgamma(3,1),b_phi3=rgamma(3,1))})

fit_5_to_17 <- brm(bf(y ~ phi1/(1+exp((phi2-DAP)/phi3)),
                      sigma~s(DAP,by=treatment),
                      phi1 + phi2 + phi3 ~ 0+treatment,
                      autocor = ~arma(~DAP|sample:treatment,1,1),nl = TRUE),
                   family = student, prior = prior1, data = df_5_17, iter = 1000,
                   cores = 2, chains = 2, backend = "cmdstanr", #threads = threading(4),
                   control = list(adapt_delta = 0.999,max_treedepth = 20),
                   inits = function(){list(b_phi1=rgamma(3,1),b_phi2=rgamma(3,1),b_phi3=rgamma(3,1))})

fit_5_to_19 <- brm(bf(y ~ phi1/(1+exp((phi2-DAP)/phi3)),
                      sigma~s(DAP,by=treatment),
                      phi1 + phi2 + phi3 ~ 0+treatment,
                      autocor = ~arma(~DAP|sample:treatment,1,1),nl = TRUE),
                   family = student, prior = prior1, data = df_5_19, iter = 1000,
                   cores = 2, chains = 2, backend = "cmdstanr", #threads = threading(4),
                   control = list(adapt_delta = 0.999,max_treedepth = 20),
                   inits = function(){list(b_phi1=rgamma(3,1),b_phi2=rgamma(3,1),b_phi3=rgamma(3,1))})

fit_5_to_21 <- brm(bf(y ~ phi1/(1+exp((phi2-DAP)/phi3)),
                      sigma~s(DAP,by=treatment),
                      phi1 + phi2 + phi3 ~ 0+treatment,
                      autocor = ~arma(~DAP|sample:treatment,1,1),nl = TRUE),
                   family = student, prior = prior1, data = df_5_21, iter = 1000,
                   cores = 2, chains = 2, backend = "cmdstanr", #threads = threading(4),
                   control = list(adapt_delta = 0.999,max_treedepth = 20),
                   inits = function(){list(b_phi1=rgamma(3,1),b_phi2=rgamma(3,1),b_phi3=rgamma(3,1))})

fit_5_to_23 <- brm(bf(y ~ phi1/(1+exp((phi2-DAP)/phi3)),
                      sigma~s(DAP,by=treatment),
                      phi1 + phi2 + phi3 ~ 0+treatment,
                      autocor = ~arma(~DAP|sample:treatment,1,1),nl = TRUE),
                   family = student, prior = prior1, data = df_5_23, iter = 1000,
                   cores = 2, chains = 2, backend = "cmdstanr", #threads = threading(4),
                   control = list(adapt_delta = 0.999,max_treedepth = 20),
                   inits = function(){list(b_phi1=rgamma(3,1),b_phi2=rgamma(3,1),b_phi3=rgamma(3,1))})

fit_5_to_25 <- brm(bf(y ~ phi1/(1+exp((phi2-DAP)/phi3)),
                      sigma~s(DAP,by=treatment),
                      phi1 + phi2 + phi3 ~ 0+treatment,
                      autocor = ~arma(~DAP|sample:treatment,1,1),nl = TRUE),
                   family = student, prior = prior1, data = df_5_25, iter = 2000,
                   cores = 2, chains = 2, backend = "cmdstanr", #threads = threading(4),
                   control = list(adapt_delta = 0.999,max_treedepth = 20),
                   inits = function(){list(b_phi1=rgamma(3,1),b_phi2=rgamma(3,1),b_phi3=rgamma(3,1))})

save(df,prior1, fit_5_to_9,fit_5_to_11,fit_5_to_13,fit_5_to_15,fit_5_to_17,
     fit_5_to_19,fit_5_to_21,fit_5_to_23, fit_5_to_25,
     file = "realData_scAvscB_fit_models_and_data_every_2_days_drought.rdata")

2+2


#*********************************************
# Gompertz Early stopping
#*********************************************
# idea here is that the real data is at least as well represented by gompertz as it is by logistic
# do our early stopping results change if we reparameterize on a gompertz model?

#condor_submit -i request_cpus=2 'requirements=(TARGET.Machine=!="cylon.datasci.danforthcenter.org" && TARGET.Machine=!="pallas.datasci.danforthcenter.org")'

library(brms)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(bayesplot)
library(viridis)

area_convert <- 13.2*3.7/46856
sv_shapes <- read.table("../data.tsv",sep="\t",header=T,stringsAsFactors = F)
sv_shapes <- sv_shapes[sv_shapes$Drought == "ABB" & sv_shapes$Microbes %in% c("SynCom B","SynCom A"),] #add "SynCom B"
sv_shapes$area_c <- sv_shapes$area * area_convert
df <- aggregate(data=sv_shapes[sv_shapes$DAP > 4,], area_c ~ Microbes+Barcodes+DAP, FUN=function(i)mean(i))
df <- df[(df$Barcodes %in% df$Barcodes[df$DAP == 25]) & (df$Barcodes %in% df$Barcodes[df$DAP == 5]),]
df$Microbes <- gsub(" ","",df$Microbes)
df<-df%>%rename(treatment=Microbes, y=area_c, sample=Barcodes, time=DAP)%>%
  mutate(treatment=ifelse(treatment=="SynComA", "a", "b"))

priorGomp <- prior(lognormal(log(130), .25),nlpar = "a") +
  prior(lognormal(log(12), .25), nlpar = "b") + 
  prior(lognormal(log(1.2), .25), nlpar = "c") + 
  prior(student_t(3,0,5), dpar="sigma") +
  prior(gamma(2,0.1), class="nu")

# for (i in seq(11,13,by=2)){
#   dataText<-paste0("df_",i,"<-df%>%filter(time>=5 & time<=",i,")")
#   eval(parse(text=dataText))
#   if(i<=13){
#     sigmaText<-paste0('sigma~s(time, by=treatment, k=',i-5,'),')
#   } else{
#     sigmaText<-paste0('sigma~s(time, by=treatment),')}
#   brmText<-paste0('fit_',i,' <- brm(bf(y ~ a*exp(-b*exp(-c*time)),',
#                   sigmaText,
#                   'a + b + c ~ 0+treatment, autocor = ~arma(~time|sample:treatment,1,1),nl = TRUE),family = student, prior = priorGomp, data = df_',i,', iter = 1000, cores = 2, chains = 2, backend = "cmdstanr",control = list(adapt_delta = 0.999,max_treedepth = 20),inits = function(){list(b_a=rgamma(2,1),b_b=rgamma(2,1),b_c=rgamma(2,1))})')
#   print(i)
#   eval(parse(text=brmText))
#}
# save(df_7, df_9, df_11, df_13, fit_7, fit_9, fit_11, fit_13, file="gompertzModel_earlyStopping_7_to_13.rdata")



# for (i in seq(15,19,by=2)){
#   dataText<-paste0("df_",i,"<-df%>%filter(time>=5 & time<=",i,")")
#   eval(parse(text=dataText))
#   if(i<13){
#     sigmaText<-paste0('sigma~s(time, by=treatment, k=',i-5,'),')
#   } else{
#     sigmaText<-paste0('sigma~s(time, by=treatment),')}
#   brmText<-paste0('fit_',i,' <- brm(bf(y ~ a*exp(-b*exp(-c*time)),',
#                   sigmaText,
#                   'a + b + c ~ 0+treatment, autocor = ~arma(~time|sample:treatment,1,1),nl = TRUE),family = student, prior = priorGomp, data = df_',i,', iter = 1000, cores = 2, chains = 2, backend = "cmdstanr",control = list(adapt_delta = 0.999,max_treedepth = 20),inits = function(){list(b_a=rgamma(2,1),b_b=rgamma(2,1),b_c=rgamma(2,1))})')
#   print(i)
#   eval(parse(text=brmText))
# }
# save(df_15, df_17, df_19, fit_15, fit_17, fit_19, file="gompertzModel_earlyStopping_15_to_19.rdata")


# for (i in seq(21,25,by=2)){
#   dataText<-paste0("df_",i,"<-df%>%filter(time>=5 & time<=",i,")")
#   eval(parse(text=dataText))
#   if(i<13){
#     sigmaText<-paste0('sigma~s(time, by=treatment, k=',i-5,'),')
#   } else{
#     sigmaText<-paste0('sigma~s(time, by=treatment),')}
#   brmText<-paste0('fit_',i,' <- brm(bf(y ~ a*exp(-b*exp(-c*time)),',
#                   sigmaText,
#                   'a + b + c ~ 0+treatment, autocor = ~arma(~time|sample:treatment,1,1),nl = TRUE),family = student, prior = priorGomp, data = df_',i,', iter = 1000, cores = 2, chains = 2, backend = "cmdstanr",control = list(adapt_delta = 0.999,max_treedepth = 20),inits = function(){list(b_a=rgamma(2,1),b_b=rgamma(2,1),b_c=rgamma(2,1))})')
#   print(i)
#   eval(parse(text=brmText))
# }
# save(df_21, df_23, df_25, fit_21, fit_23, fit_25, file="gompertzModel_earlyStopping_21_to_25.rdata")










## Looking at cluster results

print(load("realDataOutputs/realData_outputs.rdata"))

growthRate_real

realGrowthModel_Output

ggplot(realGrowthModel_Output[[1]])+
  geom_col(aes(x=model, y=loo_IC, fill=model))+
  geom_segment(aes(x=model, xend=model, y=loo_IC, yend = loo_IC+loo_IC_se), size = 1.5)+
  geom_segment(aes(x=model, xend=model, y=loo_IC, yend = loo_IC-loo_IC_se), size = 1.5)+
  geom_text(aes(x=model, y=loo_IC*1.05, label = paste0(round(loo_IC, digits=0))))+
  labs(title = "LOO IC to estimate best model", x="", y="", fill="Model")+
  scale_fill_viridis(discrete = T)+
  theme_minimal() +
  theme(#legend.position="right",
    axis.line.y.left = element_line(),
    axis.line.x.bottom = element_line(),
    axis.text.x.bottom = element_blank())


print(load("realDataOutputs/realData_fit_models_and_data_every_2_days.rdata")) #A vs Control

h_15<-hypothesis(fit_5_to_15, "phi1_treatmenta/phi1_treatmentb > 1")
h_17<-hypothesis(fit_5_to_17, "phi1_treatmenta/phi1_treatmentb > 1")
h_19<-hypothesis(fit_5_to_19, "phi1_treatmenta/phi1_treatmentb > 1")
h_21<-hypothesis(fit_5_to_21, "phi1_treatmenta/phi1_treatmentb > 1")
h_23<-hypothesis(fit_5_to_23, "phi1_treatmenta/phi1_treatmentb > 1")
h_25<-hypothesis(fit_5_to_25, "phi1_treatmenta/phi1_treatmentb > 1")
h_15
h_17
h_19
h_21
h_23
h_25

df<-df%>%rename(treatment=Microbes, y=area_c, sample=Barcodes)%>%
  mutate(treatment=ifelse(treatment=="SynComA", "a", "b")) #super important encoding line
head(df)

earlyStopping<-df%>%
  mutate(Microbes=ifelse(treatment=="a", "SynComA", "Control"))%>%
  rename(area_c=y, Barcodes=sample)%>%
  ggplot(aes(DAP,area_c))+
  #facet_wrap(~Microbes)+
  geom_smooth(aes(color=Microbes,group=Barcodes), se=F, size=0.3)+
  #geom_smooth(color="gray20", se=F)+
  lapply(seq(15,25, 2), function(i) geom_segment(aes(x=i, xend=i, 
                                                 y=min(df$y[df$DAP==i]), yend=max(df$y[df$DAP==i])),
                                                 size=0.35, linetype=5))+
  lapply(seq(15,25, 2), function(i) annotate("text", x=i-0.6, y=max(df$y[df$DAP==i])*1.05,
         label=paste0("P=",round(eval(parse(text=paste0("h_",i,"$hypothesis$Post.Prob"))), 2)), size=3))+
  #annotate("text", x=12, y=160, label="Testing Hypothesis\n'phi1 SynComA/phi1 Control > 1'", size=3)+
  labs(x="Days After Planting", y="Area (cm)", title="Testing Hypothesis\nphi1 SynComA/phi1 Control > 1")+
  lims(x=c(5,25.5))+
  theme_light()+
  theme(legend.position="bottom")
earlyStopping
ggsave("Fig_3/earlyStopping.png",earlyStopping, width = 5.04,height=4.04,dpi=300)

darker.col = function(color, by = 40){
  colorRampPalette(c(color, "black"))(100)[by]
}

fancyEarlyStopping<-df%>%
  mutate(Microbes=ifelse(treatment=="a", "SynComA", "Control"))%>%
  rename(area_c=y, Barcodes=sample)%>%
  ggplot(aes(DAP,area_c))+
  geom_smooth(aes(color=Microbes,group=Barcodes), se=F, size=0.3, alpha=0.75)+
  lapply(seq(15,25, 2), function(i) geom_segment(aes(x=i+0.05, xend=i+0.05, 
                                                     y=min(df$y[df$DAP==i & df$treatment=="a"]), 
                                                     yend=max(df$y[df$DAP==i & df$treatment=="a"])),
                                                 size=0.6, linetype=5, color=darker.col("#00BFC4", by=20)))+
  lapply(seq(15,25, 2), function(i) geom_segment(aes(x=i-0.05, xend=i-0.05, 
                                                     y=min(df$y[df$DAP==i & df$treatment=="b"]), 
                                                     yend=max(df$y[df$DAP==i & df$treatment=="b"])),
                                                 size=0.6, linetype=5, color=darker.col("#F8766D", by=20)))+
  lapply(seq(15,25, 2), function(i) geom_point(aes(x=i, y=mean(df$y[df$DAP==i & df$treatment=="a"], na.rm=T)),
                                               size=3, shape="+", color=darker.col("#00BFC4", by=50)))+
  lapply(seq(15,25, 2), function(i) geom_point(aes(x=i, y=mean(df$y[df$DAP==i & df$treatment=="b"], na.rm=T)),
                                                 size=3, shape="+", color=darker.col("#F8766D", by=50)))+
  lapply(seq(15,25, 2), function(i) annotate("text", x=i-0.6, y=max(df$y[df$DAP==i])*1.05,
                                             label=paste0("P=",round(eval(parse(text=paste0("h_",i,"$hypothesis$Post.Prob"))), 2)), size=3))+
  labs(x="Days After Planting", y="Area (cm)", title="Testing Hypothesis\nphi1 SynComA/phi1 Control > 1")+
  lims(x=c(5,25.5))+
  theme_light()+
  theme(legend.position="bottom")
fancyEarlyStopping
ggsave("Fig_3/fancyEarlyStopping.png",fancyEarlyStopping, width = 5.04,height=4.04,dpi=300)



dfFull<-df%>%
  mutate(Microbes=ifelse(treatment=="a", "SynComA", "Control"))%>%
  rename(area_c=y, Barcodes=sample)
reallyFancyEarlyStopping<-df%>%
  filter(DAP<=23)%>%
  mutate(Microbes=ifelse(treatment=="a", "SynComA", "Control"))%>%
  rename(area_c=y, Barcodes=sample)%>%
  ggplot(aes(DAP,area_c))+
  geom_smooth(aes(color=Microbes,group=Barcodes), se=F, size=0.3, alpha=0.75)+
  geom_smooth(data=dfFull, aes(color=Microbes,group=Barcodes), se=F, size=0.1, alpha=0.25)+
  lapply(seq(15,23, 2), function(i) geom_segment(aes(x=i+0.05, xend=i+0.05, 
                                                     y=min(df$y[df$DAP==i & df$treatment=="a"]), 
                                                     yend=max(df$y[df$DAP==i & df$treatment=="a"])),
                                                 size=0.6, linetype=5, color=darker.col("#00BFC4", by=20)))+
  lapply(seq(15,23, 2), function(i) geom_segment(aes(x=i-0.05, xend=i-0.05, 
                                                     y=min(df$y[df$DAP==i & df$treatment=="b"]), 
                                                     yend=max(df$y[df$DAP==i & df$treatment=="b"])),
                                                 size=0.6, linetype=5, color=darker.col("#F8766D", by=20)))+
  lapply(seq(15,23, 2), function(i) geom_point(aes(x=i, y=mean(df$y[df$DAP==i & df$treatment=="a"], na.rm=T)),
                                               size=2, shape="+", color=darker.col("#00BFC4", by=50)))+
  lapply(seq(15,23, 2), function(i) geom_point(aes(x=i, y=mean(df$y[df$DAP==i & df$treatment=="b"], na.rm=T)),
                                               size=2, shape="+", color=darker.col("#F8766D", by=50)))+
  lapply(seq(15,23, 2), function(i) annotate("text", x=i-0.6, y=max(df$y[df$DAP==i])*1.05,
                                             label=paste0("P=",round(eval(parse(text=paste0("h_",i,"$hypothesis$Post.Prob"))), 2)), size=3))+
  labs(x="Days After Planting", y="Area (cm)", title="Testing Hypothesis\nphi1 SynComA/phi1 Control > 1")+
  lims(x=c(5,25.5))+
  theme_light()+
  theme(legend.position="bottom")
reallyFancyEarlyStopping
ggsave("Fig_3/reallyFancyEarlyStopping.png",reallyFancyEarlyStopping, width = 5.04,height=4.04,dpi=300)

#*****************************************************
# SynCom A vs SynCom B
#*****************************************************

print(load("realDataOutputs/realData_scAvscB_fit_models_and_data_every_2_days.rdata")) #A vs B



HYP<-"abs(phi1_treatmenta-phi1_treatmentb) > 1"

HYP<-"phi1_treatmenta/phi1_treatmentb > 1"

h_15<-hypothesis(fit_5_to_15, HYP)
h_17<-hypothesis(fit_5_to_17,HYP)
h_19<-hypothesis(fit_5_to_19, HYP)
h_21<-hypothesis(fit_5_to_21, HYP)
h_23<-hypothesis(fit_5_to_23, HYP)
h_25<-hypothesis(fit_5_to_25, HYP)

earlyStopping<-df%>%
  mutate(Microbes=ifelse(treatment=="a", "SynComA", "SynComB"))%>%
  rename(area_c=y, Barcodes=sample)%>%
  ggplot(aes(DAP,area_c))+
  geom_smooth(aes(color=Microbes,group=Barcodes), se=F, size=0.3)+
  lapply(seq(15,25, 2), function(i) geom_segment(aes(x=i, xend=i, 
                                                     y=min(df$y[df$DAP==i]), yend=max(df$y[df$DAP==i])),
                                                 size=0.35, linetype=5))+
  lapply(seq(15,25, 2), function(i) annotate("text", x=i-0.6, y=max(df$y[df$DAP==i])*1.05,
                                             label=paste0("P=",round(eval(parse(text=paste0("h_",i,"$hypothesis$Post.Prob"))), 3)), size=3))+
  #annotate("text", x=12, y=160, label="Testing Hypothesis\n'phi1 SynComA/phi1 Control > 1'", size=3)+
  labs(x="Days After Planting", y="Area (cm)")+
  lims(x=c(5,25.5))+
  theme_light()+
  theme(legend.position="bottom")
earlyStopping
ggsave("Fig_3/earlyStopping_AvB_h2.png",earlyStopping, width = 5.04,height=4.04,dpi=300)

darker.col = function(color, by = 40){
  colorRampPalette(c(color, "black"))(100)[by]
}

fancyEarlyStopping<-df%>%
  mutate(Microbes=ifelse(treatment=="a", "SynComA", "SynComB"))%>%
  rename(area_c=y, Barcodes=sample)%>%
  ggplot(aes(DAP,area_c))+
  geom_smooth(aes(color=Microbes,group=Barcodes), se=F, size=0.3, alpha=0.75)+
  #geom_smooth(aes(color=Microbes), se=F, size=0.3, alpha=0.75)+
  lapply(seq(15,25, 2), function(i) geom_segment(aes(x=i+0.05, xend=i+0.05, 
                                                     y=min(df$y[df$DAP==i & df$treatment=="a"]), 
                                                     yend=max(df$y[df$DAP==i & df$treatment=="a"])),
                                                 size=0.6, linetype=5, color=darker.col("#F8766D", by=20)))+
  lapply(seq(15,25, 2), function(i) geom_segment(aes(x=i-0.05, xend=i-0.05, 
                                                     y=min(df$y[df$DAP==i & df$treatment=="b"]), 
                                                     yend=max(df$y[df$DAP==i & df$treatment=="b"])),
                                                 size=0.6, linetype=5, color=darker.col("#00BFC4", by=20)))+
  lapply(seq(15,25, 2), function(i) geom_point(aes(x=i, y=mean(df$y[df$DAP==i & df$treatment=="a"], na.rm=T)),
                                               size=3, shape="+", color=darker.col("#F8766D", by=50)))+
  lapply(seq(15,25, 2), function(i) geom_point(aes(x=i, y=mean(df$y[df$DAP==i & df$treatment=="b"], na.rm=T)),
                                               size=3, shape="+", color=darker.col("#00BFC4", by=50)))+
  lapply(seq(15,25, 2), function(i) annotate("text", x=i-0.6, y=max(df$y[df$DAP==i])*1.05,
                                             label=paste0("P=",round(eval(parse(text=paste0("h_",i,"$hypothesis$Post.Prob"))), 2)), size=3))+
  labs(x="Days After Planting", y=~~Area~(cm^2))+
  lims(x=c(5,25.5))+
  theme_light()+
  theme(legend.position="bottom")
fancyEarlyStopping
ggsave("Fig_3/fancyEarlyStopping_AvB.png",fancyEarlyStopping, width = 5.04,height=4.04,dpi=300)



dfFull<-df%>%
  mutate(Microbes=ifelse(treatment=="a", "SynComA", "SynComB"))%>%
  rename(area_c=y, Barcodes=sample)
reallyFancyEarlyStopping<-df%>%
  filter(DAP<=23)%>%
  mutate(Microbes=ifelse(treatment=="a", "SynComA", "SynComB"))%>%
  rename(area_c=y, Barcodes=sample)%>%
  ggplot(aes(DAP,area_c))+
  geom_smooth(aes(color=Microbes,group=Barcodes), se=F, size=0.3, alpha=0.75)+
  geom_smooth(data=dfFull, aes(color=Microbes,group=Barcodes), se=F, size=0.1, alpha=0.25)+
  lapply(seq(15,23, 2), function(i) geom_segment(aes(x=i+0.05, xend=i+0.05, 
                                                     y=min(df$y[df$DAP==i & df$treatment=="a"]), 
                                                     yend=max(df$y[df$DAP==i & df$treatment=="a"])),
                                                 size=0.6, linetype=5, color=darker.col("#F8766D", by=20)))+
  lapply(seq(15,23, 2), function(i) geom_segment(aes(x=i-0.05, xend=i-0.05, 
                                                     y=min(df$y[df$DAP==i & df$treatment=="b"]), 
                                                     yend=max(df$y[df$DAP==i & df$treatment=="b"])),
                                                 size=0.6, linetype=5, color=darker.col("#00BFC4", by=20)))+
  lapply(seq(15,23, 2), function(i) geom_point(aes(x=i, y=mean(df$y[df$DAP==i & df$treatment=="a"], na.rm=T)),
                                               size=2, shape="+", color=darker.col("#F8766D", by=50)))+
  lapply(seq(15,23, 2), function(i) geom_point(aes(x=i, y=mean(df$y[df$DAP==i & df$treatment=="b"], na.rm=T)),
                                               size=2, shape="+", color=darker.col("#00BFC4", by=50)))+
  lapply(seq(15,23, 2), function(i) annotate("text", x=i-0.6, y=max(df$y[df$DAP==i])*1.05,
                                             label=paste0("P=",round(eval(parse(text=paste0("h_",i,"$hypothesis$Post.Prob"))), 2)), size=3))+
  labs(x="Days After Planting", y="Area (cm)", title="Testing Hypothesis\n|phi1 SynComA - phi1 SynComB| > 0")+
  lims(x=c(5,25.5))+
  theme_light()+
  theme(legend.position="bottom")
reallyFancyEarlyStopping
ggsave("Fig_3/reallyFancyEarlyStopping_AvB.png",reallyFancyEarlyStopping, width = 5.04,height=4.04,dpi=300)


# Cluster results on drought data

print(load("realDataOutputs/realData_scAvscB_fit_models_and_data_every_2_days_drought.rdata"))
# realData_scAvscB_fit_models_and_data_every_2_days_drought.rdata


print(load('realDataOutputs/mods21to25.rdata'))
print(load('realDataOutputs/mods13to17.rdata'))
print(load('realDataOutputs/mods17to19.rdata'))
print(load('realDataOutputs/mods9to11.rdata'))

HYP="phi1_treatmenta/phi1_treatmentb > 1"
for (i in seq(9,25,by=2)){
  hypText<-paste0("h_",i,"<-hypothesis(fit_5_to_",i,", HYP)")
  eval(parse(text=hypText))
}

head(df)

earlyStopping<-df%>%
  mutate(Microbes=ifelse(treatment=="a", "SynComA", "SynComB"))%>%
  rename(area_c=y, Barcodes=sample)%>%
  ggplot(aes(DAP,area_c))+
  #facet_wrap(~Microbes)+
  geom_smooth(aes(color=Microbes,group=Barcodes), se=F, size=0.3)+
  lapply(seq(9,25, 2), function(i) geom_segment(aes(x=i, xend=i, 
                                                     y=min(df$y[df$DAP==i]), yend=max(df$y[df$DAP==i])),
                                                 size=0.35, linetype=5))+
  lapply(seq(9,25, 2), function(i) annotate("text", x=i-0.6, y=max(df$y[df$DAP==i])*1.05,
                                             label=paste0("P=",round(eval(parse(text=paste0("h_",i,"$hypothesis$Post.Prob"))), 2)), size=3))+
  #annotate("text", x=12, y=160, label="Testing Hypothesis\n'phi1 SynComA/phi1 Control > 1'", size=3)+
  labs(x="Days After Planting", y=~~Area~(cm^2))+
  lims(x=c(5,25.5))+
  theme_light()+
  theme(legend.position="bottom")
earlyStopping
ggsave("Fig_3/earlyStopping_drought.png",earlyStopping, width = 5.04,height=4.04,dpi=300)

darker.col = function(color, by = 40){
  colorRampPalette(c(color, "black"))(100)[by]
}

#1(#27ccdb ,#ea39ed)
library(RColorBrewer)
dcol<-brewer.pal(n = 8, name = "Dark2")
colA<-dcol[8]#"#46327e"#dcol[5]#"#addc30"#"#f89540"#dcol[6]#"#440154"#27ccdb"#"#3522e3"
colB<-dcol[4]#"#fde725"#dcol[3]#"#472d7b"#"#7e03a8"#dcol[1]#"#fde725"#ea39ed"#ed6d2d"
fancyEarlyStopping<-df%>%
  mutate(Microbes=ifelse(treatment=="a", "SynComA", "SynComB"))%>%
  rename(area_c=y, Barcodes=sample)%>%
  ggplot(aes(DAP,area_c))+
  geom_smooth(aes(color=Microbes,group=Barcodes), se=F, size=0.3, alpha=0.75)+
  scale_color_manual(values = c(colB, colA))+ #1(#27ccdb ,#ea39ed)
  lapply(seq(13,25, 2), function(i) geom_segment(aes(x=i+0.05, xend=i+0.05, 
                                                     y=min(df$y[df$DAP==i & df$treatment=="b"]), 
                                                     yend=max(df$y[df$DAP==i & df$treatment=="b"])),
                                                 size=0.6, linetype=5, color=darker.col(colA, by=20)))+
  lapply(seq(13,25, 2), function(i) geom_segment(aes(x=i-0.05, xend=i-0.05, 
                                                     y=min(df$y[df$DAP==i & df$treatment=="a"]), 
                                                     yend=max(df$y[df$DAP==i & df$treatment=="a"])),
                                                 size=0.6, linetype=5, color=darker.col(colB, by=20)))+
  lapply(seq(13,25, 2), function(i) geom_point(aes(x=i, y=mean(df$y[df$DAP==i & df$treatment=="b"], na.rm=T)),
                                               size=3, shape="+", color=darker.col(colA, by=50)))+
  lapply(seq(13,25, 2), function(i) geom_point(aes(x=i, y=mean(df$y[df$DAP==i & df$treatment=="a"], na.rm=T)),
                                               size=3, shape="+", color=darker.col(colB, by=50)))+
  lapply(seq(13,25, 2), function(i) annotate("text", x=i-0.6, y=max(df$y[df$DAP==i])*1.05,
                                             label=paste0("P=",round(eval(parse(text=paste0("h_",i,"$hypothesis$Post.Prob"))), 2)), size=3))+
  labs(x="Days After Planting", y=~~Area~(cm^2))+
  lims(x=c(5,25.5))+
  theme_light()+
  theme(legend.position="bottom")
fancyEarlyStopping
ggsave("Fig_3/fancyEarlyStopping_AvB_drought_13-25_f.png",fancyEarlyStopping, width = 5.04,height=4.04,dpi=300)



#******************************************************************************************
# getGrowthRate 2 and getGrowthRatePlus on Drought SynComA vs SyncomB
#******************************************************************************************

library(brms)
library(ggplot2)
library(stringr)
library(patchwork)
library(plyr)
library(tidyverse)
library(bayesplot)
library(viridis)
library(viridisLite)

source('sourcingFunctions.R')
area_convert <- 13.2*3.7/46856
sv_shapes <- read.table("data.tsv",sep="\t",header=T,stringsAsFactors = F)
sv_shapes <- sv_shapes[sv_shapes$Drought == "ABB" & sv_shapes$Microbes %in% c("SynCom B","SynCom A"),] #add "SynCom B"
sv_shapes$area_c <- sv_shapes$area * area_convert

df <- aggregate(data=sv_shapes[sv_shapes$DAP > 4,], area_c ~ Microbes+Barcodes+DAP, FUN=function(i)mean(i))
df <- df[(df$Barcodes %in% df$Barcodes[df$DAP == 25]) & (df$Barcodes %in% df$Barcodes[df$DAP == 5]),]
df$Microbes <- gsub(" ","",df$Microbes)

df<-df%>%rename(treatment=Microbes, y=area_c, sample=Barcodes, time=DAP)%>%
  mutate(treatment=ifelse(treatment=="SynComA", "a", "b"))
head(df)

# getGrowthModel_real_drought_AvB<-getGrowthModel(dat=df, possibleModels = c("logistic", "gompertz", "monomolecular", "linear"))
# save(getGrowthModel_real_drought_AvB, "getGrowthModel_real_drought_AvB.rdata")

getGrowthModelPlus_real_drought_AvB_3_iteration<-getGrowthModelPlus(dat=df, possibleModels = c("logistic", "gompertz", "linear", "spline"),
                                                                                               comparisonIterations = 3, comparisonModelIterations = 1000, 
                                                                                               bestModelIterations=2000, comparisonSigma="spline", sigma = "spline")
save(getGrowthModelPlus_real_drought_AvB_3_iteration, "getGrowthModelPlus_real_drought_AvB_3_iteration.rdata")

# getGrowthModelPlus_real_drought_AvB<-getGrowthModelPlus(dat=df, possibleModels = c("logistic", "gompertz", "monomolecular", "linear", "spline"))
# save(getGrowthModelPlus_real_drought_AvB, "getGrowthModelPlus_real_drought_AvB.rdata")

#******************************************************************************************
# getGrowthRate Results
#******************************************************************************************


print(load("realDataOutputs/getGrowthModelPlus_real_drought_AvB.rdata"))
getGrowthModelPlus_real_drought_AvB[2][[1]]
getGrowthModelPlus_real_drought_AvB[5][[1]]

getGrowthModelPlus_real_drought_AvB[4][[1]]


print(load("realDataOutputs/getGrowthModel_real_drought_AvB.rdata"))

getGrowthModel_real_drought_AvB[2][[1]]
getGrowthModel_real_drought_AvB[5][[1]]

getGrowthModel_real_drought_AvB[4][[1]]


#getGrowthModelPlus_real_drought_AvB_3_iteration

print(load("realDataOutputs/getGrowthModelPlus_real_drought_AvB_3_iteration.rdata"))

getGrowthModelPlus_real_drought_AvB_3_iteration[2][[1]]
output_df<-getGrowthModelPlus_real_drought_AvB_3_iteration[1][[1]]

output_df$model <- with(output_df, reorder(model, -loo_IC))

looPlot<-ggplot(output_df)+
  geom_col(aes(x=model, y=loo_IC, fill=model))+
  geom_segment(aes(x=model, xend=model, y=loo_IC, yend = loo_IC+loo_IC_se), size = 1.5)+
  geom_segment(aes(x=model, xend=model, y=loo_IC, yend = loo_IC-loo_IC_se), size = 1.5)+
  geom_text(aes(x=model, y=loo_IC*1.05, label = paste0(round(loo_IC, digits=0))))+
  labs(title = "LOO IC", x="", y="", fill="Model")+
  scale_fill_viridis(discrete = T)+
  theme_minimal() +
  theme(axis.line.y.left = element_line(),
        axis.line.x.bottom = element_line(),
        axis.text.x.bottom = element_blank())
looPlot
getGrowthModelPlus_real_drought_AvB_3_iteration[5][[1]]

getGrowthModelPlus_real_drought_AvB_3_iteration[4][[1]]

#**************************************************
# Early stopping in gompertz model of real data
#**************************************************

print(load("realDataOutputs/gompertzEarlyStopping/gompertzModel_earlyStopping_7_to_13.rdata"))
print(load("realDataOutputs/gompertzEarlyStopping/gompertzModel_earlyStopping_15_to_19.rdata"))
print(load("realDataOutputs/gompertzEarlyStopping/gompertzModel_earlyStopping_21_to_25.rdata"))

area_convert <- 13.2*3.7/46856
sv_shapes <- read.table("b51632508ce74bdd9ceffcd073b871e6-shapes_outliersRM.tsv",sep="\t",header=T,stringsAsFactors = F)
sv_shapes <- sv_shapes[sv_shapes$Drought == "ABB" & sv_shapes$Microbes %in% c("SynCom B","SynCom A"),] #add "SynCom B"
sv_shapes$area_c <- sv_shapes$area * area_convert
df <- aggregate(data=sv_shapes[sv_shapes$DAP > 4,], area_c ~ Microbes+Barcodes+DAP, FUN=function(i)mean(i))
df <- df[(df$Barcodes %in% df$Barcodes[df$DAP == 25]) & (df$Barcodes %in% df$Barcodes[df$DAP == 5]),]
df$Microbes <- gsub(" ","",df$Microbes)
df<-df%>%rename(treatment=Microbes, y=area_c, sample=Barcodes)%>%
  mutate(treatment=ifelse(treatment=="SynComA", "a", "b"))
head(df)

HYP="a_treatmenta/a_treatmentb > 1"
for (i in seq(7,25,by=2)){
  hypText<-paste0("h_",i,"<-hypothesis(fit_",i,", HYP)")
  eval(parse(text=hypText))
}

head(df)

darker.col = function(color, by = 40){
  colorRampPalette(c(color, "black"))(100)[by]
}

#1(#27ccdb ,#ea39ed)
library(RColorBrewer)
dcol<-brewer.pal(n = 8, name = "Dark2")
colA<-dcol[6]#"#46327e"#dcol[5]#"#addc30"#"#f89540"#dcol[6]#"#440154"#27ccdb"#"#3522e3"
colB<-dcol[3]#"#fde725"#dcol[3]#"#472d7b"#"#7e03a8"#dcol[1]#"#fde725"#ea39ed"#ed6d2d"
fancyEarlyStopping<-df%>%
  mutate(Microbes=ifelse(treatment=="a", "SynComA", "SynComB"))%>%
  rename(area_c=y, Barcodes=sample)%>%
  ggplot(aes(DAP,area_c))+
  geom_smooth(aes(color=Microbes,group=Barcodes), se=F, size=0.3, alpha=0.75)+
  scale_color_manual(values = c(colB, colA))+ #1(#27ccdb ,#ea39ed)
  lapply(seq(7,25, 2), function(i) geom_segment(aes(x=i+0.05, xend=i+0.05, 
                                                     y=min(df$y[df$DAP==i & df$treatment=="b"]), 
                                                     yend=max(df$y[df$DAP==i & df$treatment=="b"])),
                                                 size=0.6, linetype=5, color=darker.col(colA, by=20)))+
  lapply(seq(7,25, 2), function(i) geom_segment(aes(x=i-0.05, xend=i-0.05, 
                                                     y=min(df$y[df$DAP==i & df$treatment=="a"]), 
                                                     yend=max(df$y[df$DAP==i & df$treatment=="a"])),
                                                 size=0.6, linetype=5, color=darker.col(colB, by=20)))+
  lapply(seq(7,25, 2), function(i) geom_point(aes(x=i, y=mean(df$y[df$DAP==i & df$treatment=="b"], na.rm=T)),
                                               size=3, shape="+", color=darker.col(colA, by=50)))+
  lapply(seq(7,25, 2), function(i) geom_point(aes(x=i, y=mean(df$y[df$DAP==i & df$treatment=="a"], na.rm=T)),
                                               size=3, shape="+", color=darker.col(colB, by=50)))+
  lapply(seq(7,25, 2), function(i) annotate("text", x=i-0.6, y=max(df$y[df$DAP==i])*1.05,
                                             label=paste0("P=",round(eval(parse(text=paste0("h_",i,"$hypothesis$Post.Prob"))), 2)), size=3))+
  labs(x="Days After Planting", y=~~Area~(cm^2))+
  lims(x=c(5,25.5))+
  theme_light()+
  theme(legend.position="bottom")
fancyEarlyStopping

ggsave("Fig_3/fancyEarlyStopping_AvB_drought_7-25_gompertzMod.png",fancyEarlyStopping, width = 5.04,height=4.04,dpi=300)



