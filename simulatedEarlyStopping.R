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

growthSimMM <- function(x,a,b){
  a_r <- a+rnorm(1,mean = 0,sd=10)
  b_r <- b+rnorm(1,mean=0,sd=0.025)
  return(a_r-a_r*exp(-b_r*x))
}

x <- 1:25
df <- rbind(
  do.call(rbind,lapply(1:20,function(i) data.frame("sample"=paste0("sample_",i),"treatment"="a","time"=x,
                                                   "y"=growthSimMM(x,200,0.1),stringsAsFactors = F))),#from .33
  do.call(rbind,lapply(1:20,function(i) data.frame("sample"=paste0("sample_",i),"treatment"="b","time"=x,
                                                   "y"=growthSimMM(x,160,.12),stringsAsFactors = F)))#from .42
)

ggplot(df,aes(time,y,group=interaction(treatment,sample)))+
  geom_line(aes(color=treatment))+
  ylab(~~Area~(cm^2))+
  xlab("Time")+
  #scale_y_continuous(limits = c(0,300))+
  theme_light()+
  theme(axis.ticks.length=unit(0.2,"cm"))+
  theme(strip.background=element_rect(fill="gray50",color="gray20"),
        strip.text.x=element_text(size=14,color="white"),
        strip.text.y=element_text(size=14,color="white"))+
  theme(axis.title= element_text(size = 18))+
  theme(axis.text = element_text(size = 14))+
  theme(legend.position='top')

df_5<-df%>%filter(time<=5)
df_7<-df%>%filter(time<=7)
df_9<-df%>%filter(time<=9)
df_11<-df%>%filter(time<=11)
df_13<-df%>%filter(time<=13)
df_15<-df%>%filter(time<=15)
df_17<-df%>%filter(time<=17)
df_19<-df%>%filter(time<=19)
df_21<-df%>%filter(time<=21)
df_23<-df%>%filter(time<=23)
df_25<-df%>%filter(time<=25)

priorMono <- prior(lognormal(log(130), .25),nlpar = "a") +
  prior(lognormal(log(2), .25), nlpar = "b") + 
  prior(student_t(3,0,5), dpar="sigma") +
  prior(gamma(2,0.1), class="nu")

fitMono_5 <- brm(bf(y ~ a-a*exp(-b*time), 
                     sigma~s(time, by=treatment, k=5), 
                     a + b ~ 0+treatment, 
                     autocor = ~arma(~time|sample:treatment,1,1),nl = TRUE),
                  family = student, prior = priorMono, data = df_5, iter = 2000, 
                  cores = 2, chains = 2, backend = "cmdstanr",
                  control = list(adapt_delta = 0.999,max_treedepth = 20),
                  inits = function(){list(b_a=rgamma(2,1),b_b=rgamma(2,1),b_c=rgamma(2,1))})

fitMono_7 <- brm(bf(y ~ a-a*exp(-b*time), 
                     sigma~s(time, by=treatment, k=7), 
                     a + b ~ 0+treatment, 
                     autocor = ~arma(~time|sample:treatment,1,1),nl = TRUE),
                  family = student, prior = priorMono, data = df_7, iter = 2000, 
                  cores = 2, chains = 2, backend = "cmdstanr",
                  control = list(adapt_delta = 0.999,max_treedepth = 20),
                  inits = function(){list(b_a=rgamma(2,1),b_b=rgamma(2,1),b_c=rgamma(2,1))})

fitMono_9 <- brm(bf(y ~ a-a*exp(-b*time), 
                     sigma~s(time, by=treatment, k=9), 
                     a + b ~ 0+treatment, 
                     autocor = ~arma(~time|sample:treatment,1,1),nl = TRUE),
                  family = student, prior = priorMono, data = df_9, iter = 2000, 
                  cores = 2, chains = 2, backend = "cmdstanr",
                  control = list(adapt_delta = 0.999,max_treedepth = 20),
                  inits = function(){list(b_a=rgamma(2,1),b_b=rgamma(2,1),b_c=rgamma(2,1))})

fitMono_11 <- brm(bf(y ~ a-a*exp(-b*time), 
                     sigma~s(time, by=treatment), 
                     a + b ~ 0+treatment, 
                     autocor = ~arma(~time|sample:treatment,1,1),nl = TRUE),
                  family = student, prior = priorMono, data = df_11, iter = 2000, 
                  cores = 2, chains = 2, backend = "cmdstanr",
                  control = list(adapt_delta = 0.999,max_treedepth = 20),
                  inits = function(){list(b_a=rgamma(2,1),b_b=rgamma(2,1),b_c=rgamma(2,1))})

fitMono_13 <- brm(bf(y ~ a-a*exp(-b*time), 
                     sigma~s(time, by=treatment), 
                     a + b ~ 0+treatment, 
                     autocor = ~arma(~time|sample:treatment,1,1),nl = TRUE),
                  family = student, prior = priorMono, data = df_13, iter = 2000, 
                  cores = 2, chains = 2, backend = "cmdstanr",
                  control = list(adapt_delta = 0.999,max_treedepth = 20),
                  inits = function(){list(b_a=rgamma(2,1),b_b=rgamma(2,1),b_c=rgamma(2,1))})

fitMono_15 <- brm(bf(y ~ a-a*exp(-b*time), 
                   sigma~s(time, by=treatment), 
                   a + b ~ 0+treatment, 
                   autocor = ~arma(~time|sample:treatment,1,1),nl = TRUE),
                family = student, prior = priorMono, data = df_15, iter = 2000, 
                cores = 2, chains = 2, backend = "cmdstanr",
                control = list(adapt_delta = 0.999,max_treedepth = 20),
                inits = function(){list(b_a=rgamma(2,1),b_b=rgamma(2,1),b_c=rgamma(2,1))})

fitMono_17 <- brm(bf(y ~ a-a*exp(-b*time), 
                     sigma~s(time, by=treatment), 
                     a + b ~ 0+treatment, 
                     autocor = ~arma(~time|sample:treatment,1,1),nl = TRUE),
                  family = student, prior = priorMono, data = df_17, iter = 2000, 
                  cores = 2, chains = 2, backend = "cmdstanr", 
                  control = list(adapt_delta = 0.999,max_treedepth = 20),
                  inits = function(){list(b_a=rgamma(2,1),b_b=rgamma(2,1),b_c=rgamma(2,1))})

fitMono_19 <- brm(bf(y ~ a-a*exp(-b*time), 
                     sigma~s(time, by=treatment), 
                     a + b ~ 0+treatment, 
                     autocor = ~arma(~time|sample:treatment,1,1),nl = TRUE),
                  family = student, prior = priorMono, data = df_19, iter = 2000, 
                  cores = 2, chains = 2, backend = "cmdstanr", 
                  control = list(adapt_delta = 0.999,max_treedepth = 20),
                  inits = function(){list(b_a=rgamma(2,1),b_b=rgamma(2,1),b_c=rgamma(2,1))})

fitMono_21 <- brm(bf(y ~ a-a*exp(-b*time), 
                     sigma~s(time, by=treatment), 
                     a + b ~ 0+treatment, 
                     autocor = ~arma(~time|sample:treatment,1,1),nl = TRUE),
                  family = student, prior = priorMono, data = df_21, iter = 2000, 
                  cores = 2, chains = 2, backend = "cmdstanr", 
                  control = list(adapt_delta = 0.999,max_treedepth = 20),
                  inits = function(){list(b_a=rgamma(2,1),b_b=rgamma(2,1),b_c=rgamma(2,1))})

fitMono_23 <- brm(bf(y ~ a-a*exp(-b*time), 
                     sigma~s(time, by=treatment), 
                     a + b ~ 0+treatment, 
                     autocor = ~arma(~time|sample:treatment,1,1),nl = TRUE),
                  family = student, prior = priorMono, data = df_23, iter = 2000, 
                  cores = 2, chains = 2, backend = "cmdstanr", 
                  control = list(adapt_delta = 0.999,max_treedepth = 20),
                  inits = function(){list(b_a=rgamma(2,1),b_b=rgamma(2,1),b_c=rgamma(2,1))})

fitMono_25 <- brm(bf(y ~ a-a*exp(-b*time), 
                     sigma~s(time, by=treatment), 
                     a + b ~ 0+treatment, 
                     autocor = ~arma(~time|sample:treatment,1,1),nl = TRUE),
                  family = student, prior = priorMono, data = df_25, iter = 2000, 
                  cores = 2, chains = 2, backend = "cmdstanr", 
                  control = list(adapt_delta = 0.999,max_treedepth = 20),
                  inits = function(){list(b_a=rgamma(2,1),b_b=rgamma(2,1),b_c=rgamma(2,1))})


save(df,priorMono,fitMono_5, fitMono_7,fitMono_9,fitMono_11,fitMono_13, fitMono_15,fitMono_17,
     fitMono_19,fitMono_21,fitMono_23, fitMono_25,
     file = "monomolecular_early_stopping.rdata")

# done on cluster 1

print(load("simulatedEarlyStopping/monomolecular_early_stopping.rdata"))

HYP="a_treatmenta/a_treatmentb > 1"
h_13<-hypothesis(fitMono_13, HYP)
h_15<-hypothesis(fitMono_15, HYP)
h_17<-hypothesis(fitMono_17, HYP)
h_19<-hypothesis(fitMono_19, HYP)
h_21<-hypothesis(fitMono_21, HYP)
h_23<-hypothesis(fitMono_23, HYP)
h_25<-hypothesis(fitMono_25, HYP)

darker.col = function(color, by = 40){
  colorRampPalette(c(color, "black"))(100)[by]
}

fancyEarlyStopping<-df%>%
  ggplot(aes(time,y,group=interaction(treatment,sample)))+
  geom_line(aes(color=treatment))+
  lapply(seq(13,25, 2), function(i) geom_segment(aes(x=i+0.05, xend=i+0.05, 
                                                     y=min(df$y[df$time==i & df$treatment=="b"]), 
                                                     yend=max(df$y[df$time==i & df$treatment=="b"])),
                                                 size=0.6, linetype=5, color=darker.col("#00BFC4", by=20)))+
  lapply(seq(13,25, 2), function(i) geom_segment(aes(x=i-0.05, xend=i-0.05, 
                                                     y=min(df$y[df$time==i & df$treatment=="a"]), 
                                                     yend=max(df$y[df$time==i & df$treatment=="a"])),
                                                 size=0.6, linetype=5, color=darker.col("#F8766D", by=20)))+
  lapply(seq(13,25, 2), function(i) geom_point(aes(x=i, y=mean(df$y[df$time==i & df$treatment=="b"], na.rm=T)),
                                               size=3, shape="+", color=darker.col("#00BFC4", by=50)))+
  lapply(seq(13,25, 2), function(i) geom_point(aes(x=i, y=mean(df$y[df$time==i & df$treatment=="a"], na.rm=T)),
                                               size=3, shape="+", color=darker.col("#F8766D", by=50)))+
  lapply(seq(13,25, 2), function(i) annotate("text", x=i-0.6, y=max(df$y[df$time==i])*1.05,
                                             label=paste0("P=",round(eval(parse(text=paste0("h_",i,"$hypothesis$Post.Prob"))), 3)), size=3))+
  labs(x="Days After Planting", y="Area (cm)", title="Testing Hypothesis\n A/B > 1")+
  lims(x=c(1,25.5))+
  theme_light()+
  theme(legend.position="bottom")
fancyEarlyStopping
ggsave("simulatedEarlyStopping/monoEarlyStopping.png",fancyEarlyStopping, width = 5.04,height=4.04,dpi=300)


# done on cluster 2

print(load("simulatedEarlyStopping/monomolecular_early_stopping_2.rdata"))

HYP="a_treatmenta/a_treatmentb > 1"
h_5<-hypothesis(fitMono_5, HYP)
h_7<-hypothesis(fitMono_7, HYP)
h_9<-hypothesis(fitMono_9, HYP)
h_11<-hypothesis(fitMono_11, HYP)
h_13<-hypothesis(fitMono_13, HYP)
h_15<-hypothesis(fitMono_15, HYP)
h_17<-hypothesis(fitMono_17, HYP)
h_19<-hypothesis(fitMono_19, HYP)
h_21<-hypothesis(fitMono_21, HYP)
h_23<-hypothesis(fitMono_23, HYP)
h_25<-hypothesis(fitMono_25, HYP)

darker.col = function(color, by = 40){
  colorRampPalette(c(color, "black"))(100)[by]
}

fancyEarlyStopping<-df%>%
  ggplot(aes(time,y,group=interaction(treatment,sample)))+
  geom_line(aes(color=treatment))+
  lapply(seq(5,25, 2), function(i) geom_segment(aes(x=i+0.05, xend=i+0.05, 
                                                     y=min(df$y[df$time==i & df$treatment=="b"]), 
                                                     yend=max(df$y[df$time==i & df$treatment=="b"])),
                                                 size=0.6, linetype=5, color=darker.col("#00BFC4", by=20)))+
  lapply(seq(5,25, 2), function(i) geom_segment(aes(x=i-0.05, xend=i-0.05, 
                                                     y=min(df$y[df$time==i & df$treatment=="a"]), 
                                                     yend=max(df$y[df$time==i & df$treatment=="a"])),
                                                 size=0.6, linetype=5, color=darker.col("#F8766D", by=20)))+
  lapply(seq(5,25, 2), function(i) geom_point(aes(x=i, y=mean(df$y[df$time==i & df$treatment=="b"], na.rm=T)),
                                               size=3, shape="+", color=darker.col("#00BFC4", by=50)))+
  lapply(seq(5,25, 2), function(i) geom_point(aes(x=i, y=mean(df$y[df$time==i & df$treatment=="a"], na.rm=T)),
                                               size=3, shape="+", color=darker.col("#F8766D", by=50)))+
  lapply(seq(5,25, 2), function(i) annotate("text", x=i-0.6, y=max(df$y[df$time==i])*1.05,
                                             label=paste0("P=",
                                                          round(eval(parse(text=paste0("h_",i,"$hypothesis$Post.Prob"))),
                                                                2)), size=3))+
  labs(x="Days After Planting", y="Area (cm)", title="Testing Hypothesis A/B > 1")+
  lims(x=c(1,25.5))+
  theme_light()+
  theme(legend.position="bottom", title = element_text(size=10))
fancyEarlyStopping
ggsave("simulatedEarlyStopping/monoEarlyStopping_2.png",fancyEarlyStopping, width = 5.04,height=4.04,dpi=300)

#**************************************************
# Gompertz
#*************************************************

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

growthSim <- function(x,a,b,c){
  a_r <- a+rnorm(1,mean = 0,sd=10)
  b_r <- b+rnorm(1,mean=0,sd=2)
  c_r <- c+rnorm(1,mean=0,sd=.035)
  return(a_r*exp(-b_r*exp(-c_r*x)))}

x <- 1:25
set.seed(234)
df <- rbind(
  do.call(rbind,lapply(1:20,function(i) data.frame("sample"=paste0("sample_",i),"treatment"="a","time"=x,"y"=growthSim(x,200,13,.2),stringsAsFactors = F))),
  do.call(rbind,lapply(1:20,function(i) data.frame("sample"=paste0("sample_",i),"treatment"="b","time"=x,"y"=growthSim(x,180,13,.25),stringsAsFactors = F))))

priorGomp <- prior(lognormal(log(130), .25),nlpar = "a") +
  prior(lognormal(log(12), .25), nlpar = "b") + 
  prior(lognormal(log(1.2), .25), nlpar = "c") + 
  prior(student_t(3,0,5), dpar="sigma") +
  prior(gamma(2,0.1), class="nu")

for (i in seq(7,25,by=2)){
  dataText<-paste0("df_",i,"<-df%>%filter(time<=",i,")")
  eval(parse(text=dataText))
  if(i<13){
    sigmaText<-paste0('sigma~s(time, by=treatment, k=',i,'),')
  } else{
    sigmaText<-paste0('sigma~s(time, by=treatment),')}
  brmText<-paste0('fit_',i,' <- brm(bf(y ~ a*exp(-b*exp(-c*time)),',
                  sigmaText,
                  'a + b + c ~ 0+treatment, autocor = ~arma(~time|sample:treatment,1,1),nl = TRUE),family = student, prior = priorGomp, data = df_',i,', iter = 1000, cores = 2, chains = 2, backend = "cmdstanr",control = list(adapt_delta = 0.999,max_treedepth = 20),inits = function(){list(b_a=rgamma(2,1),b_b=rgamma(2,1),b_c=rgamma(2,1))})')
  print(i)
  eval(parse(text=brmText))
}

save(df,priorGomp,fit_7,fit_9,fit_11,fit_13,
     fit_15,fit_17,fit_19,fit_21,fit_23,fit_25,
     file = "gompertz_early_stopping.rdata")

print(load("simulatedEarlyStopping/gompertz_early_stopping.rdata"))

HYP="a_treatmenta/a_treatmentb > 1"
for (i in seq(7,25,by=2)){
  hypText<-paste0("h_",i,"<-hypothesis(fit_",i,", HYP)")
  eval(parse(text=hypText))
}

darker.col = function(color, by = 40){
  colorRampPalette(c(color, "black"))(100)[by]
}

fancyEarlyStopping<-df%>%
  ggplot(aes(time,y,group=interaction(treatment,sample)))+
  geom_line(aes(color=treatment))+
  lapply(seq(7,25, 2), function(i) geom_segment(aes(x=i+0.05, xend=i+0.05, 
                                                    y=min(df$y[df$time==i & df$treatment=="b"]), 
                                                    yend=max(df$y[df$time==i & df$treatment=="b"])),
                                                size=0.6, linetype=5, color=darker.col("#00BFC4", by=20)))+
  lapply(seq(7,25, 2), function(i) geom_segment(aes(x=i-0.05, xend=i-0.05, 
                                                    y=min(df$y[df$time==i & df$treatment=="a"]), 
                                                    yend=max(df$y[df$time==i & df$treatment=="a"])),
                                                size=0.6, linetype=5, color=darker.col("#F8766D", by=20)))+
  lapply(seq(7,25, 2), function(i) geom_point(aes(x=i, y=mean(df$y[df$time==i & df$treatment=="b"], na.rm=T)),
                                              size=3, shape="+", color=darker.col("#00BFC4", by=50)))+
  lapply(seq(7,25, 2), function(i) geom_point(aes(x=i, y=mean(df$y[df$time==i & df$treatment=="a"], na.rm=T)),
                                              size=3, shape="+", color=darker.col("#F8766D", by=50)))+
  lapply(seq(7,25, 2), function(i) annotate("text", x=i-0.6, y=max(df$y[df$time==i])*1.05,
                                            label=paste0("P=",
                                                         round(eval(parse(text=paste0("h_",i,"$hypothesis$Post.Prob"))),
                                                               2)), size=3))+
  labs(x="Days After Planting", y="Area (cm)", title="Testing Hypothesis A/B > 1")+
  lims(x=c(1,25.5))+
  theme_light()+
  theme(legend.position="bottom", title = element_text(size=10))
fancyEarlyStopping
ggsave("simulatedEarlyStopping/gompertzEarlyStopping.png",fancyEarlyStopping, width = 5.04,height=4.04,dpi=300)


#**************************************************
# Logistic
#**************************************************

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

growthSim <- function(x,phi1,phi2,phi3){ # make function to do a bunch of growthSims and check how well they 
  phi1_r <- phi1+rnorm(1,mean = 0,sd=25) #can be modeled/how well their phi params can be recovered.
  phi2_r <- phi2+rnorm(1,mean=0,sd=1)
  phi3_r <- phi3+rnorm(1,mean=0,sd=.2)
  return(phi1_r / (1 + exp( (phi2_r-x)/phi3_r) ))
}
set.seed(123)
x <- 1:25
df <- rbind(
  do.call(rbind, lapply(1:20, function(i) data.frame("sample"=paste0("sample_",i),"treatment"="a","time"=x,"y"=growthSim(x,200,13,3),stringsAsFactors = F))),
  do.call(rbind, lapply(1:20,function(i) data.frame("sample"=paste0("sample_",i),"treatment"="b","time"=x,"y"=growthSim(x,160,13,3.5),stringsAsFactors = F))))

priorLogistic <- prior(lognormal(log(130), .25),nlpar = "phi1",coef="treatmenta") +
  prior(lognormal(log(130), .25),nlpar = "phi1",coef="treatmentb") +
  prior(lognormal(log(12), .25), nlpar = "phi2",coef="treatmenta") +
  prior(lognormal(log(12), .25), nlpar = "phi2",coef="treatmentb") +
  prior(lognormal(log(3), .25), nlpar = "phi3",coef="treatmenta") +
  prior(lognormal(log(3), .25), nlpar = "phi3",coef="treatmentb") +
  prior(student_t(3,0,5), dpar="sigma") +
  prior(gamma(2,0.1), class="nu")

for (i in seq(7,25,by=2)){
  dataText<-paste0("df_",i,"<-df%>%filter(time<=",i,")")
  eval(parse(text=dataText))
  if(i<13){
    sigmaText<-paste0('sigma~s(time, by=treatment, k=',i,'),')
  } else{
    sigmaText<-paste0('sigma~s(time, by=treatment),')}
  brmText<-paste0('fit_',i,' <- brm(bf(y ~ phi1/(1+exp((phi2-time)/phi3)),',
                  sigmaText,
                  'phi1 + phi2 + phi3 ~ 0+treatment, autocor = ~arma(~time|sample:treatment,1,1),nl = TRUE),family = student, prior = priorLogistic, data = df_',i,', iter = 1000, cores = 2, chains = 2, backend = "cmdstanr",control = list(adapt_delta = 0.999,max_treedepth = 20),inits = function(){list(b_phi1=rgamma(2,1),b_phi2=rgamma(2,1),b_phi3=rgamma(2,1))})')
  print(i)
  eval(parse(text=brmText))
}

save(df,priorLogistic,fit_7,fit_9,fit_11,fit_13,
     fit_15,fit_17,fit_19,fit_21,fit_23,fit_25,
     file = "logistic_early_stopping.rdata")

print(load("simulatedEarlyStopping/logistic_early_stopping.rdata"))

HYP="phi1_treatmenta/phi1_treatmentb > 1"
for (i in seq(7,25,by=2)){
  hypText<-paste0("h_",i,"<-hypothesis(fit_",i,", HYP)")
  eval(parse(text=hypText))
}

darker.col = function(color, by = 40){
  colorRampPalette(c(color, "black"))(100)[by]
}

fancyEarlyStopping<-df%>%
  ggplot(aes(time,y,group=interaction(treatment,sample)))+
  geom_line(aes(color=treatment))+
  lapply(seq(7,25, 2), function(i) geom_segment(aes(x=i+0.05, xend=i+0.05, 
                                                    y=min(df$y[df$time==i & df$treatment=="b"]), 
                                                    yend=max(df$y[df$time==i & df$treatment=="b"])),
                                                size=0.6, linetype=5, color=darker.col("#00BFC4", by=20)))+
  lapply(seq(7,25, 2), function(i) geom_segment(aes(x=i-0.05, xend=i-0.05, 
                                                    y=min(df$y[df$time==i & df$treatment=="a"]), 
                                                    yend=max(df$y[df$time==i & df$treatment=="a"])),
                                                size=0.6, linetype=5, color=darker.col("#F8766D", by=20)))+
  lapply(seq(7,25, 2), function(i) geom_point(aes(x=i, y=mean(df$y[df$time==i & df$treatment=="b"], na.rm=T)),
                                              size=3, shape="+", color=darker.col("#00BFC4", by=50)))+
  lapply(seq(7,25, 2), function(i) geom_point(aes(x=i, y=mean(df$y[df$time==i & df$treatment=="a"], na.rm=T)),
                                              size=3, shape="+", color=darker.col("#F8766D", by=50)))+
  lapply(seq(7,25, 2), function(i) annotate("text", x=i-0.6, y=max(df$y[df$time==i])*1.05,
                                            label=paste0("P=",
                                                         round(eval(parse(text=paste0("h_",i,"$hypothesis$Post.Prob"))),
                                                               2)), size=3))+
  labs(x="Days After Planting", y="Area (cm)", title="Testing Hypothesis phi1a/phi1b > 1")+
  lims(x=c(1,25.5))+
  theme_light()+
  theme(legend.position="bottom", title = element_text(size=8))
fancyEarlyStopping
ggsave("simulatedEarlyStopping/logisticEarlyStopping.png",fancyEarlyStopping, width = 5.04,height=4.04,dpi=300)





