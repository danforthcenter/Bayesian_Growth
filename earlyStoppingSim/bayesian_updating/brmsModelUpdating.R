# Goal in this file is to follow issue 51 (https://github.com/joshqsumner/bayesian_plant-science/issues/51)
#   Want to find a way to efficiently update a brms model using new data and recycled posterior distributions
#   Currently I'm imagining fitting a model on day 7, updating it on days 10, 13, 16, 19, 22, and 25.
#   Performance would be judged against a model fit anew with the entire data for each day.
#   Foreseen problems: 
#     Divergent transitions from update()
#     Difficulty recycling posteriors as priors
#     Very small time improvements for significant added complexity
#     Diminishing performance on already established parameters (inflection for instance) on later days
#       Could be partially resolved by monitoring change in HPDI for each parameter?
#       AKA if change in hpdi for inflection approaches 0 by day 17 then 
#         picks up when updating on day 25 that seems like a red flag, also hard to monitor.

library(brms)
library(tidyverse)
library(viridis)
library(patchwork)
library(rethinking)
library(scales)

# ********************
# Set up data
# ********************

growthSim <- function(x,phi1,phi2,phi3){
  phi1_r <- phi1+rnorm(1,mean = 0,sd=25)
  phi2_r <- phi2+rnorm(1,mean=0,sd=1)
  phi3_r <- phi3+rnorm(1,mean=0,sd=.2)
  return(phi1_r / (1 + exp( (phi2_r-x)/phi3_r)))
}
set.seed(123)
dfLogistic <- rbind(
    do.call(rbind, lapply(1:20, function(i) data.frame("sample"=paste0("sample_",i),
                                                       "treatment"="a","time"=1:67,
                                                       "y"=growthSim(1:67, 200, 13, 3),
                                                       stringsAsFactors = F))),
    do.call(rbind, lapply(1:20,function(i) data.frame("sample"=paste0("sample_",i),
                                                      "treatment"="b","time"=1:67,
                                                      "y"=growthSim(1:67, 160, 13, 3.5),
                                                      stringsAsFactors = F))))

ggplot(dfLogistic,aes(time,y,group=interaction(treatment,sample)))+
  geom_line(aes(color=treatment))+
  ylab(~~Area~(cm^2))+
  xlab("Time")+
  theme_light()+
  theme(axis.ticks.length=unit(0.2,"cm"))+
  theme(strip.background=element_rect(fill="gray50",color="gray20"),
        strip.text.x=element_text(size=14,color="white"),
        strip.text.y=element_text(size=14,color="white"))+
  theme(axis.title= element_text(size = 18))+
  theme(axis.text = element_text(size = 14))+
  theme(legend.position='top')

# ********************
# Fit Full Model
# ********************

priorList <- prior(lognormal(log(130), .25),nlpar = "phi1") +
  prior(lognormal(log(12), .25), nlpar = "phi2") +
  prior(lognormal(log(3), .25), nlpar = "phi3") +
  prior(student_t(3,0,5), dpar="sigma") +
  prior(gamma(2,0.1), class="nu")

fullModel<-brm(bf(y ~ phi1/(1+exp((phi2-time)/phi3)),
                 sigma~s(time,by=treatment), 
                 phi1 + phi2 + phi3 ~ 0+treatment,
                 autocor = ~arma(~time|sample:treatment,1,1), 
                 nl = TRUE), 
              family = student, prior = priorList, data = dfLogistic, iter = 2000, 
              cores = 4, chains = 4, backend = "cmdstanr",
              control = list(adapt_delta = 0.999, max_treedepth = 20),
              inits = function(){list(b_phi1=rgamma(2,1),b_phi2=rgamma(2,1),b_phi3=rgamma(2,1))})


# ********************
# Fit Models 7:25, by 3 with updating priors
# ********************
o=1
for (i in seq(7, 67, by=3)){ #change to 25 #changed from 6,10,4
  if(o!=1){
  eval(parse(text=paste0("previousFit<-fit_",i-3)))
  
  recycledPosteriorText<-paste0('recycledPosterior <- prior(lognormal(',mean(log(as.data.frame(previousFit$fit)$b_phi1_treatmenta)),',',
                                       sd(log(as.data.frame(previousFit$fit)$b_phi1_treatmenta)),'),',
                             'nlpar = "phi1",coef="treatmenta") +',
    'prior(lognormal(',mean(log(as.data.frame(previousFit$fit)$b_phi1_treatmentb)),',',
    sd(log(as.data.frame(previousFit$fit)$b_phi1_treatmentb)),'),',
          'nlpar = "phi1",coef="treatmentb") +',
    'prior(lognormal(',mean(log(as.data.frame(previousFit$fit)$b_phi2_treatmenta)),',',
    sd(log(as.data.frame(previousFit$fit)$b_phi2_treatmenta)),'),',
          'nlpar = "phi2",coef="treatmenta") +',
  'prior(lognormal(',mean(log(as.data.frame(previousFit$fit)$b_phi2_treatmentb)),',',
  sd(log(as.data.frame(previousFit$fit)$b_phi2_treatmentb)),'),',
                    'nlpar = "phi2",coef="treatmentb") +',
    'prior(lognormal(',mean(log(as.data.frame(previousFit$fit)$b_phi3_treatmenta)),',',
  sd(log(as.data.frame(previousFit$fit)$b_phi3_treatmenta)),'),',
          'nlpar = "phi3",coef="treatmenta") +',
  'prior(lognormal(',mean(log(as.data.frame(previousFit$fit)$b_phi3_treatmentb)),',',
  sd(log(as.data.frame(previousFit$fit)$b_phi3_treatmentb)),'),',
                    'nlpar = "phi3",coef="treatmentb") +',
  'prior(student_t(3,0,5), dpar="sigma") +prior(gamma(2,0.1), class="nu")')
  eval(parse(text=recycledPosteriorText))
  }
  
  dataText<-paste0("df_",i,"<-dfLogistic%>%filter(time<=",i,", time>",i-3,")") # from 4
  eval(parse(text=dataText))
  if(o==1){
    sigmaText<-paste0('sigma~s(time, by=treatment, k=3),')# from 4
    brmText<-paste0('fit_',i,' <- brm(bf(y ~ phi1/(1+exp((phi2-time)/phi3)),',
                    sigmaText,
                    'phi1 + phi2 + phi3 ~ 0+treatment, autocor = ~arma(~time|sample:treatment,1,1),',
                    'nl = TRUE),family = student, prior = priorList, data = df_',i,
                    ', iter = 2000, cores = 4, chains = 4, backend = "cmdstanr",',
                    'control = list(adapt_delta = 0.999,max_treedepth = 20),',
                    'inits = function(){list(b_phi1=rgamma(2,1),b_phi2=rgamma(2,1),b_phi3=rgamma(2,1))})')
  } else{
    sigmaText<-paste0('sigma~s(time, by=treatment, k=3),')# from 4
  brmText<-paste0('fit_',i,' <- brm(bf(y ~ phi1/(1+exp((phi2-time)/phi3)),',
                  sigmaText,
                  'phi1 + phi2 + phi3 ~ 0+treatment, autocor = ~arma(~time|sample:treatment,1,1),',
                  'nl = TRUE),family = student, prior = recycledPosterior, data = df_',i,
                  ', iter = 2000, cores = 4, chains = 4, backend = "cmdstanr",',
                  'control = list(adapt_delta = 0.999,max_treedepth = 20),',
                  'inits = function(){list(b_phi1=rgamma(2,1),b_phi2=rgamma(2,1),b_phi3=rgamma(2,1))})')
  }
  print(i)
  eval(parse(text=brmText))
  o=o+1
}

save(df, fit_7, fit_10, fit_13, fit_16, fit_19, fit_22, 
     fit_25, fit_28, fit_31, fit_34, fit_37, fit_40, fit_43, fit_46,
     fit_49, fit_52, fit_55, fit_58, fit_61, fit_64, fit_67, fullModel,
     file="earlyStoppingIterations/updating/modelAndPriorUpdating_DataAndModels_3.rdata")

# How am I comparing these? what matters that could get fucked up?
# Delta HPDI/information density should match with fullModel when built out of updated models
# Parameter estimates and intervals should be trivially different or the same
# Look for parameters changing at the wrong times (compare informationDensity/time for all parameters, should show it)
#     Phi3 is shrinking at later times (lower RGR)
# Build density plots and hpdi line plots for the updating models
# 
#

# ********************
# Read in data
# ********************

print(load("earlyStoppingSim/bayesian_updating/modelAndPriorUpdating_DataAndModels_3.rdata"))

# ********************
# Make dataframe of parameter estimates
# ********************

Parameters<-c("phi1","phi2","phi3")
for (param in Parameters){
  PriorMean=ifelse(grepl("phi1", param, ignore.case=T), 130,
                   ifelse(grepl("phi2", param, ignore.case=T), 12, 3))
  priorCol = rlnorm(3000, log(PriorMean), .25)
  
  for (i in seq(7,67, 3)){
    eval(parse(text=paste0("fit<-fit_",i)))
    t<-rlnorm(3000, mean(log(as.data.frame(fit$fit)$b_phi1_treatmentb)), sd(log(as.data.frame(fit$fit)$b_phi1_treatmentb)))
    
    eval(parse(text=paste0("B<-rlnorm(3000, mean(log(as.data.frame(fit$fit)$b_",
                           str_to_lower(param),
                           "_treatmentb)), sd(log(as.data.frame(fit$fit)$b_",
                           str_to_lower(param),
                           "_treatmentb)))")))
    eval(parse(text=paste0("A<-rlnorm(3000, mean(log(as.data.frame(fit$fit)$b_",
                           str_to_lower(param),
                           "_treatmenta)), sd(log(as.data.frame(fit$fit)$b_",
                           str_to_lower(param),
                           "_treatmenta)))")))
    eval(parse(text=paste0("t",i,param,"<-cbind(A, B)")))
    
    if (i==67){
      eval(parse(text=paste0("df_",param,"<-data.frame(t7=t7",param,",t10=t10",param,",t13=t13",param,
                             ",t16=t16",param,",t19=t19",param,",t22=t22",param,",t25=t25",param,
                             ",t28=t28",param,",t31=t31",param,",t34=t34",param,
                             
                             ",t37=t37",param,",t40=t40",param,",t43=t43",param,
                             ",t46=t46",param,",t49=t49",param,",t52=t52",param,
                             ",t55=t55",param,",t58=t58",param,",t61=t61",param,",t64=t64",param,",t67=t67",param,

                             ", param = '",param,"', prior=priorCol)")))
      
    }
  }
}
df<-rbind(df_phi1, df_phi2, df_phi3)

# ********************
# Plot densities of parameter estimates on each day
# ********************

probs <- seq(1,67,1)
avg_pal <- plasma(n=length(probs),  alpha=0.8)#begin=0.1, end=0.9,
avg_pal2 <- mako(n=length(probs),  alpha=0.8)#begin=0.1, end=0.9,

phi1AUpdatingPlot<-df%>%
  filter(param=="phi1")%>%
  ggplot()+
  geom_density(aes(x=prior), fill="black", alpha=0.55, color="black")+
  lapply(seq(67,7,-3), function(i) geom_density(aes_string(x=paste0("t",i, ".A")), 
                                               fill=avg_pal[i],
                                               color=avg_pal[i]))+
  geom_density(data=as.data.frame(fullModel$fit), aes(x=b_phi1_treatmenta), color="black", linetype=5)+
  annotation_custom(grid::textGrob(label = expression(Phi["1A"]~"= 200"),  
                                 x = unit(0.25, "npc"), 
                                 y = unit(0.9, "npc")))+
  xlab(~~Phi["1A"]~Distribution)+
  ylab("Density")+
  coord_cartesian(xlim=c(75,250), ylim=c(0,.75))+
  theme_minimal()

phi1BUpdatingPlot<-df%>%
  filter(param=="phi1")%>%
  ggplot()+
  geom_density(aes(x=prior), fill="black", alpha=0.55, color="black")+
  lapply(seq(67,7,-3), function(i) geom_density(aes_string(x=paste0("t",i, ".B")), 
                                               fill=avg_pal2[i],
                                               color=avg_pal2[i]))+
  geom_density(data=as.data.frame(fullModel$fit), aes(x=b_phi1_treatmentb), color="black", linetype=5)+
  annotation_custom(grid::textGrob(label = expression(Phi["1B"]~"= 160"),  
                                   x = unit(0.25, "npc"), 
                                   y = unit(0.9, "npc")))+
  xlab(~~Phi["1B"]~Distribution)+
  ylab("Density")+
  coord_cartesian(xlim=c(75,250), ylim=c(0,.75))+
  theme_minimal()


phi2AUpdatingPlot<-df%>%
  filter(param=="phi2")%>%
  ggplot()+
  geom_density(aes(x=prior), fill="black", alpha=0.55, color="black")+
  lapply(seq(67,7,-3), function(i) geom_density(aes_string(x=paste0("t",i, ".A")), 
                                               fill=avg_pal[i],
                                               color=avg_pal[i]))+
  geom_density(data=as.data.frame(fullModel$fit), aes(x=b_phi2_treatmenta), color="black", linetype=5)+
  annotation_custom(grid::textGrob(label = expression(Phi["2A"]~"= 13"),  
                                   x = unit(0.25, "npc"), 
                                   y = unit(0.9, "npc")))+
  xlab(~~Phi["2A"]~Distribution)+
  ylab("Density")+
  coord_cartesian(xlim=c(5,20), ylim=c(0,2.5))+
  theme_minimal()

phi2BUpdatingPlot<-df%>%
  filter(param=="phi2")%>%
  ggplot()+
  geom_density(aes(x=prior), fill="black", alpha=0.55, color="black")+
  lapply(seq(67,7,-3), function(i) geom_density(aes_string(x=paste0("t",i, ".B")), 
                                               fill=avg_pal2[i],
                                               color=avg_pal2[i]))+
  geom_density(data=as.data.frame(fullModel$fit), aes(x=b_phi2_treatmentb), color="black", linetype=5)+
  annotation_custom(grid::textGrob(label = expression(Phi["2B"]~"= 13"),  
                                   x = unit(0.25, "npc"), 
                                   y = unit(0.9, "npc")))+
  xlab(~~Phi["2B"]~Distribution)+
  ylab("Density")+
  coord_cartesian(xlim=c(5,20), ylim=c(0,2.5))+
  theme_minimal()

phi3AUpdatingPlot<-df%>%
  filter(param=="phi3")%>%
  ggplot()+
  geom_density(aes(x=prior), fill="black", alpha=0.55, color="black")+
  lapply(seq(67,7,-3), function(i) geom_density(aes_string(x=paste0("t",i, ".A")), 
                                                fill=avg_pal[i],
                                                color=avg_pal[i]))+
  geom_density(data=as.data.frame(fullModel$fit), aes(x=b_phi3_treatmenta), color="black", linetype=5)+
  annotation_custom(grid::textGrob(label = expression(Phi["3A"]~"= 3"),  
                                   x = unit(0.25, "npc"), 
                                   y = unit(0.9, "npc")))+
  xlab(~~Phi["3A"]~Distribution)+
  ylab("Density")+
  coord_cartesian(xlim=c(1,5), ylim=c(0,9))+
  theme_minimal()

phi3BUpdatingPlot<-df%>%
  filter(param=="phi3")%>%
  ggplot()+
  geom_density(aes(x=prior), fill="black", alpha=0.55, color="black")+
  lapply(seq(67,7,-3), function(i) geom_density(aes_string(x=paste0("t",i, ".B")), 
                                                fill=avg_pal2[i],
                                                color=avg_pal2[i]))+
  geom_density(data=as.data.frame(fullModel$fit), aes(x=b_phi3_treatmentb), color="black", linetype=5)+
  annotation_custom(grid::textGrob(label = expression(Phi["3B"]~"= 3.5"),  
                                   x = unit(0.25, "npc"), 
                                   y = unit(0.9, "npc")))+
  xlab(~~Phi["3B"]~Distribution)+
  ylab("Density")+
  coord_cartesian(xlim=c(1,5), ylim=c(0,9))+
  theme_minimal()


PlasmaCol<-plasma(1, begin=0.65)
MakoCol<-mako(1, begin=0.65)
dataA<-dfLogistic%>%
  filter(treatment=="a")%>%
  ggplot(aes(time,y,group=interaction(treatment,sample)))+
  geom_line(aes(color=time))+
  lapply(seq(7,67,3), 
         function(i) geom_point(
           data = dfLogistic[dfLogistic$time==i & dfLogistic$treatment=="a" ,],
                                             aes(x=i, y=y), color="white", 
           size=1, shape="square"))+
  scale_color_viridis(option="plasma", begin=0.05, end=0.9)+
  ylab(~~Area~(cm^2))+
  xlab("Time")+
  theme_light()+
  coord_cartesian(ylim=c(0,250), xlim=c(0,35))+
  theme(axis.ticks.length=unit(0.2,"cm"))+
  theme(strip.background=element_rect(fill="gray50",color="gray20"),
        strip.text.x=element_text(size=14,color="white"),
        strip.text.y=element_text(size=14,color="white"))+
  theme(axis.title= element_text(size = 10))+
  theme(axis.text = element_text(size = 8))+
  theme(legend.position='none')
dataA
dataB<-dfLogistic%>%
  filter(treatment=="b")%>%
  ggplot(aes(time,y,group=interaction(treatment,sample)))+
  geom_line(aes(color=time))+
  lapply(seq(7,67,3), 
         function(i) geom_point(
           data = dfLogistic[dfLogistic$time==i & dfLogistic$treatment=="b" ,],
           aes(x=i, y=y), color="white", 
           size=1, shape="square"))+
  scale_color_viridis(option="mako", begin=0.05, end=0.9)+
  ylab(~~Area~(cm^2))+
  xlab("Time")+
  theme_light()+
  coord_cartesian(ylim=c(0,250), xlim=c(0,35))+
  theme(axis.ticks.length=unit(0.2,"cm"))+
  theme(strip.background=element_rect(fill="gray50",color="gray20"),
        strip.text.x=element_text(size=14,color="white"),
        strip.text.y=element_text(size=14,color="white"))+
  theme(axis.title= element_text(size = 10))+
  theme(axis.text = element_text(size = 8))+
  theme(legend.position='none')

dataB

patch<- (dataA|phi1AUpdatingPlot|phi2AUpdatingPlot|phi3AUpdatingPlot)/(dataB|phi1BUpdatingPlot|phi2BUpdatingPlot|phi3BUpdatingPlot)
patch
ggsave("modelUpdating_parameterDistributions_withDataRev_67.png", patch,width=16,height=8, dpi=300, bg="#ffffff")



# ********************
# Get HPDIs
# ********************

test<-df%>%
  group_by(param)%>%
  summarize(
    across(
      .cols = contains("t"),
      .fns = ~HPDI(., prob=0.95),
      .names = "{.col}"
    ))%>%ungroup()%>%
  t()%>%
  as.data.frame()%>%
  rownames_to_column(var='param')%>%
  rename(Phi1_lower=V1, Phi1_upper=V2,
         Phi2_lower=V3, Phi2_upper=V4,
         Phi3_lower=V5, Phi3_upper=V6)%>%
  filter(param!="param")%>%
  mutate(treatment = substr(param, nchar(param), nchar(param)),
         time = str_remove_all(param, "\\."),
         time = str_remove_all(time, "A|B"),
         time= str_remove_all(time, "t"))%>%
  select(-param)

head(test)

colname<-colnames(figureData_phi1A)
hpdis<-rbind(
  data.frame(do.call(rbind, as.data.frame(lapply(seq(1:ncol(figureData_phi1A)), 
                                                 function(i) HPDI(figureData_phi1A[,i], 
                                                                  prob=0.95)))), variable="phi1", treatment="A",
             row.names=colname)%>%rename(lower=X1, upper=X2)%>%rownames_to_column(var="time")%>%rowid_to_column(var="t")%>%
    mutate(length=upper-lower, 
           diff=lag(length,default = NA)-length,
           diffPct=100*(1-(length/lag(length)))),
  data.frame(do.call(rbind, as.data.frame(lapply(seq(1:ncol(figureData_phi1B)), 
                                                 function(i) HPDI(figureData_phi1B[,i], 
                                                                  prob=0.95)))), variable="phi1", treatment="B",
             row.names=colname)%>%rename(lower=X1, upper=X2)%>%rownames_to_column(var="time")%>%rowid_to_column(var="t")%>%
    mutate(length=upper-lower, 
           diff=lag(length,default = NA)-length,
           diffPct=100*(1-(length/lag(length)))),
  data.frame(do.call(rbind, as.data.frame(lapply(seq(1:ncol(figureData_phi2B)), 
                                                 function(i) HPDI(figureData_phi2B[,i], 
                                                                  prob=0.95)))), variable="phi2", treatment="B",
             row.names=colname)%>%rename(lower=X1, upper=X2)%>%rownames_to_column(var="time")%>%rowid_to_column(var="t")%>%
    mutate(length=upper-lower, 
           diff=lag(length,default = NA)-length,
           diffPct=100*(1-(length/lag(length)))),
  data.frame(do.call(rbind, as.data.frame(lapply(seq(1:ncol(figureData_phi2A)), 
                                                 function(i) HPDI(figureData_phi2A[,i], 
                                                                  prob=0.95)))), variable="phi2", treatment="A",
             row.names=colname)%>%rename(lower=X1, upper=X2)%>%rownames_to_column(var="time")%>%rowid_to_column(var="t")%>%
    mutate(length=upper-lower, 
           diff=lag(length,default = NA)-length,
           diffPct=100*(1-(length/lag(length)))),
  data.frame(do.call(rbind, as.data.frame(lapply(seq(1:ncol(figureData_phi3B)), 
                                                 function(i) HPDI(figureData_phi3B[,i], 
                                                                  prob=0.95)))), variable="phi3", treatment="B",
             row.names=colname)%>%rename(lower=X1, upper=X2)%>%rownames_to_column(var="time")%>%rowid_to_column(var="t")%>%
    mutate(length=upper-lower, 
           diff=lag(length,default = NA)-length,
           diffPct=100*(1-(length/lag(length)))),
  data.frame(do.call(rbind, as.data.frame(lapply(seq(1:ncol(figureData_phi3A)), 
                                                 function(i) HPDI(figureData_phi3A[,i], 
                                                                  prob=0.95)))), variable="phi3", treatment="A",
             row.names=colname)%>%rename(lower=X1, upper=X2)%>%rownames_to_column(var="time")%>%rowid_to_column(var="t")%>%
    mutate(length=upper-lower, 
           diff=lag(length,default = NA)-length,
           diffPct=100*(1-(length/lag(length)))))

hpdis<-hpdis%>%
  mutate(diff=ifelse(is.na(diff), 0, diff),
         t = ifelse(time=="prior", "t0", time),
         t = as.numeric(str_remove_all(t, "t")))






