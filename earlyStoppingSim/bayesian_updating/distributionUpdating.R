library(brms)
library(tidyverse)
library(viridis)
library(patchwork)
library(rethinking)
library(scales)
# growthSim <- function(x,phi1,phi2,phi3){
#   phi1_r <- phi1+rnorm(1,mean = 0,sd=25)
#   phi2_r <- phi2+rnorm(1,mean=0,sd=1)
#   phi3_r <- phi3+rnorm(1,mean=0,sd=.2)
#   return(phi1_r / (1 + exp( (phi2_r-x)/phi3_r)))
# }
# set.seed(123)
# df <- rbind(
#     do.call(rbind, lapply(1:20, function(i) data.frame("sample"=paste0("sample_",i),
#                                                        "treatment"="a","time"=1:25,
#                                                        "y"=growthSim(1:25, 200, 13, 3),
#                                                        stringsAsFactors = F))),
#     do.call(rbind, lapply(1:20,function(i) data.frame("sample"=paste0("sample_",i),
#                                                       "treatment"="b","time"=1:25,
#                                                       "y"=growthSim(1:25, 160, 13, 3.5),
#                                                       stringsAsFactors = F))))
# priorList <- prior(lognormal(log(130), .25),nlpar = "phi1") +
#   prior(lognormal(log(12), .25), nlpar = "phi2") +
#   prior(lognormal(log(3), .25), nlpar = "phi3") +
#   prior(student_t(3,0,5), dpar="sigma") +
#   prior(gamma(2,0.1), class="nu")
# 
# flatPriors<-set_prior("uniform(0,1000)", lb = 0, ub = 1000, nlpar = "phi1") +
#   set_prior("uniform(0,1000)", lb = 0, ub = 1000, nlpar = "phi2") +
#   set_prior("uniform(0,1000)", lb = 0, ub = 1000, nlpar = "phi3") +
#   prior(student_t(3,0,5), dpar="sigma") +
#   prior(gamma(2,0.1), class="nu")
# 
# wrongPriors<- prior(lognormal(log(50), .15),nlpar = "phi1") +
#    prior(lognormal(log(2), .1), nlpar = "phi2") +
#    prior(student_t(3, 30, 5), nlpar = "phi3") +
#    prior(student_t(3,0,5), dpar="sigma") +
#    prior(gamma(2,0.1), class="nu")
# plot(density(rstudent_t(1000,3, 30, 5)))
# for (i in seq(3,5,by=2)){
#   dataText<-paste0("df_",i,"<-df%>%filter(time<=",i,")")
#   eval(parse(text=dataText))
#   if(i<=9){
#     sigmaText<-paste0('sigma~s(time, by=treatment, k=',i,'),')
#   } else{
#     sigmaText<-paste0('sigma~s(time, by=treatment),')}
#   brmText<-paste0('fit_',i,' <- brm(bf(y ~ phi1/(1+exp((phi2-time)/phi3)),',
#                   sigmaText,
#                   'phi1 + phi2 + phi3 ~ 0+treatment, autocor = ~arma(~time|sample:treatment,1,1),',
#                   'nl = TRUE),family = student, prior = priorList, data = df_',i,
#                   ', iter = 2000, cores = 4, chains = 4, backend = "cmdstanr",',
#                   'control = list(adapt_delta = 0.999,max_treedepth = 20)',
#                   'inits = function(){list(b_phi1=rgamma(2,1),b_phi2=rgamma(2,1),b_phi3=rgamma(2,1))})')
#   print(i)
#   eval(parse(text=brmText))
# }
# save(df, fit_3, fit_5, fit_7, fit_9, fit_11, fit_13, fit_15, fit_17, fit_19, fit_21, fit_23, fit_25,
#      file="earlyStoppingIterations/updating/parameterUpdating_DataAndModels_wrongPriors.rdata")


#print(load("earlyStoppingSim/bayesian_updating/parameterUpdating_DataAndModels.rdata"))
# 

print(load("earlyStoppingSim/bayesian_updating/parameterUpdating_DataAndModels_3to25.rdata"))



# make a figure to show how all the parameters change at each timepoint
# Want to have 8-12 panels?
# run models on days 7:25 and make a panel for each?
# I want to show phi1, 2, and 3. I'll start with just the phi1 plots.


# make a version with all the plots together, viridis to show days
set.seed(123)
figureData_phi1B<-data.frame(
  prior = rlnorm(3000, log(130), .25),
  t3 = rlnorm(3000,
              mean(log(as.data.frame(fit_3$fit)$b_phi1_treatmentb)),
              sd(log(as.data.frame(fit_3$fit)$b_phi1_treatmentb))),
  t5 = rlnorm(3000,
              mean(log(as.data.frame(fit_3$fit)$b_phi1_treatmentb)),
              sd(log(as.data.frame(fit_3$fit)$b_phi1_treatmentb))),
  t7 = rlnorm(3000,
              mean(log(as.data.frame(fit_7$fit)$b_phi1_treatmentb)),
              sd(log(as.data.frame(fit_7$fit)$b_phi1_treatmentb))),
  t9 = rlnorm(3000,
              mean(log(as.data.frame(fit_9$fit)$b_phi1_treatmentb)),
              sd(log(as.data.frame(fit_9$fit)$b_phi1_treatmentb))),
  t11 = rlnorm(3000,
               mean(log(as.data.frame(fit_11$fit)$b_phi1_treatmentb)),
               sd(log(as.data.frame(fit_11$fit)$b_phi1_treatmentb))),
  t13 = rlnorm(3000,
               mean(log(as.data.frame(fit_13$fit)$b_phi1_treatmentb)),
               sd(log(as.data.frame(fit_13$fit)$b_phi1_treatmentb))),
  t15 = rlnorm(3000,
               mean(log(as.data.frame(fit_15$fit)$b_phi1_treatmentb)),
               sd(log(as.data.frame(fit_15$fit)$b_phi1_treatmentb))),
  t17 = rlnorm(3000,
               mean(log(as.data.frame(fit_17$fit)$b_phi1_treatmentb)),
               sd(log(as.data.frame(fit_17$fit)$b_phi1_treatmentb))),
  t19 = rlnorm(3000,
               mean(log(as.data.frame(fit_19$fit)$b_phi1_treatmentb)),
               sd(log(as.data.frame(fit_19$fit)$b_phi1_treatmentb))),
  t21 = rlnorm(3000,
               mean(log(as.data.frame(fit_21$fit)$b_phi1_treatmentb)),
               sd(log(as.data.frame(fit_21$fit)$b_phi1_treatmentb))),
  t23 = rlnorm(3000,
               mean(log(as.data.frame(fit_23$fit)$b_phi1_treatmentb)),
               sd(log(as.data.frame(fit_23$fit)$b_phi1_treatmentb))),
  t25 = rlnorm(3000,
               mean(log(as.data.frame(fit_25$fit)$b_phi1_treatmentb)),
               sd(log(as.data.frame(fit_25$fit)$b_phi1_treatmentb))))

figureData_phi1A<-data.frame(
  prior = rlnorm(3000, log(130), .25),
  t3 = rlnorm(3000,
              mean(log(as.data.frame(fit_3$fit)$b_phi1_treatmenta)),
              sd(log(as.data.frame(fit_3$fit)$b_phi1_treatmenta))),
  t5 = rlnorm(3000,
              mean(log(as.data.frame(fit_3$fit)$b_phi1_treatmenta)),
              sd(log(as.data.frame(fit_3$fit)$b_phi1_treatmenta))),
  t7 = rlnorm(3000,
              mean(log(as.data.frame(fit_7$fit)$b_phi1_treatmenta)),
              sd(log(as.data.frame(fit_7$fit)$b_phi1_treatmenta))),
  t9 = rlnorm(3000,
              mean(log(as.data.frame(fit_9$fit)$b_phi1_treatmenta)),
              sd(log(as.data.frame(fit_9$fit)$b_phi1_treatmenta))),
  t11 = rlnorm(3000,
               mean(log(as.data.frame(fit_11$fit)$b_phi1_treatmenta)),
               sd(log(as.data.frame(fit_11$fit)$b_phi1_treatmenta))),
  t13 = rlnorm(3000,
               mean(log(as.data.frame(fit_13$fit)$b_phi1_treatmenta)),
               sd(log(as.data.frame(fit_13$fit)$b_phi1_treatmenta))),
  t15 = rlnorm(3000,
               mean(log(as.data.frame(fit_15$fit)$b_phi1_treatmenta)),
               sd(log(as.data.frame(fit_15$fit)$b_phi1_treatmenta))),
  t17 = rlnorm(3000,
               mean(log(as.data.frame(fit_17$fit)$b_phi1_treatmenta)),
               sd(log(as.data.frame(fit_17$fit)$b_phi1_treatmenta))),
  t19 = rlnorm(3000,
               mean(log(as.data.frame(fit_19$fit)$b_phi1_treatmenta)),
               sd(log(as.data.frame(fit_19$fit)$b_phi1_treatmenta))),
  t21 = rlnorm(3000,
               mean(log(as.data.frame(fit_21$fit)$b_phi1_treatmenta)),
               sd(log(as.data.frame(fit_21$fit)$b_phi1_treatmenta))),
  t23 = rlnorm(3000,
               mean(log(as.data.frame(fit_23$fit)$b_phi1_treatmenta)),
               sd(log(as.data.frame(fit_23$fit)$b_phi1_treatmenta))),
  t25 = rlnorm(3000,
               mean(log(as.data.frame(fit_25$fit)$b_phi1_treatmenta)),
               sd(log(as.data.frame(fit_25$fit)$b_phi1_treatmenta))))

figureData_phi2A<-data.frame(
  prior = rlnorm(3000, log(12), .25),
  t3 = rlnorm(3000,
              mean(log(as.data.frame(fit_3$fit)$b_phi2_treatmenta)),
              sd(log(as.data.frame(fit_3$fit)$b_phi2_treatmenta))),
  t5 = rlnorm(3000,
              mean(log(as.data.frame(fit_3$fit)$b_phi2_treatmenta)),
              sd(log(as.data.frame(fit_3$fit)$b_phi2_treatmenta))),
  t7 = rlnorm(3000,
              mean(log(as.data.frame(fit_7$fit)$b_phi2_treatmenta)),
              sd(log(as.data.frame(fit_7$fit)$b_phi2_treatmenta))),
  t9 = rlnorm(3000,
              mean(log(as.data.frame(fit_9$fit)$b_phi2_treatmenta)),
              sd(log(as.data.frame(fit_9$fit)$b_phi2_treatmenta))),
  t11 = rlnorm(3000,
               mean(log(as.data.frame(fit_11$fit)$b_phi2_treatmenta)),
               sd(log(as.data.frame(fit_11$fit)$b_phi2_treatmenta))),
  t13 = rlnorm(3000,
               mean(log(as.data.frame(fit_13$fit)$b_phi2_treatmenta)),
               sd(log(as.data.frame(fit_13$fit)$b_phi2_treatmenta))),
  t15 = rlnorm(3000,
               mean(log(as.data.frame(fit_15$fit)$b_phi2_treatmenta)),
               sd(log(as.data.frame(fit_15$fit)$b_phi2_treatmenta))),
  t17 = rlnorm(3000,
               mean(log(as.data.frame(fit_17$fit)$b_phi2_treatmenta)),
               sd(log(as.data.frame(fit_17$fit)$b_phi2_treatmenta))),
  t19 = rlnorm(3000,
               mean(log(as.data.frame(fit_19$fit)$b_phi2_treatmenta)),
               sd(log(as.data.frame(fit_19$fit)$b_phi2_treatmenta))),
  t21 = rlnorm(3000,
               mean(log(as.data.frame(fit_21$fit)$b_phi2_treatmenta)),
               sd(log(as.data.frame(fit_21$fit)$b_phi2_treatmenta))),
  t23 = rlnorm(3000,
               mean(log(as.data.frame(fit_23$fit)$b_phi2_treatmenta)),
               sd(log(as.data.frame(fit_23$fit)$b_phi2_treatmenta))),
  t25 = rlnorm(3000,
               mean(log(as.data.frame(fit_25$fit)$b_phi2_treatmenta)),
               sd(log(as.data.frame(fit_25$fit)$b_phi2_treatmenta))))

figureData_phi2B<-data.frame(
  prior = rlnorm(3000, log(12), .25),
  t3 = rlnorm(3000,
              mean(log(as.data.frame(fit_3$fit)$b_phi2_treatmentb)),
              sd(log(as.data.frame(fit_3$fit)$b_phi2_treatmentb))),
  t5 = rlnorm(3000,
              mean(log(as.data.frame(fit_3$fit)$b_phi2_treatmentb)),
              sd(log(as.data.frame(fit_3$fit)$b_phi2_treatmentb))),
  t7 = rlnorm(3000,
              mean(log(as.data.frame(fit_7$fit)$b_phi2_treatmentb)),
              sd(log(as.data.frame(fit_7$fit)$b_phi2_treatmentb))),
  t9 = rlnorm(3000,
              mean(log(as.data.frame(fit_9$fit)$b_phi2_treatmentb)),
              sd(log(as.data.frame(fit_9$fit)$b_phi2_treatmentb))),
  t11 = rlnorm(3000,
               mean(log(as.data.frame(fit_11$fit)$b_phi2_treatmentb)),
               sd(log(as.data.frame(fit_11$fit)$b_phi2_treatmentb))),
  t13 = rlnorm(3000,
               mean(log(as.data.frame(fit_13$fit)$b_phi2_treatmentb)),
               sd(log(as.data.frame(fit_13$fit)$b_phi2_treatmentb))),
  t15 = rlnorm(3000,
               mean(log(as.data.frame(fit_15$fit)$b_phi2_treatmentb)),
               sd(log(as.data.frame(fit_15$fit)$b_phi2_treatmentb))),
  t17 = rlnorm(3000,
               mean(log(as.data.frame(fit_17$fit)$b_phi2_treatmentb)),
               sd(log(as.data.frame(fit_17$fit)$b_phi2_treatmentb))),
  t19 = rlnorm(3000,
               mean(log(as.data.frame(fit_19$fit)$b_phi2_treatmentb)),
               sd(log(as.data.frame(fit_19$fit)$b_phi2_treatmentb))),
  t21 = rlnorm(3000,
               mean(log(as.data.frame(fit_21$fit)$b_phi2_treatmentb)),
               sd(log(as.data.frame(fit_21$fit)$b_phi2_treatmentb))),
  t23 = rlnorm(3000,
               mean(log(as.data.frame(fit_23$fit)$b_phi2_treatmentb)),
               sd(log(as.data.frame(fit_23$fit)$b_phi2_treatmentb))),
  t25 = rlnorm(3000,
               mean(log(as.data.frame(fit_25$fit)$b_phi2_treatmentb)),
               sd(log(as.data.frame(fit_25$fit)$b_phi2_treatmentb))))

figureData_phi3A<-data.frame(
  prior = rlnorm(3000, log(3), .25),
  t3 = rlnorm(3000,
              mean(log(as.data.frame(fit_3$fit)$b_phi3_treatmenta)),
              sd(log(as.data.frame(fit_3$fit)$b_phi3_treatmenta))),
  t5 = rlnorm(3000,
              mean(log(as.data.frame(fit_3$fit)$b_phi3_treatmenta)),
              sd(log(as.data.frame(fit_3$fit)$b_phi3_treatmenta))),
  t7 = rlnorm(3000,
              mean(log(as.data.frame(fit_7$fit)$b_phi3_treatmenta)),
              sd(log(as.data.frame(fit_7$fit)$b_phi3_treatmenta))),
  t9 = rlnorm(3000,
              mean(log(as.data.frame(fit_9$fit)$b_phi3_treatmenta)),
              sd(log(as.data.frame(fit_9$fit)$b_phi3_treatmenta))),
  t11 = rlnorm(3000,
               mean(log(as.data.frame(fit_11$fit)$b_phi3_treatmenta)),
               sd(log(as.data.frame(fit_11$fit)$b_phi3_treatmenta))),
  t13 = rlnorm(3000,
               mean(log(as.data.frame(fit_13$fit)$b_phi3_treatmenta)),
               sd(log(as.data.frame(fit_13$fit)$b_phi3_treatmenta))),
  t15 = rlnorm(3000,
               mean(log(as.data.frame(fit_15$fit)$b_phi3_treatmenta)),
               sd(log(as.data.frame(fit_15$fit)$b_phi3_treatmenta))),
  t17 = rlnorm(3000,
               mean(log(as.data.frame(fit_17$fit)$b_phi3_treatmenta)),
               sd(log(as.data.frame(fit_17$fit)$b_phi3_treatmenta))),
  t19 = rlnorm(3000,
               mean(log(as.data.frame(fit_19$fit)$b_phi3_treatmenta)),
               sd(log(as.data.frame(fit_19$fit)$b_phi3_treatmenta))),
  t21 = rlnorm(3000,
               mean(log(as.data.frame(fit_21$fit)$b_phi3_treatmenta)),
               sd(log(as.data.frame(fit_21$fit)$b_phi3_treatmenta))),
  t23 = rlnorm(3000,
               mean(log(as.data.frame(fit_23$fit)$b_phi3_treatmenta)),
               sd(log(as.data.frame(fit_23$fit)$b_phi3_treatmenta))),
  t25 = rlnorm(3000,
               mean(log(as.data.frame(fit_25$fit)$b_phi3_treatmenta)),
               sd(log(as.data.frame(fit_25$fit)$b_phi3_treatmenta))))

figureData_phi3B<-data.frame(
  prior = rlnorm(3000, log(3), .25),
  t3 = rlnorm(3000,
              mean(log(as.data.frame(fit_3$fit)$b_phi3_treatmentb)),
              sd(log(as.data.frame(fit_3$fit)$b_phi3_treatmentb))),
  t5 = rlnorm(3000,
              mean(log(as.data.frame(fit_3$fit)$b_phi3_treatmentb)),
              sd(log(as.data.frame(fit_3$fit)$b_phi3_treatmentb))),
  t7 = rlnorm(3000,
              mean(log(as.data.frame(fit_7$fit)$b_phi3_treatmentb)),
              sd(log(as.data.frame(fit_7$fit)$b_phi3_treatmentb))),
  t9 = rlnorm(3000,
              mean(log(as.data.frame(fit_9$fit)$b_phi3_treatmentb)),
              sd(log(as.data.frame(fit_9$fit)$b_phi3_treatmentb))),
  t11 = rlnorm(3000,
               mean(log(as.data.frame(fit_11$fit)$b_phi3_treatmentb)),
               sd(log(as.data.frame(fit_11$fit)$b_phi3_treatmentb))),
  t13 = rlnorm(3000,
               mean(log(as.data.frame(fit_13$fit)$b_phi3_treatmentb)),
               sd(log(as.data.frame(fit_13$fit)$b_phi3_treatmentb))),
  t15 = rlnorm(3000,
               mean(log(as.data.frame(fit_15$fit)$b_phi3_treatmentb)),
               sd(log(as.data.frame(fit_15$fit)$b_phi3_treatmentb))),
  t17 = rlnorm(3000,
               mean(log(as.data.frame(fit_17$fit)$b_phi3_treatmentb)),
               sd(log(as.data.frame(fit_17$fit)$b_phi3_treatmentb))),
  t19 = rlnorm(3000,
               mean(log(as.data.frame(fit_19$fit)$b_phi3_treatmentb)),
               sd(log(as.data.frame(fit_19$fit)$b_phi3_treatmentb))),
  t21 = rlnorm(3000,
               mean(log(as.data.frame(fit_21$fit)$b_phi3_treatmentb)),
               sd(log(as.data.frame(fit_21$fit)$b_phi3_treatmentb))),
  t23 = rlnorm(3000,
               mean(log(as.data.frame(fit_23$fit)$b_phi3_treatmentb)),
               sd(log(as.data.frame(fit_23$fit)$b_phi3_treatmentb))),
  t25 = rlnorm(3000,
               mean(log(as.data.frame(fit_25$fit)$b_phi3_treatmentb)),
               sd(log(as.data.frame(fit_25$fit)$b_phi3_treatmentb))))


## Example panel, show comparison of two timepoints
figDat<-figureData_phi3B
colnames(figDat)<-paste0("t",c(1:ncol(figDat)))

for (i in 1:ncol(figDat)){
  if(i > 1){
eval(parse(text=paste0("hdpiPrev<-HPDI(figDat$t",i-1,", prob=0.95)")))
eval(parse(text=paste0("hdpi<-HPDI(figDat$t",i,", prob=0.95)")))
Plot<-ggplot(figDat)+
  geom_density(aes_string(x=paste0("t",i-1)), fill="lightgray", linetype=5, color="gray20")+
  geom_density(aes_string(x=paste0("t",i)), fill="navy", alpha=0.5)+
  geom_segment(aes(x= hdpi[1], xend=hdpi[2], y=9.8, yend=9.8))+
  lapply(c(1:2), function(i) annotate("point",x= hdpi[[i]], y=9.8, size=2))+
  annotate("point",x= eval(parse(text=paste0("mean(figDat$t",i,")"))), y=9.8, size=2, color="red")+
  geom_segment(aes(x= hdpiPrev[1], xend=hdpiPrev[2], y=9.6, yend=9.6),
               linetype=5, color="gray50")+
  lapply(c(1:2), function(i) annotate("point",x= hdpiPrev[[i]], y=9.6, size=2, color="gray50"))+
  annotate("point",x= eval(parse(text=paste0("mean(figDat$t",i-1,")"))), y=9.6, size=1, color="darkred")+
  scale_x_continuous(breaks=seq(1,6,1))+
  coord_cartesian(ylim=c(0,10), xlim=c(1,6))+
  labs(x=expression(Phi['3B']~Estimate), y="", title=paste0("Day ",
       str_remove_all(colnames(figureData_phi3B)[i], "t")
       ))+
  theme_light()+
  theme(axis.text.y = element_blank())
print(Plot)
}
}

2+2

probs <- seq(1,25,1)
avg_pal <- plasma(n=length(probs), begin=0.1, end=0.9, alpha=0.8)
avg_pal2 <- mako(n=length(probs), begin=0.1, end=0.9, alpha=0.8)
phi1B_vir<-ggplot(figureData_phi1B)+
  #geom_density(aes(x=prior), fill="black", alpha=0.5)+
  lapply(seq(25,3,-2), function(i) geom_density(aes_string(x=paste0("t",i)), 
                                               fill=avg_pal2[i],
                                               color=avg_pal2[i]))+
  geom_density(aes(x=prior), fill="black", alpha=0.55, color="black")+
  annotation_custom(grid::textGrob(label = expression(Phi["1B"]~"= 160"),  
                                   x = unit(0.25, "npc"), 
                                   y = unit(0.9, "npc")))+
  xlab(~~Phi["1B"]~Distribution)+
  ylab("Density")+
  coord_cartesian(xlim=c(75,250), ylim=c(0,.2))+
  theme_minimal()
phi1B_vir
phi1A_vir<-ggplot(figureData_phi1A)+
  #geom_density(aes(x=prior), fill="black", alpha=0.5)+
  lapply(seq(25,3,-2), function(i) geom_density(aes_string(x=paste0("t",i)), 
                                                fill=avg_pal[i],
                                                color=avg_pal[i]))+
  geom_density(aes(x=prior), fill="black", alpha=0.55, color="black")+
  annotation_custom(grid::textGrob(label = expression(Phi["1A"]~"= 200"),  
                                   x = unit(0.25, "npc"), 
                                   y = unit(0.9, "npc")))+
  xlab(~~Phi["1A"]~Distribution)+
  ylab("Density")+
  coord_cartesian(xlim=c(75,250), ylim=c(0,.2))+
  theme_minimal()
phi1A_vir
phi2B_vir<-ggplot(figureData_phi2B)+
  #geom_density(aes(x=prior), fill="black", alpha=0.5)+
  lapply(seq(25,3,-2), function(i) geom_density(aes_string(x=paste0("t",i)), 
                                                fill=avg_pal2[i],
                                                color=avg_pal2[i]))+
  geom_density(aes(x=prior), fill="black", alpha=0.55, color="black")+
  annotation_custom(grid::textGrob(label = expression(Phi["2B"]~"= 13"),  
                                   x = unit(0.25, "npc"), 
                                   y = unit(0.9, "npc")))+
  xlab(~~Phi["2B"]~Distribution)+
  ylab("Density")+
  coord_cartesian(xlim=c(5,20), ylim=c(0,4.25))+
  theme_minimal()

phi2A_vir<-ggplot(figureData_phi2A)+
  #geom_density(aes(x=prior), fill="black", alpha=0.5)+
  lapply(seq(25,3,-2), function(i) geom_density(aes_string(x=paste0("t",i)), 
                                                fill=avg_pal[i],
                                                color=avg_pal[i]))+
  geom_density(aes(x=prior), fill="black", alpha=0.55, color="black")+
  annotation_custom(grid::textGrob(label = expression(Phi["2A"]~"= 13"),  
                                   x = unit(0.25, "npc"), 
                                   y = unit(0.9, "npc")))+
  xlab(~~Phi["2A"]~Distribution)+
  ylab("Density")+
  coord_cartesian(xlim=c(5,20), ylim=c(0,4.25))+
  theme_minimal()

phi3B_vir<-ggplot(figureData_phi3B)+
  #geom_density(aes(x=prior), fill="black", alpha=0.5)+
  lapply(seq(25,3,-2), function(i) geom_density(aes_string(x=paste0("t",i)), 
                                                fill=avg_pal2[i],
                                                color=avg_pal2[i]))+
  geom_density(aes(x=prior), fill="black", alpha=0.55, color="black")+
  annotation_custom(grid::textGrob(label = expression(Phi["3B"]~"= 3.5"),  
                                   x = unit(0.25, "npc"), 
                                   y = unit(0.9, "npc")))+
  xlab(~~Phi["3B"]~Distribution)+
  ylab("Density")+
  coord_cartesian(xlim=c(1,5), ylim=c(0,9))+
  theme_minimal()

phi3A_vir<-ggplot(figureData_phi3A)+
  #geom_density(aes(x=prior), fill="black", alpha=0.5)+
  lapply(seq(25,3,-2), function(i) geom_density(aes_string(x=paste0("t",i)), 
                                                fill=avg_pal[i],
                                                color=avg_pal[i]))+
  geom_density(aes(x=prior), fill="black", alpha=0.55, color="black")+
  annotation_custom(grid::textGrob(label = expression(Phi["3A"]~"= 3"),  
                                   x = unit(0.25, "npc"), 
                                   y = unit(0.9, "npc")))+
  xlab(~~Phi["3A"]~Distribution)+
  ylab("Density")+
  coord_cartesian(xlim=c(1,5), ylim=c(0,9))+
  theme_minimal()

patch<- (phi1A_vir+ phi2A_vir +phi3A_vir)/(phi1B_vir+ phi2B_vir +phi3B_vir)
patch
ggsave("logisticPosteriorUpdating_2tone.png",patch,width=12,height=8, dpi=300, bg="#ffffff")

head(df) #plot df, add to patch...
PlasmaCol<-plasma(1, begin=0.65)
MakoCol<-mako(1, begin=0.65)
dataA<-df%>%
  filter(treatment=="a")%>%
  ggplot(aes(time,y,group=interaction(treatment,sample)))+
  geom_line(aes(color=time))+
  scale_color_viridis(option="plasma", begin=0.05, end=0.9)+
  ylab(~~Area~(cm^2))+
  xlab("Time")+
  theme_light()+
  coord_cartesian(ylim=c(0,250))+
  theme(axis.ticks.length=unit(0.2,"cm"))+
  theme(strip.background=element_rect(fill="gray50",color="gray20"),
        strip.text.x=element_text(size=14,color="white"),
        strip.text.y=element_text(size=14,color="white"))+
  theme(axis.title= element_text(size = 10))+
  theme(axis.text = element_text(size = 8))+
  theme(legend.position='none')
dataA
dataB<-df%>%
  filter(treatment=="b")%>%
  ggplot(aes(time,y,group=interaction(treatment,sample)))+
  geom_line(aes(color=time))+
  scale_color_viridis(option="mako", begin=0.05, end=0.9)+
  ylab(~~Area~(cm^2))+
  xlab("Time")+
  theme_light()+
  coord_cartesian(ylim=c(0,250))+
  theme(axis.ticks.length=unit(0.2,"cm"))+
  theme(strip.background=element_rect(fill="gray50",color="gray20"),
        strip.text.x=element_text(size=14,color="white"),
        strip.text.y=element_text(size=14,color="white"))+
  theme(axis.title= element_text(size = 10))+
  theme(axis.text = element_text(size = 8))+
  theme(legend.position='none')

dataB

patch<- (dataA|phi1A_vir|phi2A_vir|phi3A_vir)/(dataB|phi1B_vir|phi2B_vir|phi3B_vir)
patch
ggsave("logisticPosteriorUpdating_withDataRev.png",patch,width=16,height=8, dpi=300, bg="#ffffff")

dataB+
  geom_vline(xintercept=7)


# separate figure, want to show the rate of change in hpdi length by day for each parameter...

library(rethinking)
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

save(hpdis, file="HPDIS_3_to_25.rdata")
densities_df<-rbind(figureData_phi1A%>%mutate(param = "phi1A"),
                    figureData_phi1B%>%mutate(param = "phi1B"),
                    figureData_phi2A%>%mutate(param = "phi2A"),
                    figureData_phi2B%>%mutate(param = "phi2B"),
                    figureData_phi3A%>%mutate(param = "phi3A"),
                    figureData_phi3B%>%mutate(param = "phi3B"))
save(densities_df, file="densities_df_3_to_25.rdata")

labStrings <- as_labeller(c('A' = "Treatment A",
                            "B" = "Treatment B",
                            "phi1" = "~~Phi[1]",
                            "phi2" = "Phi 2",
                            "phi3" = "Phi 3"
                            ))
PlasmaCol<-plasma(1, begin=0.75)
MakoCol<-mako(1, begin=0.75)


hpdis$trtVar<-paste0(hpdis$treatment, hpdis$variable)
hpdis$TreatmentVar<-factor(hpdis$trtVar, labels = c("Phi[1]~Treatment~A", "Phi[2]~Treatment~A", "Phi[3]~Treatment~A",
                                                    "Phi[1]~Treatment~B", "Phi[2]~Treatment~B", "Phi[3]~Treatment~B"))

hpdiDeltaLPlot<-ggplot(hpdis)+
  facet_wrap(~TreatmentVar, labeller=label_parsed)+ 
  geom_smooth(aes(x=t, y=diff, color=treatment), se=F)+
  scale_x_continuous(breaks=seq(0,25, 5))+
  scale_color_manual(values = c(PlasmaCol, MakoCol))+
  labs(x="Time")+
  ylab("Change in HPDI Length")+
  theme_minimal()+
  theme(strip.background = element_rect(colour="white",
                                        fill="white"),
        plot.background=element_rect(fill="white", color="white"),
        axis.title.x=element_text(size=14,color="black"),
        axis.title.y=element_text(size=14,color="black"),
        panel.background = element_blank(),
        legend.position = "none")
hpdiDeltaLPlot

hpdiDeltaPctPlot<-ggplot(hpdis)+
  facet_wrap(~TreatmentVar, labeller=label_parsed)+ 
  geom_smooth(aes(x=t, y=diffPct/100, color=treatment), se=F)+
  scale_x_continuous(breaks=seq(0,25, 5))+
  scale_color_manual(values = c(PlasmaCol, MakoCol))+
  scale_y_continuous(labels=percent_format())+
  labs(x="")+
  ylab("Change in HPDI Length (%)")+
  theme_minimal()+
  theme(strip.background = element_rect(colour="white",
                                        fill="white"),
        plot.background=element_rect(fill="white", color="white"),
        axis.title.x=element_text(size=14,color="black"),
        axis.title.y=element_text(size=14,color="black"),
        axis.text.x.bottom = element_blank(),
        panel.background = element_blank(),
        legend.position = "none")
hpdiDeltaPctPlot

patch<-hpdiDeltaPctPlot/hpdiDeltaLPlot
ggsave("DeltaHPDI.png",patch,width=6,height=10, dpi=300, bg="#ffffff")


hpdis$Treatment<-factor(hpdis$treatment, labels = c("Phi[1]~Treatment~A", "Phi[1]~Treatment~B"))
p1Line<-hpdis%>%
  filter(variable=="phi1")%>%
  ggplot()+
  facet_wrap(~Treatment, labeller = label_parsed)+
  geom_line(aes(x=t, y=length, color=treatment), size=1, linetype=5)+
  geom_smooth(aes(x=t, y=length, color=treatment), se=F, size=1)+
  scale_color_manual(values = c(PlasmaCol, MakoCol))+
  scale_x_continuous(breaks=seq(0,25, 5))+
  scale_y_continuous(breaks=seq(0,150, 20))+
  labs(x=" ")+
  ylab(" ")+
  theme_minimal()+
  theme(strip.background = element_rect(colour="white",
                                        fill="white"),
        plot.background=element_rect(fill="white", color="white"),
        axis.title.x=element_text(size=14,color="black"),
        axis.title.y=element_text(size=14,color="black"),
        panel.background = element_blank(),
        legend.position = "none")
hpdis$Treatment<-factor(hpdis$treatment, labels = c("Phi[2]~Treatment~A", "Phi[2]~Treatment~B"))
p2Line<-hpdis%>%
  filter(variable=="phi2")%>%
  ggplot()+
  facet_wrap(~Treatment, labeller = label_parsed)+
  geom_line(aes(x=t, y=length, color=treatment), size=1, linetype=5)+
  geom_smooth(aes(x=t, y=length, color=treatment), se=F, size=1)+
  scale_color_manual(values = c(PlasmaCol, MakoCol))+
  scale_x_continuous(breaks=seq(0,25, 5))+
  scale_y_continuous(breaks=seq(0,13, 5))+
  labs(x="")+
  ylab("HPDI Length")+
  theme_minimal()+
  theme(strip.background = element_rect(colour="white",
                                        fill="white"),
        plot.background=element_rect(fill="white", color="white"),
        axis.title.x=element_text(size=14,color="black"),
        axis.title.y=element_text(size=14,color="black"),
        panel.background = element_blank(),
        legend.position = "none")


hpdis$Treatment<-factor(hpdis$treatment, labels = c("Phi[3]~Treatment~A", "Phi[3]~Treatment~B"))
p3Line<-hpdis%>%
  filter(variable=="phi3")%>%
  ggplot()+
  facet_wrap(~Treatment, labeller = label_parsed)+
  geom_line(aes(x=t, y=length, color=treatment), size=1, linetype=5)+
  geom_smooth(aes(x=t, y=length, color=treatment), se=F, size=1)+
  scale_color_manual(values = c(PlasmaCol, MakoCol))+
  scale_x_continuous(breaks=seq(0,25, 5))+
  scale_y_continuous(breaks=seq(0,5, 1))+
  labs(x="Time")+
  ylab(" ")+
  theme_minimal()+
  theme(strip.background = element_rect(colour="white",
                                        fill="white"),
        plot.background=element_rect(fill="white", color="white"),
        axis.title.x=element_text(size=14,color="black"),
        axis.title.y=element_text(size=14,color="black"),
        panel.background = element_blank(),
        legend.position = "none")

patch<-(p1Line)/(p2Line)/(p3Line)
patch
# try these with geom_line, see how much smoothing is happening
ggsave("logisticHDPI_lines_smooths.png",patch,width=6,height=10, dpi=300, bg="#ffffff")


#
# Information Density 
#
# 

figDat<-figureData_phi1B
colnames(figDat)<-paste0("t",c(1:ncol(figDat)))

for (i in 1:ncol(figDat)){
  if (i==1){
    infoDensityDf<-data.frame(variable="Phi1", treatment="A",
                              deltaInfoDensity = 0, time = 0)
  }
  if(i > 1){
    print(colnames(figureData_phi1A)[i])
    eval(parse(text=paste0("prevHpdi<-HPDI(figDat$t",i-1,", prob=0.95)")))
    eval(parse(text=paste0("newHpdi<-HPDI(figDat$t",i,", prob=0.95)")))
    
    eval(parse(text=paste0("dprev = density(figDat$t",i-1,")")))
    eval(parse(text=paste0("dnew = density(figDat$t",i,")")))
    
    fprev = approxfun(dprev$x, dprev$y)
    fnew = approxfun(dnew$x, dnew$y)
    
    eval(parse(text=paste0('newMean<-quantile(density(figDat$t',i,')$x, probs=0.5)')))
    
    ovrng <- c(max(min(dprev$x), min(dnew$x)), min(max(dprev$x), max(dnew$x)))
    e <- seq(min(ovrng), max(ovrng), length.out=500)
    h <- fnew(e)-fprev(e)
    area<-sum( (h[-1]+h[-length(h)]) /2 *diff(e) *(h[-1]>=0+0))
    # 
    # plot(dnew, main=paste0(colnames(figureData_phi1A)[i-1],"=green, ", 
    #                        colnames(figureData_phi1A)[i],"=black"),
    #      xlim=c(100,250), ylim=c(0,.22))
    # lines(dprev, col="green")
    # jj<-which(h>0 & seq_along(h) %% 5==0)
    # j<-e[jj]
    # segments(j, fprev(j), j, fprev(j)+h[jj])
    
    iterDf<-data.frame(variable=infoDensityDf$variable[1], 
                       treatment=infoDensityDf$treatment[1],
                       deltaInfoDensity = area,
                       time = as.numeric(str_remove_all(colnames(figureData_phi1A),"t")[i]))
    
    infoDensityDf<-rbind(infoDensityDf, iterDf)
  }
}
infoDensityDf

infoDensityDf%>%
  ggplot()+
  geom_line(aes(x=factor(time), y=deltaInfoDensity, group=variable))+
  labs(x="Time", y="Change in Area", title=expression(~~Phi["1B"]))+
  theme_minimal()






# # # Flat Priors

print(load("earlyStoppingSim/bayesian_updating/parameterUpdating_DataAndModels_flatPriors.rdata"))
# flat is non-conjugate...
plot(density(as.data.frame(fit_5$fit)$b_phi1_treatmentb), main="Day 5")
plot(density(as.data.frame(fit_7$fit)$b_phi1_treatmentb), main="Day 7")
plot(density(as.data.frame(fit_11$fit)$b_phi1_treatmentb), main="Day 11")
plot(density(as.data.frame(fit_17$fit)$b_phi1_treatmentb), main="Day 17")

dprev<-density(as.data.frame(fit_3$fit)$b_phi1_treatmentb)
fprev<-approxfun(dprev$x, dprev$y)
plot(dprev)
probs<-fprev(1:1000)
plot(probs)

probs<-probs/sum(probs)

rCustom<-function(x){
  runif(1,0,1000)
}

den<-dnorm #replace with your own density

#calculates the cdf by numerical integration
cdf<-function(x) integrate(den,0,x)[[1]]

#inverts the cdf
inverse.cdf<-function(x,cdf,starting.value=0.1){
  lower.found<-FALSE
  lower<-starting.value
  while(!lower.found){
    if(cdf(lower)>=(x-.000001))
      lower<-lower-(lower-starting.value)^2-1
    else
      lower.found<-TRUE
  }
  upper.found<-FALSE
  upper<-starting.value
  while(!upper.found){
    if(cdf(upper)<=(x+.000001))
      upper<-upper+(upper-starting.value)^2+1
    else
      upper.found<-TRUE
  }
  uniroot(function(y) cdf(y)-x,c(lower,upper))$root
}

#generates 1000 random variables of distribution 'den'
vars<-apply(matrix(runif(1000)),1,function(x) inverse.cdf(x,cdf))
hist(vars)


set.seed(123)
figureData_phi1B<-data.frame(
  prior = runif(3000, 0, 1000),
  t3 = rlnorm(3000,
              mean(log(as.data.frame(fit_3$fit)$b_phi1_treatmentb)),
              sd(log(as.data.frame(fit_3$fit)$b_phi1_treatmentb))),
  t5 = rlnorm(3000,
              mean(log(as.data.frame(fit_3$fit)$b_phi1_treatmentb)),
              sd(log(as.data.frame(fit_3$fit)$b_phi1_treatmentb))),
  t7 = rlnorm(3000,
              mean(log(as.data.frame(fit_7$fit)$b_phi1_treatmentb)),
              sd(log(as.data.frame(fit_7$fit)$b_phi1_treatmentb))),
  t9 = rlnorm(3000,
              mean(log(as.data.frame(fit_9$fit)$b_phi1_treatmentb)),
              sd(log(as.data.frame(fit_9$fit)$b_phi1_treatmentb))),
  t11 = rlnorm(3000,
               mean(log(as.data.frame(fit_11$fit)$b_phi1_treatmentb)),
               sd(log(as.data.frame(fit_11$fit)$b_phi1_treatmentb))),
  t13 = rlnorm(3000,
               mean(log(as.data.frame(fit_13$fit)$b_phi1_treatmentb)),
               sd(log(as.data.frame(fit_13$fit)$b_phi1_treatmentb))),
  t15 = rlnorm(3000,
               mean(log(as.data.frame(fit_15$fit)$b_phi1_treatmentb)),
               sd(log(as.data.frame(fit_15$fit)$b_phi1_treatmentb))),
  t17 = rlnorm(3000,
               mean(log(as.data.frame(fit_17$fit)$b_phi1_treatmentb)),
               sd(log(as.data.frame(fit_17$fit)$b_phi1_treatmentb))),
  t19 = rlnorm(3000,
               mean(log(as.data.frame(fit_19$fit)$b_phi1_treatmentb)),
               sd(log(as.data.frame(fit_19$fit)$b_phi1_treatmentb))),
  t21 = rlnorm(3000,
               mean(log(as.data.frame(fit_21$fit)$b_phi1_treatmentb)),
               sd(log(as.data.frame(fit_21$fit)$b_phi1_treatmentb))),
  t23 = rlnorm(3000,
               mean(log(as.data.frame(fit_23$fit)$b_phi1_treatmentb)),
               sd(log(as.data.frame(fit_23$fit)$b_phi1_treatmentb))),
  t25 = rlnorm(3000,
               mean(log(as.data.frame(fit_25$fit)$b_phi1_treatmentb)),
               sd(log(as.data.frame(fit_25$fit)$b_phi1_treatmentb))))

figureData_phi1A<-data.frame(
  prior = runif(3000, 0, 1000),
  t3 = rlnorm(3000,
              mean(log(as.data.frame(fit_3$fit)$b_phi1_treatmenta)),
              sd(log(as.data.frame(fit_3$fit)$b_phi1_treatmenta))),
  t5 = rlnorm(3000,
              mean(log(as.data.frame(fit_3$fit)$b_phi1_treatmenta)),
              sd(log(as.data.frame(fit_3$fit)$b_phi1_treatmenta))),
  t7 = rlnorm(3000,
              mean(log(as.data.frame(fit_7$fit)$b_phi1_treatmenta)),
              sd(log(as.data.frame(fit_7$fit)$b_phi1_treatmenta))),
  t9 = rlnorm(3000,
              mean(log(as.data.frame(fit_9$fit)$b_phi1_treatmenta)),
              sd(log(as.data.frame(fit_9$fit)$b_phi1_treatmenta))),
  t11 = rlnorm(3000,
               mean(log(as.data.frame(fit_11$fit)$b_phi1_treatmenta)),
               sd(log(as.data.frame(fit_11$fit)$b_phi1_treatmenta))),
  t13 = rlnorm(3000,
               mean(log(as.data.frame(fit_13$fit)$b_phi1_treatmenta)),
               sd(log(as.data.frame(fit_13$fit)$b_phi1_treatmenta))),
  t15 = rlnorm(3000,
               mean(log(as.data.frame(fit_15$fit)$b_phi1_treatmenta)),
               sd(log(as.data.frame(fit_15$fit)$b_phi1_treatmenta))),
  t17 = rlnorm(3000,
               mean(log(as.data.frame(fit_17$fit)$b_phi1_treatmenta)),
               sd(log(as.data.frame(fit_17$fit)$b_phi1_treatmenta))),
  t19 = rlnorm(3000,
               mean(log(as.data.frame(fit_19$fit)$b_phi1_treatmenta)),
               sd(log(as.data.frame(fit_19$fit)$b_phi1_treatmenta))),
  t21 = rlnorm(3000,
               mean(log(as.data.frame(fit_21$fit)$b_phi1_treatmenta)),
               sd(log(as.data.frame(fit_21$fit)$b_phi1_treatmenta))),
  t23 = rlnorm(3000,
               mean(log(as.data.frame(fit_23$fit)$b_phi1_treatmenta)),
               sd(log(as.data.frame(fit_23$fit)$b_phi1_treatmenta))),
  t25 = rlnorm(3000,
               mean(log(as.data.frame(fit_25$fit)$b_phi1_treatmenta)),
               sd(log(as.data.frame(fit_25$fit)$b_phi1_treatmenta))))

figureData_phi2A<-data.frame(
  prior = runif(3000, 0, 1000),
  t3 = rlnorm(3000,
              mean(log(as.data.frame(fit_3$fit)$b_phi2_treatmenta)),
              sd(log(as.data.frame(fit_3$fit)$b_phi2_treatmenta))),
  t5 = rlnorm(3000,
              mean(log(as.data.frame(fit_3$fit)$b_phi2_treatmenta)),
              sd(log(as.data.frame(fit_3$fit)$b_phi2_treatmenta))),
  t7 = rlnorm(3000,
              mean(log(as.data.frame(fit_7$fit)$b_phi2_treatmenta)),
              sd(log(as.data.frame(fit_7$fit)$b_phi2_treatmenta))),
  t9 = rlnorm(3000,
              mean(log(as.data.frame(fit_9$fit)$b_phi2_treatmenta)),
              sd(log(as.data.frame(fit_9$fit)$b_phi2_treatmenta))),
  t11 = rlnorm(3000,
               mean(log(as.data.frame(fit_11$fit)$b_phi2_treatmenta)),
               sd(log(as.data.frame(fit_11$fit)$b_phi2_treatmenta))),
  t13 = rlnorm(3000,
               mean(log(as.data.frame(fit_13$fit)$b_phi2_treatmenta)),
               sd(log(as.data.frame(fit_13$fit)$b_phi2_treatmenta))),
  t15 = rlnorm(3000,
               mean(log(as.data.frame(fit_15$fit)$b_phi2_treatmenta)),
               sd(log(as.data.frame(fit_15$fit)$b_phi2_treatmenta))),
  t17 = rlnorm(3000,
               mean(log(as.data.frame(fit_17$fit)$b_phi2_treatmenta)),
               sd(log(as.data.frame(fit_17$fit)$b_phi2_treatmenta))),
  t19 = rlnorm(3000,
               mean(log(as.data.frame(fit_19$fit)$b_phi2_treatmenta)),
               sd(log(as.data.frame(fit_19$fit)$b_phi2_treatmenta))),
  t21 = rlnorm(3000,
               mean(log(as.data.frame(fit_21$fit)$b_phi2_treatmenta)),
               sd(log(as.data.frame(fit_21$fit)$b_phi2_treatmenta))),
  t23 = rlnorm(3000,
               mean(log(as.data.frame(fit_23$fit)$b_phi2_treatmenta)),
               sd(log(as.data.frame(fit_23$fit)$b_phi2_treatmenta))),
  t25 = rlnorm(3000,
               mean(log(as.data.frame(fit_25$fit)$b_phi2_treatmenta)),
               sd(log(as.data.frame(fit_25$fit)$b_phi2_treatmenta))))

figureData_phi2B<-data.frame(
  prior = runif(3000, 0, 1000),
  t3 = rlnorm(3000,
              mean(log(as.data.frame(fit_3$fit)$b_phi2_treatmentb)),
              sd(log(as.data.frame(fit_3$fit)$b_phi2_treatmentb))),
  t5 = rlnorm(3000,
              mean(log(as.data.frame(fit_3$fit)$b_phi2_treatmentb)),
              sd(log(as.data.frame(fit_3$fit)$b_phi2_treatmentb))),
  t7 = rlnorm(3000,
              mean(log(as.data.frame(fit_7$fit)$b_phi2_treatmentb)),
              sd(log(as.data.frame(fit_7$fit)$b_phi2_treatmentb))),
  t9 = rlnorm(3000,
              mean(log(as.data.frame(fit_9$fit)$b_phi2_treatmentb)),
              sd(log(as.data.frame(fit_9$fit)$b_phi2_treatmentb))),
  t11 = rlnorm(3000,
               mean(log(as.data.frame(fit_11$fit)$b_phi2_treatmentb)),
               sd(log(as.data.frame(fit_11$fit)$b_phi2_treatmentb))),
  t13 = rlnorm(3000,
               mean(log(as.data.frame(fit_13$fit)$b_phi2_treatmentb)),
               sd(log(as.data.frame(fit_13$fit)$b_phi2_treatmentb))),
  t15 = rlnorm(3000,
               mean(log(as.data.frame(fit_15$fit)$b_phi2_treatmentb)),
               sd(log(as.data.frame(fit_15$fit)$b_phi2_treatmentb))),
  t17 = rlnorm(3000,
               mean(log(as.data.frame(fit_17$fit)$b_phi2_treatmentb)),
               sd(log(as.data.frame(fit_17$fit)$b_phi2_treatmentb))),
  t19 = rlnorm(3000,
               mean(log(as.data.frame(fit_19$fit)$b_phi2_treatmentb)),
               sd(log(as.data.frame(fit_19$fit)$b_phi2_treatmentb))),
  t21 = rlnorm(3000,
               mean(log(as.data.frame(fit_21$fit)$b_phi2_treatmentb)),
               sd(log(as.data.frame(fit_21$fit)$b_phi2_treatmentb))),
  t23 = rlnorm(3000,
               mean(log(as.data.frame(fit_23$fit)$b_phi2_treatmentb)),
               sd(log(as.data.frame(fit_23$fit)$b_phi2_treatmentb))),
  t25 = rlnorm(3000,
               mean(log(as.data.frame(fit_25$fit)$b_phi2_treatmentb)),
               sd(log(as.data.frame(fit_25$fit)$b_phi2_treatmentb))))

figureData_phi3A<-data.frame(
  prior = runif(3000, 0, 1000),
  t3 = rlnorm(3000,
              mean(log(as.data.frame(fit_3$fit)$b_phi3_treatmenta)),
              sd(log(as.data.frame(fit_3$fit)$b_phi3_treatmenta))),
  t5 = rlnorm(3000,
              mean(log(as.data.frame(fit_3$fit)$b_phi3_treatmenta)),
              sd(log(as.data.frame(fit_3$fit)$b_phi3_treatmenta))),
  t7 = rlnorm(3000,
              mean(log(as.data.frame(fit_7$fit)$b_phi3_treatmenta)),
              sd(log(as.data.frame(fit_7$fit)$b_phi3_treatmenta))),
  t9 = rlnorm(3000,
              mean(log(as.data.frame(fit_9$fit)$b_phi3_treatmenta)),
              sd(log(as.data.frame(fit_9$fit)$b_phi3_treatmenta))),
  t11 = rlnorm(3000,
               mean(log(as.data.frame(fit_11$fit)$b_phi3_treatmenta)),
               sd(log(as.data.frame(fit_11$fit)$b_phi3_treatmenta))),
  t13 = rlnorm(3000,
               mean(log(as.data.frame(fit_13$fit)$b_phi3_treatmenta)),
               sd(log(as.data.frame(fit_13$fit)$b_phi3_treatmenta))),
  t15 = rlnorm(3000,
               mean(log(as.data.frame(fit_15$fit)$b_phi3_treatmenta)),
               sd(log(as.data.frame(fit_15$fit)$b_phi3_treatmenta))),
  t17 = rlnorm(3000,
               mean(log(as.data.frame(fit_17$fit)$b_phi3_treatmenta)),
               sd(log(as.data.frame(fit_17$fit)$b_phi3_treatmenta))),
  t19 = rlnorm(3000,
               mean(log(as.data.frame(fit_19$fit)$b_phi3_treatmenta)),
               sd(log(as.data.frame(fit_19$fit)$b_phi3_treatmenta))),
  t21 = rlnorm(3000,
               mean(log(as.data.frame(fit_21$fit)$b_phi3_treatmenta)),
               sd(log(as.data.frame(fit_21$fit)$b_phi3_treatmenta))),
  t23 = rlnorm(3000,
               mean(log(as.data.frame(fit_23$fit)$b_phi3_treatmenta)),
               sd(log(as.data.frame(fit_23$fit)$b_phi3_treatmenta))),
  t25 = rlnorm(3000,
               mean(log(as.data.frame(fit_25$fit)$b_phi3_treatmenta)),
               sd(log(as.data.frame(fit_25$fit)$b_phi3_treatmenta))))

figureData_phi3B<-data.frame(
  prior = runif(3000, 0, 1000),
  t3 = rlnorm(3000,
              mean(log(as.data.frame(fit_3$fit)$b_phi3_treatmentb)),
              sd(log(as.data.frame(fit_3$fit)$b_phi3_treatmentb))),
  t5 = rlnorm(3000,
              mean(log(as.data.frame(fit_3$fit)$b_phi3_treatmentb)),
              sd(log(as.data.frame(fit_3$fit)$b_phi3_treatmentb))),
  t7 = rlnorm(3000,
              mean(log(as.data.frame(fit_7$fit)$b_phi3_treatmentb)),
              sd(log(as.data.frame(fit_7$fit)$b_phi3_treatmentb))),
  t9 = rlnorm(3000,
              mean(log(as.data.frame(fit_9$fit)$b_phi3_treatmentb)),
              sd(log(as.data.frame(fit_9$fit)$b_phi3_treatmentb))),
  t11 = rlnorm(3000,
               mean(log(as.data.frame(fit_11$fit)$b_phi3_treatmentb)),
               sd(log(as.data.frame(fit_11$fit)$b_phi3_treatmentb))),
  t13 = rlnorm(3000,
               mean(log(as.data.frame(fit_13$fit)$b_phi3_treatmentb)),
               sd(log(as.data.frame(fit_13$fit)$b_phi3_treatmentb))),
  t15 = rlnorm(3000,
               mean(log(as.data.frame(fit_15$fit)$b_phi3_treatmentb)),
               sd(log(as.data.frame(fit_15$fit)$b_phi3_treatmentb))),
  t17 = rlnorm(3000,
               mean(log(as.data.frame(fit_17$fit)$b_phi3_treatmentb)),
               sd(log(as.data.frame(fit_17$fit)$b_phi3_treatmentb))),
  t19 = rlnorm(3000,
               mean(log(as.data.frame(fit_19$fit)$b_phi3_treatmentb)),
               sd(log(as.data.frame(fit_19$fit)$b_phi3_treatmentb))),
  t21 = rlnorm(3000,
               mean(log(as.data.frame(fit_21$fit)$b_phi3_treatmentb)),
               sd(log(as.data.frame(fit_21$fit)$b_phi3_treatmentb))),
  t23 = rlnorm(3000,
               mean(log(as.data.frame(fit_23$fit)$b_phi3_treatmentb)),
               sd(log(as.data.frame(fit_23$fit)$b_phi3_treatmentb))),
  t25 = rlnorm(3000,
               mean(log(as.data.frame(fit_25$fit)$b_phi3_treatmentb)),
               sd(log(as.data.frame(fit_25$fit)$b_phi3_treatmentb))))

probs <- seq(1,25,1)
avg_pal <- plasma(n=length(probs), begin=0.1, end=0.9, alpha=0.8)
avg_pal2 <- mako(n=length(probs), begin=0.1, end=0.9, alpha=0.8)
phi1B_vir<-ggplot(figureData_phi1B)+
  #geom_density(aes(x=prior), fill="black", alpha=0.5)+
  lapply(seq(25,3,-2), function(i) geom_density(aes_string(x=paste0("t",i)), 
                                                fill=avg_pal2[i],
                                                color=avg_pal2[i]))+
  geom_density(aes(x=prior), fill="black", alpha=0.55, color="black")+
  annotation_custom(grid::textGrob(label = expression(Phi["1B"]~"= 160"),  
                                   x = unit(0.25, "npc"), 
                                   y = unit(0.9, "npc")))+
  xlab(~~Phi["1B"]~Distribution)+
  ylab("Density")+
  coord_cartesian(xlim=c(75,250), ylim=c(0,1))+
  theme_minimal()
phi1B_vir
phi1A_vir<-ggplot(figureData_phi1A)+
  #geom_density(aes(x=prior), fill="black", alpha=0.5)+
  lapply(seq(25,3,-2), function(i) geom_density(aes_string(x=paste0("t",i)), 
                                                fill=avg_pal[i],
                                                color=avg_pal[i]))+
  geom_density(aes(x=prior), fill="black", alpha=0.55, color="black")+
  annotation_custom(grid::textGrob(label = expression(Phi["1A"]~"= 200"),  
                                   x = unit(0.25, "npc"), 
                                   y = unit(0.9, "npc")))+
  xlab(~~Phi["1A"]~Distribution)+
  ylab("Density")+
  coord_cartesian(xlim=c(75,250), ylim=c(0,1))+
  theme_minimal()
phi1A_vir
phi2B_vir<-ggplot(figureData_phi2B)+
  #geom_density(aes(x=prior), fill="black", alpha=0.5)+
  lapply(seq(25,3,-2), function(i) geom_density(aes_string(x=paste0("t",i)), 
                                                fill=avg_pal2[i],
                                                color=avg_pal2[i]))+
  geom_density(aes(x=prior), fill="black", alpha=0.55, color="black")+
  annotation_custom(grid::textGrob(label = expression(Phi["2B"]~"= 13"),  
                                   x = unit(0.25, "npc"), 
                                   y = unit(0.9, "npc")))+
  xlab(~~Phi["2B"]~Distribution)+
  ylab("Density")+
  coord_cartesian(xlim=c(5,25), ylim=c(0,13))+
  theme_minimal()

phi2A_vir<-ggplot(figureData_phi2A)+
  #geom_density(aes(x=prior), fill="black", alpha=0.5)+
  lapply(seq(25,3,-2), function(i) geom_density(aes_string(x=paste0("t",i)), 
                                                fill=avg_pal[i],
                                                color=avg_pal[i]))+
  geom_density(aes(x=prior), fill="black", alpha=0.55, color="black")+
  annotation_custom(grid::textGrob(label = expression(Phi["2A"]~"= 13"),  
                                   x = unit(0.25, "npc"), 
                                   y = unit(0.9, "npc")))+
  xlab(~~Phi["2A"]~Distribution)+
  ylab("Density")+
  coord_cartesian(xlim=c(5,25), ylim=c(0,13))+
  theme_minimal()

phi3B_vir<-ggplot(figureData_phi3B)+
  #geom_density(aes(x=prior), fill="black", alpha=0.5)+
  lapply(seq(25,3,-2), function(i) geom_density(aes_string(x=paste0("t",i)), 
                                                fill=avg_pal2[i],
                                                color=avg_pal2[i]))+
  geom_density(aes(x=prior), fill="black", alpha=0.55, color="black")+
  annotation_custom(grid::textGrob(label = expression(Phi["3B"]~"= 3.5"),  
                                   x = unit(0.25, "npc"), 
                                   y = unit(0.9, "npc")))+
  xlab(~~Phi["3B"]~Distribution)+
  ylab("Density")+
  coord_cartesian(xlim=c(0,7.5), ylim=c(0,30))+
  theme_minimal()

phi3A_vir<-ggplot(figureData_phi3A)+
  #geom_density(aes(x=prior), fill="black", alpha=0.5)+
  lapply(seq(25,3,-2), function(i) geom_density(aes_string(x=paste0("t",i)), 
                                                fill=avg_pal[i],
                                                color=avg_pal[i]))+
  geom_density(aes(x=prior), fill="black", alpha=0.55, color="black")+
  annotation_custom(grid::textGrob(label = expression(Phi["3A"]~"= 3"),  
                                   x = unit(0.25, "npc"), 
                                   y = unit(0.9, "npc")))+
  xlab(~~Phi["3A"]~Distribution)+
  ylab("Density")+
  coord_cartesian(xlim=c(0,7.5), ylim=c(0,30))+
  theme_minimal()

patch<- (phi1A_vir+ phi2A_vir +phi3A_vir)/(phi1B_vir+ phi2B_vir +phi3B_vir)
patch
ggsave("logisticPosteriorUpdating_flatPriors.png",patch,width=12,height=8, dpi=300, bg="#ffffff")

head(df) #plot df, add to patch...
PlasmaCol<-plasma(1, begin=0.65)
MakoCol<-mako(1, begin=0.65)
dataA<-df%>%
  filter(treatment=="a")%>%
  ggplot(aes(time,y,group=interaction(treatment,sample)))+
  geom_line(aes(color=time))+
  scale_color_viridis(option="plasma", begin=0.05, end=0.9)+
  ylab(~~Area~(cm^2))+
  xlab("Time")+
  theme_light()+
  coord_cartesian(ylim=c(0,250))+
  theme(axis.ticks.length=unit(0.2,"cm"))+
  theme(strip.background=element_rect(fill="gray50",color="gray20"),
        strip.text.x=element_text(size=14,color="white"),
        strip.text.y=element_text(size=14,color="white"))+
  theme(axis.title= element_text(size = 10))+
  theme(axis.text = element_text(size = 8))+
  theme(legend.position='none')
dataA
dataB<-df%>%
  filter(treatment=="b")%>%
  ggplot(aes(time,y,group=interaction(treatment,sample)))+
  geom_line(aes(color=time))+
  scale_color_viridis(option="mako", begin=0.05, end=0.9)+
  ylab(~~Area~(cm^2))+
  xlab("Time")+
  theme_light()+
  coord_cartesian(ylim=c(0,250))+
  theme(axis.ticks.length=unit(0.2,"cm"))+
  theme(strip.background=element_rect(fill="gray50",color="gray20"),
        strip.text.x=element_text(size=14,color="white"),
        strip.text.y=element_text(size=14,color="white"))+
  theme(axis.title= element_text(size = 10))+
  theme(axis.text = element_text(size = 8))+
  theme(legend.position='none')

dataB

patch<- (dataA|phi1A_vir|phi2A_vir|phi3A_vir)/(dataB|phi1B_vir|phi2B_vir|phi3B_vir)
patch
ggsave("logisticPosteriorUpdating_flatPriors_withData.png",patch,width=16,height=8, dpi=300, bg="#ffffff")





  

