# This file is starting with issue 55 (https://github.com/joshqsumner/bayesian_plant-science/issues/55)
# Goal is to find a way to weight hypothesis to make appropriate early stopping rules
# First ideas on how to do this revolve around the change from prior to posterior for each day
#   For example, in logistic data on day 11 we expect to see a big change 
#                      in the posterior distribution for inflection point
#   But on day 11 we should be seeing almost no change in posterior distribution for growth rate, and a little for asymptote.
#   By tracking the change in HPDI length we can see which days are most informative.
#   Once the most informative days have passed (f''(t)  = 0), we should be in a better position to test hypotheses.
#   Using that information about which days are most informative we should be able to make good stopping criteria
# Initial steps steps
#     Review how hypothesis() works with MCMC
#         hypothesis() can use a dataframe, so if I can get posterior summaries I might be able to manipulate the spread
#         based on the information density at the current time, larger spread making it less likely to hit significance
#     Write a formula for what I want to have happen
#     Review HPDI changes over time
#     Combine above to make a weighted hypothesis test

library(brms)
library(tidyverse)
library(viridis)
library(patchwork)
library(rethinking)
library(scales)


?brms::hypothesis

# ****************************************
# Load data
# ****************************************

print(load("earlyStoppingSim/bayesian_updating/parameterUpdating_DataAndModels_3to25.rdata"))
print(load("earlyStoppingSim/hypothesisTestWeighting/HPDIS_3_to_25.rdata"))
print(load("earlyStoppingSim/hypothesisTestWeighting/densities_df_3_to_25.rdata"))

# ****************************************
# Look at plots of HPDIS
# ****************************************

labStrings <- as_labeller(c('A' = "Treatment A",
                            "B" = "Treatment B",
                            "phi1" = "~~Phi[1]",
                            "phi2" = "Phi 2",
                            "phi3" = "Phi 3"))
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
p1Line

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

# ****************************************
# Visualize hypotheses and HPDIS
# ****************************************

hyp_list<-lapply(seq(3,25,2), function(i) 
                               eval(parse(text=paste0('as.data.frame(hypothesis(fit_',
                                                                i,
                                                                ', "phi1_treatmenta/phi1_treatmentb > 1.05")$hypothesis)'))))

for (i in 1:length(hyp_list)){
  timeSeq<-seq(3,25,2)
  if (i==1){
    hyp_df = as.data.frame(hyp_list[i])%>%mutate(time=timeSeq[i])
  } else {
    hyp_df<-rbind(hyp_df, as.data.frame(hyp_list[i])%>%mutate(time=timeSeq[i]))
  }
}

PlasmaCol<-plasma(1, begin=0.75)
MakoCol<-mako(1, begin=0.75)

hyp_df%>%
  mutate(colorCol = ifelse(as.numeric(Post.Prob)>=0.95, "1","0"))%>%
  ggplot()+
  geom_line(aes(x=time, y=Post.Prob), color=plasma(1, begin=0.2))+
  geom_point(aes(x=time, y=Post.Prob, color=colorCol))+
  scale_color_viridis(discrete=T, option="plasma", begin=0.2, end=0.4)+
  #geom_line(data=hpdis[hpdis$variable=="phi1" & hpdis$treatment=="A",],
  #          aes(x=t, y=length/max(length), color=treatment), size=1, linetype=5)+
  geom_hline(yintercept=0.95, linetype=5, color=plasma(1, begin=0.4))+
  geom_smooth(data=hpdis[hpdis$variable=="phi1" & hpdis$treatment=="A",],
              aes(x=t, y=diffPct/100),color=PlasmaCol, se=F)+
  
  geom_smooth(data=hpdis[hpdis$variable=="phi1" & hpdis$treatment=="B",],
              aes(x=t, y=diffPct/100),color=MakoCol, se=F)+
  
  scale_x_continuous(breaks=seq(3,25, 2), limits=c(2,25))+
  scale_y_continuous(labels=percent_format(),breaks=seq(0,1, .2), limits=c(0, 1))+
  labs(x="Time", y="Posterior Probability", title=expression(Phi["1A"]/Phi["1B"]>"1.05"),
       subtitle="Percentage change in HPDI shown")+
  theme_minimal()+
  theme(title=element_text(size=14,color="black"),
    axis.title=element_text(size=12,color="black"),
        legend.position = "none")

head(hpdis)


drv <- function(x, y) c(NA, (diff(y) /diff(x))) 
middle_pts <- function(x) c(NA, (x[-1] - diff(x) /2))
print(load("earlyStoppingSim/hypothesisTestWeighting/HPDIS_3_to_25.rdata"))
hpdis<-hpdis%>%
  mutate(diff = length-lag(length),
         lengthPct = length/max(length),
         diffPct=100*(1-(length/lag(length))),
         second_d = drv(middle_pts(t), drv(t, length)))

hyp_df%>%
  mutate(colorCol = ifelse(as.numeric(Post.Prob)>=0.95, "1","0"))%>%
  ggplot()+
  geom_line(aes(x=time, y=Post.Prob), color=plasma(1, begin=0.2))+
  geom_point(aes(x=time, y=Post.Prob, color=colorCol))+
  scale_color_viridis(discrete=T, option="plasma", begin=0.2, end=0.4)+
  
  geom_smooth(data=hpdis[hpdis$variable=="phi1" & hpdis$treatment=="A",],
           aes(x=t, y=length/max(length)),se=F, color=plasma(1, begin=0.4), size=0.5, linetype=1)+
  
  geom_smooth(data=hpdis[hpdis$variable=="phi1" & hpdis$treatment=="A",],
            aes(x=t, y=diffPct/100), se=F, color=plasma(1, begin=0.4), size=0.5, linetype=5)+
  
  geom_smooth(data=hpdis[hpdis$variable=="phi1" & hpdis$treatment=="A",],
              aes(x=t, y=(second_d)),
              se=F, color=plasma(1, begin=0.4), size=0.5, linetype=3)+
  
  geom_hline(yintercept=0.95, linetype=5, color=plasma(1, begin=0.4))+
  scale_x_continuous(breaks=seq(3,25, 2), limits=c(2,25))+
  scale_y_continuous(labels=percent_format(),breaks=seq(0,1, .2), limits=c(0, 1))+
  
  labs(x="Time", y="Posterior Probability", title=expression(Phi["1A"]/Phi["1B"]>"1.05"),
       subtitle=expression(paste(sigma, ", d", sigma,"/dt, d"[2], sigma,"/d"[2],"t"))
       )+
  theme_minimal()+
  theme(title=element_text(size=14,color="black"),
        axis.title=element_text(size=12,color="black"),
        legend.position = "none")

hyp_list<-lapply(seq(3,25,2), function(i) 
  eval(parse(text=paste0('as.data.frame(hypothesis(fit_',
                         i,
                         ', "phi2_treatmenta/phi2_treatmentb > 1.05")$hypothesis)'))))

for (i in 1:length(hyp_list)){
  timeSeq<-seq(3,25,2)
  if (i==1){
    hyp_df_2 = as.data.frame(hyp_list[i])%>%mutate(time=timeSeq[i])
  } else {
    hyp_df_2<-rbind(hyp_df_2, as.data.frame(hyp_list[i])%>%mutate(time=timeSeq[i]))
  }}

hyp_df_2%>%
  mutate(colorCol = ifelse(as.numeric(Post.Prob)>=0.95, "1","0"))%>%
  ggplot()+
  geom_line(aes(x=time, y=Post.Prob), color=plasma(1, begin=0.2))+
  geom_point(aes(x=time, y=Post.Prob, color=colorCol))+
  scale_color_viridis(discrete=T, option="plasma", begin=0.2, end=0.4)+
  #geom_line(data=hpdis[hpdis$variable=="phi1" & hpdis$treatment=="A",],
  #          aes(x=t, y=length/max(length), color=treatment), size=1, linetype=5)+
  geom_hline(yintercept=0.95, linetype=5, color=plasma(1, begin=0.4))+
  geom_smooth(data=hpdis[hpdis$variable=="phi2" & hpdis$treatment=="A",],
              aes(x=t, y=diffPct/100),color=PlasmaCol, se=F)+
  
  geom_smooth(data=hpdis[hpdis$variable=="phi2" & hpdis$treatment=="B",],
              aes(x=t, y=diffPct/100),color=MakoCol, se=F)+
  
  scale_x_continuous(breaks=seq(3,25, 2), limits=c(2,25))+
  scale_y_continuous(labels=percent_format(),breaks=seq(0,1, .2), limits=c(0, 1))+
  labs(x="Time", y="Posterior Probability", title=expression(Phi["2A"]/Phi["2B"]>"1.05"),
       subtitle="Percentage change in HPDI shown")+
  theme_minimal()+
  theme(title=element_text(size=14,color="black"),
        axis.title=element_text(size=12,color="black"),
        legend.position = "none")


hyp_list<-lapply(seq(3,25,2), function(i) 
  eval(parse(text=paste0('as.data.frame(hypothesis(fit_',
                         i,
                         ', "phi3_treatmenta/phi3_treatmentb > 1.05")$hypothesis)'))))

for (i in 1:length(hyp_list)){
  timeSeq<-seq(3,25,2)
  if (i==1){
    hyp_df_3 = as.data.frame(hyp_list[i])%>%mutate(time=timeSeq[i])
  } else {
    hyp_df_3<-rbind(hyp_df_3, as.data.frame(hyp_list[i])%>%mutate(time=timeSeq[i]))
  }}

hyp_df_3%>%
  mutate(colorCol = ifelse(as.numeric(Post.Prob)>=0.95, "1","0"))%>%
  ggplot()+
  geom_line(aes(x=time, y=Post.Prob), color=plasma(1, begin=0.2))+
  geom_point(aes(x=time, y=Post.Prob, color=colorCol))+
  scale_color_viridis(discrete=T, option="plasma", begin=0.2, end=0.4)+
  geom_line(data=hpdis[hpdis$variable=="phi3" & hpdis$treatment=="A",],
            aes(x=t, y=length/max(length)),color=PlasmaCol, size=1, linetype=5)+
  geom_line(data=hpdis[hpdis$variable=="phi3" & hpdis$treatment=="B",],
            aes(x=t, y=length/max(length)), color=MakoCol, size=1, linetype=5)+
  geom_hline(yintercept=0.95, linetype=5, color=plasma(1, begin=0.4))+
  #geom_smooth(data=hpdis[hpdis$variable=="phi3" & hpdis$treatment=="A",],
  #            aes(x=t, y=diffPct/100),color=PlasmaCol, se=F)+
  
  #geom_smooth(data=hpdis[hpdis$variable=="phi3" & hpdis$treatment=="B",],
  #            aes(x=t, y=diffPct/100),color=MakoCol, se=F)+
  
  scale_x_continuous(breaks=seq(3,25, 2), limits=c(2,25))+
  scale_y_continuous(labels=percent_format(),breaks=seq(0,1, .2), limits=c(0, 1))+
  labs(x="Time", y="Posterior Probability", title=expression(Phi["3A"]/Phi["3B"]>"1.05"),
       subtitle="Percentage change in HPDI shown")+
  theme_minimal()+
  theme(title=element_text(size=14,color="black"),
        axis.title=element_text(size=12,color="black"),
        legend.position = "none")


# ****************************************
# Weighting posterior draws
# ****************************************
# Two options that I see (assuming conjugate priors), one is more robust (doesn't require conjugate priors)
# 1: So if hypothesis would be sampling from two posteriors that look like: 
#     log(N(199,.1) vs log(N(180,.15) but we are penalizing tests then make data that is 
#     log(N(199,.15) vs log(N(180,.2)  or something and sample from that instead.
# That option requires a parameterized distribution to increase the variance of.
#
# 2: But the goal is really to make the tails blurry so maybe there is a way to weight 
#     the sampling from the posterior instead of changing the posterior? Similar to how PSIS LOO works.
#     Something where draws are weighted by how close they are to the mean, 
#     with max amount of weight changing by the criteria we've been thinking about?
#     There is a possible route for this already in brms...
#     https://mc-stan.org/posterior/articles/posterior.html#using-custom-functions-1
#     This also opens up the option of positive weights, making a hypothesis EASIER to fulfill,
#     if we believe that should happen

library(posterior)
library(rlang)

getWeightedDraws<-function(mod, param){
  eval(parse(text=paste0("draws_df<-as.data.frame(prepare_predictions(",mod,")$nlpars$",param,"$fe$b)")))
  
  trtA<-draws_df%>%
    rowid_to_column(var="draw")%>%
    select(contains("treatmenta"), draw)%>%
    arrange(desc(b_phi1_treatmenta))%>%
    mutate(rankA = 1:nrow(draws_df),
           distA = abs(b_phi1_treatmenta - mean(b_phi1_treatmenta)),
           distNormA = (distA-mean(distA))/sd(distA),
           distNorm2A = (distA-min(distA))/(max(distA)-min(distA)),
           weightExpA = .27/(exp(distNormA)),
           weightHarshA = 1/(abs(distA)+1),
           weightNormA = (-1*distNorm2A)+1,
           weightMomentA = ((-1*distA/sd(distA))+max(distA))/max(distA))
  trtB<-draws_df%>%
    rowid_to_column(var="draw")%>%
    select(contains("treatmentb"), draw)%>%
    arrange(desc(b_phi1_treatmentb))%>%
    mutate(rankB = 1:nrow(draws_df),
           distB = abs(b_phi1_treatmentb - mean(b_phi1_treatmentb)),
           distNormB = (distB-mean(distB))/sd(distB),
           distNorm2B = (distB-min(distB))/(max(distB)-min(distB)),
           weightExpB = .27/(exp(distNormB)),
           weightHarshB = 1/(abs(distB)+1),
           weightNormB = (-1*distNorm2B)+1,
           weightMomentB = ((-1*distB/sd(distB))+max(distB))/max(distB))
  dat<-left_join(trtA, trtB, by="draw")
  return(dat)
  }

phi1Draws<-getWeightedDraws("fit_25", "phi1")
head(phi1Draws)

phi1Draws%>%
  select(ends_with("A", ignore.case=T))%>%
  rename_with(.fn = function(i) substr(i, 1, nchar(i)-1), .cols = everything())%>%
  ggplot()+
  geom_point(aes(x=rank, y=weightHarsh, color="Harsh Penalty"))+
  geom_point(aes(x=rank, y=weightExp, color="Exp Penalty"))+
  geom_point(aes(x=rank, y=weightNorm, color="Normalized Penalty"))+
  geom_point(aes(x=rank, y=weightMoment, color="1/rank Penalty"))+
  scale_color_viridis(option="plasma", discrete=T)+
  labs(color="Penalty", x="Rank", y="Weight")+
  theme_minimal()+
  theme(title=element_text(size=14,color="black"),
        axis.title=element_text(size=12,color="black"),
        legend.position = "bottom")

# I just totally smoothbrained this whole idea below, 
# I probably need to go to sleep.

unweighted<-hypothesis(phi1Draws, "b_phi1_treatmenta/b_phi1_treatmentb>1.05")
exp<-hypothesis(phi1Draws, "weightExpA/weightExpB>1.05")
harsh<-hypothesis(phi1Draws, "weightHarshA/weightHarshB>1.05")
norm<-hypothesis(phi1Draws, "weightNormA/weightNormB>1.05")
moment<-hypothesis(phi1Draws, "weightMomentA/weightMomentB>1.05")


compareHypotheses<-function(dat, params=c()){
  o=1
  
  dat<-dat%>%
    rename_with(.fn = function(i) str_replace_all(i, "a", "A"), .cols = everything())%>%
    rename_with(.fn = function(i) str_replace_all(i, "b", "B"), .cols = everything())
  
 for (param in params){
   param<-str_replace_all(param, "a", "A")
   param<-str_replace_all(param, "b", "B")
   eval(parse(text=paste0(param,'Hyp<-hypothesis(',quote(dat),
                          ', "',substr(param, 1,nchar(param)-1),'A/',
                          substr(param, 1,nchar(param)-1),'B>1.05")$hypothesis')))
   if(o==1){
     eval(parse(text=paste0("hypDf<-",param,"Hyp")))
     o=o+1
   } else{
     eval(parse(text=paste0("hypDf<-rbind(hypDf, ",param,"Hyp)")))
   }
   }
  return(hypDf)
}

compareHypotheses(dat=phi1Draws, params=c("b_phi1_treatmenta", "weightExpA", "weightNormA"))











data.frame(
  
)


