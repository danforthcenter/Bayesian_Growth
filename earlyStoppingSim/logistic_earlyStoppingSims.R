 # This file is starting with github issue 50 (https://github.com/joshqsumner/bayesian_plant-science/issues/50)
 # This seems like it will be a lot like bayesSSD1.Rmd that I was using for sample size estimates in the RCT class
 # I'll start this with a more fleshed out verison of that function
 # 
 # 
 # 
 # 
 # 
 # 
 # 
 # 
 # 
 # 

# *********************************************
# Load packages
# *********************************************
#condor_submit -i request_cpus=4 request_memory=2GB 'requirements=(TARGET.Machine=!="cylon.datasci.danforthcenter.org" && TARGET.Machine=!="pallas.datasci.danforthcenter.org")'

library(brms)
library(patchwork)
library(plyr)
library(tidyverse)
library(bayesplot)
library(viridis)
library(scales)
library(broom)

# *********************************************
# Define growthSim(), makeData(), and runSims()
# *********************************************

growthSim <- function(x,phi1,phi2,phi3){ # make function to do a bunch of growthSims and check how well they 
  phi1_r <- phi1+rnorm(1,mean = 0,sd=25) #can be modeled/how well their phi params can be recovered.
  phi2_r <- phi2+rnorm(1,mean=0,sd=1)
  phi3_r <- phi3+rnorm(1,mean=0,sd=.2)
  return(phi1_r / (1 + exp( (phi2_r-x)/phi3_r)))
}

# *********************************************
# Define makeData()
# *********************************************

makeData<-function(asymptote1=200, 
                   inflectionPoint1=13, 
                   growthRate1=3,
                   asymptote2=160,
                   inflectionPoint2=13,
                   growthRate2=3, 
                   lastTimepoint=25){
  df <- rbind(
    do.call(rbind, lapply(1:20, function(i) data.frame("sample"=paste0("sample_",i),
                                                       "treatment"="a","time"=1:lastTimepoint,
                                                       "y"=growthSim(1:lastTimepoint, asymptote1, inflectionPoint1, growthRate1),
                                                       stringsAsFactors = F))),
    do.call(rbind, lapply(1:20,function(i) data.frame("sample"=paste0("sample_",i),
                                                      "treatment"="b","time"=1:lastTimepoint,
                                                      "y"=growthSim(1:lastTimepoint, asymptote2, inflectionPoint2, growthRate2),
                                                      stringsAsFactors = F))))
  
}


# *********************************************
# Define runSims()
# *********************************************

runSims<-function(asymptote1 = 200,
                  asymptote2 = 160,
                  iterations = 100,
                  firstTestDay = 11,
                  lastTestDay=NULL){

  priorList <- prior(lognormal(log(130), .25),nlpar = "phi1") +
    prior(lognormal(log(12), .25), nlpar = "phi2") +
    prior(lognormal(log(3), .25), nlpar = "phi3") +
    prior(student_t(3,0,5), dpar="sigma") +
    prior(gamma(2,0.1), class="nu")

  HYPOTHESISTOTESTA<-"phi1_treatmenta/phi1_treatmentb > 1"
  HYPOTHESISTOTESTB<-"phi2_treatmenta/phi2_treatmentb > 1"
  HYPOTHESISTOTESTC<-"phi3_treatmenta/phi3_treatmentb > 1"
  
  for (i in seq(firstTestDay, lastTestDay%||%firstTestDay, by=2)){

    for (iter in 1:iterations){
      print(paste0("day: ",i))
      print(paste0("iteration: ", iter, "/", iterations))
      if (iter==1){
        iterDf<-makeData(asymptote1 = asymptote1, asymptote2 = asymptote2, lastTimepoint = i)
        test<-t.test(iterDf[iterDf$treatment=="a" & iterDf$time==i,]$y, iterDf[iterDf$treatment=="b" & iterDf$time==i,]$y, alternative="g")
        testDf<-broom::tidy(test)
        
        firstFit<-brm(bf(y ~ phi1/(1+exp((phi2-time)/phi3)),
                       sigma~s(time,by=treatment), 
                       phi1 + phi2 + phi3 ~ 0+treatment,
                       autocor = ~arma(~time|sample:treatment,1,1), 
                       nl = TRUE), 
                    family = student, prior = priorList, data = iterDf, iter = 2000, 
                    cores = 4, chains = 4, backend = "cmdstanr",
                    control = list(adapt_delta = 0.999, max_treedepth = 20),
                    inits = function(){list(b_phi1=rgamma(2,1),b_phi2=rgamma(2,1),b_phi3=rgamma(2,1))})
#hypothesis(firstFit, c(HYPOTHESISTOTESTA, HYPOTHESISTOTESTB, HYPOTHESISTOTESTC)) # do a speed test on these
        h<-hypothesis(firstFit, HYPOTHESISTOTESTA)
        HYP1<-h$hypothesis%>%mutate(iterationNumber=iter, stoppingDay = i, trueAsymptote1=asymptote1, trueAsymptote2=asymptote2,tTestP=testDf$p.value, meanA=mean(iterDf[iterDf$treatment=="a" & iterDf$time==i,]$y), meanB = mean(iterDf[iterDf$treatment=="b" & iterDf$time==i,]$y))
        h<-hypothesis(firstFit, HYPOTHESISTOTESTB)
        HYP2<-h$hypothesis%>%mutate(iterationNumber=iter, stoppingDay = i, trueAsymptote1=asymptote1, trueAsymptote2=asymptote2,tTestP=testDf$p.value, meanA=mean(iterDf[iterDf$treatment=="a" & iterDf$time==i,]$y), meanB = mean(iterDf[iterDf$treatment=="b" & iterDf$time==i,]$y))
        h<-hypothesis(firstFit, HYPOTHESISTOTESTC)
        HYP3<-h$hypothesis%>%mutate(iterationNumber=iter, stoppingDay = i, trueAsymptote1=asymptote1, trueAsymptote2=asymptote2,tTestP=testDf$p.value, meanA=mean(iterDf[iterDf$treatment=="a" & iterDf$time==i,]$y), meanB = mean(iterDf[iterDf$treatment=="b" & iterDf$time==i,]$y))
        allH<-rbind(HYP1, HYP2, HYP3)
        } else{
        iterDf<-makeData(asymptote1 = asymptote1, asymptote2 = asymptote2, lastTimepoint = i)
        test<-t.test(iterDf[iterDf$treatment=="a" & iterDf$time==i,]$y, iterDf[iterDf$treatment=="b" & iterDf$time==i,]$y, alternative="l")
        testDf<-broom::tidy(test)
        fitIter<-brm(bf(y ~ phi1/(1+exp((phi2-time)/phi3)),
                        sigma~s(time,by=treatment), 
                        phi1 + phi2 + phi3 ~ 0+treatment,
                        autocor = ~arma(~time|sample:treatment,1,1), 
                        nl = TRUE), 
                     family = student, prior = priorList, data = iterDf, iter = 2000, 
                     cores = 4, chains = 4, backend = "cmdstanr", #threads = threading(4), 
                     control = list(adapt_delta = 0.999, max_treedepth = 20),
                     inits = function(){list(b_phi1=rgamma(2,1),b_phi2=rgamma(2,1),b_phi3=rgamma(2,1))})
       h<-hypothesis(fitIter, HYPOTHESISTOTESTA)
       HYP1<-h$hypothesis%>%mutate(iterationNumber=iter, stoppingDay = i, trueAsymptote1=asymptote1, trueAsymptote2=asymptote2,tTestP=testDf$p.value, meanA=mean(iterDf[iterDf$treatment=="a" & iterDf$time==i,]$y), meanB = mean(iterDf[iterDf$treatment=="b" & iterDf$time==i,]$y))
       h<-hypothesis(fitIter, HYPOTHESISTOTESTB)
       HYP2<-h$hypothesis%>%mutate(iterationNumber=iter, stoppingDay = i, trueAsymptote1=asymptote1, trueAsymptote2=asymptote2,tTestP=testDf$p.value, meanA=mean(iterDf[iterDf$treatment=="a" & iterDf$time==i,]$y), meanB = mean(iterDf[iterDf$treatment=="b" & iterDf$time==i,]$y))
       h<-hypothesis(fitIter, HYPOTHESISTOTESTC)
       HYP3<-h$hypothesis%>%mutate(iterationNumber=iter, stoppingDay = i, trueAsymptote1=asymptote1, trueAsymptote2=asymptote2,tTestP=testDf$p.value, meanA=mean(iterDf[iterDf$treatment=="a" & iterDf$time==i,]$y), meanB = mean(iterDf[iterDf$treatment=="b" & iterDf$time==i,]$y))
       allH<-rbind(allH, HYP1, HYP2, HYP3)
    }
    } # close for loop on 1:iterations
    renameAllHText<-paste0("allH_", i, "<-allH")
    eval(parse(text=renameAllHText))
    
    if(i == firstTestDay){
      combineResultsSoFarText<-paste0("outputDf<-rbind(allH_",i,")")
      eval(parse(text=combineResultsSoFarText))
    } else if(i > firstTestDay){
      combineResultsSoFarText<-paste0("outputDf<-rbind(outputDf, allH_",i,")")
      eval(parse(text=combineResultsSoFarText))
        }
    }#close for loop on days to try stopping

  return(outputDf)
}

# *********************************************
# Use runSims()
# *********************************************

# runSims_160v200_100iter<-runSims(asymptote1 = 200,
#                                  asymptote2 = 160,
#                                  iterations = 100,
#                                  firstTestDay = 13,
#                                  lastTestDay=25)
# save(runSims_160v200_100iter, file="runSims_160v200_100iter_1Hyp_2.rdata")

# runSims_170v200_100iter<-runSims(asymptote1 = 200,
#                                  asymptote2 = 170,
#                                  iterations = 100,
#                                  firstTestDay = 13,
#                                  lastTestDay=25)
# save(runSims_170v200_100iter, file="runSims_170v200_100iter_1Hyp_2.rdata")

# runSims_180v200_100iter<-runSims(asymptote1 = 200,
#                                  asymptote2 = 180,
#                                  iterations = 100,
#                                  firstTestDay = 13,
#                                  lastTestDay=25)
# save(runSims_180v200_100iter, file="runSims_180v200_100iter_1Hyp_2.rdata")

# runSims_190v200_100iter<-runSims(asymptote1 = 200,
#                                  asymptote2 = 190,
#                                  iterations = 100,
#                                  firstTestDay = 13,
#                                  lastTestDay=25)
# save(runSims_190v200_100iter, file="runSims_190v200_100iter_1Hyp_2.rdata")

runSims_200v200_100iter<-runSims(asymptote1 = 200,
                                 asymptote2 = 200,
                                 iterations = 100,
                                 firstTestDay = 13,
                                 lastTestDay=25)
save(runSims_200v200_100iter, file="runSims_200v200_100iter_1Hyp_2.rdata")

# *********************************************
# Check Original runSims() Results 
# *********************************************

dir("earlyStoppingSim/earlyStoppingIterations/", pattern = paste0('*.rdata'))

print(load("earlyStoppingSim/earlyStoppingIterations/runSims_160v200_100iter.rdata"))
print(load("earlyStoppingSim/earlyStoppingIterations/runSims_170v200_100iter.rdata"))
print(load("earlyStoppingSim/earlyStoppingIterations/runSims_180v200_100iter.rdata"))
print(load("earlyStoppingSim/earlyStoppingIterations/runSims_190v200_100iter.rdata"))
print(load("earlyStoppingSim/earlyStoppingIterations/runSims_200v200_100iter.rdata"))

estimatesPlot<-rbind(runSims_160v200_100iter,runSims_170v200_100iter, runSims_180v200_100iter,runSims_190v200_100iter,runSims_200v200_100iter)%>%
  group_by(stoppingDay,trueAsymptote2)%>%
  summarize(estimate = mean(Estimate), trueDiff = (200/trueAsymptote2)-1)%>%distinct(.keep_all = T)%>%
  ggplot()+
  geom_col(aes(x=trueAsymptote2, y=estimate, group=factor(stoppingDay),fill=factor(stoppingDay)), position=position_dodge2())+
  geom_segment(aes(x = trueAsymptote2-5,xend = trueAsymptote2+5 ,y = trueDiff, yend = trueDiff), linetype=5)+
  scale_fill_viridis(discrete=T,option="C")+
  labs(x="Phi1 B", y="Average Estimate", title="Phi1A/Phi1B-1 Estimates", subtitle="Phi1A=200", fill="Stopping Day")+
  theme_minimal() +
  theme(legend.position="bottom")+
  theme(axis.line.y.left = element_line(),
        axis.line.x.bottom = element_line())

p160<-runSims_160v200_100iter%>%
  mutate(foundDiff = ifelse(Post.Prob>=0.95, 1,0))%>%
  group_by(stoppingDay)%>%summarize(foundDiff = mean(foundDiff))%>%
  ggplot()+
  geom_col(aes(x=factor(stoppingDay), y=foundDiff, fill=factor(stoppingDay)))+
  geom_jitter(data=runSims_160v200_100iter, aes(x=factor(stoppingDay), y=Post.Prob), alpha=0.6, width=0.1)+
  geom_jitter(data=runSims_160v200_100iter, aes(x=factor(stoppingDay), y=Post.Prob), alpha=0.6, width=0.1)+
  scale_fill_viridis(discrete=T)+
  scale_y_continuous(labels = percent_format(), limits=c(0,1))+
  geom_hline(yintercept = 0.95, linetype=5)+
  labs(x=" ", y="Percent Detecting Difference", title="160 vs 200 Asymptote")+
  theme_minimal() +
  theme(legend.position="none")+
  theme(axis.line.y.left = element_line(),
        axis.line.x.bottom = element_line())


p170<-runSims_170v200_100iter%>%
  mutate(foundDiff = ifelse(Post.Prob>=0.95, 1,0))%>%
  group_by(stoppingDay)%>%summarize(foundDiff = mean(foundDiff))%>%
  ggplot()+
  geom_col(aes(x=factor(stoppingDay), y=foundDiff, fill=factor(stoppingDay)))+
  geom_jitter(data=runSims_170v200_100iter, aes(x=factor(stoppingDay), y=Post.Prob), alpha=0.6, width=0.1)+  scale_fill_viridis(discrete=T)+
  scale_y_continuous(labels = percent_format(), limits=c(0,1))+
  geom_hline(yintercept = 0.95, linetype=5)+
  labs(x=" ", y=" ", title="170 vs 200 Asymptote")+
  theme_minimal() +
  theme(legend.position="none")+
  theme(axis.line.y.left = element_line(),
        axis.line.x.bottom = element_line())

p180<-runSims_180v200_100iter%>%
  mutate(foundDiff = ifelse(Post.Prob>=0.95, 1,0))%>%
  group_by(stoppingDay)%>%summarize(foundDiff = mean(foundDiff))%>%
  ggplot()+
  geom_col(aes(x=factor(stoppingDay), y=foundDiff, fill=factor(stoppingDay)))+
  geom_jitter(data=runSims_180v200_100iter, aes(x=factor(stoppingDay), y=Post.Prob), alpha=0.6, width=0.1)+  scale_fill_viridis(discrete=T)+
  scale_y_continuous(labels = percent_format(), limits=c(0,1))+
  geom_hline(yintercept = 0.95, linetype=5)+
  labs(x=" ", y=" ", title="180 vs 200 Asymptote")+
  theme_minimal() +
  theme(legend.position="none")+
  theme(axis.line.y.left = element_line(),
        axis.line.x.bottom = element_line())

p190<-runSims_190v200_100iter%>%
  mutate(foundDiff = ifelse(Post.Prob>=0.95, 1,0))%>%
  group_by(stoppingDay)%>%summarize(foundDiff = mean(foundDiff))%>%
  ggplot()+
  geom_col(aes(x=factor(stoppingDay), y=foundDiff, fill=factor(stoppingDay)))+
  geom_jitter(data=runSims_190v200_100iter, aes(x=factor(stoppingDay), y=Post.Prob), alpha=0.6, width=0.1)+  scale_fill_viridis(discrete=T)+
  scale_y_continuous(labels = percent_format(), limits=c(0,1))+
  geom_hline(yintercept = 0.95, linetype=5)+
  labs(x="Stopping Day", y="Percent Detecting Difference", title="190 vs 200 Asymptote")+
  theme_minimal() +
  theme(legend.position="none")+
  theme(axis.line.y.left = element_line(),
        axis.line.x.bottom = element_line())

p200<-runSims_200v200_100iter%>%
  mutate(foundDiff = ifelse(Post.Prob>=0.95, 1,0))%>%
  group_by(stoppingDay)%>%summarize(foundDiff = mean(foundDiff))%>%
  ggplot()+
  geom_col(aes(x=factor(stoppingDay), y=foundDiff, fill=factor(stoppingDay)))+
  geom_jitter(data=runSims_200v200_100iter, aes(x=factor(stoppingDay), y=Post.Prob), alpha=0.6, width=0.1)+  scale_fill_viridis(discrete=T)+
  scale_y_continuous(labels = percent_format(), limits=c(0,1))+
  geom_hline(yintercept = 0.95, linetype=5)+
  labs(x="Stopping Day", y=" ", title="200 vs 200 Asymptote")+
  theme_minimal() +
  theme(legend.position="none")+
  theme(axis.line.y.left = element_line(),
        axis.line.x.bottom = element_line())

patch<-(p160+p170+p180)/(p190+p200+estimatesPlot)
ggsave("160through200vs200_Hyp1_v2.png",patch,width=12,height=8, dpi=300, bg="#ffffff")

# *********************************************
# Check runSims() Results from HypRatio >1.025
# *********************************************

dir("earlyStoppingSim/1025Hyp_take2/", pattern = paste0('*.rdata'))

print(load("earlyStoppingSim/1025Hyp_take2/runSims_160v200_100iter_1025Hyp_2.rdata"))
print(load("earlyStoppingSim/1025Hyp_take2/runSims_170v200_100iter_1025Hyp_2.rdata"))
print(load("earlyStoppingSim/1025Hyp_take2/runSims_180v200_100iter_1025Hyp_2.rdata"))
print(load("earlyStoppingSim/1025Hyp_take2/runSims_190v200_100iter_1025Hyp_2.rdata"))
print(load("earlyStoppingSim/1025Hyp_take2/runSims_200v200_100iter_1025Hyp_2.rdata"))


estimatesPlot<-rbind(runSims_160v200_100iter,
                     runSims_170v200_100iter, 
                     runSims_180v200_100iter,
                     runSims_190v200_100iter,
                     runSims_200v200_100iter)%>%
  filter(grepl("phi1", Hypothesis))%>%
  group_by(stoppingDay,trueAsymptote2)%>%
  summarize(estimate = mean(Estimate), trueDiff = (200/trueAsymptote2)-1.025)%>%
  distinct(.keep_all = T)%>%
  ggplot()+
  geom_col(aes(x=trueAsymptote2, y=estimate, group=factor(stoppingDay),
               fill=factor(stoppingDay)), position=position_dodge2())+
  geom_segment(aes(x = trueAsymptote2-5,xend = trueAsymptote2+5 ,
                   y = trueDiff, yend = trueDiff), linetype=5)+
  scale_fill_viridis(discrete=T,option="C")+
  labs(x="Phi1 B", y="Average Estimate", title="Phi1A/Phi1B-1.025 Estimates", subtitle="Phi1A=200", fill="Stopping Day")+
  theme_minimal() +
  theme(legend.position="bottom")+
  theme(axis.line.y.left = element_line(),
        axis.line.x.bottom = element_line())

p160<-runSims_160v200_100iter%>%
  filter(grepl("phi1", Hypothesis))%>%
  mutate(foundDiff = ifelse(Post.Prob>=0.95, 1,0))%>%
  group_by(stoppingDay)%>%summarize(foundDiff = mean(foundDiff))%>%
  ggplot()+
  geom_col(aes(x=factor(stoppingDay), y=foundDiff, fill=factor(stoppingDay)))+
  geom_jitter(data=filter(runSims_160v200_100iter, grepl("phi1", Hypothesis)),
              aes(x=factor(stoppingDay),
                  y=Post.Prob), alpha=0.5, width=0.1, color="black")+ 
  scale_fill_viridis(discrete=T)+
  scale_y_continuous(labels = percent_format(), limits=c(-0.01,1.05))+
  geom_hline(yintercept = 0.95, linetype=5)+
  geom_text(aes(x=factor(stoppingDay),
                y=1.04,
                label=paste0(round(foundDiff,2))), size=3)+
  labs(x=" ", y="Percent Detecting Difference", title="160 vs 200 Asymptote")+
  theme_minimal() +
  theme(legend.position="none")+
  theme(axis.line.y.left = element_line(),
        axis.line.x.bottom = element_line())
p160
p170<-runSims_170v200_100iter%>%
  filter(grepl("phi1", Hypothesis))%>%
  mutate(foundDiff = ifelse(Post.Prob>=0.95, 1,0))%>%
  group_by(stoppingDay)%>%summarize(foundDiff = mean(foundDiff))%>%
  ggplot()+
  geom_col(aes(x=factor(stoppingDay), y=foundDiff, fill=factor(stoppingDay)))+
  geom_jitter(data=filter(runSims_170v200_100iter, grepl("phi1", Hypothesis)),
              aes(x=factor(stoppingDay),
                  y=Post.Prob), alpha=0.5, width=0.1, color="black")+ 
  scale_fill_viridis(discrete=T)+
  scale_y_continuous(labels = percent_format(), limits=c(-0.01,1.05))+
  geom_hline(yintercept = 0.95, linetype=5)+
  geom_text(aes(x=factor(stoppingDay),
                y=1.04,
                label=paste0(round(foundDiff,2))), size=3)+
  labs(x=" ", y=" ", title="170 vs 200 Asymptote")+
  theme_minimal() +
  theme(legend.position="none")+
  theme(axis.line.y.left = element_line(),
        axis.line.x.bottom = element_line())

p180<-runSims_180v200_100iter%>%
  filter(grepl("phi1", Hypothesis))%>%
  mutate(foundDiff = ifelse(Post.Prob>=0.95, 1,0))%>%
  group_by(stoppingDay)%>%summarize(foundDiff = mean(foundDiff))%>%
  ggplot()+
  geom_col(aes(x=factor(stoppingDay), y=foundDiff, fill=factor(stoppingDay)))+
  geom_jitter(data=filter(runSims_180v200_100iter, grepl("phi1", Hypothesis)),
              aes(x=factor(stoppingDay),
                  y=Post.Prob), alpha=0.5, width=0.1, color="black")+ 
  scale_fill_viridis(discrete=T)+
  scale_y_continuous(labels = percent_format(), limits=c(-0.01,1.05))+
  geom_hline(yintercept = 0.95, linetype=5)+
  geom_text(aes(x=factor(stoppingDay),
                y=1.04,
                label=paste0(round(foundDiff,2))), size=3)+
  labs(x=" ", y=" ", title="180 vs 200 Asymptote")+
  theme_minimal() +
  theme(legend.position="none")+
  theme(axis.line.y.left = element_line(),
        axis.line.x.bottom = element_line())

p190<-runSims_190v200_100iter%>%
  filter(grepl("phi1", Hypothesis))%>%
  mutate(foundDiff = ifelse(Post.Prob>=0.95, 1,0))%>%
  group_by(stoppingDay)%>%summarize(foundDiff = mean(foundDiff))%>%
  ggplot()+
  geom_col(aes(x=factor(stoppingDay), y=foundDiff, fill=factor(stoppingDay)))+
  geom_jitter(data=filter(runSims_190v200_100iter, grepl("phi1", Hypothesis)),
              aes(x=factor(stoppingDay),
                  y=Post.Prob), alpha=0.5, width=0.1, color="black")+ 
  scale_fill_viridis(discrete=T)+
  scale_y_continuous(labels = percent_format(), limits=c(-0.01,1.05))+
  geom_hline(yintercept = 0.95, linetype=5)+
  geom_text(aes(x=factor(stoppingDay),
                y=1.04,
                label=paste0(round(foundDiff,2))), size=3)+
  labs(x="Stopping Day", y="Percent Detecting Difference", title="190 vs 200 Asymptote")+
  theme_minimal() +
  theme(legend.position="none")+
  theme(axis.line.y.left = element_line(),
        axis.line.x.bottom = element_line())
p190
p200<-runSims_200v200_100iter%>%
  filter(grepl("phi1", Hypothesis))%>%#,
  #       tTestP<0.95)%>%
  mutate(foundDiff = ifelse(Post.Prob>=0.95, 1,0))%>%
  group_by(stoppingDay)%>%summarize(foundDiff = mean(foundDiff))%>%
  ggplot()+
  geom_col(aes(x=factor(stoppingDay), y=foundDiff, fill=factor(stoppingDay)))+
  geom_jitter(data=filter(runSims_200v200_100iter, grepl("phi1", Hypothesis)),
              aes(x=factor(stoppingDay),
                  y=ifelse(tTestP>=0.95, NA, Post.Prob)
                  ), alpha=0.5, width=0.1, color="black")+ 
  geom_jitter(data=runSims_200v200_100iter,
              aes(x=factor(stoppingDay),
                  y=ifelse(tTestP>=0.95, Post.Prob, NA)
                  ), alpha=0.65, width=0.1,  color="red")+ 
  scale_fill_viridis(discrete=T)+
  scale_y_continuous(labels = percent_format(), limits=c(-0.01,1.05))+
  geom_hline(yintercept = 0.95, linetype=5)+
  geom_text(aes(x=factor(stoppingDay),
                y=1.04,
                label=paste0(round(foundDiff,2))), size=3)+
  labs(x="Stopping Day", y=" ", title="200 vs 200 Asymptote")+
  theme_minimal() +
  theme(legend.position="none")+
  theme(axis.line.y.left = element_line(),
        axis.line.x.bottom = element_line())
p200
patch<-(p160+p170+p180)/(p190+p200+estimatesPlot)
ggsave("160through200vs200_Hyp1025_detailed.png",patch,width=12,height=8, dpi=300, bg="#ffffff")


# *********************************************
# Check runSims() Estimates
# *********************************************

print(load("earlyStoppingSim/105Hyp/runSims_160v200_100iter_105Hyp.rdata"))
print(load("earlyStoppingSim/105Hyp/runSims_170v200_100iter_105Hyp.rdata"))
print(load("earlyStoppingSim/105Hyp/runSims_180v200_100iter_105Hyp.rdata"))
print(load("earlyStoppingSim/105Hyp/runSims_190v200_100iter_105Hyp.rdata"))
print(load("earlyStoppingSim/105Hyp/runSims_200v200_100iter_105Hyp.rdata"))

laps<-lapply(seq(13,25,2), 
             function(i) 
               round(summary(200/(runSims_200v200_100iter[runSims_200v200_100iter$stoppingDay==i & 
                                                            runSims_200v200_100iter$Hypothesis=='(phi1_treatmenta/phi1_treatmentb)-(1.05) > 0',]$Estimate+1.05))[[5]]-
                     summary(200/(runSims_200v200_100iter[runSims_200v200_100iter$stoppingDay==i  &
                                                            runSims_200v200_100iter$Hypothesis=='(phi1_treatmenta/phi1_treatmentb)-(1.05) > 0',]$Estimate+1.05))[[2]],1))
runSims_200v200_100iter%>%
  filter(grepl("phi1", Hypothesis))%>%
  mutate(color = ifelse(Post.Prob>=0.95, 1, 0))%>%
  mutate(adjEstimate = 200/(Estimate+1.05))%>%
  ggplot()+
  geom_point(aes(x=factor(stoppingDay), y=adjEstimate, group=factor(iterationNumber),
                 color=factor(stoppingDay)), position=position_dodge2(), alpha=0.5)+
  geom_point(aes(x=factor(stoppingDay), y=ifelse(abs(color*adjEstimate)>0,color*adjEstimate,NA), 
                 group=factor(iterationNumber)),color="black", position=position_dodge2(), size=0.5)+
  geom_point(aes(x=factor(stoppingDay), y=ifelse(abs(adjEstimate-200)>20,adjEstimate,NA), 
                 group=factor(iterationNumber)),
             color="black", position=position_dodge2(), size=5, shape="O")+
  lapply(c(1:7), function(i ) annotate("text", x=i, y=235, label=laps[i]))+
  scale_color_viridis(discrete=T,option="C")+
  labs(x="Stopping Day", y="Estimate", title="Phi1B Estimates", 
       subtitle="(True Phi1B = 160) IQR Range Labelled", fill="Stopping Day")+
  theme_minimal() +
  theme(legend.position="none")+
  theme(axis.line.y.left = element_line(),
        axis.line.x.bottom = element_line())

# *********************************************
# Check data iterations
# *********************************************

df1<-makeData(asymptote1=200, 
              inflectionPoint1=13, 
              growthRate1=3,
              asymptote2=200,
              inflectionPoint2=13,
              growthRate2=3, 
              lastTimepoint=25)

ggplot(df1,aes(time,y,group=interaction(treatment,sample)))+
  facet_wrap(vars(treatment))+
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


N<-10000
set.seed(123)
for (i in 1:N){
  df1<-makeData(asymptote1=200, 
                inflectionPoint1=13, 
                growthRate1=3,
                asymptote2=200,
                inflectionPoint2=13,
                growthRate2=3, 
                lastTimepoint=25)
  
  maxA<-df1%>%filter(treatment=="a")%>%slice_max(order_by=y, n=1, with_ties=F)%>%pull(y)
  maxB<-df1%>%filter(treatment=="b")%>%slice_max(order_by=y, n=1, with_ties=F)%>%pull(y)
  minA<-df1%>%filter(treatment=="a")%>%slice_min(order_by=y, n=1, with_ties=F)%>%pull(y)
  minB<-df1%>%filter(treatment=="b")%>%slice_min(order_by=y, n=1, with_ties=F)%>%pull(y)
  
  test<-t.test(df1[df1$treatment=="a" & df1$time==25,]$y, df1[df1$treatment=="b" & df1$time==25,]$y, alternative="l")
  testDf<-broom::tidy(test)
  if(i==1){
    out_df<-data.frame(maxA = maxA, maxB = maxB, minA=minA, minB=minB, iter=i, 
                       tTestP=testDf$p.value, meanA_t25 = testDf$estimate1, 
                       meanB_t25 = testDf$estimate2, 
                       significantDiff = ifelse(testDf$p.value<=0.05, 1, 0))
  } else{
    iter_df<-data.frame(maxA = maxA, maxB = maxB, minA=minA, minB=minB, iter=i, 
                        tTestP=testDf$p.value, meanA_t25 = testDf$estimate1, 
                        meanB_t25 = testDf$estimate2, 
                        significantDiff = ifelse(testDf$p.value<=0.05, 1, 0))
    out_df<-rbind(out_df, iter_df)
  }
  
}
head(out_df)

out_df%>%
  mutate(significantDiff=ifelse(significantDiff==1, "Yes", "No"))%>%
  group_by(significantDiff)%>%
  summarize(pct=n()/10000)%>%
  ggplot()+
  geom_col(aes(x=factor(significantDiff), y=pct, fill=factor(significantDiff)))+
  scale_fill_viridis(discrete=T)+
  scale_y_continuous(labels = percent_format(), limits=c(0,1))+
  geom_text(aes(x=factor(significantDiff), y=pct+0.025, label=paste0(round(pct, 3))))+
  labs(x="Significant T Test", y="")+
  theme(legend.position="none")
  
out_df%>%
  mutate(ratio = meanA_t25/meanB_t25,
         ratioB = ifelse(ratio >=1.05, "Yes", "No"))%>%
  group_by(ratioB)%>%
  summarize(n=n())%>%
  ggplot()+
  geom_col(aes(x=factor(ratioB), y=n, fill=factor(ratioB)))+
  scale_fill_viridis(discrete=T)+
  geom_text(aes(x=factor(ratioB), y=n+125, label=n))+
  labs(x="mean(A)/mean(B) > 1.05", y="")+
  theme(legend.position="none")

out_df%>%
  mutate(absoluteDiff = maxA-maxB,
         absoluteDiffB = ifelse(absoluteDiff >=25, "Yes", "No"))%>%
  group_by(absoluteDiffB)%>%
  summarize(n=n())%>%
  ggplot()+
  geom_col(aes(x=factor(absoluteDiffB), y=n, fill=factor(absoluteDiffB)))+
  scale_fill_viridis(discrete=T)+
  geom_text(aes(x=factor(absoluteDiffB), y=n+125, label=n))+
  labs(x="max(A)-max(B) > 25", y="")+
  theme(legend.position="none")


