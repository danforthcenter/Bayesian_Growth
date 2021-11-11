# This file is following issue 53 (https://github.com/joshqsumner/bayesian_plant-science/issues/53)
# Plan is to use gompertz data and do the same early Stopping Simulations as have been done for logistic data
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

gompertzSim <- function(x,a,b,c){
  a_r <- a+rnorm(1,mean = 0,sd=10)
  b_r <- b+rnorm(1,mean=0,sd=2)
  c_r <- c+rnorm(1,mean=0,sd=.035)
  return(a_r*exp(-b_r*exp(-c_r*x)))
}
makeData<-function(asymptote1=200, 
                   inflectionPoint1=13, 
                   growthRate1=.2,
                   asymptote2=160,
                   inflectionPoint2=13,
                   growthRate2=.2, 
                   lastTimepoint=25){
  df <- rbind(
    do.call(rbind, lapply(1:20, function(i) data.frame("sample"=paste0("sample_",i),
                                                       "treatment"="a","time"=1:lastTimepoint,
                                                       "y"=gompertzSim(1:lastTimepoint, asymptote1, inflectionPoint1, growthRate1),
                                                       stringsAsFactors = F))),
    do.call(rbind, lapply(1:20,function(i) data.frame("sample"=paste0("sample_",i),
                                                      "treatment"="b","time"=1:lastTimepoint,
                                                      "y"=gompertzSim(1:lastTimepoint, asymptote2, inflectionPoint2, growthRate2),
                                                      stringsAsFactors = F))))
  
}


# *********************************************
# Define runSims()
# *********************************************

runSims<-function(asymptote1 = 200,
                  asymptote2 = 160,
                  iterations = 100,
                  firstTestDay = 9,
                  lastTestDay=NULL){
  
  priorList <- prior(lognormal(log(130), .25),nlpar = "a") +
    prior(lognormal(log(12), .25), nlpar = "b") + 
    prior(lognormal(log(1.2), .25), nlpar = "c") + 
    prior(student_t(3,0,5), dpar="sigma") +
    prior(gamma(2,0.1), class="nu")
  
  HYPOTHESISTOTESTA<-"a_treatmenta/a_treatmentb > 1.05"
  HYPOTHESISTOTESTB<-"b_treatmenta/b_treatmentb > 1.05"
  HYPOTHESISTOTESTC<-"c_treatmenta/c_treatmentb > 1.05"
  
  for (i in seq(firstTestDay, lastTestDay%||%firstTestDay, by=2)){
    
    for (iter in 1:iterations){
      print(paste0("day: ",i))
      print(paste0("iteration: ", iter, "/", iterations))
      if (iter==1){
        iterDf<-makeData(asymptote1 = asymptote1, asymptote2 = asymptote2, lastTimepoint = i)
        test<-t.test(iterDf[iterDf$treatment=="a" & iterDf$time==i,]$y, iterDf[iterDf$treatment=="b" & iterDf$time==i,]$y, alternative="g")
        testDf<-broom::tidy(test)
        
        if(i<=9){
          firstFit<-brm(bf(y ~ a*exp(-b*exp(-c*time)),
                           sigma~s(time,by=treatment, k=9), 
                           a + b + c ~ 0+treatment,
                           autocor = ~arma(~time|sample:treatment,1,1), nl = TRUE), 
                        family = student, prior = priorList, data = iterDf, iter = 2000, 
                        cores = 4, chains = 4, backend = "cmdstanr",
                        control = list(adapt_delta = 0.999, max_treedepth = 20),
                        inits = function(){list(b_a=rgamma(2,1),b_b=rgamma(2,1),b_c=rgamma(2,1))})
          
        } else{

        firstFit<-brm(bf(y ~ a*exp(-b*exp(-c*time)),
                           sigma~s(time,by=treatment), 
                           a + b + c ~ 0+treatment,
                           autocor = ~arma(~time|sample:treatment,1,1), 
                           nl = TRUE), 
                        family = student, prior = priorList, data = iterDf, iter = 2000, 
                        cores = 4, chains = 4, backend = "cmdstanr",
                        control = list(adapt_delta = 0.999, max_treedepth = 20),
                        inits = function(){list(b_a=rgamma(2,1),b_b=rgamma(2,1),b_c=rgamma(2,1))})}
        
        #hypothesis(firstFit, c(HYPOTHESISTOTESTA, HYPOTHESISTOTESTB, HYPOTHESISTOTESTC)) # do a speed test on these
        h<-hypothesis(firstFit, HYPOTHESISTOTESTA)
        HYP1<-h$hypothesis%>%mutate(iterationNumber=iter, stoppingDay = i, trueAsymptote1=asymptote1, trueAsymptote2=asymptote2,tTestP=testDf$p.value, meanA=mean(iterDf[iterDf$treatment=="a" & iterDf$time==i,]$y), meanB = mean(iterDf[iterDf$treatment=="b" & iterDf$time==i,]$y))
        h<-hypothesis(firstFit, HYPOTHESISTOTESTB)
        HYP2<-h$hypothesis%>%mutate(iterationNumber=iter, stoppingDay = i, trueAsymptote1=asymptote1, trueAsymptote2=asymptote2,tTestP=testDf$p.value, meanA=mean(iterDf[iterDf$treatment=="a" & iterDf$time==i,]$y), meanB = mean(iterDf[iterDf$treatment=="b" & iterDf$time==i,]$y))
        h<-hypothesis(firstFit, HYPOTHESISTOTESTC)
        HYP3<-h$hypothesis%>%mutate(iterationNumber=iter, stoppingDay = i, trueAsymptote1=asymptote1, trueAsymptote2=asymptote2,tTestP=testDf$p.value, meanA=mean(iterDf[iterDf$treatment=="a" & iterDf$time==i,]$y), meanB = mean(iterDf[iterDf$treatment=="b" & iterDf$time==i,]$y))
        allH<-rbind(HYP1, HYP2, HYP3)
        } else{ # if not first iteration on given day
        iterDf<-makeData(asymptote1 = asymptote1, asymptote2 = asymptote2, lastTimepoint = i)
        test<-t.test(iterDf[iterDf$treatment=="a" & iterDf$time==i,]$y, iterDf[iterDf$treatment=="b" & iterDf$time==i,]$y, alternative="l")
        testDf<-broom::tidy(test)
        
        if(i<=9){
          fitIter<-brm(bf(y ~ a*exp(-b*exp(-c*time)),
                           sigma~s(time,by=treatment, k=9), 
                           a + b + c ~ 0+treatment,
                           autocor = ~arma(~time|sample:treatment,1,1), 
                           nl = TRUE), 
                        family = student, prior = priorList, data = iterDf, iter = 2000, 
                        cores = 4, chains = 4, backend = "cmdstanr",
                        control = list(adapt_delta = 0.999, max_treedepth = 20),
                        inits = function(){list(b_a=rgamma(2,1),b_b=rgamma(2,1),b_c=rgamma(2,1))})
          
        } else{
          
          fitIter<-brm(bf(y ~ a*exp(-b*exp(-c*time)),
                           sigma~s(time,by=treatment), 
                           a + b + c ~ 0+treatment,
                           autocor = ~arma(~time|sample:treatment,1,1), 
                           nl = TRUE), 
                        family = student, prior = priorList, data = iterDf, iter = 2000, 
                        cores = 4, chains = 4, backend = "cmdstanr",
                        control = list(adapt_delta = 0.999, max_treedepth = 20),
                        inits = function(){list(b_a=rgamma(2,1),b_b=rgamma(2,1),b_c=rgamma(2,1))})}
        
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
#                                  firstTestDay = 9,
#                                  lastTestDay=25)
# save(runSims_160v200_100iter, file="earlyStoppingIterations/gompertz105Hyp/gompertzRunSims_160v200_100iter_105Hyp.rdata")

# runSims_170v200_100iter<-runSims(asymptote1 = 200,
#                                  asymptote2 = 170,
#                                  iterations = 100,
#                                  firstTestDay = 9,
#                                  lastTestDay=25)
# save(runSims_170v200_100iter, file="earlyStoppingIterations/gompertz105Hyp/gompertzRunSims_170v200_100iter_105Hyp.rdata")

# runSims_180v200_100iter<-runSims(asymptote1 = 200,
#                                  asymptote2 = 180,
#                                  iterations = 100,
#                                  firstTestDay = 9,
#                                  lastTestDay=25)
# save(runSims_180v200_100iter, file="earlyStoppingIterations/gompertz105Hyp/gompertzRunSims_180v200_100iter_105Hyp.rdata")

# runSims_190v200_100iter<-runSims(asymptote1 = 200,
#                                  asymptote2 = 190,
#                                  iterations = 100,
#                                  firstTestDay = 9,
#                                  lastTestDay=25)
# save(runSims_190v200_100iter, file="earlyStoppingIterations/gompertz105Hyp/gompertzRunSims_190v200_100iter_105Hyp.rdata")

# runSims_200v200_100iter<-runSims(asymptote1 = 200,
#                                  asymptote2 = 200,
#                                  iterations = 100,
#                                  firstTestDay = 9,
#                                  lastTestDay=25)
# save(runSims_200v200_100iter, file="earlyStoppingIterations/gompertz105Hyp/gompertzRunSims_200v200_100iter_105Hyp.rdata")
