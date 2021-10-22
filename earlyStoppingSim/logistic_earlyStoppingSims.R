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

library(brms)
library(patchwork)
library(plyr)
library(tidyverse)
library(bayesplot)
library(viridis)

# *********************************************
# Define growthSim(), makeData(), and runSims()
# *********************************************

growthSim <- function(x,phi1,phi2,phi3){ # make function to do a bunch of growthSims and check how well they 
  phi1_r <- phi1+rnorm(1,mean = 0,sd=25) #can be modeled/how well their phi params can be recovered.
  phi2_r <- phi2+rnorm(1,mean=0,sd=1)
  phi3_r <- phi3+rnorm(1,mean=0,sd=.2)
  return(phi1_r / (1 + exp( (phi2_r-x)/phi3_r) ))
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
                                                      "y"=growthSim(1:lastTimepoint, asymptote2, inflectionPoint2, growthRate1),
                                                      stringsAsFactors = F))))
  
}


# *********************************************
# Define runSims()
# *********************************************


runSims<-function(asymptote1 = 200,
                  asymptote2 = 160,
                  iterations = 100,
                  firstTestDay = 15,
                  lastTestDay=NULL){
  
  priorList <- prior(lognormal(log(130), .25),nlpar = "phi1") +
    prior(lognormal(log(12), .25), nlpar = "phi2") +
    prior(lognormal(log(3), .25), nlpar = "phi3") +
    prior(student_t(3,0,5), dpar="sigma") +
    prior(gamma(2,0.1), class="nu")

  for (i in seq(firstTestDay, lastTestDay%||%firstTestDay, by=2)){

    for (iter in 1:iterations){
      print(paste0("day: ",i))
      print(paste0("iteration: ", iter, "/", iterations))
      if (iter==1){
        iterDf<-makeData(asymptote1 = asymptote1, asymptote2 = asymptote2, lastTimepoint = i)
        
        firstFit<-brm(bf(y ~ phi1/(1+exp((phi2-time)/phi3)),
                       sigma~s(time,by=treatment), 
                       phi1 + phi2 + phi3 ~ 0+treatment,
                       autocor = ~arma(~time|sample:treatment,1,1), 
                       nl = TRUE), 
                    family = student, prior = priorList, data = iterDf, iter = 2000, 
                    cores = 4, chains = 4, backend = "cmdstanr", #threads = threading(4), 
                    control = list(adapt_delta = 0.999, max_treedepth = 20),
                    inits = function(){list(b_phi1=rgamma(2,1),b_phi2=rgamma(2,1),b_phi3=rgamma(2,1))})
      
        h<-hypothesis(firstFit, "phi1_treatmenta/phi1_treatmentb > 1")
        allH<-h$hypothesis%>%mutate(iterationNumber=iter, stoppingDay = i, trueAsymptote1=asymptote1, trueAsymptote2=asymptote2)
      } else{
        iterDf<-makeData(asymptote1 = asymptote1, asymptote2 = asymptote2, lastTimepoint = i)
        fitIter<-brm(bf(y ~ phi1/(1+exp((phi2-time)/phi3)),
                        sigma~s(time,by=treatment), 
                        phi1 + phi2 + phi3 ~ 0+treatment,
                        autocor = ~arma(~time|sample:treatment,1,1), 
                        nl = TRUE), 
                     family = student, prior = priorList, data = iterDf, iter = 2000, 
                     cores = 4, chains = 4, backend = "cmdstanr", #threads = threading(4), 
                     control = list(adapt_delta = 0.999, max_treedepth = 20),
                     inits = function(){list(b_phi1=rgamma(2,1),b_phi2=rgamma(2,1),b_phi3=rgamma(2,1))})
       h<-hypothesis(fitIter, "phi1_treatmenta/phi1_treatmentb > 1")
       HYP<-h$hypothesis%>%mutate(iterationNumber=iter, stoppingDay = i, trueAsymptote1=asymptote1, trueAsymptote2=asymptote2)
       allH<-rbind(allH, HYP)
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
#         asymptote2 = 160,
#         iterations = 100,
#         firstTestDay = 15,
#         lastTestDay=25)
# save(runSims_160v200_100iter, "runSims_160v200_100iter.rdata")
# 
# runSims_170v200_100iter<-runSims(asymptote1 = 200,
#                                  asymptote2 = 170,
#                                  iterations = 100,
#                                  firstTestDay = 15,
#                                  lastTestDay=25)
# save(runSims_170v200_100iter, "runSims_170v200_100iter.rdata")

# runSims_180v200_100iter<-runSims(asymptote1 = 200,
#                                  asymptote2 = 180,
#                                  iterations = 100,
#                                  firstTestDay = 15,
#                                  lastTestDay=25)
# save(runSims_180v200_100iter, "runSims_180v200_100iter.rdata")

# runSims_190v200_100iter<-runSims(asymptote1 = 200,
#                                  asymptote2 = 190,
#                                  iterations = 100,
#                                  firstTestDay = 15,
#                                  lastTestDay=25)
# save(runSims_190v200_100iter, "runSims_190v200_100iter.rdata")

runSims_200v200_100iter<-runSims(asymptote1 = 200,
                                 asymptote2 = 200,
                                 iterations = 100,
                                 firstTestDay = 15,
                                 lastTestDay=25)
save(runSims_200v200_100iter, "runSims_200v200_100iter.rdata")

