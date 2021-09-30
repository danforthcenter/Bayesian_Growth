library(brms)
library(ggplot2)
library(stringr)
library(patchwork)
library(cowplot)
library(plyr)
library(tidyverse) 
library(bayesplot)
library(viridis)

# Will have this and a similar function called getHeteroskedasticity (need to shorten that)
# This function will be getGrowthRate
# this will retrieve whatever multiplier is being put on time from a given brms model.

# Make some example data/models to use
growthSimL <- function(x,a){
  a_r <- a+rnorm(1,mean = 0,sd=0.2)
  return(a_r*x) 
}
set.seed(123)
x <- 1:25
df <- rbind(
  do.call(rbind,lapply(1:20,function(i) data.frame("sample"=paste0("sample_",i),"treatment"="a","time"=x,"y"=growthSimL(x,1.10),stringsAsFactors = F))),
  do.call(rbind,lapply(1:20,function(i) data.frame("sample"=paste0("sample_",i),"treatment"="b","time"=x,"y"=growthSimL(x,0.95),stringsAsFactors = F)))
)

# Fit Model
priorLinear <- prior(lognormal(log(1), .25),nlpar = "a") +
  prior(student_t(3,0,5), dpar="sigma") +
  prior(gamma(2,0.1), class="nu")

fitLinear1 <- brm(bf(y ~ a*time, 
                     sigma~time:treatment, 
                     a ~ 0+treatment, 
                     autocor = ~arma(~time|sample:treatment,1,1),nl = T),
                  family = student, prior = priorLinear, data = df, iter = 1000, 
                  cores = 2, chains = 2, backend = "cmdstanr", #threads = threading(4),
                  control = list(adapt_delta = 0.999,max_treedepth = 20),
                  inits = function(){list(b_a=rgamma(2,1),b_b=rgamma(2,1),b_c=rgamma(2,1))})

# Logistic
growthSim <- function(x,phi1,phi2,phi3){
  phi1_r <- phi1+rnorm(1,mean = 0,sd=25)
  phi2_r <- phi2+rnorm(1,mean=0,sd=1)
  phi3_r <- phi3+rnorm(1,mean=0,sd=.2)
  return(phi1_r/(1+exp((phi2_r-x)/phi3_r)))
}

set.seed(123)
x <- 1:25
df <- rbind(
  do.call(rbind,
          lapply(1:20,
                 function(i) data.frame("sample"=paste0("sample_",i),"treatment"="a","time"=x,"y"=growthSim(x,200,13,3),stringsAsFactors = F))),
  do.call(rbind,
          lapply(1:20,
                 function(i) data.frame("sample"=paste0("sample_",i),"treatment"="b","time"=x,"y"=growthSim(x,160,13,3.5),stringsAsFactors = F)))
)

prior_linear <- prior(lognormal(log(130), .25),nlpar = "phi1",coef="treatmenta") +
  prior(lognormal(log(130), .25),nlpar = "phi1",coef="treatmentb") +
  prior(lognormal(log(12), .25), nlpar = "phi2",coef="treatmenta") +
  prior(lognormal(log(12), .25), nlpar = "phi2",coef="treatmentb") +
  prior(lognormal(log(3), .25), nlpar = "phi3",coef="treatmenta") +
  prior(lognormal(log(3), .25), nlpar = "phi3",coef="treatmentb") +
  prior(student_t(3,0,5), dpar="sigma") +
  prior(gamma(2,0.1), class="nu")

fit_logistic <- brm(bf(y ~ phi1/(1+exp((phi2-time)/phi3)),
                       sigma~time+time:treatment, #linear [no s()] variation
                       phi1 + phi2 + phi3 ~ 0+treatment,
                       autocor = ~arma(~time|sample:treatment,1,1),nl = TRUE),
                    family = student, prior = prior_linear, data = df, iter = 2000,
                    cores = 2, chains = 4, backend = "cmdstanr", threads = threading(4),
                    control = list(adapt_delta = 0.999,max_treedepth = 20),
                    inits = function(){list(b_phi1=rgamma(2,1),b_phi2=rgamma(2,1),b_phi3=rgamma(2,1))})

# MonoMolecular

growthSimMM <- function(x,a,b,c){
  a_r <- a+rnorm(1,mean = 0,sd=10)
  b_r <- b+rnorm(1,mean=0,sd=0.025)
  return(a_r-a_r*exp(-b_r*x)) # a_r-a_r*exp(-(b_r*x/c_r))
}

x <- 1:25
df <- rbind(
  do.call(rbind,lapply(1:20,function(i) data.frame("sample"=paste0("sample_",i),"treatment"="a","time"=x,"y"=growthSimMM(x,200,.33),stringsAsFactors = F))),
  do.call(rbind,lapply(1:20,function(i) data.frame("sample"=paste0("sample_",i),"treatment"="b","time"=x,"y"=growthSimMM(x,160,.42),stringsAsFactors = F)))
)

priorMono <- prior(lognormal(log(130), .25),nlpar = "a") +
  prior(lognormal(log(2), .25), nlpar = "b") + 
  prior(student_t(3,0,5), dpar="sigma") +
  prior(gamma(2,0.1), class="nu")

fit_mono <- brm(bf(y ~ a-a*exp(-b*time), 
                   sigma~time:treatment, 
                   a + b ~ 0+treatment, 
                   autocor = ~arma(~time|sample:treatment,1,1),nl = TRUE),
                family = student, prior = priorMono, data = df, iter = 1000, 
                cores = 2, chains = 2, backend = "cmdstanr", #threads = threading(4),
                control = list(adapt_delta = 0.999,max_treedepth = 20),
                inits = function(){list(b_a=rgamma(2,1),b_b=rgamma(2,1),b_c=rgamma(2,1))})

# Review model objects
fitLinear1
fit_logistic
fit_mono

load("fit_quad_cluster.rdata")
fit_quad_cluster
load("spline_model_and_priors.rdata")
fit_spline


growthRate<-function(model, timeVar="time"){
  string<-model$formula[[1]] # pull out formula object
  formula<-as.character(string[3]) # pick the right hand side of formula
  formulaNoPar<-str_replace_all(formula, "(\\(|\\))", " ") # get rid of parenthesis
  formulaTrimmed<-str_extract(formulaNoPar, # pull out time/stuff or time*stuff or stuff*time or stuff/time
                              paste0("(",timeVar, "\\s*(\\/\\s*\\S*|\\*\\s*\\S*))|(((\\S*\\s*\\*\\s*)|(\\S*\\s*\\/\\s?))",timeVar,")"))
  notEvenMyFinalForm<-str_remove_all(formulaTrimmed, "\\s") # get rid of spaces
  finalForm<-str_remove(notEvenMyFinalForm, paste0("(",timeVar,".)|(.",timeVar,")"))  #remove time variable label
  
  p13<-model[13] #get the estimates
  dat<-as.data.frame(p13$fit)%>%select(contains(paste0("b_",str_remove_all(finalForm, "\\W"),"_")))%>%
    summarize( #grab the relevant columns, take the mean of each (should be one per treatment group)
      across(.cols=everything(),
             .fns = mean)
    )
  formulaOut<-str_remove_all(formula, "\\s") #return the formula too for good measure.
  
  fullSummary<-summary(model)
  
  output<-list(finalForm, formulaOut, dat, fullSummary)
  return(output)
}

growthRate(fitLinear1)
growthRate(fit_mono)
growthRate(fit_logistic)

growthRate(fit_spline)














