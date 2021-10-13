library(brms)
library(ggplot2)
library(stringr)
library(patchwork)
library(cowplot)
library(plyr)
library(tidyverse) 
library(bayesplot)
library(viridis)

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

fitLinear_spline <- brm(bf(y ~ a*time, 
                     sigma~s(time,by=treatment),
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

prior_logistic <- prior(lognormal(log(130), .25),nlpar = "phi1",coef="treatmenta") +
  prior(lognormal(log(130), .25),nlpar = "phi1",coef="treatmentb") +
  prior(lognormal(log(12), .25), nlpar = "phi2",coef="treatmenta") +
  prior(lognormal(log(12), .25), nlpar = "phi2",coef="treatmentb") +
  prior(lognormal(log(3), .25), nlpar = "phi3",coef="treatmenta") +
  prior(lognormal(log(3), .25), nlpar = "phi3",coef="treatmentb") +
  prior(student_t(3,0,5), dpar="sigma") +
  prior(gamma(2,0.1), class="nu")

fit_logistic_spline <- brm(bf(y ~ phi1/(1+exp((phi2-time)/phi3)),
                       sigma~s(time,by=treatment),
                       phi1 + phi2 + phi3 ~ 0+treatment,
                       autocor = ~arma(~time|sample:treatment,1,1),nl = TRUE),
                    family = student, prior = prior_logistic, data = df, iter = 1000,
                    cores = 2, chains = 2, backend = "cmdstanr", #threads = threading(4),
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

fit_mono_spline <- brm(bf(y ~ a-a*exp(-b*time), 
                   sigma~s(time,by=treatment), 
                   a + b ~ 0+treatment, 
                   autocor = ~arma(~time|sample:treatment,1,1),nl = TRUE),
                family = student, prior = priorMono, data = df, iter = 1000, 
                cores = 2, chains = 2, backend = "cmdstanr", #threads = threading(4),
                control = list(adapt_delta = 0.999,max_treedepth = 20),
                inits = function(){list(b_a=rgamma(2,1),b_b=rgamma(2,1),b_c=rgamma(2,1))})

# Review model objects

fitLinear_spline
fit_logistic_spline
fit_mono_spline
fit_quad_cluster




load("fit_quad_cluster.rdata")
fit_quad_cluster
load("spline_model_and_priors.rdata")
fit_spline


# make function to pull out sigma components

getHetSked<-function(model){
  
  string<-model$formula[[2]]$sigma # pull out sigma formula object
  formula<-paste0("sigma ~ ", as.character(string[3])) # pick the right hand side of formula
  
  p13<-model[13] #get the estimates
  dat<-as.data.frame(p13$fit)%>%
    select(contains("sigma"), -contains("Intercept"))%>% # splines will return a bunch of columns
    summarize( #grab the relevant columns, take the mean of each (should be one per treatment group)
      across(.cols=everything(),
             .fns = mean)
    )%>%t()%>%as.data.frame()%>%rownames_to_column() #flip it around to be more legible
  colnames(dat)<-c("Parameter", "Value")
  
  if (any(grepl("^zs_", colnames(as.data.frame(p13$fit))))){ #warn if there are splines
    warning("Spline components found, interpretation is complex: ", formula)
    # warning("Spline components found, Rethink your life")
  }
  outList<-list(formula, dat)
  return(outList)
}

splineHetSked<-getHetSked(fit_spline)
splineHetSked[[1]]
head(splineHetSked[[2]])
getHetSked(fit_quad_cluster)


# making splines to model heteroskedasticity using s() from brms, which calls s() from mgcv.
# the zs_ columns are from a matrix made by mgcv functions smooth2random( smoothCon( SMOOTH TERMS ) )
# but I have not been able to understand how to combine the smoothed terms to make some function (linear or GAM)
#       to return a measure of the heteroskedasticity at timepoint x











## Small Reprex

library(brms)
library(tidyverse)

logisticData <- function(x,phi1,phi2,phi3){
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
                 function(i) data.frame("sample"=paste0("sample_",i),"treatment"="a","time"=x,"y"=logisticData(x,200,13,3), stringsAsFactors = F))),
  do.call(rbind,
          lapply(1:20,
                 function(i) data.frame("sample"=paste0("sample_",i),"treatment"="b","time"=x,"y"=logisticData(x,160,13,3.5), stringsAsFactors = F)))
)

prior_logistic <- prior(lognormal(log(130), .25),nlpar = "phi1",coef="treatmenta") +
  prior(lognormal(log(130), .25),nlpar = "phi1",coef="treatmentb") +
  prior(lognormal(log(12), .25), nlpar = "phi2",coef="treatmenta") +
  prior(lognormal(log(12), .25), nlpar = "phi2",coef="treatmentb") +
  prior(lognormal(log(3), .25), nlpar = "phi3",coef="treatmenta") +
  prior(lognormal(log(3), .25), nlpar = "phi3",coef="treatmentb") +
  prior(student_t(3,0,5), dpar="sigma") +
  prior(gamma(2,0.1), class="nu")

fit_logistic_spline <- brm(bf(y ~ phi1/(1+exp((phi2-time)/phi3)),
                              sigma~s(time,by=treatment),
                              phi1 + phi2 + phi3 ~ 0+treatment,
                              autocor = ~arma(~time|sample:treatment,1,1),nl = TRUE),
                           family = student, prior = prior_logistic, data = df, iter = 1000,
                           cores = 2, chains = 2, backend = "cmdstanr",
                           control = list(adapt_delta = 0.999,max_treedepth = 20),
                           inits = function(){list(b_phi1=rgamma(2,1),b_phi2=rgamma(2,1),b_phi3=rgamma(2,1))})

p13<-fit_logistic_spline[13] 
as.data.frame(p13$fit)%>%
  select(contains("sigma"), -contains("Intercept"))%>% # splines return a bunch of columns
  summarize(across(.cols=everything(), .fns = mean))

p13$fit #8 rows for s_sigma_stimetreatmenta_1
        # 8 rows for s_sigma_stimetreatmentb_1
        # 8 rows for zs_sigma_1_1
        # 8 rows for zs_sigma_2_1
# all these have to be coming from the matrices out of mgcv, but getting 8 columns is not clearw

