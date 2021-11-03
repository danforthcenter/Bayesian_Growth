library(brms)
library(stringr)
library(plyr)
library(tidyverse)


# make up some logistic data
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

# fit a model to our fake data (please skip this and load the data below, model fitting takes time)

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

# Load data

print(load("fit_quad_cluster.rdata"))
fit_quad_cluster
print(load("spline_model_and_priors.rdata"))
fit_spline


# Goal of this function is to pull out variance components

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
  colnames(dat)<-c("Parameter", "Value") # rename columns for ease of use
  
  if (any(grepl("^zs_", colnames(as.data.frame(p13$fit))))){ #warn if there are splines
    warning("Spline components found, interpretation is complex: ", formula) 
  }
  outList<-list(formula, dat) #return formula and dataframe
  return(outList)
}

# Looking at how the function does:

quadHetSked<-getHetSked(fit_quad_cluster)
quadHetSked[[1]] # easy to read formula: "sigma ~ time + timeSQ + time:treatment + timeSQ:treatment"
head(quadHetSked[[2]]) #data with estimates that are easy to interpret
# Parameter         Value
# 1              b_sigma_time  0.3071404958
# 2            b_sigma_timeSQ -0.0080115424
# 3   b_sigma_time:treatmentb -0.0058986760
# 4 b_sigma_timeSQ:treatmentb  0.0004361837



splineHetSked<-getHetSked(fit_spline)
splineHetSked[[1]] # only describes what was done for R's setup : "sigma ~ s(time, by = treatment)"
head(splineHetSked[[2]], n=10) #tons of stuff that I don't see a way to put back together
tail(splineHetSked[[2]], n=10)
# the square bracket[] numbers seem to correspond to columns of a matrix at some background point in mgcv.

# Parameter      Value
# 1   bs_sigma_stime:treatmenta_1  5.5918562
# 2   bs_sigma_stime:treatmentb_1  4.5650439
# 3   sds_sigma_stimetreatmenta_1  1.5118055
# 4   sds_sigma_stimetreatmentb_1  1.2089812
# 5  s_sigma_stimetreatmenta_1[1]  0.2669645
# 6  s_sigma_stimetreatmenta_1[2] -1.0377094
# 7  s_sigma_stimetreatmenta_1[3]  0.8254088
# 8  s_sigma_stimetreatmenta_1[4]  1.5527680
# 9  s_sigma_stimetreatmenta_1[5]  0.1404215
# 10 s_sigma_stimetreatmenta_1[6] -0.6004418

# Parameter       Value
# 27 zs_sigma_1_1[7] -0.01212021
# 28 zs_sigma_1_1[8]  0.24052911
# 29 zs_sigma_2_1[1]  0.11675898
# 30 zs_sigma_2_1[2] -0.84749670
# 31 zs_sigma_2_1[3]  0.28676063
# 32 zs_sigma_2_1[4]  1.02742871
# 33 zs_sigma_2_1[5] -0.10475658
# 34 zs_sigma_2_1[6] -0.32349652
# 35 zs_sigma_2_1[7]  0.03384145
# 36 zs_sigma_2_1[8]  0.20181178



# To recap, I'm making splines to model heteroskedasticity using s() from brms, which calls s() from mgcv.
# the zs_ columns are from a matrix made by mgcv functions smooth2random( smoothCon( SMOOTH_TERMS ) )
# but I have not been able to understand how to combine the smoothed terms to make some function (linear or GAM)
#       to return a measure of the heteroskedasticity at timepoint x
# relevant brms source code is on line 145 https://github.com/paul-buerkner/brms/blob/master/R/data-predictor.R


# What I would like to be able to do is put together a usable formula for the smooth curve that
#     the thin plate regression splines ended up making by using these estimates 

#example outside of the function(){} for ease.

p13<-fit_spline[13] 
as.data.frame(p13$fit)%>%
  select(contains("sigma"), -contains("Intercept"))%>% # splines return a bunch of columns
  summarize(across(.cols=everything(), .fns = mean))%>%knitr::kable()

p13$fit #8 rows for s_sigma_stimetreatmenta_1
        # 8 rows for s_sigma_stimetreatmentb_1
        # 8 rows for zs_sigma_1_1
        # 8 rows for zs_sigma_2_1
# all these have to be coming from the matrices out of mgcv, but getting 8 columns is not clear.

