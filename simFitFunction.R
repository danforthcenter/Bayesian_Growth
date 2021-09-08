library(brms)
library(ggplot2)
library(stringr)
library(patchwork)
library(cowplot)
library(tidyverse) #sorry, I'll use it sparingly

growthSim <- function(x,phi1,phi2,phi3){ # make function to do a bunch of growthSims and check how well they 
  phi1_r <- phi1+rnorm(1,mean = 0,sd=25) #can be modeled/how well their phi params can be recovered.
  phi2_r <- phi2+rnorm(1,mean=0,sd=1)
  phi3_r <- phi3+rnorm(1,mean=0,sd=.2)
  return(phi1_r/(1+exp((phi2_r-x)/phi3_r)))
}

#function to make Data (only to make the next part easier to read)
makeData<-function(x.=x,nSamples.=nSamples, phi1_1.=phi1_a, phi2_1.=phi2_a, phi3_1.=phi3_a, phi1_2.=phi1_b, phi2_2.=phi2_b, phi3_2.=phi3_b){
  df<-rbind(
    do.call(rbind,
            lapply(1:nSamples.,
                   function(i) data.frame("sample"=paste0("sample_",i),
                                          "treatment"="a",
                                          "time"=x.,
                                          "y"=growthSim(x., phi1_1., phi2_1., phi3_1.),
                                          stringsAsFactors = F))),
    do.call(rbind,
            lapply(1:nSamples.,
                   function(i) data.frame("sample"=paste0("sample_",i),
                                          "treatment"="b",
                                          "time"=x.,
                                          "y"=growthSim(x.,phi1_2., phi2_2., phi3_2.),
                                          stringsAsFactors = F)))
  )
  return(df)
}


# define function that makes (however many) iterations of df, using provided parameters.
modelSims<-function(iterations = 10, xTime=25, nSamples = 20, phi1_a=200, phi2_a=13, phi3_a=3, phi1_b=160, phi2_b=13, phi3_b=3.5){
  #loop through iterations making data, model, and gathering results
  for (i in 1:iterations){
    cat("\nStarting Iteration ", i, "/",iterations,"\n")
  iteration_row<-data.frame(iteration = i, elpd_loo = NA,elpd_loo_se=NA, p_loo=NA, p_loo_se=NA, loo_IC=NA,loo_IC_se=NA) #store iteration number
  x<-1:xTime
  dat <- makeData(phi1_1.=phi1_a, phi2_1.=phi2_a, phi3_1.=phi3_a, phi1_2.=phi1_b, phi2_2.=phi2_b, phi3_2.=phi3_b) #make the data
        # model the data
  prior_none <- prior(lognormal(log(130), .25),nlpar = "phi1",coef="treatmenta") +
    prior(lognormal(log(130), .25),nlpar = "phi1",coef="treatmentb") +
    prior(lognormal(log(12), .25), nlpar = "phi2",coef="treatmenta") +
    prior(lognormal(log(12), .25), nlpar = "phi2",coef="treatmentb") +
    prior(lognormal(log(3), .25), nlpar = "phi3",coef="treatmenta") +
    prior(lognormal(log(3), .25), nlpar = "phi3",coef="treatmentb") +
    prior(gamma(2,0.1), class="nu")
  
  fit_none <- brm(bf(y ~ phi1/(1+exp((phi2-time)/phi3)),
                     phi1 + phi2 + phi3 ~ 0+treatment, # if we don't model sigma then homoscedasticity is assumed
                     autocor = ~arma(~time|sample:treatment,1,1),nl = TRUE),
                  family = student, prior = prior_none, data = df, iter = 1000,
                  cores = 2, chains = 2, backend = "cmdstanr", #threads = threading(4),
                  control = list(adapt_delta = 0.999,max_treedepth = 20),
                  inits = function(){list(b_phi1=rgamma(2,1),b_phi2=rgamma(2,1),b_phi3=rgamma(2,1))})
        #get model metrics
  fit_none<-add_criterion(fit_none, "loo")
        #store iteration info in dataframe row
  iteration_row$elpd_loo[1]<-fit_none$criteria$loo$estimates[1,1]
  iteration_row$elpd_loo_se[1]<-fit_none$criteria$loo$estimates[1,2]
  iteration_row$p_loo[1]<-fit_none$criteria$loo$estimates[2,1]
  iteration_row$p_loo_se[1]<-fit_none$criteria$loo$estimates[2,2]
  iteration_row$loo_IC[1]<-fit_none$criteria$loo$estimates[3,1]
  iteration_row$loo_IC_se[1]<-fit_none$criteria$loo$estimates[3,2]
  fitSum<-summary(fit_none)
  fixed<-fitSum$fixed
  fixed <- cbind(rownames(fixed), data.frame(fixed, row.names=NULL))
  colnames(fixed)[1]="label"
  fitSumWide<-pivot_wider(fixed, names_from = label, values_from = colnames(fixed)[2:8])
  fitSumWide<-select(fitSumWide, contains("Estimate"), contains("Est.Error"))
  iteration_row<-cbind(iteration_row, fitSumWide)
        # append rows
  if(i==1){
    metrics_df<-iteration_row
  } else{
    metrics_df<-rbind(metrics_df, iteration_row)}
  }#close iterations for loop
  return(metrics_df)
}

start<-Sys.time()
metrics_df<-modelSims(iterations = 10)
Sys.time()-start

summ<-summarize(metrics_df, across(.cols=everything(), .fns = list(mean, sd), .names = "{.col}_{.fn}"))
colnames(summ)<-str_replace_all(str_replace_all(colnames(summ), "_1", "_mean"), "_2", "_sd")
View(summ)

r1<-select(summ, contains("_mean"))
colnames(r1)<-str_remove_all(colnames(r1), "_mean")
r2<- select(summ, contains("_sd"))
colnames(r2)<-str_remove_all(colnames(r2), "_sd")
summ2<-rbind(r1, r2)
row.names(summ2)<-c("mean", "sd")
summ2<-summ2%>%select(contains("looIC"), contains("phi"), contains("Estimate"))
View(summ2)




#would it be valuable to store the values we originally input and display them for reference/show the actual differences?




