library(brms)
library(ggplot2)
library(stringr)
library(patchwork)
library(cowplot)
library(plyr)
library(tidyverse)
library(bayesplot)
library(viridisLite)

#******************************************************************************************
# Could have various growthSim functions to make simulated growth data of different types
#******************************************************************************************

growthSimLogistic <- function(x,phi1,phi2,phi3){
  phi1_r <- phi1+rnorm(1,mean = 0,sd=25)
  phi2_r <- phi2+rnorm(1,mean=0,sd=1)
  phi3_r <- phi3+rnorm(1,mean=0,sd=.2)
  return(phi1_r/(1+exp((phi2_r-x)/phi3_r)))
}

growthSimGompertz <- function(x,a,b,c){
  a_r <- a+rnorm(1,mean = 0,sd=10)
  b_r <- b+rnorm(1,mean=0,sd=2)
  c_r <- c+rnorm(1,mean=0,sd=.035)
  return(a_r*exp(-b_r*exp(-c_r*x)))
}

growthSimMM <- function(x,a,b){
  a_r <- a+rnorm(1,mean = 0,sd=10)
  b_r <- b+rnorm(1,mean=0,sd=0.025)
  return(a_r-a_r*exp(-b_r*x))
}

growthSimExp <- function(x,a,b){
  a_r <- a+rnorm(1,mean = 0,sd=1.5)
  b_r <- b+rnorm(1,mean=0,sd=0.01)
  return(a_r*exp(b_r * x))
}

growthSimPL <- function(x,a,b){
  a_r <- a+rnorm(1,mean = 0,sd=2)
  b_r <- b+rnorm(1,mean=0,sd=0.075)
  return(a_r * x^(b_r)) # b_r/c_r < 1
}

growthSimL <- function(x,a){
  a_r <- a+rnorm(1,mean = 0,sd=0.2)
  return(a_r*x)
}

#******************************************************************************************
# Make a set of data to test 
#******************************************************************************************

x<-1:25
df<-rbind(
  do.call(rbind,
          lapply(1:20,
                 function(i) data.frame("sample"=paste0("sample_",i),
                                        "treatment"="a",
                                        "time"=x,
                                        "y"=growthSimLogistic(x, 200, 13, 3),
                                        stringsAsFactors = F))),
  do.call(rbind,
          lapply(1:20,
                 function(i) data.frame("sample"=paste0("sample_",i),
                                        "treatment"="b",
                                        "time"=x,
                                        "y"=growthSimLogistic(x,160, 13, 3.5),
                                        stringsAsFactors = F))))


#******************************************************************************************
# All this marked with smaller dividers would end up being the growthModelApproximation function
#******************************************************************************************

#**************************************
# Try fitting models based on each kind of growth to whichever data is being used
#**************************************
# Logistic
# Gompertz
# Power Law
# Monomolecular
# Exponential
# Linear
# Polynomial (picking degree could be problematic, I am having trouble with how this would work)

#possibleModels<-c("logistic", "gompertz", "powerlaw", "monomolecular", "exponential", "linear")

getGrowthModel<-function(dat=df, possibleModels = c("logistic", "gompertz", "powerlaw", "monomolecular", "exponential", "linear"), 
         comparisonIterations = 1, iterations = 1, sigma = "linear", timeVar = "time"){

if(!any(grepl(timeVar, colnames(df)))){
  renameText<-paste0("dat <- dat%>%rename(time=",timeVar,")")
  eval(parse(text=renameText))
}
  
for (model in possibleModels){

  if(model=="logistic"){
    specificFormula = "y ~ a/(1+exp((b-time)/c)),"
    priorList <- prior(lognormal(log(130), .25),nlpar = "a",coef="treatmenta") +
      prior(lognormal(log(130), .25),nlpar = "a",coef="treatmentb") + # need to be reallly thick priors if we have no idea?
      prior(lognormal(log(12), .25), nlpar = "b",coef="treatmenta") +
      prior(lognormal(log(12), .25), nlpar = "b",coef="treatmentb") +
      prior(lognormal(log(3), .25), nlpar = "c",coef="treatmenta") +
      prior(lognormal(log(3), .25), nlpar = "c",coef="treatmentb") +
      prior(student_t(3,0,5), dpar="sigma") +
      prior(gamma(2,0.1), class="nu")
    paramsToEstimate = "a + b + c"
    INITS_LIST<-list(b_a=rgamma(2,1),b_b=rgamma(2,1),b_c=rgamma(2,1))
  } else if(model=="gompertz"){
    specificFormula = "y ~ a*exp(-b*exp(-c*time)),"
    priorList <- prior(lognormal(log(130), .25),nlpar = "a") +
      prior(lognormal(log(12), .25), nlpar = "b") + 
      prior(lognormal(log(1.2), .25), nlpar = "c") + 
      prior(student_t(3,0,5), dpar="sigma") +
      prior(gamma(2,0.1), class="nu")
    paramsToEstimate = "a + b + c"
    INITS_LIST<-list(b_a=rgamma(2,1),b_b=rgamma(2,1),b_c=rgamma(2,1))
  } else if(model=="powerlaw"){
    specificFormula = "y ~ a* time^b,"
    priorList <- prior(lognormal(log(13), .25),nlpar = "a") +
      prior(lognormal(log(2), .25), nlpar = "b") + 
      prior(student_t(3,0,5), dpar="sigma") +
      prior(gamma(2,0.1), class="nu")
    paramsToEstimate = "a + b"
    INITS_LIST<-list(b_a=rgamma(2,1),b_b=rgamma(2,1))
  } else if(model=="monomolecular"){
    specificFormula = "y ~ a-a*exp(-b*time),"
    priorList <- prior(lognormal(log(130), .25),nlpar = "a") +
      prior(lognormal(log(2), .25), nlpar = "b") + 
      prior(student_t(3,0,5), dpar="sigma") +
      prior(gamma(2,0.1), class="nu")
    paramsToEstimate = "a + b"
    INITS_LIST<-list(b_a=rgamma(2,1),b_b=rgamma(2,1))
  } else if(model=="exponential"){
    specificFormula = "y ~ a*exp(b*time),"
    priorList <- prior(lognormal(log(15), .25),nlpar = "a") +
      prior(lognormal(log(0.1), .1), nlpar = "b") + 
      prior(student_t(3,0,5), dpar="sigma") +
      prior(gamma(2,0.1), class="nu")
    paramsToEstimate = "a + b"
    INITS_LIST<-list(b_a=rgamma(2,1),b_b=rgamma(2,1))
  } else if(model=="linear"){
    specificFormula = "y ~ a*time,"
    priorList <- prior(lognormal(log(1), .25),nlpar = "a") +
      prior(student_t(3,0,5), dpar="sigma") +
      prior(gamma(2,0.1), class="nu")
    paramsToEstimate ="a"
    INITS_LIST<-list(b_a=rgamma(2,1))
  } else {
    warning('Warning: "', model,'" is not recognized.\n Model options must be from: "logistic", "gompertz", "powerlaw", "monomolecular", "exponential", "linear"')
    break
    }
  
  bfText<-paste0("bf(", specificFormula, "sigma~time:treatment, ",paramsToEstimate," ~ 0+treatment,autocor = ~arma(~time|sample:treatment,1,1),nl = TRUE)")
  bayesFormula<-eval(parse(text = bfText))
  for (i in 1:comparisonIterations){
    cat("\nStarting Iteration ", i, "/",comparisonIterations," for ", model, "\n")
    iteration_row<-data.frame(iteration = i, elpd_loo = NA,elpd_loo_se=NA, p_loo=NA, p_loo_se=NA, loo_IC=NA,loo_IC_se=NA) #store iteration number
    
    fit <- brm(bayesFormula,
                    family = student, prior = priorList, data = dat, iter = 1000,
                    cores = 2, chains = 2, backend = "cmdstanr",
                    control = list(adapt_delta = 0.999,max_treedepth = 20),
                    inits = function(){INITS_LIST})
    
    fit<-add_criterion(fit, "loo")
    iteration_row$elpd_loo[1]<-fit$criteria$loo$estimates[1,1]
    iteration_row$elpd_loo_se[1]<-fit$criteria$loo$estimates[1,2]
    iteration_row$p_loo[1]<-fit$criteria$loo$estimates[2,1]
    iteration_row$p_loo_se[1]<-fit$criteria$loo$estimates[2,2]
    iteration_row$loo_IC[1]<-fit$criteria$loo$estimates[3,1]
    iteration_row$loo_IC_se[1]<-fit$criteria$loo$estimates[3,2]
    fitSum<-summary(fit)
    fixed<-fitSum$fixed
    fixed <- cbind(rownames(fixed), data.frame(fixed, row.names=NULL))
    colnames(fixed)[1]="label"
    fitSumWide<-pivot_wider(fixed, names_from = label, values_from = colnames(fixed)[2:ncol(fixed)])
    fitSumWide<-dplyr::select(fitSumWide, contains("Estimate"), contains("Est.Error"))
    iteration_row<-cbind(iteration_row, fitSumWide)
    if(i==1){
      metrics_df<-iteration_row
    } else{
      metrics_df<-rbind.fill(metrics_df, iteration_row)
    }
    metrics_df_means<-metrics_df%>%select(contains("loo"))%>%summarize(across(.cols=everything(), .fns = mean))
    saveDataText<-paste0(model,"Metrics_df<-metrics_df_means")
    eval(parse(text = saveDataText))
  } #end loop on iterations
} # end loop on models

  # now need to do the comparison from each *Metrics_df
  # new for loop that goes through all the models and does stuff to the *Metrics_df object,
  # probably appending them all 
  
  for (i in 1:length(possibleModels)){
    if (i ==1){
      # make a combined_df object out of the first *Metrics_df
      startCombinedDataText<-paste0("combined_df<-",
                                    possibleModels[i],
                                    "Metrics_df%>%select(contains('loo'))%>%mutate(model='",possibleModels[i],"')")
      eval(parse(text = startCombinedDataText))
      # combined_df<-
    } else{
      # rbind next *Metrics_dfs to the combined_df object (only LOO things)
      mutateText<-paste0("append_df<-",possibleModels[i],"Metrics_df%>%select(contains('loo'))%>%mutate(model='",possibleModels[i],"')")
      eval(parse(text = mutateText))
      combined_df<-rbind(combined_df, append_df)
    }
  }
  combined_df$model <- with(combined_df, reorder(model, -loo_IC))
  
  looPlot<-ggplot(combined_df)+
    geom_col(aes(x=model, y=loo_IC, fill=model))+
    geom_segment(aes(x=model, xend=model, y=loo_IC, yend = loo_IC+loo_IC_se), size = 1.5)+
    geom_segment(aes(x=model, xend=model, y=loo_IC, yend = loo_IC-loo_IC_se), size = 1.5)+
    geom_text(aes(x=model, y=loo_IC*1.05, label = paste0(round(loo_IC, digits=0))))+
    labs(title = "LOO IC to estimate best model", x="", y="", fill="Model")+
    scale_fill_viridis(discrete = T)+
    theme_minimal() +
    theme(axis.line.y.left = element_line(),
          axis.line.x.bottom = element_line(),
          axis.text.x.bottom = element_blank())
  
  bestModel<-combined_df%>%
    slice_min(order_by=loo_IC, n=1)%>%
    pull(model)
  print(bestModel)
  bestModel<-bestModel[1]
  print(bestModel)
  # run a new brms model for "bestModel"
  # these priors should probably be really thick since we have very little idea what is happening
  #     unless I make a smarter function that grabs the parameters from the model that was from the 'bestModel' test...
  if(bestModel=="logistic"){
    specificFormula = "y ~ a/(1+exp((b-time)/c)),"
    priorList <- prior(lognormal(log(130), .25),nlpar = "a",coef="treatmenta") +
      prior(lognormal(log(130), .25),nlpar = "a",coef="treatmentb") +
      prior(lognormal(log(12), .25), nlpar = "b",coef="treatmenta") +
      prior(lognormal(log(12), .25), nlpar = "b",coef="treatmentb") +
      prior(lognormal(log(3), .25), nlpar = "c",coef="treatmenta") +
      prior(lognormal(log(3), .25), nlpar = "c",coef="treatmentb") +
      prior(student_t(3,0,5), dpar="sigma") +
      prior(gamma(2,0.1), class="nu")
    paramsToEstimate ="a + b + c"
    INITS_LIST<-list(b_a=rgamma(2,1),b_b=rgamma(2,1),b_c=rgamma(2,1))
  } else if(bestModel=="gompertz"){
    specificFormula = "y ~ a*exp(-b*exp(-c*time)),"
    priorList <- prior(lognormal(log(130), .25),nlpar = "a") +
      prior(lognormal(log(12), .25), nlpar = "b") + 
      prior(lognormal(log(1.2), .25), nlpar = "c") + 
      prior(student_t(3,0,5), dpar="sigma") +
      prior(gamma(2,0.1), class="nu")
    paramsToEstimate ="a + b + c"
    INITS_LIST<-list(b_a=rgamma(2,1),b_b=rgamma(2,1),b_c=rgamma(2,1))
  } else if(bestModel=="powerlaw"){
    specificFormula = "y ~ a* time^b,"
    priorList <- prior(lognormal(log(13), .25),nlpar = "a") +
      prior(lognormal(log(2), .25), nlpar = "b") + 
      prior(student_t(3,0,5), dpar="sigma") +
      prior(gamma(2,0.1), class="nu")
    paramsToEstimate ="a + b"
    INITS_LIST<-list(b_a=rgamma(2,1),b_b=rgamma(2,1))
  } else if(bestModel=="monomolecular"){
    specificFormula = "y ~ a-a*exp(-b*time),"
    priorList <- prior(lognormal(log(130), .25),nlpar = "a") +
      prior(lognormal(log(2), .25), nlpar = "b") + 
      prior(student_t(3,0,5), dpar="sigma") +
      prior(gamma(2,0.1), class="nu")
    paramsToEstimate ="a + b"
    INITS_LIST<-list(b_a=rgamma(2,1),b_b=rgamma(2,1))
  } else if(bestModel=="exponential"){
    specificFormula = "y ~ a*exp(b*time),"
    priorList <- prior(lognormal(log(15), .25),nlpar = "a") +
      prior(lognormal(log(0.1), .1), nlpar = "b") + 
      prior(student_t(3,0,5), dpar="sigma") +
      prior(gamma(2,0.1), class="nu")
    paramsToEstimate ="a + b"
    INITS_LIST<-list(b_a=rgamma(2,1),b_b=rgamma(2,1))
  } else if(bestModel=="linear"){
    specificFormula = "y ~ a*time,"
    priorList <- prior(lognormal(log(1), .25),nlpar = "a") +
      prior(student_t(3,0,5), dpar="sigma") +
      prior(gamma(2,0.1), class="nu")
    paramsToEstimate ="a"
    INITS_LIST<-list(b_a=rgamma(2,1))
  }
  
  print(specificFormula)
  print(INITS_LIST)
  print("Line 298")
  
  # if the time variable is not named time then this whole thing will just die.
  # Might want to rename some variable to time early on
  
  sigma_<-ifelse(sigma=="linear", "sigma~time+time:treatment,", 
                 ifelse(sigma=="spline", "sigma~s(time,by=treatment),",
                 ifelse(sigma=="quad", "lf(sigma~ time + timeSQ + time:treatment + timeSQ:treatment),",
                 ifelse(sigma=="exp", "sigma~time+time:treatment,", #Not evaluated, doesn't matter
                 ifelse(sigma=="none", "", paste0(sigma))))))
  bfText<-paste0("bf(", specificFormula, 
                 sigma_, paramsToEstimate, 
                 " ~ 0+treatment,autocor = ~arma(~time|sample:treatment,1,1), nl = TRUE)")
  
  print(bfText)
  print("Line 313")
  
  bayesFormula<-eval(parse(text = bfText))
  
  if (sigma=="quad"){
    prior_list<-prior_list+
      prior(student_t(3,0,5), dpar="sigma", coef="time") + 
      prior(student_t(3,0,5), dpar="sigma", coef="time:treatmentb") + 
      prior(student_t(3,0,5), dpar="sigma", coef="timeSQ:treatmentb") + 
      prior(student_t(3,0,5), dpar="sigma", coef="timeSQ")
  } else if(sigma=="exp"){
    prior_list <-prior_list +
      prior(student_t(1, 0, 5), nlpar="d") + 
      prior(student_t(1, 0, 5), nlpar="r")
    
    bfText<-paste0("bf(", specificFormula, ", ",paramsToEstimate, 
" ~ 0+treatment,autocor = ~arma(~time|sample:treatment,1,1), nl = TRUE)+nlf(sigma ~ d * exp(r * time))+lf(d + r ~ 0+treatment)")
    bayesFormula<-eval(parse(text = bfText))
  }
  
  bestFit <- brm(bayesFormula,
         family = student, prior = priorList, data = dat, iter = 4000,
         cores = 2, chains = 4, backend = "cmdstanr",
         control = list(adapt_delta = 0.999,max_treedepth = 20),
         inits = function(){INITS_LIST})
  
  outputList<-list(combined_df, looPlot, bestModel, bestFit)
  return(outputList)
  
} #end function

# dat=df, possibleModels = c("logistic", "gompertz", "powerlaw", "monomolecular", "exponential", "linear"), 
# comparisonIterations = 1, iterations = 1, sigma = "linear", timeVar = "time"

df<-rbind(
  do.call(rbind,
          lapply(1:20,
                 function(i) data.frame("sample"=paste0("sample_",i),
                                        "treatment"="a",
                                        "time"=x,
                                        "y"=growthSimLogistic(x, 200, 13, 3),
                                        stringsAsFactors = F))),
  do.call(rbind,
          lapply(1:20,
                 function(i) data.frame("sample"=paste0("sample_",i),
                                        "treatment"="b",
                                        "time"=x,
                                        "y"=growthSimLogistic(x,160, 13, 3.5),
                                        stringsAsFactors = F))))

df_dap<-rbind(
  do.call(rbind,
          lapply(1:20,
                 function(i) data.frame("sample"=paste0("sample_",i),
                                        "treatment"="a",
                                        "DAP"=x,
                                        "y"=growthSimLogistic(x, 200, 13, 3),
                                        stringsAsFactors = F))),
  do.call(rbind,
          lapply(1:20,
                 function(i) data.frame("sample"=paste0("sample_",i),
                                        "treatment"="b",
                                        "DAP"=x,
                                        "y"=growthSimLogistic(x,160, 13, 3.5),
                                        stringsAsFactors = F))))
#list(combined_df, looPlot, bestModel, bestFit)
growthModelOutput_from_df<-getGrowthModel(dat=df, 
               possibleModels = c("logistic", "gompertz","linear"),
               comparisonIterations = 1, iterations = 1, 
               sigma = "linear", timeVar = "time")

growthModelOutput_from_df[1]
growthModelOutput_from_df[4]


#**************************************
# Compare all fit models, rank them by LOO IC or WAIC 
#**************************************

#**************************************
# Return best model, comparisons of models, label for best growth rate formula estimate,
# Estimates of coefficients for that formula, and how about a nice ribbon plot too
#**************************************





