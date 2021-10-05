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

getGrowthModel<-function(dat=df, possibleModels = c("logistic", "gompertz", "powerlaw", "monomolecular", "exponential", "linear"), 
         comparisonIterations = 1, comparisonModelIterations = 1000, bestModelIterations=4000, comparisonSigma="linear", sigma = "spline", timeVar = "time"){

if(!any(grepl(timeVar, colnames(df)))){
  renameText<-paste0("dat <- dat%>%rename(time=",timeVar,")")
  eval(parse(text=renameText))
}
  
for (model in possibleModels){
  #**************************************
  # Define priors and formulas for this iteration's model
  #**************************************
  if(model=="logistic"){ # need to all be reallly thick priors if we have no idea?
    specificFormula = "y ~ a/(1+exp((b-time)/c)),"
    priorList <- prior(lognormal(log(130), .25),nlpar = "a",coef="treatmenta") + # lognormal(log(100), .75)
      prior(lognormal(log(130), .25),nlpar = "a",coef="treatmentb") +  # lognormal(log(100), .75)
      prior(lognormal(log(12), .25), nlpar = "b",coef="treatmenta") + # lognormal(log(20), .75)
      prior(lognormal(log(12), .25), nlpar = "b",coef="treatmentb") + # lognormal(log(20), .75)
      prior(lognormal(log(3), .25), nlpar = "c",coef="treatmenta") + # lognormal(log(20), .75)
      prior(lognormal(log(3), .25), nlpar = "c",coef="treatmentb") + # lognormal(log(20), .75)
      prior(student_t(3,0,5), dpar="sigma") +
      prior(gamma(2,0.1), class="nu")
    paramsToEstimate = "a + b + c"
    INITS_LIST<-list(b_a=rgamma(2,1),b_b=rgamma(2,1),b_c=rgamma(2,1))
  } else if(model=="gompertz"){
    specificFormula = "y ~ a*exp(-b*exp(-c*time)),"
    priorList <- prior(lognormal(log(130), .25),nlpar = "a") + # lognormal(log(100), .75)
      prior(lognormal(log(12), .25), nlpar = "b") + # lognormal(log(20), .75)
      prior(lognormal(log(1.2), .25), nlpar = "c") + # lognormal(log(10), 1)
      prior(student_t(3,0,5), dpar="sigma") +
      prior(gamma(2,0.1), class="nu")
    paramsToEstimate = "a + b + c"
    INITS_LIST<-list(b_a=rgamma(2,1),b_b=rgamma(2,1),b_c=rgamma(2,1))
  } else if(model=="powerlaw"){
    specificFormula = "y ~ a* time^b,"
    priorList <- prior(lognormal(log(13), .25),nlpar = "a") + # lognormal(log(20), .75)
      prior(lognormal(log(2), .25), nlpar = "b") + # lognormal(log(10), 1)
      prior(student_t(3,0,5), dpar="sigma") +
      prior(gamma(2,0.1), class="nu")
    paramsToEstimate = "a + b"
    INITS_LIST<-list(b_a=rgamma(2,1),b_b=rgamma(2,1))
  } else if(model=="monomolecular"){
    specificFormula = "y ~ a-a*exp(-b*time),"
    priorList <- prior(lognormal(log(130), .25),nlpar = "a") + # lognormal(log(100), .75)
      prior(lognormal(log(2), .25), nlpar = "b") + # lognormal(log(10), 1)
      prior(student_t(3,0,5), dpar="sigma") +
      prior(gamma(2,0.1), class="nu")
    paramsToEstimate = "a + b"
    INITS_LIST<-list(b_a=rgamma(2,1),b_b=rgamma(2,1))
  } else if(model=="exponential"){
    specificFormula = "y ~ a*exp(b*time),"
    priorList <- prior(lognormal(log(15), .25),nlpar = "a") + # lognormal(log(20), .75)
      prior(lognormal(log(0.1), .1), nlpar = "b") + # lognormal(log(10), 1)
      prior(student_t(3,0,5), dpar="sigma") +
      prior(gamma(2,0.1), class="nu")
    paramsToEstimate = "a + b"
    INITS_LIST<-list(b_a=rgamma(2,1),b_b=rgamma(2,1))
  } else if(model=="linear"){
    specificFormula = "y ~ a*time,"
    priorList <- prior(lognormal(log(1), .25),nlpar = "a") + # lognormal(log(10), 1)
      prior(student_t(3,0,5), dpar="sigma") +
      prior(gamma(2,0.1), class="nu")
    paramsToEstimate ="a"
    INITS_LIST<-list(b_a=rgamma(2,1))
  } else {
    warning('Warning: "', model,'" is not recognized.\n Model options are: "logistic", "gompertz", "powerlaw", "monomolecular", "exponential", "linear"')
    break
  }
  #**************************************
  # Apply comparisonSigma
  #**************************************
  
  sigma_<-ifelse(comparisonSigma=="linear", "sigma~time+time:treatment,", 
                 ifelse(comparisonSigma=="spline", "sigma~s(time,by=treatment),",
                        ifelse(comparisonSigma=="quad", "lf(sigma~ time + timeSQ + time:treatment + timeSQ:treatment),",
                               ifelse(comparisonSigma=="exp", "sigma~time+time:treatment,", 
                                      ifelse(comparisonSigma=="none", "", paste0(comparisonSigma))))))
  bfText<-paste0("bf(", specificFormula, 
                 sigma_, paramsToEstimate, 
                 " ~ 0+treatment,autocor = ~arma(~time|sample:treatment,1,1), nl = TRUE)")
  
  bayesFormula<-eval(parse(text = bfText))
  
  if (comparisonSigma=="quad"){
    prior_list<-prior_list+
      prior(student_t(3,0,5), dpar="sigma", coef="time") + 
      prior(student_t(3,0,5), dpar="sigma", coef="time:treatmentb") + 
      prior(student_t(3,0,5), dpar="sigma", coef="timeSQ:treatmentb") + 
      prior(student_t(3,0,5), dpar="sigma", coef="timeSQ")
  } else if(comparisonSigma=="exp"){
    prior_list <-prior_list +
      prior(student_t(1, 0, 5), nlpar="d") + 
      prior(student_t(1, 0, 5), nlpar="r")
    
    bfText<-paste0("bf(", specificFormula, ", ",paramsToEstimate, 
                   " ~ 0+treatment,autocor = ~arma(~time|sample:treatment,1,1), nl = TRUE)+nlf(sigma ~ d * exp(r * time))+lf(d + r ~ 0+treatment)")
    bayesFormula<-eval(parse(text = bfText))
  }
  
  #**************************************
  # Try fitting models based on each kind of growth to whichever data is being used
  #**************************************
  bfText<-paste0("bf(", specificFormula, "sigma~time:treatment, ",paramsToEstimate," ~ 0+treatment,autocor = ~arma(~time|sample:treatment,1,1),nl = TRUE)")
  bayesFormula<-eval(parse(text = bfText))
  for (i in 1:comparisonIterations){
    cat("\nStarting Iteration ", i, "/",comparisonIterations," for ", model, "\n")
    iteration_row<-data.frame(iteration = i, elpd_loo = NA,elpd_loo_se=NA, p_loo=NA, p_loo_se=NA, loo_IC=NA,loo_IC_se=NA) #store iteration number
    
    fit <- brm(bayesFormula,
                    family = student, prior = priorList, data = dat, iter = comparisonModelIterations,
                    cores = 2, chains = 2, backend = "cmdstanr",
                    control = list(adapt_delta = 0.999,max_treedepth = 20),
                    inits = function(){INITS_LIST})
    #**************************************
    # Get model metrics
    #**************************************
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
  } #end loop on iterations
    #**************************************
    # Save metrics for this model under a new name
    #**************************************
    metrics_df_means<-metrics_df%>%select(contains("loo"))%>%summarize(across(.cols=everything(), .fns = mean))
    saveDataText<-paste0(model,"Metrics_df<-metrics_df_means")
    eval(parse(text = saveDataText))
    
    # saveModelText<-paste0(model,"_mod<-fit")
    # eval(parse(text = saveModelText))
    # 
    # retrieveBestModelText<-paste0("oldBestFit<-",bestModel,"_mod")
    # eval(parse(text=retrieveBestModelText))
    modelEstimates <- as.data.frame(fit$fit)%>%
      select(contains("b_"))%>%
      summarize(across(.cols=everything(), .fns=mean))
    
    storeModelEstimatesText<-paste0(model,"ModelEstimates<-modelEstimates")
    eval(parse(text=storeModelEstimatesText))
    
} # end loop on models

  #**************************************
  # combined metrics from all models into one dataframe
  #**************************************
  
  for (i in 1:length(possibleModels)){
    if (i ==1){
      startCombinedDataText<-paste0("combined_df<-",
                                    possibleModels[i],
                                    "Metrics_df%>%select(contains('loo'))%>%mutate(model='",possibleModels[i],"')")
      eval(parse(text = startCombinedDataText))
    } else{
      mutateText<-paste0("append_df<-",possibleModels[i],"Metrics_df%>%select(contains('loo'))%>%mutate(model='",possibleModels[i],"')")
      eval(parse(text = mutateText))
      combined_df<-rbind(combined_df, append_df)
    }
  }
  combined_df$model <- with(combined_df, reorder(model, -loo_IC))
  #**************************************
  # Plot metrics from all models
  #**************************************
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
  #**************************************
  # Select best model
  #**************************************
  bestModel<-combined_df%>%
    slice_min(order_by=loo_IC, n=1)%>%
    pull(model)

  bestModel<-bestModel[1]
  
  #**************************************
  # Run new model for the best type of model
  #**************************************  
  # retrieveBestModelText<-paste0("oldBestFit<-",bestModel,"_mod")
  # eval(parse(text=retrieveBestModelText))
  # oldEstimates <- as.data.frame(oldBestFit$fit)%>%select(contains("b_"))%>%summarize(across(.cols=everything(), .fns=mean))
  
  if(bestModel=="logistic"){
    
    # A<-logisticModelEstimates%>%select(contains("b_a"))%>%t()%>%mean() #need to either paste these together
    # B<-logisticModelEstimates%>%select(contains("b_b"))%>%t()%>%mean() # or just use the numbers I had earlier
    # C<-logisticModelEstimates%>%select(contains("b_c"))%>%t()%>%mean() #these are not evaluated here, so it 
    specificFormula = "y ~ a/(1+exp((b-time)/c))," #                    gets messed up (A isn't in C++?)
    
    #priorText<-paste0('priorList<-prior(lognormal(log(',A,'), .25),nlpar="a" + prior(lognormal(log(',B,'), .25), nlpar = "b") + prior(lognormal(log(',C,'), .25), nlpar = "c") + prior(student_t(3,0,5), dpar="sigma") +prior(gamma(2,0.1), class="nu")')
    #eval(parse(text=priorText))
    priorList <- prior(lognormal(log(130), .25),nlpar = "a") +
      prior(lognormal(log(12), .25), nlpar = "b") +
      prior(lognormal(log(3), .25), nlpar = "c") +
      prior(student_t(3,0,5), dpar="sigma") +
      prior(gamma(2,0.1), class="nu")
    paramsToEstimate ="a + b + c"
    INITS_LIST<-list(b_a=rgamma(2,1),b_b=rgamma(2,1),b_c=rgamma(2,1))
  } else if(bestModel=="gompertz"){
    # A<-gompertzModelEstimates%>%select(contains("b_a"))%>%t()%>%mean()
    # B<-gompertzModelEstimates%>%select(contains("b_b"))%>%t()%>%mean()
    # C<-gompertzModelEstimates%>%select(contains("b_c"))%>%t()%>%mean()
    specificFormula = "y ~ a*exp(-b*exp(-c*time)),"
    
    # priorText<-paste0('priorList<-prior(lognormal(log(',A,'), .25),nlpar="a" +
    #                   prior(lognormal(log(',B,'), .25), nlpar = "b") + 
    #                   prior(lognormal(log(',C,'), .25), nlpar = "c") + 
    #                   prior(student_t(3,0,5), dpar="sigma") +
    #                   prior(gamma(2,0.1), class="nu")')
    # eval(parse(text=priorText))
    
    priorList <- prior(lognormal(log(130), .25),nlpar = "a") +
      prior(lognormal(log(12), .25), nlpar = "b") + 
      prior(lognormal(log(1.2), .25), nlpar = "c") + 
      prior(student_t(3,0,5), dpar="sigma") +
      prior(gamma(2,0.1), class="nu")
    paramsToEstimate ="a + b + c"
    INITS_LIST<-list(b_a=rgamma(2,1),b_b=rgamma(2,1),b_c=rgamma(2,1))
  } else if(bestModel=="powerlaw"){
    # A<-powerlawModelEstimates%>%select(contains("b_a"))%>%t()%>%mean()
    # B<-powerlawModelEstimates%>%select(contains("b_b"))%>%t()%>%mean()
    specificFormula = "y ~ a* time^b,"
    
    # priorText<-paste0('priorList<-prior(lognormal(log(',A,'), .25),nlpar="a" +
    #                   prior(lognormal(log(',B,'), .25), nlpar = "b") + 
    #                   prior(student_t(3,0,5), dpar="sigma") +
    #                   prior(gamma(2,0.1), class="nu")')
    # eval(parse(text=priorText))
    
    priorList <- prior(lognormal(log(13), .25),nlpar = "a") +
      prior(lognormal(log(2), .25), nlpar = "b") + 
      prior(student_t(3,0,5), dpar="sigma") +
      prior(gamma(2,0.1), class="nu")
    paramsToEstimate ="a + b"
    INITS_LIST<-list(b_a=rgamma(2,1),b_b=rgamma(2,1))
  } else if(bestModel=="monomolecular"){
    # A<-monomolecularModelEstimates%>%select(contains("b_a"))%>%t()%>%mean()
    # B<-monomolecularModelEstimates%>%select(contains("b_b"))%>%t()%>%mean()
    specificFormula = "y ~ a-a*exp(-b*time),"
    
    # priorText<-paste0('priorList<-prior(lognormal(log(',A,'), .25),nlpar="a" +
    #                   prior(lognormal(log(',B,'), .25), nlpar = "b") + 
    #                   prior(student_t(3,0,5), dpar="sigma") +
    #                   prior(gamma(2,0.1), class="nu")')
    # eval(parse(text=priorText))
    
    priorList <- prior(lognormal(log(130), .25),nlpar = "a") +
      prior(lognormal(log(2), .25), nlpar = "b") + 
      prior(student_t(3,0,5), dpar="sigma") +
      prior(gamma(2,0.1), class="nu")
    paramsToEstimate ="a + b"
    INITS_LIST<-list(b_a=rgamma(2,1),b_b=rgamma(2,1))
  } else if(bestModel=="exponential"){
    # A<-exponentialModelEstimates%>%select(contains("b_a"))%>%t()%>%mean()
    # B<-exponentialModelEstimates%>%select(contains("b_b"))%>%t()%>%mean()
    specificFormula = "y ~ a*exp(b*time),"
    
    # priorText<-paste0('priorList<-prior(lognormal(log(',A,'), .25),nlpar="a" +
    #                   prior(lognormal(log(',B,'), .25), nlpar = "b") + 
    #                   prior(student_t(3,0,5), dpar="sigma") +
    #                   prior(gamma(2,0.1), class="nu"')
    # eval(parse(text=priorText))
    
    priorList <- prior(lognormal(log(15), .25),nlpar = "a") +
      prior(lognormal(log(.1), .1), nlpar = "b") + 
      prior(student_t(3,0,5), dpar="sigma") +
      prior(gamma(2,0.1), class="nu")
    paramsToEstimate ="a + b"
    INITS_LIST<-list(b_a=rgamma(2,1),b_b=rgamma(2,1))
  } else if(bestModel=="linear"){
    # A<-linearModelEstimates%>%select(contains("b_a"))%>%t()%>%mean()
    specificFormula = "y ~ a*time,"
    
    # priorText<-paste0('priorList<-prior(lognormal(log(',A,'), .25),nlpar="a" + prior(student_t(3,0,5), dpar="sigma") + prior(gamma(2,0.1), class="nu")')
    # eval(parse(text=priorText))
    
    priorList <- prior(lognormal(log(1), .25),nlpar = "a") +
      prior(student_t(3,0,5), dpar="sigma") +
      prior(gamma(2,0.1), class="nu")
    paramsToEstimate ="a"
    INITS_LIST<-list(b_a=rgamma(2,1))
  }
  
  sigma_<-ifelse(sigma=="linear", "sigma~time+time:treatment,", 
                 ifelse(sigma=="spline", "sigma~s(time,by=treatment),",
                 ifelse(sigma=="quad", "lf(sigma~ time + timeSQ + time:treatment + timeSQ:treatment),",
                 ifelse(sigma=="exp", "sigma~time+time:treatment,", #Not evaluated, doesn't matter
                 ifelse(sigma=="none", "", paste0(sigma))))))
  bfText<-paste0("bf(", specificFormula, 
                 sigma_, paramsToEstimate, 
                 " ~ 0+treatment,autocor = ~arma(~time|sample:treatment,1,1), nl = TRUE)")

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
  print("Compiling Best Fit Model")
  bestFit <- brm(bayesFormula,
         family = student, prior = priorList, data = dat, iter = bestModelIterations,
         cores = 2, chains = 4, backend = "cmdstanr",
         control = list(adapt_delta = 0.999,max_treedepth = 20),
         inits = function(){INITS_LIST})
  
  #**************************************
  # Ribbon graph of that model
  #**************************************  
  
  probs <- seq(from=99, to=1, by=-2)/100
  avg_pal <- plasma(n=length(probs))
  df_test <- rbind(data.frame("treatment"="a",time=1:25,sample="new1"),data.frame("treatment"="b",time=1:25,sample="new2"))
  df_pred <- predict(bestFit, df_test, probs=probs) # add fit_spline in place of fit1 for first checks.
  test <- cbind(df_test,df_pred)
  
  ribbonPlot <- ggplot(test,aes(time,Estimate))+
    facet_wrap(~treatment)+
    lapply(seq(1,49,2),function(i) geom_ribbon(aes_string(ymin=paste("Q",i,sep = ""),ymax=paste("Q",100-i,sep = "")),fill=avg_pal[i],alpha=0.5))+
    geom_line(data=dat,aes(time,y,group=interaction(treatment,sample)),color="gray20", size=0.3)+
    ylab("Area (cm, simulated)")+
    xlab("Time")+
    theme_light()+
    theme(axis.ticks.length=unit(0.2,"cm"))+
    theme(strip.background=element_rect(fill="gray50",color="gray20"),
          strip.text.x=element_text(size=14,color="white"),
          strip.text.y=element_text(size=14,color="white"))+
    theme(axis.title= element_text(size = 18))+
    theme(axis.text = element_text(size = 14))+
    theme(legend.position='top')
  
  outputList<-list(combined_df, looPlot, bestModel, bestFit, ribbonPlot)
  return(outputList)
  
} #end function

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
# list(combined_df, looPlot, bestModel, bestFit)
# c("logistic", "gompertz", "powerlaw", "monomolecular", "exponential", "linear")

growthModelOutput_from_df<-getGrowthModel(dat=df, 
               possibleModels = c("logistic", "gompertz","linear", "monomolecular"),
               comparisonIterations = 1, iterations = 1, 
               sigma = "linear", timeVar = "time")

growthModelOutput_from_df[1]
growthModelOutput_from_df[2]
growthModelOutput_from_df[3]
growthModelOutput_from_df[4]


growthModelOutput_from_df_all<-getGrowthModel(dat=df, 
                                          possibleModels = c("logistic", "gompertz", "monomolecular", "exponential", "linear"),
                                          comparisonIterations = 1, iterations = 1, 
                                          sigma = "linear", timeVar = "time")

growthModelOutput_from_df_all[1]
growthModelOutput_from_df_all[2]
growthModelOutput_from_df_all[3]
growthModelOutput_from_df_all[4]

growthModelOutput_from_df_logistic<-getGrowthModel(dat=df, 
                                              possibleModels = c("logistic", "gompertz", "monomolecular", "linear"),
                                              comparisonIterations = 1, iterations = 1, 
                                              sigma = "linear", timeVar = "time")

growthModelOutput_from_df_all_logistic[1]
growthModelOutput_from_df_all_logistic[2]
growthModelOutput_from_df_all_logistic[3]
growthModelOutput_from_df_all_logistic[4]
growthModelOutput_from_df_all_logistic[5]

#dat=df, possibleModels = c("logistic", "gompertz", "powerlaw", "monomolecular", "exponential", "linear"), 
#comparisonIterations = 1, iterations = 1, 
# comparisonModelIterations = 1000, bestModelIterations=4000, 
# comparisonSigma="linear", sigma = "spline", timeVar = "time"
  
growthModelOutput_MoreArgs<-getGrowthModel(dat=df, 
                                           possibleModels= c("logistic", "gompertz", "monomolecular", "linear"),
                                           comparisonIterations = 1,
                                           comparisonModelIterations = 1000, bestModelIterations=4000, 
                                           comparisonSigma="linear", sigma = "spline", timeVar = "time")

save(growthModelOutput_MoreArgs, file="growthModelOutput.rdata")

#**************************************
# Look at the data and see what you think
#**************************************


ggplot(df,aes(time,y,group=interaction(treatment,sample)))+
  geom_line(aes(color=treatment))+
  ylab("Area (cm, simulated)")+
  xlab("Time")+
  labs(title="Input Data")+
  theme_light()+
  theme(axis.ticks.length=unit(0.2,"cm"))+
  theme(strip.background=element_rect(fill="gray50",color="gray20"),
        strip.text.x=element_text(size=14,color="white"),
        strip.text.y=element_text(size=14,color="white"))+
  theme(axis.title= element_text(size = 18))+
  theme(title= element_text(size = 20))+
  theme(axis.text = element_text(size = 14))+
  theme(legend.position='top')


#**************************************
# playing with how to recapitulate priors
#**************************************


A = 180.145722
B = 13.0723891
C = 3.2576259
priorText<-paste0('prior_string(lognormal(log(',round(A),
                  '), .25),nlpar="a" + prior_string(lognormal(log(',round(B),
                  '), .25), nlpar = "b") + prior_string(lognormal(log(',round(C),
                  '), .25), nlpar="c") + prior_string(student_t(3,0,5), dpar="sigma") +prior(gamma(2,0.1), class="nu")')
priorText
priorList<-eval(parse(text=priorText))

fit <- brm(bf(y ~ a/(1+exp((b-time)/c)),
                     sigma~time+time:treatment, #linear [no s()] variation
                     a + b + c ~ 0+treatment,
                     autocor = ~arma(~time|sample:treatment,1,1),nl = TRUE),
                  family = student, prior = priorList, data = df, iter = 1000,
                  cores = 2, chains = 2, backend = "cmdstanr", 
                  control = list(adapt_delta = 0.999,max_treedepth = 20),
                  inits = function(){list(b_a=rgamma(2,1),b_b=rgamma(2,1),b_c=rgamma(2,1))})



