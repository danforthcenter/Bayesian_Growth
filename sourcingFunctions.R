library(brms)
library(ggplot2)
library(stringr)
library(patchwork)
library(cowplot)
library(plyr)
library(tidyverse)
library(bayesplot)
library(viridis)
library(viridisLite)

## Get Growth Model Local
getGrowthModelLocal<-function(dat=df, possibleModels = c("logistic", "gompertz", "powerlaw", "monomolecular", "exponential", "linear"), 
                         comparisonIterations = 1, comparisonModelIterations = 1000, bestModelIterations=4000, comparisonSigma="linear", sigma = "spline", timeVar = "time"){
  
  if(timeVar!="time"){
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
    theme(#legend.position="right",
      axis.line.y.left = element_line(),
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
  
  outputList<-list(combined_df, looPlot, bestModel, bestFit, ribbonPlot, dat)
  return(outputList)
  
} #end function

## Get Growth Model

getGrowthModel<-function(dat=df, possibleModels = c("logistic", "gompertz", "powerlaw", "monomolecular", "exponential", "linear"), 
                              comparisonIterations = 1, comparisonModelIterations = 1000, 
                         bestModelIterations=4000, comparisonSigma="linear", sigma = "spline", 
                         timeVar = "time", treatmentVar="treatment", identityVar = "sample", outcomeVar="y",
                         controlTreatment = "b"){

  renameText<-paste0("dat <- dat%>%rename(time=",timeVar,", treatment = ",treatmentVar,", sample=",identityVar,", y=",outcomeVar,")")
  eval(parse(text=renameText))

  mutateText<-paste0('dat <- dat%>%mutate(treatment=ifelse(treatment=="',controlTreatment,'", "b", "a"))')
  eval(parse(text=renameText))
  
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
    theme(#legend.position="right",
      axis.line.y.left = element_line(),
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
                 cores = 4, chains = 4, backend = "cmdstanr",
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
  
  outputList<-list(combined_df, looPlot, bestModel, bestFit, ribbonPlot, dat)
  return(outputList)
  
} #end function


## Get Het Sked

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

## Get Growth Rate

growthRate<-function(model, timeVar="time"){
  string<-model$formula[[1]] # pull out formula object
  formula<-as.character(string[3]) # pick the right hand side of formula
  formulaNoPar<-str_replace_all(formula, "(\\(|\\))", " ") # get rid of parenthesis
  expForm<-str_extract(formulaNoPar, 
                       paste0("(", timeVar,"\\s*\\^\\s*\\S*)"))
  formulaTrimmed<-str_extract(formulaNoPar, # pull out time/stuff or time*stuff or stuff*time or stuff/time
                              paste0("(",timeVar, "\\s*(\\/\\s*\\S*|\\*\\s*\\S*))|(((\\S*\\s*\\*\\s*)|(\\S*\\s*\\/\\s?))",timeVar,")"))
  if(!is.na(expForm)){
    formulaTrimmed<-expForm
  }
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
