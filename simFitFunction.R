library(brms)
library(ggplot2)
library(stringr)
library(patchwork)
library(cowplot)
library(plyr)
library(tidyverse) #sorry, I'll use it sparingly

growthSim <- function(x,phi1,phi2,phi3){
  phi1_r <- phi1+rnorm(1,mean = 0,sd=25)
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
modelSims<-function(iterations = 10, sigma = "none", xTime=25, nSamples = 20, phi1_a=200, phi2_a=13, phi3_a=3, phi1_b=160, phi2_b=13, phi3_b=3.5){
  sigma<-ifelse(sigma=="linear", "sigma~time+time:treatment,", 
                ifelse(sigma=="spline", "sigma~s(time,by=treatment),", 
                       ifelse(sigma=="none", "", 
                              paste0(sigma, ","))))
  bfText<-paste0("bf(y ~ phi1/(1+exp((phi2-time)/phi3)),", 
         sigma, 
         "phi1 + phi2 + phi3 ~ 0+treatment,autocor = ~arma(~time|sample:treatment,1,1),nl = TRUE)")
  prior_list <- prior(lognormal(log(130), .25),nlpar = "phi1",coef="treatmenta") +
    prior(lognormal(log(130), .25),nlpar = "phi1",coef="treatmentb") +
    prior(lognormal(log(12), .25), nlpar = "phi2",coef="treatmenta") +
    prior(lognormal(log(12), .25), nlpar = "phi2",coef="treatmentb") +
    prior(lognormal(log(3), .25), nlpar = "phi3",coef="treatmenta") +
    prior(lognormal(log(3), .25), nlpar = "phi3",coef="treatmentb") +
    prior(gamma(2,0.1), class="nu")
  if(sigma!="none"){
    prior_list<-prior_list+
      prior(student_t(3,0,5), dpar="sigma")}
  bayesFormula<-eval(parse(text = bfText))
  for (i in 1:iterations){ #loop through iterations making data, model, and gathering results
    cat("\nStarting Iteration ", i, "/",iterations,"\n") # print progress
  iteration_row<-data.frame(iteration = i, elpd_loo = NA,elpd_loo_se=NA, p_loo=NA, p_loo_se=NA, loo_IC=NA,loo_IC_se=NA) #store iteration number
  x<-1:xTime
  dat <- makeData(x. = x, nSamples.=nSamples, phi1_1.=phi1_a, phi2_1.=phi2_a, phi3_1.=phi3_a, phi1_2.=phi1_b, phi2_2.=phi2_b, phi3_2.=phi3_b) #make the data
        # model the data
  fit_none <- brm(bayesFormula,
                  family = student, prior = prior_list, data = dat, iter = 1000,
                  cores = 2, chains = 2, backend = "cmdstanr", #threads = threading(4),
                  control = list(adapt_delta = 0.999,max_treedepth = 20),
                  inits = function(){list(b_phi1=rgamma(2,1),b_phi2=rgamma(2,1),b_phi3=rgamma(2,1))})
        # get model metrics
  fit_none<-add_criterion(fit_none, "loo")
        # store iteration info in dataframe row
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
    metrics_df<-rbind.fill(metrics_df, iteration_row)}
  }#close iterations for loop
  summ<-summarize(metrics_df, across(.cols=everything(), .fns = list(mean, sd), .names = "{.col}_{.fn}"))
  colnames(summ)<-str_replace_all(str_replace_all(colnames(summ), "_1", "_mean"), "_2", "_sd")
  r1<-select(summ, contains("_mean"))
  colnames(r1)<-str_remove_all(colnames(r1), "_mean")
  r2<- select(summ, contains("_sd"))
  colnames(r2)<-str_remove_all(colnames(r2), "_sd")
  summ2<-rbind.fill(r1, r2)
  row.names(summ2)<-c("mean", "sd (between iterations)")
  summary_df<-summ2%>%select(contains("looIC"), contains("phi"), contains("Estimate"))
  
  plot<-metrics_df%>%
    dplyr::select(contains("Estimate"))%>%
    pivot_longer(cols = starts_with("Estimate"))%>%
    dplyr::mutate(name=factor(str_remove_all(name, "Estimate_")))%>%
    ggplot()+
    geom_boxplot(aes(x=value, y=name))+
    geom_point(aes(x=phi1_a, y=1), color="red", size=0.5)+
    geom_point(aes(x=phi1_b, y=2), color="red", size=0.5)+
    geom_point(aes(x=phi2_a, y=3), color="red", size=0.5)+
    geom_point(aes(x=phi2_b, y=4), color="red", size=0.5)+
    geom_point(aes(x=phi3_a, y=5), color="red", size=0.5)+
    geom_point(aes(x=phi3_b, y=6), color="red", size=0.5)+
    labs(title = paste0("modelSims Output on ", iterations, " iterations"), x="", y="")+
    theme_minimal() +
    theme(axis.line.y.left = element_line(),
          axis.line.x.bottom = element_line())

  phiOnePlot<-metrics_df%>%
    dplyr::select(contains("Estimate_phi1"))%>%
    pivot_longer(cols = starts_with("Estimate"))%>%
    dplyr::mutate(name=factor(str_remove_all(name, "Estimate_")))%>%
    ggplot()+
    geom_boxplot(aes(x=value, y=name))+
    geom_point(aes(x=phi1_a, y=1), color="red", size=0.5)+
    geom_point(aes(x=phi1_b, y=2), color="red", size=0.5)+
    annotate("text", x=phi1_a, y=0.55, label=paste0(phi1_a), color="red", size=3)+
    annotate("text", x=phi1_b, y=1.55, label=paste0(phi1_b), color="red", size=3)+
    annotate("text", x = median(metrics_df$Estimate_phi1_treatmenta),
             y = 1.45,
             label = paste("Est:", round(median(metrics_df$Estimate_phi1_treatmenta))),
             size = 3)+
    annotate("text", x = median(metrics_df$Estimate_phi1_treatmentb),
             y = 2.45,
             label = paste("Est:", round(median(metrics_df$Estimate_phi1_treatmentb))),
             size = 3)+
    labs(title = paste0("Phi 1 Estimation on ", iterations, " iterations"), x="", y="")+
    theme_minimal() +
    theme(axis.line.y.left = element_line(),
          axis.line.x.bottom = element_line())

  phiTwoPlot<-metrics_df%>%
    dplyr::select(contains("Estimate_phi2"))%>%
    pivot_longer(cols = starts_with("Estimate"))%>%
    dplyr::mutate(name=factor(str_remove_all(name, "Estimate_")))%>%
    ggplot()+
    geom_boxplot(aes(x=value, y=name))+
    geom_point(aes(x=phi2_a, y=1), color="red", size=0.5)+
    geom_point(aes(x=phi2_b, y=2), color="red", size=0.5)+
    annotate("text", x=phi2_a, y=0.55, label=paste0(phi2_a), color="red", size=3)+
    annotate("text", x=phi2_b, y=1.55, label=paste0(phi2_b), color="red", size=3)+
    annotate("text", x = median(metrics_df$Estimate_phi2_treatmenta),
             y = 1.45,
             label = paste("Est:", round(median(metrics_df$Estimate_phi2_treatmenta), digits=2)),
             size = 3)+
    annotate("text", x = median(metrics_df$Estimate_phi2_treatmentb),
             y = 2.45,
             label = paste("Est:", round(median(metrics_df$Estimate_phi2_treatmentb), digits=2)),
             size = 3)+
    labs(title = paste0("Phi 2 Estimation on ", iterations, " iterations"), x="", y="")+
    theme_minimal() +
    theme(axis.line.y.left = element_line(),
          axis.line.x.bottom = element_line())

  phiThreePlot<-metrics_df%>%
    dplyr::select(contains("Estimate_phi3"))%>%
    pivot_longer(cols = starts_with("Estimate"))%>%
    dplyr::mutate(name=factor(str_remove_all(name, "Estimate_")))%>%
    ggplot()+
    geom_boxplot(aes(x=value, y=name))+
    geom_point(aes(x=phi3_a, y=1), color="red", size=0.5)+
    geom_point(aes(x=phi3_b, y=2), color="red", size=0.5)+
    annotate("text", x=phi3_a, y=0.55, label=paste0(phi3_a), color="red", size=3)+
    annotate("text", x=phi3_b, y=1.55, label=paste0(phi3_b), color="red", size=3)+
    annotate("text", x = median(metrics_df$Estimate_phi3_treatmenta),
             y = 1.45,
             label = paste("Est:", round(median(metrics_df$Estimate_phi3_treatmenta), digits=2)),
             size = 3)+
    annotate("text", x = median(metrics_df$Estimate_phi3_treatmentb),
             y = 2.45,
             label = paste("Est:", round(median(metrics_df$Estimate_phi3_treatmentb), digits=2)),
             size = 3)+
    labs(title = paste0("Phi 3 Estimation on ", iterations, " iterations"), x="", y="")+
    theme_minimal() +
    theme(axis.line.y.left = element_line(),
          axis.line.x.bottom = element_line())

  
  modelSimsOutput<-list(metrics_df, summary_df, plot, phiOnePlot, phiTwoPlot, phiThreePlot)
  return(modelSimsOutput)
}

start<-Sys.time()
set.seed(123)
modelSimsOutput<-modelSims(iterations = 5, sigma="none", xTime=25, nSamples = 20, phi1_a=200, phi2_a=13, phi3_a=3, phi1_b=160, phi2_b=13, phi3_b=3.5)
Sys.time()-start

# Check results
modelSimsOutput[3]
modelSimsOutput[4]
modelSimsOutput[5]
modelSimsOutput[6]

# run different sigmas to compare
start<-Sys.time()
modelSimsOutput_None<-modelSims(iterations = 5, sigma="none", xTime=25, nSamples = 20, phi1_a=200, phi2_a=13, phi3_a=3, phi1_b=160, phi2_b=13, phi3_b=3.5)
Sys.time()-start
Sys.sleep(120)
start<-Sys.time()
modelSimsOutput_linear<-modelSims(iterations = 5, sigma="linear", xTime=25, nSamples = 20, phi1_a=200, phi2_a=13, phi3_a=3, phi1_b=160, phi2_b=13, phi3_b=3.5)
Sys.time()-start
Sys.sleep(120)
start<-Sys.time()
modelSimsOutput_Spline<-modelSims(iterations = 5, sigma="spline", xTime=25, nSamples = 20, phi1_a=200, phi2_a=13, phi3_a=3, phi1_b=160, phi2_b=13, phi3_b=3.5)
Sys.time()-start



save(modelSimsOutput_None,modelSimsOutput_linear,modelSimsOutput_Spline,file ="modelSimsOutputs.rdata")








# Comparing outputs functions

#grab Loo_IC from [1] of each object entered.
m1<-modelSimsOutput
m2<-modelSimsOutput

compareSimModels<-function(...){
argnames <- sys.call()
arguments<-list(...)
i =1
for (dat in arguments){
  name = as.character(rlang::enquo(dat))
  df<-as.data.frame(dat[1])
  if(i==1){ #initialize data binding 
    loo_IC_Mean<-mean(df$loo_IC)
    loo_IC_SE<-mean(df$loo_IC_se)
    # needs to get an identifier from the model, preferably the argument name as a string,
    # but for some reason that is eluding me right now
    name<-unlist(lapply(argnames[-1], as.character))[i]
    output_df<-data.frame(name=name, loo_IC_Mean= loo_IC_Mean, loo_IC_SE=loo_IC_SE)
  } else {  #bind data to existing rows
    loo_IC_Mean<-mean(df$loo_IC)
    loo_IC_SE<-mean(df$loo_IC_se)
    name<-unlist(lapply(argnames[-1], as.character))[i]
    output_df_new_row<-data.frame(name=name, loo_IC_Mean= loo_IC_Mean, loo_IC_SE=loo_IC_SE)
    output_df<-rbind(output_df, output_df_new_row)
  }
  i=i+1
}
looPlot<-ggplot(output_df)+
  geom_col(aes(x=name, y=loo_IC_Mean, fill=name))+
  geom_point(aes(x=name, y=loo_IC_Mean, color=name))+
  geom_point(aes(x=name, y=loo_IC_Mean+loo_IC_SE, color=name))+
  geom_point(aes(x=name, y=loo_IC_Mean-loo_IC_SE, color=name))+
  geom_segment(aes(x=name, x_end=name, y=loo_IC_Mean, y_end = loo_IC_Mean+loo_IC_SE))+
  geom_segment(aes(x=name, x_end=name, y=loo_IC_Mean, y_end = loo_IC_Mean-loo_IC_SE))+
  labs(title = "LOO IC", x="", y="")+
  theme_minimal() +
  theme(axis.line.y.left = element_line(),
        axis.line.x.bottom = element_line())
outputList<-list(output_df, looPlot)
return(outputList)
}
test<-compareSimModels(m1, m2)
head(test)











## testing with stored metrics
df<-read.csv("metrics_df_100.csv")
df%>%
  select(contains("Estimate"))%>%
  pivot_longer(cols = starts_with("Estimate"))%>%
  mutate(name=factor(str_remove_all(name, "Estimate_")))%>%
  ggplot()+
  geom_boxplot(aes(x=value, y=name))+
  geom_point(aes(x=phi1_a, y=1), color="red", size=0.5)+
  geom_point(aes(x=phi1_b, y=2), color="red", size=0.5)+
  geom_point(aes(x=phi2_a, y=3), color="red", size=0.5)+
  geom_point(aes(x=phi2_b, y=4), color="red", size=0.5)+
  geom_point(aes(x=phi3_a, y=5), color="red", size=0.5)+
  geom_point(aes(x=phi3_b, y=6), color="red", size=0.5)+
  labs(title = paste0("modelSims Output on ", iterations, " iterations"), x="", y="")+
  theme_minimal() +
  theme(axis.line.y.left = element_line(),
        axis.line.x.bottom = element_line())






















