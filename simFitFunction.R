library(brms)
library(ggplot2)
library(stringr)
library(patchwork)
library(cowplot)
library(plyr)
library(tidyverse) #sorry, I'll use it sparingly
library(bayesplot)
library(viridisLite)

growthSim <- function(x,phi1,phi2,phi3){
  phi1_r <- phi1+rnorm(1,mean = 0,sd=25)
  phi2_r <- phi2+rnorm(1,mean=0,sd=1)
  phi3_r <- phi3+rnorm(1,mean=0,sd=.2)
  return(phi1_r/(1+exp((phi2_r-x)/phi3_r)))
}

################################## Make Data Function ################################## 
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

################################## define modelSims() ################################## 
modelSims<-function(iterations = 10, sigma = "none", xTime=25, nSamples = 20, phi1_a=200, phi2_a=13, phi3_a=3, phi1_b=160, phi2_b=13, phi3_b=3.5){
  sigma_<-ifelse(sigma=="linear", "sigma~time+time:treatment,", 
                ifelse(sigma=="spline", "sigma~s(time,by=treatment),",
                       ifelse(sigma=="quad", "lf(sigma~ time + timeSQ + time:treatment + timeSQ:treatment),",
                              ifelse(sigma=="exp", "sigma~time+time:treatment,", #Not evaluated, doesn't matter
                                     ifelse(sigma=="none", "", 
                                            paste0(sigma, ","))))))
  bfText<-paste0("bf(y ~ phi1/(1+exp((phi2-time)/phi3)),", 
         sigma_, 
         "phi1 + phi2 + phi3 ~ 0+treatment,autocor = ~arma(~time|sample:treatment,1,1),nl = TRUE)")
  bayesFormula<-eval(parse(text = bfText))
  
  prior_list <- prior(lognormal(log(130), .25),nlpar = "phi1",coef="treatmenta") +
    prior(lognormal(log(130), .25),nlpar = "phi1",coef="treatmentb") +
    prior(lognormal(log(12), .25), nlpar = "phi2",coef="treatmenta") +
    prior(lognormal(log(12), .25), nlpar = "phi2",coef="treatmentb") +
    prior(lognormal(log(3), .25), nlpar = "phi3",coef="treatmenta") +
    prior(lognormal(log(3), .25), nlpar = "phi3",coef="treatmentb") +
    prior(gamma(2,0.1), class="nu")
  if(sigma %in% c("linear", "spline")){
    prior_list<-prior_list+
      prior(student_t(3,0,5), dpar="sigma")
    } else if (sigma=="quad"){
        prior_list<-prior_list+
          prior(student_t(3,0,5), dpar="sigma", coef="time") + 
          prior(student_t(3,0,5), dpar="sigma", coef="time:treatmentb") + 
          prior(student_t(3,0,5), dpar="sigma", coef="timeSQ:treatmentb") + 
          prior(student_t(3,0,5), dpar="sigma", coef="timeSQ")
    } else if(sigma=="exp"){
      prior_list <-prior_list +
        prior(student_t(1, 0, 5), nlpar="a") + 
        prior(student_t(1, 0, 5), nlpar="r")
      bayesFormula<-bf(y ~ phi1/(1+exp((phi2-time)/(phi3))),
                       phi1 + phi2 + phi3 ~ 0+treatment,
                       autocor = ~arma(~time|sample:treatment,1,1), nl = TRUE)+
        nlf(sigma ~ a * exp(r * time))+
        lf(a + r ~ 0+treatment)
    }
  
  for (i in 1:iterations){
    cat("\nStarting Iteration ", i, "/",iterations,"\n") 
  iteration_row<-data.frame(iteration = i, elpd_loo = NA,elpd_loo_se=NA, p_loo=NA, p_loo_se=NA, loo_IC=NA,loo_IC_se=NA) #store iteration number
  x<-1:xTime
  dat <- makeData(x. = x, nSamples.=nSamples, phi1_1.=phi1_a, phi2_1.=phi2_a, phi3_1.=phi3_a, phi1_2.=phi1_b, phi2_2.=phi2_b, phi3_2.=phi3_b) #make the data
  if(sigma=="quad"){dat$timeSQ<-dat$time^2}
  fit_none <- brm(bayesFormula,
                  family = student, prior = prior_list, data = dat, iter = 4000,
                  cores = 2, chains = 4, backend = "cmdstanr",
                  control = list(adapt_delta = 0.999,max_treedepth = 20),
                  inits = function(){list(b_phi1=rgamma(2,1),b_phi2=rgamma(2,1),b_phi3=rgamma(2,1))})
  
  fit_none<-add_criterion(fit_none, "loo")
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
  fitSumWide<-pivot_wider(fixed, names_from = label, values_from = colnames(fixed)[2:ncol(fixed)])
  fitSumWide<-dplyr::select(fitSumWide, contains("Estimate"), contains("Est.Error"))
  iteration_row<-cbind(iteration_row, fitSumWide)
  if(i==1){
    metrics_df<-iteration_row
  } else{
    metrics_df<-rbind.fill(metrics_df, iteration_row)}
  }
  summ<-dplyr::summarize(metrics_df, across(.cols=everything(), .fns = list(mean, sd), .names = "{.col}_{.fn}"))
  colnames(summ)<-str_replace_all(str_replace_all(colnames(summ), "_1", "_mean"), "_2", "_sd")
  r1<-select(summ, contains("_mean"))
  colnames(r1)<-str_remove_all(colnames(r1), "_mean")
  r2<- dplyr::select(summ, contains("_sd"))
  colnames(r2)<-str_remove_all(colnames(r2), "_sd")
  summ2<-rbind.fill(r1, r2)
  row.names(summ2)<-c("mean", "sd (between iterations)")
  summary_df<-dplyr::select(summ2, contains("looIC"), contains("phi"), contains("Estimate"))

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

  patchPlot<-phiOnePlot+phiTwoPlot+phiThreePlot
  
  modelSimsOutput<-list(metrics_df, summary_df, patchPlot)
  return(modelSimsOutput)
}



################################## run different sigmas to compare ################################## 

# modelSimsOutput_Exp<-modelSims(iterations = 20, sigma="exp", xTime=25, nSamples = 20, phi1_a=200, phi2_a=13, phi3_a=3, phi1_b=160, phi2_b=13, phi3_b=3.5)
# 
# save(modelSimsOutput_Exp,file ="modelSimsOutputs_exp_20.rda")

# start<-Sys.time()
# modelSimsOutput_None<-modelSims(iterations = 2, sigma="none", xTime=25, nSamples = 20, phi1_a=200, phi2_a=13, phi3_a=3, phi1_b=160, phi2_b=13, phi3_b=3.5)
# Sys.time()-start
# 
# start<-Sys.time()
# modelSimsOutput_linear<-modelSims(iterations = 20, sigma="linear", xTime=25, nSamples = 20, phi1_a=200, phi2_a=13, phi3_a=3, phi1_b=160, phi2_b=13, phi3_b=3.5)
# Sys.time()-start
# 
# start<-Sys.time()
# modelSimsOutput_Spline<-modelSims(iterations = 25, sigma="spline", xTime=25, nSamples = 20, phi1_a=200, phi2_a=13, phi3_a=3, phi1_b=160, phi2_b=13, phi3_b=3.5)
# Sys.time()-start
#save(modelSimsOutput_Spline,file ="modelSimsOutputs_spline_20_v2.rda")

# start<-Sys.time()
# modelSimsOutput_quad<-modelSims(iterations = 25, sigma="quad", xTime=25, nSamples = 20, phi1_a=200, phi2_a=13, phi3_a=3, phi1_b=160, phi2_b=13, phi3_b=3.5)
# Sys.time()-start
#save(modelSimsOutput_quad,file ="modelSimsOutputs_quad_20.rda")

################################## Comparing outputs function ################################## 

compareSimModels<-function(...){
argnames <- sys.call()
arguments<-list(...)
i =1
for (dat in arguments){
  name = as.character(rlang::enquo(dat))
  df<-as.data.frame(dat[1])

  if(i==1){ 
    loo_IC_Mean<-mean(df$loo_IC)
    loo_IC_SE<-mean(df$loo_IC_se)
    phi1_estimate_a<-mean(df$Estimate_phi1_treatmenta)
    phi1_estimate_b<-mean(df$Estimate_phi1_treatmentb)
    phi2_estimate_a<-mean(df$Estimate_phi2_treatmenta)
    phi2_estimate_b<-mean(df$Estimate_phi2_treatmentb)
    phi3_estimate_a<-mean(df$Estimate_phi3_treatmenta)
    phi3_estimate_b<-mean(df$Estimate_phi3_treatmentb)
    name<-unlist(lapply(argnames[-1], as.character))[i]
    output_df<-data.frame(modelName=name, loo_IC_Mean= loo_IC_Mean, loo_IC_SE=loo_IC_SE, 
                          phi1_estimate_a=phi1_estimate_a, phi1_estimate_b=phi1_estimate_b,
                          phi2_estimate_a=phi2_estimate_a, phi2_estimate_b=phi2_estimate_b,
                          phi3_estimate_a=phi3_estimate_a,phi3_estimate_b=phi3_estimate_b)
  } else {  
    loo_IC_Mean<-mean(df$loo_IC)
    loo_IC_SE<-mean(df$loo_IC_se)
    phi1_estimate_a<-mean(df$Estimate_phi1_treatmenta)
    phi1_estimate_b<-mean(df$Estimate_phi1_treatmentb)
    phi2_estimate_a<-mean(df$Estimate_phi2_treatmenta)
    phi2_estimate_b<-mean(df$Estimate_phi2_treatmentb)
    phi3_estimate_a<-mean(df$Estimate_phi3_treatmenta)
    phi3_estimate_b<-mean(df$Estimate_phi3_treatmentb)
    name<-unlist(lapply(argnames[-1], as.character))[i]
    output_df_new_row<-data.frame(modelName=name, loo_IC_Mean= loo_IC_Mean, loo_IC_SE=loo_IC_SE, 
                                  phi1_estimate_a=phi1_estimate_a, phi1_estimate_b=phi1_estimate_b,
                                  phi2_estimate_a=phi2_estimate_a, phi2_estimate_b=phi2_estimate_b,
                                  phi3_estimate_a=phi3_estimate_a,phi3_estimate_b=phi3_estimate_b)
    output_df<-rbind(output_df, output_df_new_row)
  }
  i=i+1
}
output_df$modelName <- with(output_df, reorder(modelName, -loo_IC_Mean))

looPlot<-ggplot(output_df)+
  geom_col(aes(x=modelName, y=loo_IC_Mean, fill=modelName))+
  geom_segment(aes(x=modelName, xend=modelName, y=loo_IC_Mean, yend = loo_IC_Mean+loo_IC_SE), size = 1.5)+
  geom_segment(aes(x=modelName, xend=modelName, y=loo_IC_Mean, yend = loo_IC_Mean-loo_IC_SE), size = 1.5)+
  geom_text(aes(x=modelName, y=loo_IC_Mean*1.05, label = paste0(round(loo_IC_Mean, digits=0))))+
  labs(title = "LOO IC", x="", y="", fill="Model")+
  scale_fill_viridis(discrete = T)+
  theme_minimal() +
  theme(axis.line.y.left = element_line(),
        axis.line.x.bottom = element_line(),
        axis.text.x.bottom = element_blank())

output_df_lean<-dplyr::select(output_df, modelName, loo_IC_Mean, loo_IC_SE=loo_IC_SE)%>%arrange(desc(loo_IC_Mean))
outputList<-list(output_df_lean, output_df, looPlot)
return(outputList)
}

`Homoskedastic Model`<-modelSimsOutput_None
`Linear Heteroskedasticity`<-modelSimsOutput_linear
`Quadratic Heteroskedasticity`<-modelSimsOutput_quad
`Spline Heteroskedasticity`<-modelSimsOutput_Spline

test<-compareSimModels(`Homoskedastic Model`, `Linear Heteroskedasticity`, `Quadratic Heteroskedasticity`, `Spline Heteroskedasticity`)

looICPlot<-test[3][[1]]
looICPlot
ggsave("Fig2/LOO_IC_sim_data.png",looICPlot, width = 7.04, height=4.04, dpi=300, bg = "#ffffff")


################################## Function to Recapitulate growthSims Params ################################## 

recapNoise<-function(x, rescale=4){
  
  if (!is.data.frame(x)){
    df = x[[1]]} else{df = x}
  sd_phi1_treatmenta<-sd(df$Estimate_phi1_treatmenta)*rescale
  sd_phi1_treatmentb<-sd(df$Estimate_phi1_treatmentb)*rescale
  sd_phi2_treatmenta<-sd(df$Estimate_phi2_treatmenta)*rescale
  sd_phi2_treatmentb<-sd(df$Estimate_phi2_treatmentb)*rescale
  sd_phi3_treatmenta<-sd(df$Estimate_phi3_treatmenta)*rescale
  sd_phi3_treatmentb<-sd(df$Estimate_phi3_treatmentb)*rescale
  
  phi1_noise_sd<-mean(sd_phi1_treatmenta, sd_phi1_treatmentb)
  phi2_noise_sd<-mean(sd_phi2_treatmenta, sd_phi2_treatmentb)
  phi3_noise_sd<-mean(sd_phi3_treatmenta, sd_phi3_treatmentb)
  
  outputData = data.frame(phi1_noise_sd=phi1_noise_sd, phi2_noise_sd=phi2_noise_sd, phi3_noise_sd=phi3_noise_sd)
  return(outputData)
}

metdf<-read_csv("metrics_df_100.csv")
recapNoise(x=metdf)

################################## Make Dataframe for plots ##################################

xTime=25
nSamples = 20
phi1_1=200
phi2_1=13
phi3_1=3
phi1_2=160
phi2_2=13
phi3_2=3.5
x<-1:xTime
df<-rbind(
  do.call(rbind,
          lapply(1:nSamples,
                 function(i) data.frame("sample"=paste0("sample_",i),
                                        "treatment"="a",
                                        "time"=x,
                                        "y"=growthSim(x, phi1_1, phi2_1, phi3_1),
                                        stringsAsFactors = F))),
  do.call(rbind,
          lapply(1:nSamples,
                 function(i) data.frame("sample"=paste0("sample_",i),
                                        "treatment"="b",
                                        "time"=x,
                                        "y"=growthSim(x,phi1_2, phi2_2, phi3_2),
                                        stringsAsFactors = F))))
  

############################### Quadratic Via TimeSquared Column with NLF ##################################

df$timeSQ<- df$time^2

prior_quad <- prior(lognormal(log(130), .25),nlpar = "phi1",coef="treatmenta") +
  prior(lognormal(log(130), .25),nlpar = "phi1",coef="treatmentb") +
  prior(lognormal(log(3), .25), nlpar = "phi2",coef="treatmenta") +
  prior(lognormal(log(3), .25), nlpar = "phi2",coef="treatmentb") +
  prior(lognormal(log(3), .25), nlpar = "phi3",coef="treatmenta") +
  prior(lognormal(log(3), .25), nlpar = "phi3",coef="treatmentb") +
  prior(student_t(3,0,5), dpar="sigma", coef="time") + #works with default flat sigma prior (t for intercept)
  prior(student_t(3,0,5), dpar="sigma", coef="time:treatmentb") + #also works with specifying sigma dpar and coefs
  prior(student_t(3,0,5), dpar="sigma", coef="timeSQ:treatmentb") + #no divergence, priors are there,
  prior(student_t(3,0,5), dpar="sigma", coef="timeSQ") + #estimates are reasonable,
  prior(gamma(2,0.1), class="nu") #but it is slow, I think it's good to throw at the cluster though.
start<-Sys.time()
fit_quad_withSQ_col_init0 <- brm(bf(y ~ phi1/(1+exp((phi2-time)/phi3)),
                              lf(sigma~ time + timeSQ + time:treatment + timeSQ:treatment),
                              phi1 + phi2 + phi3 ~ 0+treatment,
                              autocor = ~arma(~time|sample:treatment,1,1),nl = TRUE),
                           family = student, prior = prior_quad, data = df, iter = 1000,
                           cores = 2, chains = 2, backend = "cmdstanr", 
                           control = list(adapt_delta = 0.999,max_treedepth = 20),
                           inits = function(){list(b_phi1=rgamma(2,1),
                                                   b_phi2=rgamma(2,1),
                                                   b_phi3=rgamma(2,1), 
                                                   sigma = 0#rgamma(3,1)
                                                   )})#,sigma=rgamma(2,0.1)
Sys.time()-start


################################## Exponential via a*exp(bx) ##################################

prior_exp <- prior(lognormal(log(130), .25),nlpar = "phi1",coef="treatmenta") +
  prior(lognormal(log(130), .25),nlpar = "phi1",coef="treatmentb") +
  prior(lognormal(log(12), .25), nlpar = "phi2",coef="treatmenta") +
  prior(lognormal(log(12), .25), nlpar = "phi2",coef="treatmentb") +
  prior(lognormal(log(3), .25), nlpar = "phi3",coef="treatmenta") +
  prior(lognormal(log(3), .25), nlpar = "phi3",coef="treatmentb") +
  prior(student_t(1, 0, 5), nlpar="a") + # tends to run with student_t(3, 0, 5) prior, but a and r should be positive here 
  prior(student_t(1, 0, 5), nlpar="r") + # failing with strictly positive priors or with inits specified
  prior(gamma(2,0.1), class="nu")

form = bf(y ~ phi1/(1+exp((phi2-time)/(phi3))),
          phi1 + phi2 + phi3 ~ 0+treatment,
          autocor = ~arma(~time|sample:treatment,1,1), nl = TRUE)+
  nlf(sigma ~ a * exp(r * time))+
  lf(a + r ~ 0+treatment)

fit_exp_2 <- brm(form,
                 family = student, prior = prior_exp, 
                 data = df, iter =10000,
                 cores = 2, chains = 4, backend = "cmdstanr",
                 control = list(adapt_delta = 0.999, max_treedepth = 20),
                 inits = function(){list(b_phi1=rgamma(2,1),
                                         b_phi2=rgamma(2,1),
                                         b_phi3=rgamma(2,1))})


################################## Graphing Ribbons Spline ################################## 

prior_spline <- prior(lognormal(log(130), .25),nlpar = "phi1",coef="treatmenta") +
  prior(lognormal(log(130), .25),nlpar = "phi1",coef="treatmentb") +
  prior(lognormal(log(12), .25), nlpar = "phi2",coef="treatmenta") +
  prior(lognormal(log(12), .25), nlpar = "phi2",coef="treatmentb") +
  prior(lognormal(log(3), .25), nlpar = "phi3",coef="treatmenta") +
  prior(lognormal(log(3), .25), nlpar = "phi3",coef="treatmentb") +
  prior(student_t(3,0,5), dpar="sigma") +
  prior(gamma(2,0.1), class="nu")

fit_spline <- brm(bf(y ~ phi1/(1+exp((phi2-time)/phi3)),
                     sigma~s(time,by=treatment), 
                     phi1 + phi2 + phi3 ~ 0+treatment,
                     autocor = ~arma(~time|sample:treatment,1,1),
                     nl = TRUE),
                  family = student, prior = prior_spline, data = df, iter = 4000, 
                  cores = 2, 
                  chains = 4, 
                  backend = "cmdstanr", 
                  control = list(adapt_delta = 0.999,
                                 max_treedepth = 20),
                  inits = function(){list(b_phi1=rgamma(2,1),b_phi2=rgamma(2,1),b_phi3=rgamma(2,1))}) 

save(fit_spline, prior_spline, file="spline_model_and_priors.rdata")

load("spline_model_and_priors.rdata")

probs <- seq(from=99, to=1, by=-2)/100
avg_pal <- turbo(n=length(probs))
df_test_spline <- rbind(data.frame("treatment"="a",time=1:25, sample="new1"),data.frame("treatment"="b",time=1:25, sample="new2"))
df_pred_spline <- predict(fit_spline, df_test_spline, probs=probs)
test_spline<- cbind(df_test_spline, df_pred_spline)

pSpline<-ggplot(test_spline,aes(time,Estimate))+
  facet_wrap(~treatment)+
  lapply(seq(1,49,2),function(i) geom_ribbon(aes_string(ymin=paste("Q",i,sep = ""),ymax=paste("Q",100-i,sep = "")),fill=avg_pal[i],alpha=0.5))+
  geom_line(data=df,aes(time, y, group=interaction(treatment,sample)),color="gray20", size=0.3)+
  ylab("Area (cm, simulated)")+
  xlab("Time")+
  theme_light()+
  theme(axis.ticks.length=unit(0.2,"cm"))+
  theme(strip.background=element_rect(fill="gray50",color="gray20"),
        strip.text.x=element_text(size=14,color="white"),
        strip.text.y=element_text(size=14,color="white"))+
  theme(axis.title= element_text(size = 18))+
  theme(title= element_text(size = 20))+
  theme(axis.text = element_text(size = 14))+
  theme(legend.position='top')
pSpline

ggsave("Fig2/spline_ribbons_with_sim_data.png",pSpline, width = 7.04, height=4.04, dpi=300)

################################## Graphing Ribbons None ################################## 


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
                family = student, prior = prior_none, data = df, iter = 4000,
                cores = 2, chains = 4, backend = "cmdstanr", #threads = threading(4),
                control = list(adapt_delta = 0.999,max_treedepth = 20),
                inits = function(){list(b_phi1=rgamma(2,1),b_phi2=rgamma(2,1),b_phi3=rgamma(2,1))})

save(fit_none, prior_none, file="none_model_and_priors.rdata")

load("none_model_and_priors.rdata")

library(viridisLite)
probs <- seq(from=99, to=1, by=-2)/100
avg_pal <- turbo(n=length(probs))
df_test_none <- rbind(data.frame("treatment"="a",time=1:25,sample="new1"),data.frame("treatment"="b",time=1:25,sample="new2"))
df_pred_none <- predict(fit_none, df_test_none, probs=probs) 
test_none <- cbind(df_test_none,df_pred_none)

pNone <- ggplot(test_none, aes(time,Estimate))+
  facet_wrap(~treatment)+
  lapply(seq(1,49,2),function(i) geom_ribbon(aes_string(ymin=paste("Q",i,sep = ""),ymax=paste("Q",100-i,sep = "")),fill=avg_pal[i],alpha=0.5))+
  geom_line(data=df,aes(time,y,group=interaction(treatment,sample)),color="gray20", size=0.3)+
  ylab("Area (cm, simulated)")+
  xlab("Time")+
  coord_cartesian(ylim=c(0,250))+
  theme_light()+
  theme(axis.ticks.length=unit(0.2,"cm"))+
  theme(strip.background=element_rect(fill="gray50",color="gray20"),
        strip.text.x=element_text(size=14,color="white"),
        strip.text.y=element_text(size=14,color="white"))+
  theme(axis.title= element_text(size = 18))+
  theme(title= element_text(size = 20))+
  theme(axis.text = element_text(size = 14))+
  theme(legend.position='top')
pNone

ggsave("Fig2/homoskedastic_ribbons_with_sim_data.png",pNone, width = 7.04, height=4.04, dpi=300)

################################## Graphing Ribbons Linear ################################## 

prior_linear <- prior(lognormal(log(130), .25),nlpar = "phi1",coef="treatmenta") +
  prior(lognormal(log(130), .25),nlpar = "phi1",coef="treatmentb") +
  prior(lognormal(log(12), .25), nlpar = "phi2",coef="treatmenta") +
  prior(lognormal(log(12), .25), nlpar = "phi2",coef="treatmentb") +
  prior(lognormal(log(3), .25), nlpar = "phi3",coef="treatmenta") +
  prior(lognormal(log(3), .25), nlpar = "phi3",coef="treatmentb") +
  prior(student_t(3,0,5), dpar="sigma") +
  prior(gamma(2,0.1), class="nu")

fit_linear <- brm(bf(y ~ phi1/(1+exp((phi2-time)/phi3)),
                     sigma~time+time:treatment,
                     phi1 + phi2 + phi3 ~ 0+treatment,
                     autocor = ~arma(~time|sample:treatment,1,1),nl = TRUE),
                  family = student, prior = prior_linear, data = df, iter = 4000,
                  cores = 2, chains = 4, backend = "cmdstanr", #threads = threading(4),
                  control = list(adapt_delta = 0.999,max_treedepth = 20),
                  inits = function(){list(b_phi1=rgamma(2,1),b_phi2=rgamma(2,1),b_phi3=rgamma(2,1))})

save(fit_linear, prior_linear, file="linear_model_and_priors.rdata")

load("linear_model_and_priors.rdata")


probs <- seq(from=99, to=1, by=-2)/100
avg_pal <- turbo(n=length(probs))
df_test_linear <- rbind(data.frame("treatment"="a",time=1:25,sample="new1"),
                 data.frame("treatment"="b",time=1:25,sample="new2"))
df_pred_linear <- predict(fit_linear, df_test_linear, probs=probs) 
test_linear <- cbind(df_test_linear, df_pred_linear)

pLinear <- ggplot(test_linear,aes(time,Estimate))+
  facet_wrap(~treatment)+
  lapply(seq(1,49,2),function(i) geom_ribbon(aes_string(ymin=paste("Q",i,sep = ""),ymax=paste("Q",100-i,sep = "")),fill=avg_pal[i],alpha=0.5))+
  geom_line(data=df,aes(time,y,group=interaction(treatment,sample)),color="gray20", size=0.3)+
  ylab("Area (cm, simulated)")+
  xlab("Time")+
  theme_light()+
  theme(axis.ticks.length=unit(0.2,"cm"))+
  theme(strip.background=element_rect(fill="gray50",color="gray20"),
        strip.text.x=element_text(size=14,color="white"),
        strip.text.y=element_text(size=14,color="white"))+
  theme(axis.title= element_text(size = 18))+
  theme(title= element_text(size = 20))+
  theme(axis.text = element_text(size = 14))+
  theme(legend.position='top')
pLinear
ggsave("Fig2/linear_ribbons_with_sim_data.png",pLinear, width = 7.04, height=4.04, dpi=300)

################################## Graphing Ribbons Quad ################################## 

df$timeSQ<- df$time^2

prior_quad <- prior(lognormal(log(130), .25),nlpar = "phi1",coef="treatmenta") +
  prior(lognormal(log(130), .25),nlpar = "phi1",coef="treatmentb") +
  prior(lognormal(log(3), .25), nlpar = "phi2",coef="treatmenta") +
  prior(lognormal(log(3), .25), nlpar = "phi2",coef="treatmentb") +
  prior(lognormal(log(3), .25), nlpar = "phi3",coef="treatmenta") +
  prior(lognormal(log(3), .25), nlpar = "phi3",coef="treatmentb") +
  prior(student_t(3,0,5), dpar="sigma", coef="time") + 
  prior(student_t(3,0,5), dpar="sigma", coef="time:treatmentb") + 
  prior(student_t(3,0,5), dpar="sigma", coef="timeSQ:treatmentb") + 
  prior(student_t(3,0,5), dpar="sigma", coef="timeSQ") + 
  prior(gamma(2,0.1), class="nu") 
start<-Sys.time()
# fit_quad_cluster <- brm(bf(y ~ phi1/(1+exp((phi2-time)/phi3)),
#                                     lf(sigma~ time + timeSQ + time:treatment + timeSQ:treatment),
#                                     phi1 + phi2 + phi3 ~ 0+treatment,
#                                     autocor = ~arma(~time|sample:treatment,1,1),nl = TRUE),
#                                  family = student, prior = prior_quad, data = df, iter = 10000,
#                                  cores = 4, chains = 4, backend = "cmdstanr", 
#                                  control = list(adapt_delta = 0.999,max_treedepth = 20),
#                                  inits = function(){list(b_phi1=rgamma(2,1),
#                                                          b_phi2=rgamma(2,1),
#                                                          b_phi3=rgamma(2,1), 
#                                                          sigma = 0
#                                  )})
# Sys.time()-start
# save(fit_quad_cluster, "fit_quad_cluster.rdata")



load("fit_quad_cluster.rdata") # loads fit_quad_cluster model

probs <- seq(from=99, to=1, by=-2)/100
avg_pal <- turbo(n=length(probs))
df_test_quad <- rbind(data.frame("treatment"="a",time=1:25, sample="new1"),
                 data.frame("treatment"="b",time=1:25, sample="new2"))
df_test_quad$timeSQ<-df_test_quad$time^2
df_pred_quad <- predict(fit_quad_cluster, df_test_quad, probs=probs) 
test_quad <- cbind(df_test_quad, df_pred_quad)

pQuad <- ggplot(test_quad,aes(time,Estimate))+
  facet_wrap(~treatment)+
  lapply(seq(1,49,2),function(i) geom_ribbon(aes_string(ymin=paste("Q",i,sep = ""),ymax=paste("Q",100-i,sep = "")),fill=avg_pal[i],alpha=0.5))+
  geom_line(data=df,aes(time,y,group=interaction(treatment,sample)),color="gray20")+
  ylab("Area (cm, simulated)")+
  xlab("Time")+
  theme_light()+
  theme(axis.ticks.length=unit(0.2,"cm"))+
  theme(strip.background=element_rect(fill="gray50",color="gray20"),
        strip.text.x=element_text(size=14,color="white"),
        strip.text.y=element_text(size=14,color="white"))+
  theme(axis.title= element_text(size = 18))+
  theme(title= element_text(size = 20))+
  theme(axis.text = element_text(size = 14))+
  theme(legend.position='top')
pQuad

ggsave("Fig2/quadratic_ribbons_with_sim_data.png",pQuad, width = 7.04, height=4.04, dpi=300)

################################### Ribbons for EXP ###############################

load("exponential_fit_2.rdata") #loads fit_exp_2

probs <- seq(from=99, to=1, by=-2)/100
avg_pal <- turbo(n=length(probs))
df_test_exp_a <- rbind(data.frame("treatment"="a",time=1:25,sample="new1"))
df_test_exp_b <- rbind(data.frame("treatment"="b",time=1:25,sample="new2"))
df_test_exp<-rbind(df_test_exp_a, df_test_exp_b)
df_pred_exp <- predict(fit_exp_2, df_test_exp, probs=probs) 
test_exp <- cbind(df_test_exp, df_pred_exp)

#summary(fit_exp_2)

pExp <- ggplot(test_exp,aes(time,Estimate))+
  facet_wrap(~treatment)+
  lapply(seq(1,49,2),function(i) geom_ribbon(aes_string(ymin=paste("Q",i,sep = ""),ymax=paste("Q",100-i,sep = "")),fill=avg_pal[i],alpha=0.5))+
  geom_line(data=df,aes(time,y,group=interaction(treatment,sample)),color="gray20")+
  ylab("Area (cm, simulated)")+
  xlab("Time")+
  theme_light()+
  theme(axis.ticks.length=unit(0.2,"cm"))+
  theme(strip.background=element_rect(fill="gray50",color="gray20"),
        strip.text.x=element_text(size=14,color="white"),
        strip.text.y=element_text(size=14,color="white"))+
  theme(axis.title= element_text(size = 18))+
  theme(title= element_text(size = 20))+
  theme(axis.text = element_text(size = 14))+
  theme(legend.position='top')
pExp

ggsave("Fig2/exponential_ribbons_with_sim_data.png",pExp, width = 7.04, height=4.04, dpi=300)

################################## Patchwork ################################## 

pSpline_titled<-pSpline+labs(title="Splines")
pLinear_titled<-pLinear+labs(title="Linear")
pExp_titled<-pExp+labs(title="Exponential")
pQuad_titled<-pQuad+labs(title="Quadratic")
pNone_titled<-pNone+labs(title="None")

patchesOhoolihan<- (pSpline_titled +pQuad_titled)/ (pNone_titled + pLinear_titled)
patchesOhoolihan

ggsave("Fig2/CombinedRibbons.png", patchesOhoolihan,width = 9.54, height=9.54, dpi=300)













