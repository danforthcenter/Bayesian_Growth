library(brms)
library(ggplot2)
library(stringr)
library(patchwork)
library(cowplot)
library(viridis)
library(viridisLite)

growthSimMM <- function(x,a,b){
  a_r <- a+rnorm(1,mean = 0,sd=10)
  b_r <- b+rnorm(1,mean=0,sd=0.02)
  return(a_r-a_r*exp(-b_r*x))
}
set.seed(321)
x <- 1:25
df <- rbind(
  do.call(rbind,lapply(1:20,function(i) data.frame("sample"=paste0("sample_",i),"treatment"="a","time"=x,"y"=growthSimMM(x,200,.08),stringsAsFactors = F))),
  do.call(rbind,lapply(1:20,function(i) data.frame("sample"=paste0("sample_",i),"treatment"="b","time"=x,"y"=growthSimMM(x,160,.1),stringsAsFactors = F)))
)

p <- ggplot(df,aes(time,y,group=interaction(treatment,sample)))+
  geom_line(aes(color=treatment))+
  ylab("Area (cm, simulated)")+
  xlab("Time")+
  labs(title="Monomolecular")+
  theme_light()+
  theme(axis.ticks.length=unit(0.2,"cm"))+
  theme(strip.background=element_rect(fill="gray50",color="gray20"),
        strip.text.x=element_text(size=14,color="white"),
        strip.text.y=element_text(size=14,color="white"))+
  theme(axis.title= element_text(size = 18))+
  theme(title= element_text(size = 20))+
  theme(axis.text = element_text(size = 14))+
  theme(legend.position='top')
p




priorMono <- prior(lognormal(log(130), .25),nlpar = "a") +
  prior(lognormal(log(2), .25), nlpar = "b") + 
  prior(student_t(3,0,5), dpar="sigma") +
  prior(gamma(2,0.1), class="nu")

fitMono1 <- brm(bf(y ~ a-a*exp(-b*time), 
               sigma~time:treatment, 
               a + b ~ 0+treatment, 
               autocor = ~arma(~time|sample:treatment,1,1),nl = TRUE),
            family = student, prior = priorMono, data = df, iter = 4000, 
            cores = 2, chains = 4, backend = "cmdstanr", #threads = threading(4),
            control = list(adapt_delta = 0.999,max_treedepth = 20),
            inits = function(){list(b_a=rgamma(2,1),b_b=rgamma(2,1),b_c=rgamma(2,1))})



h <- hypothesis(fitMono1, "a_treatmenta/a_treatmentb > 1")

plot(h)

probs <- seq(from=99, to=1, by=-2)/100
avg_pal <- plasma(n=length(probs))
df_test <- rbind(data.frame("treatment"="a",time=1:25,sample="new1"),data.frame("treatment"="b",time=1:25,sample="new2"))
df_pred <- predict(fitMono1,df_test,probs=probs) # add fit_spline in place of fit1 for first checks.
test <- cbind(df_test,df_pred)

p <- ggplot(test,aes(time,Estimate))+
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
  theme(axis.text = element_text(size = 14))+
  theme(legend.position='top')
p
ggsave("Fig2/monomolecular_linear_ribbon.png",p, width = 7.04, height=4.04, dpi=300)




#***************************************************************************************
################################## Make Data Function ################################## 
#***************************************************************************************

library(brms)
library(ggplot2)
library(stringr)
library(patchwork)
library(plyr)
library(tidyverse)
library(cowplot)
library(viridis)
library(viridisLite)

growthSimMM <- function(x,a,b){
  a_r <- a+rnorm(1,mean = 0,sd=10)
  b_r <- b+rnorm(1,mean=0,sd=0.02)
  return(a_r-a_r*exp(-b_r*x))
}


makeData<-function(x.=x,nSamples.=nSamples, a_1.=a_1, a_2.=a_2, b_1.= b_1, b_2.=b_2){
  df<-rbind(
    do.call(rbind,
            lapply(1:nSamples.,
                   function(i) data.frame("sample"=paste0("sample_",i),
                                          "treatment"="a",
                                          "time"=x.,
                                          "y"=growthSimMM(x., a_1., b_1.),
                                          stringsAsFactors = F))),
    do.call(rbind,
            lapply(1:nSamples.,
                   function(i) data.frame("sample"=paste0("sample_",i),
                                          "treatment"="b",
                                          "time"=x.,
                                          "y"=growthSimMM(x., a_2., b_2.),
                                          stringsAsFactors = F)))
  )
  return(df)
}

#*******************************************************************************************************************
################################## define modelSims() for Power Law Growth Models ################################## 200,13,.2
#*******************************************************************************************************************

modelSimsMM<-function(iterations = 5, sigma = "none", xTime=25, nSamples = 20, a_1 = 200, a_2=160, b_1 = .08, b_2=.1){
  sigma_<-ifelse(sigma=="linear", "sigma~time+time:treatment,", 
                 ifelse(sigma=="spline", "sigma~s(time,by=treatment),",
                        ifelse(sigma=="quad", "lf(sigma~ time + timeSQ + time:treatment + timeSQ:treatment),",
                               ifelse(sigma=="exp", "sigma~time+time:treatment,", #Not evaluated, doesn't matter
                                      ifelse(sigma=="none", "", 
                                             paste0(sigma, ","))))))
  bfText<-paste0("bf(y ~ a-a*exp(-b*time),", 
                 sigma_, 
                 "a + b ~ 0+treatment,autocor = ~arma(~time|sample:treatment,1,1),nl = TRUE)")
  bayesFormula<-eval(parse(text = bfText))
  
  prior_list<-prior(lognormal(log(130), .25),nlpar = "a") +
    prior(lognormal(log(2), .25), nlpar = "b") + 
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
      prior(student_t(1, 0, 5), nlpar="d") + 
      prior(student_t(1, 0, 5), nlpar="r")
    
    bayesFormula<-bf(y ~ a-a*exp(-b*time),
                     a + b ~ 0+treatment,
                     autocor = ~arma(~time|sample:treatment,1,1), nl = TRUE)+
      nlf(sigma ~ d * exp(r * time))+
      lf(d + r ~ 0+treatment)
  }
  
  for (i in 1:iterations){
    cat("\nStarting Iteration ", i, "/",iterations,"\n") 
    iteration_row<-data.frame(iteration = i, elpd_loo = NA,elpd_loo_se=NA, p_loo=NA, p_loo_se=NA, loo_IC=NA,loo_IC_se=NA) #store iteration number
    x<-1:xTime
    dat <- makeData(x.=x,nSamples.=nSamples, a_1.=a_1, a_2.=a_2, b_1.= b_1, b_2.=b_2) #make the data
    if(sigma=="quad"){dat$timeSQ<-dat$time^2}
    fit_none <- brm(bayesFormula,
                    family = student, prior = prior_list, data = dat, iter = 2000,
                    cores = 2, chains = 2, backend = "cmdstanr",
                    control = list(adapt_delta = 0.999,max_treedepth = 20),
                    inits = function(){list(b_a=rgamma(2,1),b_b=rgamma(2,1),b_c=rgamma(2,1))})
    
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
  summary_df<-dplyr::select(summ2, contains("looIC"), contains("Estimate"))
  
  aPlot<-metrics_df%>%
    dplyr::select(contains("Estimate_a"))%>%
    pivot_longer(cols = starts_with("Estimate"))%>%
    dplyr::mutate(name=factor(str_remove_all(name, "Estimate_")))%>%
    ggplot()+
    geom_boxplot(aes(x=value, y=name))+
    geom_point(aes(x=a_1, y=1), color="red", size=0.5)+
    geom_point(aes(x=a_2, y=2), color="red", size=0.5)+
    annotate("text", x=a_1, y=0.55, label=paste0(a_1), color="red", size=3)+
    annotate("text", x=a_2, y=1.55, label=paste0(a_2), color="red", size=3)+
    annotate("text", x = median(metrics_df$Estimate_a_treatmenta),
             y = 1.45,
             label = paste("Est:", round(median(metrics_df$Estimate_a_treatmenta), 2)),
             size = 3)+
    annotate("text", x = median(metrics_df$Estimate_a_treatmentb),
             y = 2.45,
             label = paste("Est:", round(median(metrics_df$Estimate_a_treatmentb), 2)),
             size = 3)+
    labs(title = paste0("'a' Estimation on ", iterations, " iterations"), x="", y="")+
    theme_minimal() +
    theme(axis.line.y.left = element_line(),
          axis.line.x.bottom = element_line())
  
  bPlot<-metrics_df%>%
    dplyr::select(contains("Estimate_b"))%>%
    pivot_longer(cols = starts_with("Estimate"))%>%
    dplyr::mutate(name=factor(str_remove_all(name, "Estimate_")))%>%
    ggplot()+
    geom_boxplot(aes(x=value, y=name))+
    geom_point(aes(x=b_1, y=1), color="red", size=0.5)+
    geom_point(aes(x=b_2, y=2), color="red", size=0.5)+
    annotate("text", x=b_1, y=0.55, label=paste0(b_1), color="red", size=3)+
    annotate("text", x=b_2, y=1.55, label=paste0(b_2), color="red", size=3)+
    annotate("text", x = median(metrics_df$Estimate_b_treatmenta),
             y = 1.45,
             label = paste("Est:", round(median(metrics_df$Estimate_b_treatmenta), 2)),
             size = 3)+
    annotate("text", x = median(metrics_df$Estimate_b_treatmentb),
             y = 2.45,
             label = paste("Est:", round(median(metrics_df$Estimate_b_treatmentb), 2)),
             size = 3)+
    labs(title = paste0("'b' Estimation on ", iterations, " iterations"), x="", y="")+
    theme_minimal() +
    theme(axis.line.y.left = element_line(),
          axis.line.x.bottom = element_line())
  
  patchPlot<-aPlot+bPlot
  
  renameModText<-paste0("fit_MM_",sigma,"<-fit_none")
  eval(parse(text=renameModText))
  outputText<-paste0("modelSimsOutput<-list(metrics_df, summary_df, patchPlot, fit_MM_",sigma,")")
  eval(parse(text=outputText))
  
  return(modelSimsOutput)
}

# modelSims_MM_homo<-modelSimsMM(sigma = "none")
# save(modelSims_MM_homo, file="modelSimsOutputFiles/monomolecularOutputs/modelSims_MM_homo_2.rdata")
# 2+2
# 
# modelSims_MM_linear<-modelSimsMM(sigma = "linear")
# save(modelSims_MM_linear, file="modelSims_MM_linear_2.rdata")
# 2+2

modelSims_MM_spline<-modelSimsMM(sigma = "spline")
save(modelSims_MM_spline, file="modelSims_MM_spline_2.rdata")
2+2

# modelSims_MM_exp<-modelSimsMM(sigma = "exp")
# save(modelSims_MM_exp, file="modelSims_MM_exp_2.rdata")
# 2+2

# modelSims_MM_quad<-modelSimsMM(sigma = "quad")
# save(modelSims_MM_quad, file="modelSims_MM_quad_2.rdata")
# 2+2


#*******************************************************************************************************************
################################## mono modelSims() Outputs ##################################
#*******************************************************************************************************************

# for gradual
growthSimMM <- function(x,a,b){
  a_r <- a+rnorm(1,mean = 0,sd=10)
  b_r <- b+rnorm(1,mean=0,sd=0.02)
  return(a_r-a_r*exp(-b_r*x))
}

set.seed(321)
x <- 1:25
df <- rbind(
  do.call(rbind,lapply(1:20,function(i) data.frame("sample"=paste0("sample_",i),"treatment"="a","time"=x,"y"=growthSimMM(x,200,.08),stringsAsFactors = F))),
  do.call(rbind,lapply(1:20,function(i) data.frame("sample"=paste0("sample_",i),"treatment"="b","time"=x,"y"=growthSimMM(x,160,.1),stringsAsFactors = F)))
)

# for sharp

growthSimMM <- function(x,a,b){
  a_r <- a+rnorm(1,mean = 0,sd=10)
  b_r <- b+rnorm(1,mean=0,sd=0.025)
  return(a_r-a_r*exp(-b_r*x))
}
set.seed(123)
x <- 1:25
df <- rbind(
  do.call(rbind,lapply(1:20,function(i) data.frame("sample"=paste0("sample_",i),"treatment"="a","time"=x,"y"=growthSimMM(x,200,.33),stringsAsFactors = F))),
  do.call(rbind,lapply(1:20,function(i) data.frame("sample"=paste0("sample_",i),"treatment"="b","time"=x,"y"=growthSimMM(x,160,.42),stringsAsFactors = F)))
)




print(load("modelSimsOutputs/monomolecularSims/modelSims_MM_exp.rdata"))

probs <- seq(from=99, to=1, by=-2)/100
avg_pal <- turbo(n=length(probs))
df_test <- rbind(data.frame("treatment"="a",time=1:25,sample="new1"),data.frame("treatment"="b",time=1:25,sample="new2"))
df_pred <- predict(modelSims_MM_exp[[4]],df_test,probs=probs) # add fit_spline in place of fit1 for first checks.
test <- cbind(df_test,df_pred)

mExpPlot <- ggplot(test,aes(time,Estimate))+
  facet_wrap(~treatment)+
  lapply(seq(1,49,2),function(i) geom_ribbon(aes_string(ymin=paste("Q",i,sep = ""),ymax=paste("Q",100-i,sep = "")),fill=avg_pal[i],alpha=0.5))+
  geom_line(data=df,aes(time,y,group=interaction(treatment,sample)),color="gray20", size=0.25)+
  coord_cartesian(ylim=c(0,250))+
  ylab(~~Area~(cm^2))+
  xlab("Time")+
  theme_light()+
  theme(axis.ticks.length=unit(0.2,"cm"))+
  theme(strip.background=element_rect(fill="gray50",color="gray20"),
        strip.text.x=element_text(size=14,color="white"),
        strip.text.y=element_text(size=14,color="white"))+
  theme(axis.title= element_text(size = 18))+
  theme(axis.text = element_text(size = 14))+
  theme(legend.position='top')
mExpPlot




print(load("modelSimsOutputs/monomolecularSims/modelSims_MM_homo.rdata"))

probs <- seq(from=99, to=1, by=-2)/100
avg_pal <- turbo(n=length(probs))
df_test <- rbind(data.frame("treatment"="a",time=1:25,sample="new1"),data.frame("treatment"="b",time=1:25,sample="new2"))
df_pred <- predict(modelSims_MM_homo[[4]],df_test,probs=probs) # add fit_spline in place of fit1 for first checks.
test <- cbind(df_test,df_pred)

mHomoPlot <- ggplot(test,aes(time,Estimate))+
  facet_wrap(~treatment)+
  lapply(seq(1,49,2),function(i) geom_ribbon(aes_string(ymin=paste("Q",i,sep = ""),ymax=paste("Q",100-i,sep = "")),fill=avg_pal[i],alpha=0.5))+
  geom_line(data=df,aes(time,y,group=interaction(treatment,sample)),color="gray20", size=0.25)+
  ylab(~~Area~(cm^2))+
  xlab("Time")+
  coord_cartesian(ylim=c(0,250))+
  theme_light()+
  theme(axis.ticks.length=unit(0.2,"cm"))+
  theme(strip.background=element_rect(fill="gray50",color="gray20"),
        strip.text.x=element_text(size=14,color="white"),
        strip.text.y=element_text(size=14,color="white"))+
  theme(axis.title= element_text(size = 18))+
  theme(axis.text = element_text(size = 14))+
  theme(legend.position='top')
mHomoPlot

print(load("modelSimsOutputs/monomolecularSims/modelSims_MM_linear.rdata"))

probs <- seq(from=99, to=1, by=-2)/100
avg_pal <- turbo(n=length(probs))
df_test <- rbind(data.frame("treatment"="a",time=1:25,sample="new1"),data.frame("treatment"="b",time=1:25,sample="new2"))
df_pred <- predict(modelSims_MM_linear[[4]],df_test,probs=probs) # add fit_spline in place of fit1 for first checks.
test <- cbind(df_test,df_pred)

mLinearPlot <- ggplot(test,aes(time,Estimate))+
  facet_wrap(~treatment)+
  lapply(seq(1,49,2),function(i) geom_ribbon(aes_string(ymin=paste("Q",i,sep = ""),ymax=paste("Q",100-i,sep = "")),fill=avg_pal[i],alpha=0.5))+
  geom_line(data=df,aes(time,y,group=interaction(treatment,sample)),color="gray20", size=0.25)+
  coord_cartesian(ylim=c(0,250))+
  ylab(~~Area~(cm^2))+
  xlab("Time")+
  theme_light()+
  theme(axis.ticks.length=unit(0.2,"cm"))+
  theme(strip.background=element_rect(fill="gray50",color="gray20"),
        strip.text.x=element_text(size=14,color="white"),
        strip.text.y=element_text(size=14,color="white"))+
  theme(axis.title= element_text(size = 18))+
  theme(axis.text = element_text(size = 14))+
  theme(legend.position='top')
mLinearPlot

print(load("modelSimsOutputs/monomolecularSims/modelSims_MM_quad.rdata"))

probs <- seq(from=99, to=1, by=-2)/100
avg_pal <- turbo(n=length(probs))
df_test <- rbind(data.frame("treatment"="a",time=1:25,sample="new1"),data.frame("treatment"="b",time=1:25,sample="new2"))%>%mutate(timeSQ=time^2)
df_pred <- predict(modelSims_MM_quad[[4]],df_test,probs=probs) # add fit_spline in place of fit1 for first checks.
test <- cbind(df_test,df_pred)

mQuadPlot <- ggplot(test,aes(time,Estimate))+
  facet_wrap(~treatment)+
  lapply(seq(1,49,2),function(i) geom_ribbon(aes_string(ymin=paste("Q",i,sep = ""),ymax=paste("Q",100-i,sep = "")),fill=avg_pal[i],alpha=0.5))+
  geom_line(data=df,aes(time,y,group=interaction(treatment,sample)),color="gray20", size=0.25)+
  coord_cartesian(ylim=c(0,250))+
  ylab(~~Area~(cm^2))+
  xlab("Time")+
  theme_light()+
  theme(axis.ticks.length=unit(0.2,"cm"))+
  theme(strip.background=element_rect(fill="gray50",color="gray20"),
        strip.text.x=element_text(size=14,color="white"),
        strip.text.y=element_text(size=14,color="white"))+
  theme(axis.title= element_text(size = 18))+
  theme(axis.text = element_text(size = 14))+
  theme(legend.position='top')
mQuadPlot



print(load("modelSimsOutputs/monomolecularSims/modelSims_MM_spline.rdata"))

probs <- seq(from=99, to=1, by=-2)/100
avg_pal <- turbo(n=length(probs))
df_test <- rbind(data.frame("treatment"="a",time=1:25,sample="new1"),data.frame("treatment"="b",time=1:25,sample="new2"))
df_pred <- predict(modelSims_MM_spline[[4]],df_test,probs=probs) # add fit_spline in place of fit1 for first checks.
test <- cbind(df_test,df_pred)

mSplinePlot <- ggplot(test,aes(time,Estimate))+
  facet_wrap(~treatment)+
  lapply(seq(1,49,2),function(i) geom_ribbon(aes_string(ymin=paste("Q",i,sep = ""),ymax=paste("Q",100-i,sep = "")),fill=avg_pal[i],alpha=0.5))+
  geom_line(data=df,aes(time,y,group=interaction(treatment,sample)),color="gray20", size=0.25)+
  coord_cartesian(ylim=c(0,250))+
  ylab(~~Area~(cm^2))+
  xlab("Time")+
  theme_light()+
  theme(axis.ticks.length=unit(0.2,"cm"))+
  theme(strip.background=element_rect(fill="gray50",color="gray20"),
        strip.text.x=element_text(size=14,color="white"),
        strip.text.y=element_text(size=14,color="white"))+
  theme(axis.title= element_text(size = 18))+
  theme(axis.text = element_text(size = 14))+
  theme(legend.position='top')
mSplinePlot

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

`Homoskedastic Model`<-modelSims_MM_homo
`Linear Heteroskedasticity`<-modelSims_MM_linear
`Quadratic Heteroskedasticity`<-modelSims_MM_quad
`Spline Heteroskedasticity`<-modelSims_MM_spline
`Exponential Heteroskedasticity`<-modelSims_MM_exp

comparisons<-compareSimModels(`Homoskedastic Model`, `Linear Heteroskedasticity`, `Quadratic Heteroskedasticity`, `Spline Heteroskedasticity`,`Exponential Heteroskedasticity`)

looICPlot<-comparisons[3][[1]]
looICPlot
ggsave("LOO_IC_monomolecular_hetSkeds.png",looICPlot, width = 7.04, height=4.04, dpi=300, bg = "#ffffff")


pSpline_titled<-mSplinePlot+labs(title="Splines")+coord_cartesian(ylim=c(0,250))+
  theme(axis.text = element_text(size = 10),axis.title= element_text(size = 8),
        title = element_text(size=10))
pLinear_titled<-mLinearPlot+labs(title="Linear")+coord_cartesian(ylim=c(0,250))+ylab("")+xlab("")+theme(axis.text = element_text(size = 8),
                                                                                                        axis.title= element_text(size = 8),
                                                                                                        title = element_text(size=10))
pExp_titled<-mExpPlot+labs(title="Exponential")+coord_cartesian(ylim=c(0,250))+ylab("")+xlab("")+theme(axis.text = element_text(size = 8),
                                                                                                       axis.title= element_text(size = 8),
                                                                                                       title = element_text(size=10))
pQuad_titled<-mQuadPlot+labs(title="Quadratic")+coord_cartesian(ylim=c(0,250))+ylab("")+theme(axis.text = element_text(size = 10),
                                                                                              axis.title= element_text(size = 8),
                                                                                              title = element_text(size=10))
pNone_titled<-mHomoPlot+labs(title="Homoskedastic")+coord_cartesian(ylim=c(0,250))+xlab("")+theme(axis.text = element_text(size = 8),
                                                                                                  axis.title= element_text(size = 8),
                                                                                                  title = element_text(size=10))

#looICPlot_to_add<-
looDf<-comparisons[[2]]%>%
  mutate(modelName = str_remove_all(modelName, "Model"))%>%
  mutate(modelName = str_remove_all(modelName, "Heteroskedasticity"))
looDf$modelName <- with(looDf, reorder(modelName, -loo_IC_Mean))

looICPlot_to_add<-looDf%>%
  ggplot()+
  geom_col(aes(x=modelName, y=loo_IC_Mean, fill=modelName))+
  geom_segment(aes(x=modelName, xend=modelName, y=loo_IC_Mean, yend = loo_IC_Mean+loo_IC_SE), size = 1.5)+
  geom_segment(aes(x=modelName, xend=modelName, y=loo_IC_Mean, yend = loo_IC_Mean-loo_IC_SE), size = 1.5)+
  geom_text(aes(x=modelName, y=loo_IC_Mean*1.08, label = paste0(round(loo_IC_Mean, digits=0))), size=3)+
  labs(title = "LOO IC", x="", y="", fill="Model")+
  scale_fill_viridis(discrete = T)+
  theme_minimal() +
  theme(axis.line.y.left = element_line(),
        axis.line.x.bottom = element_line(),
        axis.text.x.bottom = element_blank(),
        axis.text.y = element_blank(),
        title = element_text(size=10),
        legend.title = element_text(size=6.5),
        legend.text = element_text(size=6.5))
looICPlot_to_add
patchesOhoolihan<- (pNone_titled | pLinear_titled | pExp_titled)/(pSpline_titled | pQuad_titled | looICPlot_to_add)
patchesOhoolihan

ggsave("fig2/monomolecularFigure2_badInput.png", patchesOhoolihan, dpi=300, width = 10, height = 5)


