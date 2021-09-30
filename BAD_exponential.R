# a * exp(r * time)
library(brms)
library(ggplot2)
library(stringr)
library(patchwork)
library(cowplot)
library(viridisLite)

growthSimExp <- function(x,a,b,c){
  a_r <- a+rnorm(1,mean = 0,sd=1.5)
  b_r <- b+rnorm(1,mean=0,sd=0.01)
  return(a_r*exp(b_r * x))
}

# x<-1:25
# a<-1
# b<-.15
# y<-a*exp(b*x)
# plot(y)

x <- 1:25
set.seed(432)
df <- rbind(
  do.call(rbind,lapply(1:20,function(i) data.frame("sample"=paste0("sample_",i),
                                                   "treatment"="a","time"=x,
                                                   "y"=growthSimExp(x,17, 0.095),stringsAsFactors = F))),
  do.call(rbind,lapply(1:20,function(i) data.frame("sample"=paste0("sample_",i),
                                                   "treatment"="b","time"=x,
                                                   "y"=growthSimExp(x,16, .085),stringsAsFactors = F)))
)

p <- ggplot(df,aes(time,y,group=interaction(treatment,sample)))+
  geom_line(aes(color=treatment))+
  ylab("Area (cm, simulated)")+
  xlab("Time")+
  labs(title="Exponential")+
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




priorExp <- prior(lognormal(log(15), .25),nlpar = "a") +
  prior(lognormal(log(0.1), .1), nlpar = "b") + 
  prior(student_t(3,0,5), dpar="sigma") +
  prior(gamma(2,0.1), class="nu")

fitExp1 <- brm(bf(y ~ a*exp(b*time), 
                 sigma~time:treatment, 
                 a + b ~ 0+treatment, 
                 autocor = ~arma(~time|sample:treatment,1,1),nl = TRUE),
              family = student, prior = priorExp, data = df, iter = 2000, 
              cores = 2, chains = 2, backend = "cmdstanr", #threads = threading(4),
              control = list(adapt_delta = 0.999,max_treedepth = 20),
              inits = function(){list(b_a=rgamma(2,1),b_b=rgamma(2,1))})

fitExp_spline <- brm(bf(y ~ a*exp(b*time), 
                  sigma~s(time,by=treatment), 
                  a + b ~ 0+treatment, 
                  autocor = ~arma(~time|sample:treatment,1,1),nl = TRUE),
               family = student, prior = priorExp, data = df, iter = 2000, 
               cores = 2, chains = 2, backend = "cmdstanr", #threads = threading(4),
               control = list(adapt_delta = 0.999,max_treedepth = 20),
               inits = function(){list(b_a=rgamma(2,1),b_b=rgamma(2,1))})
save(fitExp_spline, "exponentialData_splineModel.rdata")

h <- hypothesis(fitExp1, "a_treatmenta/a_treatmentb > 1")
h
plot(h)

probs <- seq(from=99, to=1, by=-2)/100
avg_pal <- plasma(n=length(probs))
df_test <- rbind(data.frame("treatment"="a",time=1:25,sample="new1"),data.frame("treatment"="b",time=1:25,sample="new2"))
df_pred <- predict(fitExp1,df_test,probs=probs) # add fit_spline in place of fit1 for first checks.
test <- cbind(df_test,df_pred)

p <- ggplot(test,aes(time,Estimate))+
  facet_wrap(~treatment)+
  lapply(seq(1,49,2),function(i) geom_ribbon(aes_string(ymin=paste("Q",i,sep = ""),ymax=paste("Q",100-i,sep = "")),fill=avg_pal[i],alpha=0.5))+
  geom_line(data=df,aes(time,y,group=interaction(treatment,sample)),color="gray20", size=0.3)+
  ylab("Area (cm, simulated)")+
  xlab("Time")+
  #coord_cartesian(ylim=c(0,750))+
  #lims(y=c(0,500))+
  #scale_y_continuous(limits = c(0,300))+
  theme_light()+
  theme(axis.ticks.length=unit(0.2,"cm"))+
  theme(strip.background=element_rect(fill="gray50",color="gray20"),
        strip.text.x=element_text(size=14,color="white"),
        strip.text.y=element_text(size=14,color="white"))+
  theme(axis.title= element_text(size = 18))+
  theme(axis.text = element_text(size = 14))+
  theme(legend.position='top')
p
ggsave("Fig2/exponential_linear_ribbon.png",p, width = 7.04, height=4.04, dpi=300)


#***************************************************************************************
################################## Make Data Function ################################## 
#***************************************************************************************

makeData<-function(x.=x,nSamples.=nSamples, a_1.=a_1, a_2.=a_2, b_1.= b_1, b_2.=b_2){
  df<-rbind(
    do.call(rbind,
            lapply(1:nSamples.,
                   function(i) data.frame("sample"=paste0("sample_",i),
                                          "treatment"="a",
                                          "time"=x.,
                                          "y"=growthSimExp(x., a_1., b_1.),
                                          stringsAsFactors = F))),
    do.call(rbind,
            lapply(1:nSamples.,
                   function(i) data.frame("sample"=paste0("sample_",i),
                                          "treatment"="b",
                                          "time"=x.,
                                          "y"=growthSimExp(x., a_2., b_2.),
                                          stringsAsFactors = F)))
  )
  return(df)
}
#*******************************************************************************************************************
################################## define modelSims() for Power Law Growth Models ################################## 
#*******************************************************************************************************************

modelSimsExp<-function(iterations = 5, sigma = "none", xTime=25, nSamples = 20, a_1 = 17, a_2=0.095, b_1 = 16, b_2=0.085){
  sigma_<-ifelse(sigma=="linear", "sigma~time+time:treatment,", 
                 ifelse(sigma=="spline", "sigma~s(time,by=treatment),",
                        ifelse(sigma=="quad", "lf(sigma~ time + timeSQ + time:treatment + timeSQ:treatment),",
                               ifelse(sigma=="exp", "sigma~time+time:treatment,", #Not evaluated, doesn't matter
                                      ifelse(sigma=="none", "", 
                                             paste0(sigma, ","))))))
  bfText<-paste0("bf(y ~ a*exp(b*time), ", 
                 sigma_, 
                 "a + b  ~ 0+treatment,autocor = ~arma(~time|sample:treatment,1,1),nl = TRUE)")
  bayesFormula<-eval(parse(text = bfText))
  
  prior_list<-prior(lognormal(log(15), .25),nlpar = "a") +
    prior(lognormal(log(0.1), .1), nlpar = "b") + 
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
      prior(student_t(1, 0, 5), nlpar="c") + 
      prior(student_t(1, 0, 5), nlpar="r")
    
    bayesFormula<-bf(y ~ a*exp(b*time), 
                     a + b ~ 0+treatment,
                     autocor = ~arma(~time|sample:treatment,1,1), nl = TRUE)+
      nlf(sigma ~ c * exp(r * time))+
      lf(c + r ~ 0+treatment)
  }
  
  for (i in 1:iterations){
    cat("\nStarting Iteration ", i, "/",iterations,"\n") 
    iteration_row<-data.frame(iteration = i, elpd_loo = NA,elpd_loo_se=NA, p_loo=NA, p_loo_se=NA, loo_IC=NA,loo_IC_se=NA) #store iteration number
    x<-1:xTime
    dat <- makeData(x.=x,nSamples.=nSamples, a_1.=a_1, a_2.=a_2, b_1.= b_1, b_2.=b_2) #make the data
    if(sigma=="quad"){dat$timeSQ<-dat$time^2}
    fit_none <- brm(bayesFormula,
                    family = student, prior = prior_list, data = dat, iter = 1000,
                    cores = 2, chains = 2, backend = "cmdstanr",
                    control = list(adapt_delta = 0.999,max_treedepth = 20),
                    inits = function(){list(b_a=rgamma(2,1),b_b=rgamma(2,1))})
    
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
  
  modelSimsOutput<-list(metrics_df, summary_df, patchPlot)
  return(modelSimsOutput)
}

modelSimsExp(iterations=2)
