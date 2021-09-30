library(brms)
library(ggplot2)
library(stringr)
library(patchwork)
library(cowplot)
library(viridisLite)

growthSimL <- function(x,a){
  a_r <- a+rnorm(1,mean = 0,sd=0.2)
  return(a_r*x) #1-(1-a) * exp(-b * time)
}
set.seed(345)
x <- 1:25
df <- rbind(
  do.call(rbind,lapply(1:20,function(i) data.frame("sample"=paste0("sample_",i),"treatment"="a","time"=x,"y"=growthSimL(x,1.10),stringsAsFactors = F))),
  do.call(rbind,lapply(1:20,function(i) data.frame("sample"=paste0("sample_",i),"treatment"="b","time"=x,"y"=growthSimL(x,0.95),stringsAsFactors = F)))
)

p <- ggplot(df,aes(time,y,group=interaction(treatment,sample)))+
  geom_line(aes(color=treatment))+
  ylab("Area (cm, simulated)")+
  xlab("Time")+
  labs(title="Linear")+
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


priorLinear <- prior(lognormal(log(1), .25),nlpar = "a") +
  prior(student_t(3,0,5), dpar="sigma") +
  prior(gamma(2,0.1), class="nu")

fitLinear1 <- brm(bf(y ~ a*time, 
                   sigma~time:treatment, 
                   a ~ 0+treatment, 
                   autocor = ~arma(~time|sample:treatment,1,1),nl = T),
                family = student, prior = priorLinear, data = df, iter = 1000, 
                cores = 2, chains = 2, backend = "cmdstanr", #threads = threading(4),
                control = list(adapt_delta = 0.999,max_treedepth = 20),
                inits = function(){list(b_a=rgamma(2,1))})

h <- hypothesis(fitLinear1, "a_treatmenta/a_treatmentb > 1")

plot(h)

probs <- seq(from=99, to=1, by=-2)/100
avg_pal <- plasma(n=length(probs))
df_test <- rbind(data.frame("treatment"="a",time=1:25,sample="new1"),data.frame("treatment"="b",time=1:25,sample="new2"))
df_pred <- predict(fitLinear1,df_test,probs=probs) # add fit_spline in place of fit1 for first checks.
test <- cbind(df_test,df_pred)

p <- ggplot(test,aes(time,Estimate))+
  facet_wrap(~treatment)+
  lapply(seq(1,49,2),function(i) geom_ribbon(aes_string(ymin=paste("Q",i,sep = ""),ymax=paste("Q",100-i,sep = "")),fill=avg_pal[i],alpha=0.5))+
  geom_line(data=df,aes(time,y,group=interaction(treatment,sample)),color="gray20", size=0.3)+
  ylab("Area (cm, simulated)")+
  xlab("Time")+
  labs(title="Linear")+
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
ggsave("Fig2/linearGrowth_linear_ribbon.png",p, width = 7.04, height=4.04, dpi=300)



#***************************************************************************************
################################## Make Data Function ################################## 
#***************************************************************************************

makeData<-function(x.=x,nSamples.=nSamples, a_1.=a_1, a_2.=a_2){
  df<-rbind(
    do.call(rbind,
            lapply(1:nSamples.,
                   function(i) data.frame("sample"=paste0("sample_",i),
                                          "treatment"="a",
                                          "time"=x.,
                                          "y"=growthSimL(x., a_1.),
                                          stringsAsFactors = F))),
    do.call(rbind,
            lapply(1:nSamples.,
                   function(i) data.frame("sample"=paste0("sample_",i),
                                          "treatment"="b",
                                          "time"=x.,
                                          "y"=growthSimL(x., a_2.),
                                          stringsAsFactors = F)))
  )
  return(df)
}
#*******************************************************************************************************************
################################## define modelSims() for Power Law Growth Models ################################## 
#*******************************************************************************************************************

modelSimsL<-function(iterations = 5, sigma = "none", xTime=25, nSamples = 20, a_1 = 1.1, a_2=0.95){
  sigma_<-ifelse(sigma=="linear", "sigma~time+time:treatment,", 
                 ifelse(sigma=="spline", "sigma~s(time,by=treatment),",
                        ifelse(sigma=="quad", "lf(sigma~ time + timeSQ + time:treatment + timeSQ:treatment),",
                               ifelse(sigma=="exp", "sigma~time+time:treatment,", #Not evaluated, doesn't matter
                                      ifelse(sigma=="none", "", 
                                             paste0(sigma, ","))))))
  bfText<-paste0("bf(y ~ a*time,", 
                 sigma_, 
                 "a ~ 0+treatment,autocor = ~arma(~time|sample:treatment,1,1),nl = TRUE)")
  bayesFormula<-eval(parse(text = bfText))
  
  prior_list<-prior(lognormal(log(13), .25),nlpar = "a") +
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
    
    bayesFormula<-bf(y ~ a*time,
                     a ~ 0+treatment,
                     autocor = ~arma(~time|sample:treatment,1,1), nl = TRUE)+
      nlf(sigma ~ c * exp(r * time))+
      lf(c + r ~ 0+treatment)
  }
  
  for (i in 1:iterations){
    cat("\nStarting Iteration ", i, "/",iterations,"\n") 
    iteration_row<-data.frame(iteration = i, elpd_loo = NA,elpd_loo_se=NA, p_loo=NA, p_loo_se=NA, loo_IC=NA,loo_IC_se=NA) #store iteration number
    x<-1:xTime
    dat <- makeData(x.=x,nSamples.=nSamples, a_1.=a_1, a_2.=a_2) #make the data
    if(sigma=="quad"){dat$timeSQ<-dat$time^2}
    fit_none <- brm(bayesFormula,
                    family = student, prior = prior_list, data = dat, iter = 1000,
                    cores = 2, chains = 2, backend = "cmdstanr",
                    control = list(adapt_delta = 0.999,max_treedepth = 20),
                    inits = function(){list(b_a=rgamma(2,1))})
    
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

  modelSimsOutput<-list(metrics_df, summary_df, aPlot)
  return(modelSimsOutput)
}

linearTest<-modelSimsL()
linearTest















