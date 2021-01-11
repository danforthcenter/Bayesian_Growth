library(brms)
library(ggplot2)
library(stringr)
library(patchwork)
library(cowplot)

setwd("~/Danforth/Datasci/Bayesian_adaptive")

#*************************************************************************************************
# Simulated Data
#*************************************************************************************************
growthSim <- function(x,a,b,c){
  a_r <- a+rnorm(1,mean = 0,sd=10)
  b_r <- b+rnorm(1,mean=0,sd=2)
  c_r <- c+rnorm(1,mean=0,sd=.035)
  return(a_r*exp(-b_r*exp(-c_r*x)))
}

x <- 1:25
df <- rbind(
  do.call(rbind,lapply(1:20,function(i) data.frame("sample"=paste0("sample_",i),"treatment"="a","time"=x,"y"=growthSim(x,200,13,.2),stringsAsFactors = F))),
  do.call(rbind,lapply(1:20,function(i) data.frame("sample"=paste0("sample_",i),"treatment"="b","time"=x,"y"=growthSim(x,160,13,.25),stringsAsFactors = F)))
)

p <- ggplot(df,aes(time,y,group=interaction(treatment,sample)))+
  geom_line(aes(color=treatment))+
  ylab("Area (cm, simulated)")+
  xlab("Time")+
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




prior1 <- prior(lognormal(log(130), .25),nlpar = "a") +
  prior(lognormal(log(12), .25), nlpar = "b") + 
  prior(lognormal(log(1.2), .25), nlpar = "c") + 
  prior(student_t(3,0,5), dpar="sigma") +
  prior(gamma(2,0.1), class="nu")

fit1 <- brm(bf(y ~ a*exp(-b*exp(-c*time)), 
               sigma~time:treatment, 
               a + b + c ~ 0+treatment, 
               autocor = ~arma(~time|sample:treatment,1,1),nl = TRUE),
            family = student, prior = prior1, data = df, iter = 10000, 
            cores = 4, chains = 4, backend = "cmdstanr", threads = threading(4),
            control = list(adapt_delta = 0.999,max_treedepth = 20),
            inits = function(){list(b_a=rgamma(2,1),b_b=rgamma(2,1),b_c=rgamma(2,1))})
