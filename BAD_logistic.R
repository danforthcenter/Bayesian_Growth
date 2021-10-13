#issue 3
library(brms)
library(ggplot2)
library(stringr)
library(patchwork)
library(cowplot)


#setwd('/home/jberry/Danforth/Datasci/Bayesian_Growth/')

#*************************************************************************************************
# Simulated Data
#*************************************************************************************************

growthSim <- function(x,phi1,phi2,phi3){ # make function to do a bunch of growthSims and check how well they 
  phi1_r <- phi1+rnorm(1,mean = 0,sd=25) #can be modeled/how well their phi params can be recovered.
  phi2_r <- phi2+rnorm(1,mean=0,sd=1)
  phi3_r <- phi3+rnorm(1,mean=0,sd=.2)
  return(phi1_r / (1 + exp( (phi2_r-x)/phi3_r) ))
}
set.seed(123)
x <- 1:25
simdf <- rbind(
  do.call(rbind,
          lapply(1:20,
                 function(i) data.frame("sample"=paste0("sample_",i),"treatment"="a","time"=x,"y"=growthSim(x,200,13,3),stringsAsFactors = F))),
  do.call(rbind,
          lapply(1:20,
                 function(i) data.frame("sample"=paste0("sample_",i),"treatment"="b","time"=x,"y"=growthSim(x,160,13,3.5),stringsAsFactors = F)))
)
#simulated data is taking a vector of x values (1 to 25) 
  # and using growthSim to make a simulated y value that adds random noise to the 3 numbers we provide then uses those 
  # along with the iterating X value, the higher the X value gets the more it subtracts from the numerator of the exp()
  # So the exponent gets smaller, which makes the overall denominator smaller, thus exponential model of growth.
  # phi1_r/(1+exp((phi2_r-x)/phi3_r)) <- right there.
# As for the do.call()s, they are just making a df with 20 samples (lapplying 1:20) 
  #to the function we have to make a row for each value of x.
  #Wrapped in rbind by do.call it makes 25 rows per iteration in 1:20,
  # That is happening once with one set of growthSim params, and once with another set of params, for total of 1k rows.


head(df)
# write.csv(df,"simulated_data_v1.csv",row.names = F,quote = F)
df <- read.csv("simulated_data_v1.csv",stringsAsFactors = F)

f1a <- ggplot(df,aes(time,y,group=interaction(treatment,sample)))+
  geom_line(aes(color=treatment))+
  ylab(~~Area~(cm^2))+
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
f1a
p1 <- ggplot(data=data.frame("y"=0:300,"x"=dlnorm(0:300,log(170),.25)),aes(y,x))+
  geom_area(fill="gray80",color="gray40")+
  coord_flip()+
  scale_x_continuous(limits = c(0,300))+
  theme_void()

p2 <- ggplot(data=data.frame("y"=0:25,"x"=dlnorm(0:25,log(12),.25)),aes(y,x))+
  geom_area(fill="gray80",color="gray40")+
  scale_y_reverse()+
  theme_void()

(f1a+p1+p2+plot_spacer())+plot_layout(nrow=2,widths = c(3,1),heights = c(3,1))



f1b <- f1a + 
  facet_wrap(~treatment)+
  geom_vline(data=df[df$treatment == "a",],aes(xintercept=13),color="gray20")+
  geom_vline(data=df[df$treatment == "b",],aes(xintercept=13),color="gray20")+
  geom_hline(data=df[df$treatment == "a",],aes(yintercept=200),color="gray20")+
  geom_hline(data=df[df$treatment == "b",],aes(yintercept=160),color="gray20")+
  geom_abline(data=df[df$treatment == "a",],aes(intercept=(0.5*200)-(200*(1/3)*(1/4))*13,slope=200*(1/3)*(1/4)),color="gray20")+
  geom_abline(data=df[df$treatment == "b",],aes(intercept=(0.5*160)-(160*(1/3.5)*(1/4))*13,slope=160*(1/3.5)*(1/4)),color="gray20")
f1b
ggsave("Fig1/bayes_simulated_f1b.png",f1b,width = 7.04,height=4.04,dpi=300)

f1c<-f1a + 
  facet_wrap(~treatment)+
  theme(legend.position = "none")
f1c
ggsave("Fig1/bayes_simulated_f1c.png",f1c,width = 7.04,height=4.04,dpi=300)

#*************************************************************************************************
# Viewing/modeling heteroscedasticity
#*************************************************************************************************
var_df <- aggregate(data=df,y~treatment+time,FUN=function(i)var(i))
# equivalent to:
# df%>%group_by(treatment, time)%>%summarize(y = var(y))
  #nchar("aggregate(data=df,y~treatment+time,FUN=function(i)var(i))")
  #nchar("df%>%group_by(treatment,time)%>%summarize(y=var(y))")
    # 6 less characters with dplyr! I'm taking this as a moral victory. (this is in jest, please don't hate me)
# anyway, point is to summarize the data over treatment and timepoint, then get the variance at each time for each treatment.
aspline <- data.frame(spline(var_df$y[var_df$treatment == "a"]),"treatment"="a")
bspline <- data.frame(spline(var_df$y[var_df$treatment == "b"]),"treatment"="b")

#splines
  # Splines are a way of fitting a line to some distribution when we aren't sure where it came from.
  # They are piecewise polynomials for fitting a smooth curve to complicated data points.
  # splines try to match first and second derivatives of the points they model
ggplot(data=var_df,aes(time,y))+
  geom_line(aes(color=treatment))+
  geom_line(data=bspline,aes(x,y), linetype=5)+ #the spline is doing a very good job of covering the b treatment.
  geom_line(data=aspline,aes(x,y), linetype=5) #new

#***** Spline Var

#there is a lot going on here
# prior() function is from brms and it defines prior distributions. Each line here is adding a prior distribution into the
  # "brmsprior"  "data.frame" object `priotr_spline`.  We have 6 log normal prior distributions 
  #(corresponding to combinations of 3 PHIs and 2 treatment groups), as well as a T distribution and a gamma.
  # all these distribution functions are also from BRMS and are for use in making prior distributions like this.
prior_spline <- prior(lognormal(log(130), .25),nlpar = "phi1",coef="treatmenta") +
          prior(lognormal(log(130), .25),nlpar = "phi1",coef="treatmentb") +
          prior(lognormal(log(12), .25), nlpar = "phi2",coef="treatmenta") +
          prior(lognormal(log(12), .25), nlpar = "phi2",coef="treatmentb") +
          prior(lognormal(log(3), .25), nlpar = "phi3",coef="treatmenta") +
          prior(lognormal(log(3), .25), nlpar = "phi3",coef="treatmentb") +
          prior(student_t(3,0,5), dpar="sigma") +
          prior(gamma(2,0.1), class="nu")

# Things holding me up on this:
  # No real understanding of brm() function, bf() function, arma() function, s() function,
  # use of priors in brm, what threading does, control parameter of brm(), inits parameter.


# This keeps throwing a fatal error when I try it. Seems to have been about the cores/chains,
  # Weird behavior with fatal errors while loading brms package afterwards
fit_spline <- brm(bf(y ~ phi1/(1+exp((phi2-time)/phi3)),
               sigma~s(time,by=treatment), #spline var vs linear var vs no var changes this line,
               # splines will probably outperform everything since it is so flexible.
               # the most Heteroscedastic of the heteroscedastic options
               #drawback is time, these take a while.
               phi1 + phi2 + phi3 ~ 0+treatment,
               autocor = ~arma(~time|sample:treatment,1,1), #ar1 vs ar2 would be possible avenues to explore
               nl = TRUE),
            family = student, prior = prior_spline, data = df, iter = 1000, #10k is serious overkill, lets us be certain of global max
            cores = 4, # x cores from machine, all about the physical processor (my mac has 2 cores)
            chains = 4, # MCMC monte carlo markov chain, iterative algorithm going x times, 
            # hoping they all "converge" to the same point.
                        # If we keep this the same as cores then it runs a chain per core in parallel.
            backend = "cmdstanr",  #should break here on cmdstanr
            #threads = threading(4), # threading 4 times for each core, adds speed.
            control = list(adapt_delta = 0.999, # black box that helps convergence
                           max_treedepth = 20), #tree depth to look backwards during MCMC
            inits = function(){list(b_phi1=rgamma(2,1),b_phi2=rgamma(2,1),b_phi3=rgamma(2,1))}) 
                                # controlling our likelihood space for MCMC, 
                                # randomly sampling from a gamma makes sure it is positive, keeps it from breaking

post_spline <- data.frame(posterior_summary(fit_spline),stringsAsFactors = F)
save(fit_spline,prior_spline,post_spline,file ="logistic_model_Var.rdata")

#***** Linear Var
prior_linear <- prior(lognormal(log(130), .25),nlpar = "phi1",coef="treatmenta") +
          prior(lognormal(log(130), .25),nlpar = "phi1",coef="treatmentb") +
          prior(lognormal(log(12), .25), nlpar = "phi2",coef="treatmenta") +
          prior(lognormal(log(12), .25), nlpar = "phi2",coef="treatmentb") +
          prior(lognormal(log(3), .25), nlpar = "phi3",coef="treatmenta") +
          prior(lognormal(log(3), .25), nlpar = "phi3",coef="treatmentb") +
          prior(student_t(3,0,5), dpar="sigma") +
          prior(gamma(2,0.1), class="nu")

fit_linear <- brm(bf(y ~ phi1/(1+exp((phi2-time)/phi3)),
               sigma~time+time:treatment, #linear [no s()] variation
               phi1 + phi2 + phi3 ~ 0+treatment,
               autocor = ~arma(~time|sample:treatment,1,1),nl = TRUE),
            family = student, prior = prior_linear, data = df, iter = 10000,
            cores = 4, chains = 4, backend = "cmdstanr", threads = threading(4),
            control = list(adapt_delta = 0.999,max_treedepth = 20),
            inits = function(){list(b_phi1=rgamma(2,1),b_phi2=rgamma(2,1),b_phi3=rgamma(2,1))})

post_linear <- data.frame(posterior_summary(fit_linear),stringsAsFactors = F)
save(fit_linear,prior_linear,post_linear,file ="logistic_model_linearVar.rdata")

#***** No Var
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
            family = student, prior = prior_none, data = df, iter = 10000,
            cores = 2, chains = 2, backend = "cmdstanr", #threads = threading(4),
            control = list(adapt_delta = 0.999,max_treedepth = 20),
            inits = function(){list(b_phi1=rgamma(2,1),b_phi2=rgamma(2,1),b_phi3=rgamma(2,1))})

post_none <- data.frame(posterior_summary(fit_none),stringsAsFactors = F)
save(fit_none,prior_none,post_none,file ="logistic_model_noneVar.rdata")


#*************************************************************************************************
# Bayesian credible intervals
#*************************************************************************************************
probs <- seq(from=99, to=1, by=-2)/100
avg_pal <- turbo(n=length(probs))
df_test <- rbind(data.frame("treatment"="a",time=1:25,sample="new1"),data.frame("treatment"="b",time=1:25,sample="new2"))
df_pred <- predict(fit_spline,df_test,probs=probs) # add fit_spline in place of fit1 for first checks.
test <- cbind(df_test,df_pred)

p <- ggplot(test,aes(time,Estimate))+
  facet_wrap(~treatment)+
  lapply(seq(1,49,2),function(i) geom_ribbon(aes_string(ymin=paste("Q",i,sep = ""),ymax=paste("Q",100-i,sep = "")),fill=avg_pal[i],alpha=0.5))+
 # geom_line(data=df,aes(time,y,group=interaction(treatment,sample)),color="gray20")+
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
ggsave("Fig1/bayes_logistic_credibleIntervals.png",p,width = 7.04,height=4.04,dpi=300)



f1c <- f1b +
  geom_vline(data=df[df$treatment == "a",],aes(xintercept=post["b_phi2_treatmenta","Estimate"]),color="gray40",linetype="dashed")+
  geom_vline(data=df[df$treatment == "b",],aes(xintercept=post["b_phi2_treatmentb","Estimate"]),color="gray40",linetype="dashed")+
  geom_hline(data=df[df$treatment == "a",],aes(yintercept=post["b_phi1_treatmenta","Estimate"]),color="gray40",linetype="dashed")+
  geom_hline(data=df[df$treatment == "b",],aes(yintercept=post["b_phi1_treatmentb","Estimate"]),color="gray40",linetype="dashed")+
  geom_abline(data=df[df$treatment == "a",],aes(intercept=(0.5*post["b_phi1_treatmenta","Estimate"])-(post["b_phi1_treatmenta","Estimate"]*(1/post["b_phi3_treatmenta","Estimate"])*(1/4))*post["b_phi2_treatmenta","Estimate"],slope=post["b_phi1_treatmenta","Estimate"]*(1/post["b_phi3_treatmenta","Estimate"])*(1/4)),color="gray40",linetype="dashed")+
  geom_abline(data=df[df$treatment == "b",],aes(intercept=(0.5*post["b_phi1_treatmentb","Estimate"])-(post["b_phi1_treatmentb","Estimate"]*(1/post["b_phi3_treatmentb","Estimate"])*(1/4))*post["b_phi2_treatmentb","Estimate"],slope=post["b_phi1_treatmentb","Estimate"]*(1/post["b_phi3_treatmentb","Estimate"])*(1/4)),color="gray40",linetype="dashed")
f1c

#*************************************************************************************************
# Pull out specific day to do NHST
#*************************************************************************************************
posterior_smooths(fit1)
df_fit <- predict(fit1,newdata=df,scale="response")
test <- cbind(df,df_fit)
head(test)
plot(test$y,test$y-test$Estimate)

#*************************************************************************************************
# Bayesian updating using posteriors as priors
#*************************************************************************************************
post <- data.frame(posterior_summary(fit1),stringsAsFactors = F)
post$n_sd <- post$Est.Error*1.96
with(post,(exp(Est.Error^2)-1)*exp(2*Estimate+Est.Error^2))
new_priors <- prior1

apply(prior1,MARGIN = 1, function(i) paste0("normal(",post[gsub("__","",paste0(i[c(2,7,3)],collapse = "_")),"Estimate"],",",as.numeric(post[gsub("__","",paste0(i[c(2,7,3)],collapse = "_")),"Est.Error"])*1.96,")"))

new_priors$prior <- c(as.character(unlist(apply(post[stringr::str_detect(rownames(post),"b_"),],MARGIN = 1,function(i) paste0("normal(",i[1],",",i[5],")")))),paste0("student_t(3,",post["sigma",1],",",post["sigma",5]))
new_fit <- update(fit1,newdata = df[df$time %in% 15:17,],prior=new_priors, cores=4)


#*************************************************************************************************
# (NON-FUNCTIONAL) add step to find each statistical distribution of the priors so it's more robust
#*************************************************************************************************
uupap <- function(fit,initial_priors,newdata){
  post <- data.frame(posterior_summary(fit),stringsAsFactors = F)
  post$n_sd <- post$Est.Error*1.96
  new_priors <- initial_priors
  new_priors$prior <- apply(prior1,MARGIN = 1, function(i) paste0("normal(",post[gsub("__","",paste0(i[c(2,7,3)],collapse = "_")),"Estimate"],",",as.numeric(post[gsub("__","",paste0(i[c(2,7,3)],collapse = "_")),"Est.Error"])*1.96,")"))
  new_fit <- update(fit,newdata = newdata,prior=new_priors, cores=4)
  return(new_fit)
}
fit2 <- uupap(fit1,prior1,df[df$time %in% 15:17,])

prior_summary(fit1)

posterior_summary(fit1)

fit2 <- update(fit1,newdata = df[df$time %in% 1:17,], cores=4)
fit3 <- update(fit2,newdata = df[df$time %in% 1:19,], cores=4)
fit4 <- update(fit3,newdata = df[df$time %in% 1:21,], cores=4)
fit5 <- update(fit4,newdata = df[df$time %in% 1:23,], cores=4)
fit6 <- update(fit5,newdata = df[df$time %in% 1:25,], cores=4)

hypothesis(fit1, "phi1_treatmenta/phi1_treatmentb > 1")
hypothesis(fit2, "phi1_treatmenta/phi1_treatmentb > 1")
hypothesis(fit3, "phi1_treatmenta/phi1_treatmentb > 1")
hypothesis(fit4, "phi1_treatmenta/phi1_treatmentb > 1")
hypothesis(fit5, "phi1_treatmenta/phi1_treatmentb > 1")
hypothesis(fit6, "phi1_treatmenta/phi1_treatmentb > 1")


#*************************************************************************************************
# Trace plots of bayesian posterior estimations
#*************************************************************************************************
samples <- data.frame(brms::posterior_samples(fit1, add_chain=T))
model_vars <- names(samples)

list_of_plots <- lapply(model_vars[-((length(model_vars)-2):length(model_vars))],function(i){
  ggplot(samples, aes_string(color="chain", x = "iter", y = i)) +
    geom_line(alpha=0.5)+
    ggtitle(i)+
    ylab("Estimate")+
    theme_light()+
    theme(axis.ticks.length=unit(0.2,"cm"))+
    theme(axis.title= element_text(size = 18))+
    theme(axis.text = element_text(size = 14))
})
gridExtra::grid.arrange(grobs=list_of_plots)


#*************************************************************************************************
# Real Data
#*************************************************************************************************
setwd("~/Danforth/Datasci/Bayesian_adaptive")

area_convert <- 13.2*3.7/46856

sv_shapes <- read.table("b51632508ce74bdd9ceffcd073b871e6-shapes_outliersRM.tsv",sep="\t",header=T,stringsAsFactors = F)
sv_shapes <- sv_shapes[sv_shapes$Drought == "AAA" & sv_shapes$Microbes %in% c("Control","SynCom A"),] #add "SynCom B"
sv_shapes$area_c <- sv_shapes$area * area_convert

df <- aggregate(data=sv_shapes[sv_shapes$DAP > 4,], area_c ~ Microbes+Barcodes+DAP, FUN=function(i)mean(i))
df <- df[(df$Barcodes %in% df$Barcodes[df$DAP == 25]) & (df$Barcodes %in% df$Barcodes[df$DAP == 5]),]
df$Microbes <- gsub(" ","",df$Microbes)

ggplot(df,aes(DAP,area_c))+
  facet_wrap(~Microbes)+
  geom_smooth(aes(color=Microbes,group=Barcodes))+
  geom_smooth(color="gray20")+
  theme_light()

prior1 <- 
  prior(lognormal(log(125), .25),nlpar = "phi1",coef="MicrobesControl") +
  prior(lognormal(log(125), .25),nlpar = "phi1",coef="MicrobesSynComA") +
  prior(lognormal(log(15), .25), nlpar = "phi2",coef="MicrobesControl") +
  prior(lognormal(log(15), .25), nlpar = "phi2",coef="MicrobesSynComA") +
  prior(lognormal(log(5), .25), nlpar = "phi3",coef="MicrobesControl") +
  prior(lognormal(log(5), .25), nlpar = "phi3",coef="MicrobesSynComA") +
  prior(student_t(3,0,5), dpar="sigma") +
  prior(gamma(2,0.1), class="nu")

fit1 <- brm(bf(area_c ~ phi1/(1+exp((phi2-DAP)/phi3)),
               sigma~DAP:Microbes,
               phi1 + phi2 + phi3 ~ 0+Microbes,
               autocor = ~arma(~DAP|Barcodes:Microbes,1,1),nl = TRUE),
            family = student, prior = prior1, data = df, iter = 10000,
            cores = 2, chains = 4, backend = "cmdstanr", #threads = threading(4),
            control = list(adapt_delta = 0.999,max_treedepth = 20),
            inits = function(){list(b_phi1=rgamma(3,1),b_phi2=rgamma(3,1),b_phi3=rgamma(3,1))})

hypothesis(fit1, c("sigma_DAP:MicrobesControl > 0",
                    "sigma_DAP:MicrobesSynComA > 0",
                    #"sigma_DAP:MicrobesSynComB > 0", #this one was for "B", which is not in the subsetted data
                    #"phi1_MicrobesSynComA/phi1_MicrobesSynComB > 1",#this one was for "B", which is not in the subsetted data
                    "phi2_MicrobesControl/phi2_MicrobesSynComA > 1"#,
                    #"phi3_MicrobesSynComA/phi3_MicrobesSynComB > 1"#this one was for "B", which is not in the subsetted data
                  ))

fit2 <- update(fit1,newdata = df[df$DAP %in% 1:17,], cores=4)
fit3 <- update(fit2,newdata = df[df$DAP %in% 1:19,], cores=4)
fit4 <- update(fit3,newdata = df[df$DAP %in% 1:21,], cores=4)
fit5 <- update(fit4,newdata = df[df$DAP %in% 1:23,], cores=4)
fit6 <- update(fit5,newdata = df[df$DAP %in% 1:25,], cores=4)

hypothesis(fit1, "phi1_MicrobesSynComB/phi1_MicrobesControl < 1")
hypothesis(fit2, "phi1_MicrobesSynComB/phi1_MicrobesControl < 1")
hypothesis(fit3, "phi1_MicrobesSynComB/phi1_MicrobesControl < 1")
hypothesis(fit4, "phi1_MicrobesSynComB/phi1_MicrobesControl < 1")
hypothesis(fit5, "phi1_MicrobesSynComB/phi1_MicrobesControl < 1")
hypothesis(fit6, "phi1_MicrobesSynComB/phi1_MicrobesControl < 1")





