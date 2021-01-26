library(brms)
library(ggplot2)
library(stringr)
library(patchwork)
library(cowplot)

setwd("~/Danforth/Datasci/Bayesian_adaptive")

#*************************************************************************************************
# Simulated Data
#*************************************************************************************************
growthSim <- function(x,phi1,phi2,phi3){
  phi1_r <- phi1+rnorm(1,mean = 0,sd=25)
  phi2_r <- phi2+rnorm(1,mean=0,sd=1)
  phi3_r <- phi3+rnorm(1,mean=0,sd=.2)
  return(phi1_r/(1+exp((phi2_r-x)/phi3_r)))
}

# x <- 1:25
# df <- rbind(
#   do.call(rbind,lapply(1:20,function(i) data.frame("sample"=paste0("sample_",i),"treatment"="a","time"=x,"y"=growthSim(x,200,13,3),stringsAsFactors = F))),
#   do.call(rbind,lapply(1:20,function(i) data.frame("sample"=paste0("sample_",i),"treatment"="b","time"=x,"y"=growthSim(x,160,13,3.5),stringsAsFactors = F)))
# )
# write.csv(df,"simulated_data_v1.csv",row.names = F,quote = F)
df <- read.csv("simulated_data_v1.csv",stringsAsFactors = F)

f1a <- ggplot(df,aes(time,y,group=interaction(treatment,sample)))+
  geom_line(aes(color=treatment))+
  ylab("Area (cm, simulated)")+
  xlab("Time")+
  scale_y_continuous(limits = c(0,300))+
  theme_light()+
  theme(axis.ticks.length=unit(0.2,"cm"))+
  theme(strip.background=element_rect(fill="gray50",color="gray20"),
        strip.text.x=element_text(size=14,color="white"),
        strip.text.y=element_text(size=14,color="white"))+
  theme(axis.title= element_text(size = 18))+
  theme(axis.text = element_text(size = 14))+
  theme(legend.position='top')

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


ggsave("Fig1/bayes_simulated_f1a.png",f1a,width = 5.04,height=4.04,dpi=300)

f1b <- f1a + 
  facet_wrap(~treatment)+
  geom_vline(data=df[df$treatment == "a",],aes(xintercept=13),color="gray40",linetype="dashed")+
  geom_vline(data=df[df$treatment == "b",],aes(xintercept=13),color="gray40",linetype="dashed")+
  geom_hline(data=df[df$treatment == "a",],aes(yintercept=200),color="gray40",linetype="dashed")+
  geom_hline(data=df[df$treatment == "b",],aes(yintercept=160),color="gray40",linetype="dashed")+
  geom_abline(data=df[df$treatment == "a",],aes(intercept=(0.5*200)-(200*(1/3)*(1/4))*13,slope=200*(1/3)*(1/4)),color="gray40",linetype="dashed")+
  geom_abline(data=df[df$treatment == "b",],aes(intercept=(0.5*160)-(160*(1/3.5)*(1/4))*13,slope=160*(1/3.5)*(1/4)),color="gray40",linetype="dashed")
f1b
ggsave("Fig1/bayes_simulated_f1b.png",f1b,width = 7.04,height=4.04,dpi=300)

prior1 <- prior(lognormal(log(130), .25),nlpar = "phi1",coef="treatmenta") +
          prior(lognormal(log(130), .25),nlpar = "phi1",coef="treatmentb") +
          prior(lognormal(log(12), .25), nlpar = "phi2",coef="treatmenta") + 
          prior(lognormal(log(12), .25), nlpar = "phi2",coef="treatmentb") + 
          prior(lognormal(log(3), .25), nlpar = "phi3",coef="treatmenta") + 
          prior(lognormal(log(3), .25), nlpar = "phi3",coef="treatmentb") +
          prior(student_t(3,0,5), dpar="sigma") +
          prior(gamma(2,0.1), class="nu")

fit1 <- brm(bf(y ~ phi1/(1+exp((phi2-time)/phi3)), 
               sigma~time:treatment, 
               phi1 + phi2 + phi3 ~ 0+treatment, 
               autocor = ~arma(~time|sample:treatment,1,1),nl = TRUE),
            family = student, prior = prior1, data = df, iter = 10000, 
            cores = 4, chains = 4, backend = "cmdstanr", threads = threading(4),
            control = list(adapt_delta = 0.999,max_treedepth = 20),
            inits = function(){list(b_phi1=rgamma(2,1),b_phi2=rgamma(2,1),b_phi3=rgamma(2,1))})

plot(conditional_effects(fit1),points=T)
plot(posterior_interval(fit1))


hist(rlnorm(2000,log(130), .25),breaks = 50)
hist(rlnorm(2000,log(12), .25),breaks = 50)
hist(rlnorm(2000,log(3), .25),breaks = 50)

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
#**********************************************
#*add step to find each statistical distribution of the priors so it's more robust
#**********************************************
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


data("LakeHuron")
LakeHuron <- as.data.frame(LakeHuron)
get_prior(bf(x ~ arma(p = 2, q = 1)), data = LakeHuron)
fit <- brm(x ~ arma(p = 2, q = 1), data = LakeHuron)
summary(fit)

#*************************************************************************************************
# Trace plots of bayesian priors 
#*************************************************************************************************
samples <- data.frame(brms::posterior_samples(fit1, add_chain=T))
library(ggplot2)
# b_sigma_time x treatmenta plot 
ggplot(samples, aes(color=chain, x = iter, y = b_sigma_time.treatmenta)) + 
  geom_line()

#*************************************************************************************************
# Real Data
#*************************************************************************************************
setwd("~/Danforth/Datasci/Bayesian_adaptive")

area_convert <- 13.2*3.7/46856

sv_shapes <- read.table("b51632508ce74bdd9ceffcd073b871e6-shapes_outliersRM.tsv",sep="\t",header=T,stringsAsFactors = F)
sv_shapes <- sv_shapes[sv_shapes$Drought == "AAA" & sv_shapes$Microbes %in% c("Control","SynCom A"),]
sv_shapes$area_c <- sv_shapes$area * area_convert

df <- aggregate(data=sv_shapes[sv_shapes$DAP > 4,],area_c~Microbes+Barcodes+DAP,FUN=function(i)mean(i))
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
            cores = 4, chains = 4, backend = "cmdstanr", threads = threading(4),
            control = list(adapt_delta = 0.999,max_treedepth = 20),
            inits = function(){list(b_phi1=rgamma(3,1),b_phi2=rgamma(3,1),b_phi3=rgamma(3,1))})

hypothesis(fit1, c(
                    "sigma_DAP:MicrobesControl > 0",
                    "sigma_DAP:MicrobesSynComA > 0",
                    "sigma_DAP:MicrobesSynComB > 0",
                    "phi1_MicrobesSynComA/phi1_MicrobesSynComB > 1",
                    "phi2_MicrobesControl/phi2_MicrobesSynComA > 1",
                    "phi3_MicrobesSynComA/phi3_MicrobesSynComB > 1"
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





hist(rbeta(5000,91,11),breaks = 100)
median(rbeta(5000,91,11))

9/10





