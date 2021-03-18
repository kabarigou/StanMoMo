library(tidyr)
library(dplyr)
library(tibble)
library(matrixcalc)
library(StanMoMo)
library(tidyverse)

#This script details how to reproduce the results of our paper of 
#Section 6: Impact of Covid-type effect

#The first part of the script explains how to compute model weights when data 
#are perturbed as considered in our paper. The case without perturbations can
#be obtained similarly.

years<-1979:2018
years.fit<-years
  
  
deathFR<-load_HMD_data('FRATNP', 'Deaths_1x1', years, 50:90, "Male")$mat
exposureFR<-load_HMD_data('FRATNP', 'Exposures_1x1', years, 50:90, "Male")$mat

#Percentages of perturbations
percentage1<-0.10
percentage2<--0.05

#Perturbations applied to the data on last three years
deathFR[,"2016"]<-deathFR[,"2016"]*(1+percentage1)
deathFR[,"2017"]<-deathFR[,"2017"]*(1+percentage1)
deathFR[,"2018"]<-deathFR[,"2018"]*(1+percentage2)

# Just as Section 5 (application to real mortality data), we determine the weights
# from stacking, pseudo-BMA and BMA
#First, stacking and pseudo-BMA with training and validation sets and then BMA with the whole dataset

library(parallel)
ages.fit<-50:90
years.fit<-1979:2018
iter<-2000
validationyears<-10
forecastyears<-10

samplingfunction<-function(x){
  if (x==1) res<-lc_stan(death = deathFR,exposure=exposureFR, validation=validationyears,forecast = forecastyears, family = "nb",chains=1,cores=1,iter=iter)
  else if (x==2) res<-rh_stan(death = deathFR,exposure=exposureFR, validation=validationyears,forecast = forecastyears, family = "nb",chains=1,cores=1,iter=iter)
  else if (x==3) res<-apc_stan(death = deathFR,exposure=exposureFR, validation=validationyears,forecast = forecastyears, family = "nb",chains=1,cores=1,iter=iter)
  else if (x==4) res<-cbd_stan(death = deathFR,exposure=exposureFR, age=ages.fit, 
                               validation=validationyears,forecast = forecastyears, family = "nb",chains=1,cores=1,iter=iter)
  else if (x==5) res<-m6_stan(death = deathFR,exposure=exposureFR, age=ages.fit, 
                              validation=validationyears,forecast = forecastyears, family = "nb",chains=1,cores=1,iter=iter)
}

cl <- makeCluster(5)
clusterExport(cl,c('deathFR','exposureFR','forecastyears','validationyears','iter','ages.fit','lc_stan','rh_stan','apc_stan','cbd_stan','m6_stan'))
system.time({out <- parLapply(cl, c(1:5),samplingfunction)})

stopCluster(cl)

#Weights from stacking and Pseudo-BMA
model_weights<-mortality_weights(out)
print(model_weights)

#Weigths from the standard BMA approach via Bridgesampling
#From our empirical analysis, we found that bridgesampling requires
#a large sample. Therefore, we used 4 chains per model.

validationyears<-0
forecast<-10
nchains<-4
ncores<-4
iter<-2000

samplingfunction2<-function(x){
  if (x==1) {
    res<-lc_stan(death = deathFR,exposure=exposureFR, validation=validationyears,forecast = forecastyears, family = "nb",chains=nchains,cores=ncores,iter=iter)
    logml <- bridgesampling::bridge_sampler(res, silent = TRUE)
  }
  else if (x==2) {
    res<-rh_stan(death = deathFR,exposure=exposureFR, validation=validationyears,forecast = forecastyears, family = "nb",chains=nchains,cores=ncores,iter=iter)
    logml <- bridgesampling::bridge_sampler(res, silent = TRUE)
  }
  else if (x==3) {
    res<-apc_stan(death = deathFR,exposure=exposureFR, validation=validationyears,forecast = forecastyears, family = "nb",chains=nchains,cores=ncores,iter=iter)
    logml <- bridgesampling::bridge_sampler(res, silent = TRUE)
  }
  else if (x==4){
    res<-cbd_stan(death = deathFR,exposure=exposureFR, age=ages.fit,validation=validationyears,forecast = forecastyears, family = "nb",chains=nchains,cores=ncores,iter=iter)
    logml <- bridgesampling::bridge_sampler(res, silent = TRUE)
  } 
  else if (x==5) {
    res<-m6_stan(death = deathFR,exposure=exposureFR, age=ages.fit, 
                 validation=validationyears,forecast = forecastyears, family = "nb",chains=nchains,cores=ncores,iter=iter)
    logml <- bridgesampling::bridge_sampler(res, silent = TRUE)
  }
  return(list(stan_output = res, logml = logml))
}

cl <- makeCluster(20)
clusterExport(cl,c('deathFR','exposureFR','forecastyears','validationyears','iter','ages.fit','lc_stan','rh_stan','apc_stan','cbd_stan','m6_stan','ncores','nchains'))
system.time({out <- parLapply(cl, c(1:5),samplingfunction2)})


#The BMA weights can directly be obtained via the post_prob function of the `bridgesampling' package
post1 <- post_prob(out[[1]]$logml, out[[2]]$logml,out[[3]]$logml,out[[4]]$logml,out[[5]]$logml)
print(post1)

#In case bridgesampling still does not converge, we also implemented the Harmonic Mean Estimator as proxy.
mortality_models <- c("lc", "rh", "apc", "cbd", "m6")
post<-compute_weights_BMA(out,mortality_models)

# All posterior distributions of all parameters can be extracted with the 'extract' function of the 'rstan' package. 
out2<-lapply(out, `[[`, 1)   

params<-lapply(out2,rstan::extract)


#Hereafter are the weights obtained during our analysis (order: BMA, stacking, pseudo-BMA)
#The user should obtain similar weights.
#Here below, we show how to generate the figures of the paper.

#With perturbations
weights1<-c(0,1,0,0,0)
weights2<-c(0.153,0.140,0.698,0,0)
weights3<-c(0,0,1,0,0)

#Without perturbations
weights1<-c(0,1,0,0,0)
weights2<-c(0,0.21,0.79,0,0)
weights3<-c(0,0,1,0,0)


# Compute the forecast death rates from the model averaging approaches

pred1<-weights1[1]*params[[1]]$mufor
pred2<-weights2[1]*params[[1]]$mufor
pred3<-weights3[1]*params[[1]]$mufor

for (i in 2:5){
  pred1<-pred1+weights1[i]*params[[i]]$mufor
  pred2<-pred2+weights2[i]*params[[i]]$mufor
  pred3<-pred3+weights3[i]*params[[i]]$mufor
}


# Resize the forecast deaths rates as an array "Number of draws X Ages X Years to predict"
samplesize<-4000
years.fit<-1979:2018
years.predict<-2019:2028
pred1<-array(pred1,dim=list(samplesize,length(ages.fit),length(years.predict)),
             dimnames = list(c(1:samplesize),formatC(ages.fit),formatC(years.predict)))
pred2<-array(pred2,dim=list(samplesize,length(ages.fit),length(years.predict)),
             dimnames = list(c(1:samplesize),formatC(ages.fit),formatC(years.predict)))
pred3<-array(pred3,dim=list(samplesize,length(ages.fit),length(years.predict)),
             dimnames = list(c(1:samplesize),formatC(ages.fit),formatC(years.predict)))


# Plot which represents period survival probability at age 50 truncated age 90.

qxt <- deathFR / exposureFR
qxt<-exp(-qxt)
qxt<-apply(qxt,2,cumprod)

p5090<-colSums(qxt)

#Predictions of the period survival probability for each method
prob1<-exp(-pred1)
prob1<-apply(prob1,c(1,3),cumprod)
prob1<-t(colSums(prob1, dims=1))

prob2<-exp(-pred2)
prob2<-apply(prob2,c(1,3),cumprod)
prob2<-t(colSums(prob2, dims=1))

prob3<-exp(-pred3)
prob3<-apply(prob3,c(1,3),cumprod)
prob3<-t(colSums(prob3, dims=1))


matplot(c(2009:2018), p5090[31:40],xlim = c(2009, 2028), ylim = range(25,35), pch = 20, col = "black",
        xlab = "year", ylab=" ",log = "y",main="Life expectancy at age 50")


#Compute quantiles and median for the first method (BMA)

mxtpred2.5 <- apply(prob1, 1, quantile, probs = 0.025)
mxtpred97.5 <- apply(prob1, 1, quantile, probs = 0.975)
matlines(years.predict, mxtpred2.5, lty = 1, col = "black")
matlines(years.predict, mxtpred97.5, lty = 1, col = "black")

#Compute quantiles and median for the second method (stacking)
mxtpred2.5 <- apply(prob2, 1, quantile, probs = 0.025)
mxtpred97.5 <- apply(prob2, 1, quantile, probs = 0.975)
matlines(years.predict, mxtpred2.5, lty = 2, col = "green")
matlines(years.predict, mxtpred97.5, lty = 2, col = "green")

#Compute quantiles and median for the third method (pseudo-BMA)
mxtpred2.5 <- apply(prob3, 1, quantile, probs = 0.025)
mxtpred97.5 <- apply(prob3, 1, quantile, probs = 0.975)
matlines(years.predict, mxtpred2.5, lty = 3, col = "red")
matlines(years.predict, mxtpred97.5, lty = 3, col = "red")

#In our paper, we compared the three model averaging approaches under
#perturbed data with the LC model without perturbations. 

#The LC without perturbations was generated with the StMoMo pacakge (frequentist LC model)

#Comparison with the frequentist Poisson Lee-Carter model
library(StMoMo)
LC <- lc(link = "log")
ages.fit <- 50:90
LCfit <- fit(LC, Dxt=deathFR,Ext=exposureFR)
nsim<-10000
LCsim <- simulate(LCfit, nsim = nsim, h = 10)
qxtLC<-aperm(LCsim$rates,c(3,1,2))

prob4<- exp(-qxtLC) 
prob4<-apply(prob4,c(1,3),cumprod)
prob4<-t(colSums(prob4, dims=1))

#Compute quantiles and median for the Lee-Carter Poisson model
mxtpred2.5 <- apply(prob4, 1, quantile, probs = 0.025)
mxtpred97.5 <- apply(prob4, 1, quantile, probs = 0.975)
matlines(years.predict, mxtpred2.5, lty = 5, col = "blue")
matlines(years.predict, mxtpred97.5, lty = 5, col = "blue")

legend("topleft", c("BMA","Stacking","Pseudo-BMA","LC (no perturbations)"),
       col=c("black","green","red","blue"),lty=c(1,2,3,5),cex = 0.8,bty = "n")

# legend("topleft", c("BMA","Stacking","Pseudo-BMA"),
#        col=c("black","green","red"),lty=c(1,2,3),cex = 0.8,bty = "n")

#The median life expectancy described in the paper can be computed with 
#the following lines.

medianBMA<- apply(prob1, 1, quantile, probs = 0.5)
medianstacking<- apply(prob2, 1, quantile, probs = 0.5)
medianpseudo<- apply(prob3, 1, quantile, probs = 0.5)
medtable<-cbind(medianBMA,medianstacking,medianpseudo)

# Hereafter, we show how to compare the cohort effect in RH and APC.
# In particular, the cohort effect for APC is much more stable.

ages<-ages.fit
params<-lapply(out2,rstan::extract)

#Extract gamma from RH model
params1<-params[[2]]
gammarh<-params1$gf

#Extract gamma from APC model
params2<-params[[3]]
gammaapc<-params2$gf

library(ggfan)
library(abind)
library(dplyr)
library(ggplot2)
library(latex2exp)

cohort.fit<-(years[1]-ages[length(ages)]):(years[length(years)]-ages[1])
cohort.predict<-(cohort.fit[length(cohort.fit)]+1):(cohort.fit[length(cohort.fit)]+length(years.predict))
cohort<-c(cohort.fit,cohort.predict)

# We resize gamma appropriately in a dataframe for the use of ggplot
gamma<-abind(gammarh,gammaapc,rev.along = 0) 
dimnames(gamma)<-list(sim=formatC(c(1:1000)),cohort=formatC(cohort),model=c("rh","apc"))
gamma<-as.data.frame.table(gamma,responseName = "Obs")
gamma$cohort <-as.numeric(as.character(gamma$cohort))
levels(gamma$model) <- c("Renshaw-Haberman (RH)", "Age-Period-Cohort (APC)")

#The following ggplot compares the cohort effect of RH and APC
gammap <- ggplot(gamma, aes(x=cohort, y=Obs)) + 
  geom_interval(intervals=c(0.95),show.legend = FALSE) + 
  theme_classic()+
  facet_wrap(~model)+
  labs(title="Cohort effect",x="Year",y=TeX("$\\gamma_t$"))+
  theme(plot.title = element_text(hjust = 0.5))+
  geom_fan()

print(gammap)




