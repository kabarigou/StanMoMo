library(tidyr)
library(dplyr)
library(tibble)
library(matrixcalc)
library(StanMoMo)
library(tidyverse)
library(ggplot2)
library(magrittr)
library(ggfan)
library(abind)
library(ggpubr)

#This script details how to reproduce the results of our paper from
#Section 5: Application to real mortality data

#In a first step we compute the weights from stacking, pseudo-BMA and BMA for French data.
#The same procedure can be repeated for UK, USA and Japan given the country codes of HMD

country<- 'FRATNP'
#code UK 'GBR_NP', code USA 'USA', code Japan 'JPN'
ages.fit<-50:90
years.fit<-1979:2008

deathFR<-load_HMD_data(country, 'Deaths_1x1', years.fit, ages.fit, "Male")$mat
exposureFR<-load_HMD_data(country, 'Exposures_1x1', years.fit,ages.fit, "Male")$mat

#We run the 5 mortality models in parallel for computational efficiency
#For this, we group the mortality functions into a wrapper function and call the parLapply function
#First, we compute the weights via stacking and pseudo-BMA with training and validation sets and then BMA with the whole dataset

library(parallel)

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

#The following function gives weights according to stacking and pseudo-BMA
model_weights<-mortality_weights(out)
print(model_weights)

#Now, we fit the models and determine the weights according to the standard BMA 
# approach. This is performed via the use of the R package bridgesampling

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

# The procedure above can be repeated for each country and it delivers the model weights for BMA, stacking and Pseudo-BMA as stated in the paper.


# Here below, we show how to generate the different figures of death rates, 
#survival probabilities and Mean Absolute Error given the model weights

# All posterior distributions of all parameters can be extracted with the 'extract' function of the 'rstan' package. 
out2<-lapply(out, `[[`, 1)   
params<-lapply(out2,rstan::extract)

#Weights obtained for the different countries (order: BMA, stacking, pseudo-BMA)
countryname<-"France"
#Again, we show the case of France. Other plots are generated similarly 
#changing the country name

# countryname<-"United Kingdom"
# country<- 'GBR_NP'
# countryname<-"United States"
# country<-"USA"
# countryname<-"Japan"
# country<-"JPN"

#France
weights1<-c(0,1,0,0,0)
weights2<-c(0.093,0.750,0.157,0,0)
weights3<-c(0,1,0,0,0)

#UK
weights1<-c(0,0,0,0,1)
weights2<-c(0,0.342,0,0,0.658)
weights3<-c(0,0,0,0,1)

#USA
weights1<-c(0,1,0,0,0)
weights2<-c(0,0.548,0.452,0,0)
weights3<-c(0,0.982,0.018,0,0)

#Japan
weights1<-c(0,1,0,0,0)
weights2<-c(0.292,0.239,0.468,0,0)
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
years.predict<-2009:2018

pred1<-array(pred1,dim=list(samplesize,length(ages.fit),length(years.predict)),
             dimnames = list(sim=c(1:samplesize),age=formatC(ages.fit),year=formatC(years.predict)))

pred2<-array(pred2,dim=list(samplesize,length(ages.fit),length(years.predict)),
             dimnames = list(c(1:samplesize),formatC(ages.fit),formatC(years.predict)))
pred3<-array(pred3,dim=list(samplesize,length(ages.fit),length(years.predict)),
             dimnames = list(c(1:samplesize),formatC(ages.fit),formatC(years.predict)))

# Forecast death rates at age 65,75,85. 
pred<-abind(pred1,pred2,pred3,rev.along = 0)
pred1d<-as.data.frame.table(pred)
colnames(pred1d)<-c("Sim","Age","Year","Method","Obs")
pred1d$Age <-as.numeric(as.character(pred1d$Age))
pred1d$Year <-as.numeric(as.character(pred1d$Year))


ages.fit<-50:90
years.fit<-1979:2018
deathFR<-load_HMD_data(country, 'Deaths_1x1', years.fit, ages.fit, "Male")$mat
exposureFR<-load_HMD_data(country, 'Exposures_1x1', years.fit,ages.fit, "Male")$mat

qxt <- deathFR / exposureFR
qxt<-as.data.frame.table(qxt)
colnames(qxt)<-c("Age","Year","Obs")
qxt$Age <-as.numeric(as.character(qxt$Age))
qxt$Year <-as.numeric(as.character(qxt$Year))


#The following ggplot represents death rates of ages 65,75 and 85 
#with their 10-year forecast intervals

p <-ggplot()+
  xlim(c(1999, 2018)) + 
  geom_point(data = subset(qxt, Age %in% c(65,75,85) & Year %in% c(1999:2008)),aes(x=Year,y=Obs,shape=as.factor(Age)))+
  geom_point(data = subset(qxt, Age %in% c(65,75,85) & Year %in% c(2009:2018)),aes(x=Year,y=Obs),
             shape=4,col="steelblue")+
  scale_y_log10(limits = c(0.005, 0.2))+
  geom_interval(data=subset(pred1d, Age %in% c(65,75,85) & Year %in% c(2009:2018) & Method %in% c("A")),
                intervals=c(0.95),aes(x=Year,y=Obs,linetype=Method,color=Method,group=Age))+
  theme_bw()+ 
  geom_interval(data=subset(pred1d, Age %in% c(65,75,85) & Year %in% c(2009:2018) & Method %in% c("B")),
                intervals=c(0.95),aes(x=Year,y=Obs,linetype=Method, color=Method,group=Age))+
  geom_interval(data=subset(pred1d, Age %in% c(65,75,85) & Year %in% c(2009:2018) & Method %in% c("C")),
                intervals=c(0.95),aes(x=Year,y=Obs,linetype=Method, color=Method,group=Age))+
  labs(title=countryname, y="death rate (log scale)", x="Year",shape="Age")+
  scale_linetype_manual(name="Method",labels=c("BMA","Stacking","Pseudo-BMA"),values= c('solid', 'longdash', 'dotted'))+
  scale_colour_manual(name="Method", labels=c("BMA","Stacking","Pseudo-BMA"),values = c("black", "green","red"))+
  theme(plot.title = element_text(hjust = 0.5))

print(p)

#The following two lines are used to merge 4 figures into one.
#Assuming p1,p2,p3,p4 are the plots of the different countries

# ggarrange(p1, p2, p3, p4, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")
# ggsave("deathrates.pdf", width = 297, height = 210, units = "mm")


# Plot with the survival probability at age 50 until 90
qxt <- deathFR / exposureFR
prob<-cbind(c(1979:2018),exp(-colSums(qxt)))

prob<-as.data.frame(prob)
colnames(prob)<-c("Year","Obs")

#Predictions of the period survival probability for each method
prob1<-apply(pred1,1,colSums)
prob1<-exp(-prob1)

prob2<-apply(pred2,1,colSums)
prob2<-exp(-prob2)

prob3<-apply(pred3,1,colSums)
prob3<-exp(-prob3)

# The following ggplot represents the survival probability at age 50 truncated 
# at age 90 as represented in the paper
probpred<-abind(prob1,prob2,prob3,rev.along = 0) %>% as.data.frame.table()
colnames(probpred)<-c("Year","Sim","Method","Obs")
probpred$Year <-as.numeric(as.character(probpred$Year))

q <- ggplot()+
  geom_point(data = subset(prob, Year %in% c(1999:2008)),aes(x=Year,y=Obs,shape=19))+
  geom_point(data = subset(prob, Year %in% c(2009:2018)),aes(x=Year,y=Obs),
             shape=4,size=2,col="steelblue")+
  scale_shape_identity()+
  theme_bw()+ 
  geom_interval(data=subset(probpred, Year %in% c(2009:2018)),
                intervals=c(0.95),aes(x=Year,y=Obs,linetype=Method,color=Method,group=Method))+
  labs(title=countryname, y="survival probability (log scale)", x="Year",shape="Age")+
  scale_linetype_manual(name="Method",labels=c("BMA","Stacking","Pseudo-BMA"),values= c('solid', 'longdash', 'dotted'))+
  scale_colour_manual(name="Method", labels=c("BMA","Stacking","Pseudo-BMA"),values = c("black", "green","red"))+
  theme(plot.title = element_text(hjust = 0.5))+
  scale_y_log10()

print(q)

#The following two lines are used to merge 4 figures into one.
#Assuming p1,p2,p3,p4 are the plots of the different countries

# ggarrange(q1, q2, q3, q4, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")
# ggsave("survivalprobab.pdf", width = 297, height = 210, units = "mm")

# Mean Death Rates for each BMA approach
pred1_mean<-apply(pred1,c(2,3),mean)
pred2_mean<-apply(pred2,c(2,3),mean)
pred3_mean<-apply(pred3,c(2,3),mean)

# Determine the difference between observed deaths and forecast deaths
deathFR_2<-deathFR[,formatC(c(2009:2018))]
exposureFR_2<-exposureFR[,formatC(c(2009:2018))]

error1<-deathFR_2-exposureFR_2*pred1_mean
error2<-deathFR_2-exposureFR_2*pred2_mean
error3<-deathFR_2-exposureFR_2*pred3_mean

#Mean Absolute error per age for each BMA approach
mf1<-rowMeans(abs(error1))
mf2<-rowMeans(abs(error2))
mf3<-rowMeans(abs(error3))

#The following ggplot represents MAE per Age for each method as represented in the paper

mae<-data.frame(x=50:90,values=c(mf1,mf2,mf3),method=rep(c("1","2","3"),each=41))
r<-ggplot(mae,aes(x, values, color = method,linetype=method)) +
  geom_line()+
  labs(title=countryname, y="Mean Absolute Error", x="Age")+
  scale_linetype_manual(name="Method",labels=c("BMA","Stacking","Pseudo-BMA"),values= c('solid', 'longdash', 'dotted'))+
  scale_colour_manual(name="Method", labels=c("BMA","Stacking","Pseudo-BMA"),values = c("black", "green","red"))+
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5))
  
print(r)

#The following two lines are used to merge 4 figures into one.
#Assuming p1,p2,p3,p4 are the plots of the different countries
# ggarrange(r1, r2, r3, r4, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")
# ggsave("mae.pdf", width = 297, height = 210, units = "mm")
