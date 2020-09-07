StanMoMo
================

The StanMoMo package performs Bayesian Mortality Modeling with Stan for a wide variety of popular mortality models.The StanMoMo package performs Bayesian Mortality Modeling with Stan for a variety of popular mortality models. The current package supports the Lee-Carter (LC) model, the Renshaw-Haberman model (LC with cohort effect), the Age-Period-Cohort (APC) model, the Cairns-Blake-Dowd (CBD) model and the M6 model (CBD with cohort effect). By a simple call, the user inputs deaths and exposures and the package outputs the MCMC simulations for each parameter, the log likelihoods and predictions. Moreover, the package includes tools for model selection and Bayesian model averaging by future-out cross-validation.

Installation
------------

To install `StanMoMo` from GitHub you will need `devtools`:

``` r
install.packages("devtools")
devtools::install_github('kabarigou/StanMoMo')
```

After installing the package, you have to load the package via:

``` r
library(StanMoMo)
```

Introduction: Bayesian Lee-Carter with Stan
===========================================

We first explain how to estimate the Lee-Carter model using the `lc_stan` function in the `StanMoMo` package.

For illustration, the package already includes the object `FRMaleData` containing deaths (`FRMaleData$Dxt`) and central exposures (`FRMaleData$Ext`) for French males for the period 1816-2017 and for ages 0-100. In our example, we concentrated on ages 50-90 and the period 1970-2017. This can be obtained via:

``` r
ages.fit<-50:90
years.fit<-1970:2017
deathFR<-FRMaleData$Dxt[formatC(ages.fit),formatC(years.fit)]
exposureFR<-FRMaleData$Ext[formatC(ages.fit),formatC(years.fit)]
```

As a reminder, the Lee-Carter model assumes that mortality dynamics are given by

![ 
\\begin{aligned}
D\_{xt} &\\sim Poisson (E\_{xt}\\mu\_{xt})\\\\
\\log \\mu\_{xt}&=\\alpha\_x+\\beta\_x \\kappa\_t
\\end{aligned}
](https://latex.codecogs.com/png.latex?%20%0A%5Cbegin%7Baligned%7D%0AD_%7Bxt%7D%20%26%5Csim%20Poisson%20%28E_%7Bxt%7D%5Cmu_%7Bxt%7D%29%5C%5C%0A%5Clog%20%5Cmu_%7Bxt%7D%26%3D%5Calpha_x%2B%5Cbeta_x%20%5Ckappa_t%0A%5Cend%7Baligned%7D%0A " 
\begin{aligned}
D_{xt} &\sim Poisson (E_{xt}\mu_{xt})\\
\log \mu_{xt}&=\alpha_x+\beta_x \kappa_t
\end{aligned}
")

 To ensure identifiability of the model, we assume that

![
\\sum\_x \\beta\_x=1,\\kappa\_1=0
](https://latex.codecogs.com/png.latex?%0A%5Csum_x%20%5Cbeta_x%3D1%2C%5Ckappa_1%3D0%0A "
\sum_x \beta_x=1,\kappa_1=0
")

Moreover, we assume that the period parameter follows a first order autoregressive process (AR(1)) with linear trend:

![ 
\\kappa\_t \\sim \\mathcal{N}(c+\\delta t+\\rho \\kappa\_{t-1},\\sigma)
](https://latex.codecogs.com/png.latex?%20%0A%5Ckappa_t%20%5Csim%20%5Cmathcal%7BN%7D%28c%2B%5Cdelta%20t%2B%5Crho%20%5Ckappa_%7Bt-1%7D%2C%5Csigma%29%0A " 
\kappa_t \sim \mathcal{N}(c+\delta t+\rho \kappa_{t-1},\sigma)
")

 We note that the random walk with drift appears as a special case. The choice of priors can be found in the paper.

Bayesian Estimation
-------------------

Now, the model can be estimated and forecast for the next 10 years by a simple call to the `lc_stan` function:

``` r
fitLC=lc_stan(death = deathFR,exposure=exposureFR, forecast = 10, family = "poisson",cores=4)
```

By default, Stan samples four Markov chains with 2000 iterations with 1000 warmup iterations for each chain (hence, 4000 draws in total). The output is an object of class `stanfit` (cf. the **rstan** package) which contains the posterior draws of each parameter, the posterior log likelihoods as well as the mortality forecasts.

The user can have access to an interface for interactive MCMC diagnostics and plots and tables helpful for analyzing posterior samples through the `shinystan` package.

``` r
library(shinystan)
launch_shinystan(fitLC)
```

We can also plot the sampling for each parameter with the use of the `bayesplot` package. For instance, the posterior uncertainty intervals for ![\\alpha\_x](https://latex.codecogs.com/png.latex?%5Calpha_x "\alpha_x"), for ![x=50, 51, 52, 53](https://latex.codecogs.com/png.latex?x%3D50%2C%2051%2C%2052%2C%2053 "x=50, 51, 52, 53"), based on the 4000 draws can be obtained via

``` r
library(bayesplot)
```

    ## Warning: package 'bayesplot' was built under R version 3.6.3

    ## This is bayesplot version 1.7.1

    ## - Online documentation and vignettes at mc-stan.org/bayesplot

    ## - bayesplot theme set to bayesplot::theme_default()

    ##    * Does _not_ affect other ggplot2 plots

    ##    * See ?bayesplot_theme_set for details on theme setting

``` r
posterior <- as.array(fitLC)
mcmc_areas(posterior,pars = c("a[1]", "a[2]", "a[3]", "a[4]"),prob = 0.8, prob_outer = 0.99, point_est = "mean")
```

![](readme_files/figure-markdown_github/unnamed-chunk-5-1.png)

We can also use the output to produce fan charts depicting the uncertainty around each parameter of the model fit.

``` r
library("fanplot")
library("RColorBrewer")
library(latex2exp)
params<-rstan::extract(fitLC)
#Alpha
plot(ages.fit, colMeans(params$a), ylim=range(params$a),ylab=TeX("$\\alpha_x$"), xlab="Age: x")
fan(data=params$a, start=ages.fit[1],type = "interval", ln=NULL,probs = seq(0.01,0.99,0.01),
    fan.col = colorRampPalette(colors = rev(brewer.pal(9,"Oranges"))))
```

![](readme_files/figure-markdown_github/unnamed-chunk-6-1.png)

``` r
#Beta
plot(ages.fit, colMeans(params$b), ylim=range(params$b),ylab=TeX("$\\beta_x$"), xlab="Age: x")
fan(data=params$b, start=ages.fit[1],type = "interval", ln = NULL,probs = seq(0.01,0.99,0.01),
    fan.col = colorRampPalette(colors = rev(brewer.pal(9,"Reds"))))
```

![](readme_files/figure-markdown_github/unnamed-chunk-6-2.png)

``` r
#Kappa
plot(years.fit, colMeans(params$k), ylim=range(params$k),ylab=TeX("$\\kappa_t$"), xlab="Year: t")
fan(data=params$k, start=years.fit[1],type = "percentile",ln = NULL, probs = seq(0.01,0.99,0.01),
    fan.col = colorRampPalette(colors = rev(brewer.pal(9,"Blues"))))
```

![](readme_files/figure-markdown_github/unnamed-chunk-6-3.png)

Forecasting
-----------

Predictions of death rates for the next 10 years with confidence intervals can be obtained as follows:

``` r
# Resize the forecast deaths rates as an array "Number of draws X Ages X Years to predict"
samplesize<-4000
years.predict<-2018:2027
pred<-array(params$mufor,dim=list(samplesize,length(ages.fit),length(years.predict)),
            dimnames = list(c(1:samplesize),formatC(ages.fit),formatC(years.predict)))
#Fan plots for ages 65,75,85
probs = c(2.5, 10, 25, 50, 75, 90, 97.5)
qxt <- deathFR / exposureFR
matplot(years.fit, t(qxt[c("65", "75", "85"), ]),
        xlim = c(1970, 2027), ylim = c(0.0025, 0.2), pch = 20, col = "black",
        log = "y", xlab = "year", ylab = "death rate (log scale)")
fan(pred[,"65" , ], start = 2018, probs = probs, n.fan = 4,
    fan.col =  colorRampPalette(colors = rev(brewer.pal(9,"Reds"))), ln = NULL)
fan(pred[,"75" , ], start = 2018, probs = probs, n.fan = 4,
    fan.col =  colorRampPalette(colors = rev(brewer.pal(9,"Greens"))), ln = NULL)
fan(pred[,"85" , ], start = 2018, probs = probs, n.fan = 4,
    fan.col =  colorRampPalette(colors = rev(brewer.pal(9,"Blues"))), ln = NULL)
text(1980, qxt[c("65", "75", "85"), "2000"],
     labels = c("x = 65", "x = 75", "x = 85"))
```

![](readme_files/figure-markdown_github/unnamed-chunk-7-1.png)

Other models
============

In a similar fashion, the models RH, APC, CBD and M6 can be fitted with the following calls:

``` r
fitRH=rh_stan(death = deathFR,exposure=exposureFR, forecast = 10, family = "poisson",cores=4)
fitAPC=apc_stan(death = deathFR,exposure=exposureFR, forecast = 10, family = "poisson",cores=4)
fitCBD=cbd_stan(death = deathFR,exposure=exposureFR, age=ages.fit, forecast=10,family = "poisson",cores=4)
fitM6=m6_stan(death = deathFR,exposure=exposureFR, age=ages.fit,forecast = 10, family = "poisson",cores=4)
```

Bayesian model averaging
========================

The package also includes the `mortality_weights` function which performs Bayesian model averaging via future-out cross-validation. The weights are obtained by stacking of predictive distributions and Pseudo-BMA weighting, based on Yao et al. (2018) but adapted for mortality forecasting in our paper. Roughly speaking, the last ![K](https://latex.codecogs.com/png.latex?K "K") years of data are left out and the weights search for the optimal mortality model on these last ![K](https://latex.codecogs.com/png.latex?K "K") years. For ![K=10](https://latex.codecogs.com/png.latex?K%3D10 "K=10") and 20 years of prediction, the weights can be obtained via

``` r
fitLC=lc_stan(death = deathFR,exposure=exposureFR, forecast = 20, validation=10,family = "poisson",cores=4)
fitRH=rh_stan(death = deathFR,exposure=exposureFR, forecast = 20, validation=10,family = "poisson",cores=4)
fitAPC=apc_stan(death = deathFR,exposure=exposureFR, forecast = 20, validation=10,family = "poisson",cores=4)
fitCBD=cbd_stan(death = deathFR,exposure=exposureFR, age=ages.fit, forecast=10,validation=10,family = "poisson",cores=4)
fitM6=m6_stan(death = deathFR,exposure=exposureFR, age=ages.fit,forecast = 20, validation=10,family = "poisson",cores=4)
model_weights<-mortality_weights(list(fitLC,fitRH,fitAPC,fitCBD,fitM6))
```

Fit with a negative-binomial
----------------------------

All the models can be fitted with a Negative-Binomial distribution with the argument `family="nb"`.
