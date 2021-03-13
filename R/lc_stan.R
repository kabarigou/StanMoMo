# Save this file as `R/lc_stan.R`

#' Bayesian Lee-Carter with Stan
#'
#' Fit and Forecast Bayesian Lee-Carter model. The model can be fitted with a Poisson or
#' Negative-Binomial distribution. The function outputs posteriors distributions for each parameter,
#' predicted death rates and log-likehoods.
#'
#' The created model is either a log-Poisson or a log-Negative-Binomial version of
#' the Lee-Carter model:
#' \deqn{D_{x,t} \sim \mathcal{P}(\mu_{x,t} e_{x,t})}
#' or
#' \deqn{D_{x,t}\sim NB\left(\mu_{x,t} e_{x,t},\phi\right)}
#' with
#' \deqn{\log \mu_{xt} = \alpha_x + \beta_x\kappa_t.}
#'
#' To ensure the identifiability of th model, we impose
#' \deqn{\sum_x\beta_x = 1,\kappa_1=0.}
#'
#' For the priors, the model chooses relatively wide priors:
#' \deqn{\alpha_x \sim N(0,100),\beta_{x} \sim Dir(1,\dots,1),\frac{1}{\phi} \sim Half-N(0,1).}
#'
#' For the period term, we consider a first order autoregressive process (AR(1)) with linear trend:
#' \deqn{\kappa_{t}=c+\kappa_{t-1}+\epsilon_{t},\epsilon_{t}\sim N(0,\sigma^2)}
#' with \eqn{c \sim N(0,10),\sigma \sim Exp(0.1)}. 
#'
#' @export
#' @param death Matrix of deaths.
#' @param exposure Matrix of exposures.
#' @param forecast Number of years to forecast.
#' @param validation Number of years for validation.
#' @param family specifies the random component of the mortality model. \code{"Poisson"} assumes a
#' Poisson model with log link and \code{"nb"} assumes a negative-binomial model
#' with log link and overdispersion parameter \eqn{\phi}.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`.
#'
#'@references
#'Lee, R. D., & Carter, L. R. (1992). Modeling and forecasting U.S. mortality.
#' Journal of the American Statistical Association, 87(419), 659-671.
#'
#' @examples
#'
#' \dontrun{
#' #10-year forecasts for French data for ages 50-90 and years 1970-2017 with a log-Poisson model
#' ages.fit<-50:90
#' years.fit<-1970:2017
#' deathFR<-FRMaleData$Dxt[formatC(ages.fit),formatC(years.fit)]
#' exposureFR<-FRMaleData$Ext[formatC(ages.fit),formatC(years.fit)]
#' fitLC=lc_stan(death = deathFR,exposure=exposureFR, forecast = 10, family = "poisson")
#' }
#'
lc_stan <- function(death, exposure,forecast, validation=0, family=c("poisson","nb"), ...) {
  Tval<-0
  if (validation !=0) Tval=validation

  if (Tval==0) {
  death1<-death
  exposure1<-exposure
  death2<-vector('integer')
  exposure2<-vector('integer')
  } else {
    T<- ncol(death)-Tval
    death1<-death[,1:T]
    death2<-death[,(T+1):ncol(death)]

    exposure1<-exposure[,1:T]
    exposure2<-exposure[,(T+1):ncol(exposure)]
  }

  family<-match.arg(family)
  if (family == "poisson") {
    family <- 0
  } else if (family == "nb") {
    family <- 1
  }

  standata<- list(J=nrow(death1),T=ncol(death1),
       d=as.integer(as.vector(death1)),
       e=as.integer(as.vector(exposure1)),
       dval=as.integer(as.vector(death2)),
       eval=as.integer(as.vector(exposure2)),
       Tfor=forecast,Tval=Tval,
       family=family)
  suppressWarnings( {
    out <- rstan::sampling(stanmodels$leecarter, data = standata, ...)
      } )
  return(out)
}
