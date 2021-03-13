# Save this file as `R/rh_stan.R`

#' Bayesian M6 model with Stan
#'
#' Fit and Forecast Bayesian M6 model (CBD with cohort effect) introduced in Cairns et al (2009).
#' The model can be fitted with a Poisson or Negative-Binomial distribution. The function outputs posteriors distributions for each parameter,
#' predicted death rates and log-likehoods.
#'
#' The created model is either a log-Poisson or a log-Negative-Binomial version of
#' the M6 model:
#' \deqn{D_{x,t} \sim \mathcal{P}(\mu_{x,t} e_{x,t})}
#' or
#' \deqn{D_{x,t}\sim NB\left(\mu_{x,t} e_{x,t},\phi\right)}
#' with
#' \deqn{\log \mu_{xt} = \kappa_t^{(1)} + (x-\bar{x})\kappa_t^{(2)}+\gamma_{t-x},}
#' where \eqn{\bar{x}} is the average age in the data.
#'
#' To ensure the identifiability of th model, we impose
#' \deqn{\gamma_1=0,\gamma_C=0,}
#' where \eqn{C} represents the most recent cohort in the data.
#'
#' For the period terms, we consider a multivariate random walk with drift:
#' \deqn{\boldsymbol{\kappa}_{t}=\boldsymbol{c}+ \boldsymbol{\kappa}_{t-1}+\boldsymbol{\epsilon}_{t}^{\kappa},\quad \boldsymbol{\kappa}_{t}=\left(\begin{array}{c}\kappa_{t}^{(1)} \\\kappa_{t}^{(2)}\end{array}\right), \quad \boldsymbol{\epsilon}_{t}^{\kappa} \sim N\left(\mathbf{0}, \Sigma\right),}
#' with normal priors: \eqn{\boldsymbol{c} \sim N(0,10)}.
#'
#' The variance-covariance matrix of the error term is defined by
#' \deqn{\boldsymbol{\Sigma}=\left(\begin{array}{cc}\sigma_1^{2} & \rho_{\Sigma} \sigma_1 \sigma_2 \\\rho_{\Sigma} \sigma_1 \sigma_{Y} & \sigma_2^{2}\end{array}\right)}
#' where the variance coefficients have independent exponential priors: \eqn{\sigma_1, \sigma_2 \sim Exp(0.1)}
#' and the correlation parameter has a uniform prior: \eqn{\rho_{\Sigma} \sim U\left[-1,1\right]}.
#' As for the other models, the overdisperion parameter has a prior distribution given by
#' \deqn{\frac{1}{\phi} \sim Half-N(0,1).}
#'
#' For the cohort term, we consider a second order autoregressive process (AR(2)):
#' \deqn{\gamma_{c}=\psi_1 \gamma_{c-1}+\psi_2 \gamma_{c-2}+\epsilon^{\gamma}_{t},\quad \epsilon^{\gamma}_{t}\sim N(0,\sigma_{\gamma}).}
#'
#' To close the model specification, we impose some vague priors assumptions on the hyperparameters:
#' \deqn{\psi_1,\psi_2 \sim N(0,10),\quad \sigma_{\gamma}\sim Exp(0.1).}
#'
#' @export
#' @param death Matrix of deaths.
#' @param exposure Matrix of exposures.
#' @param forecast Number of years to forecast.
#' @param age Vector of ages.
#' @param validation Number of years for validation.
#' @param family specifies the random component of the mortality model. \code{"Poisson"} assumes a
#' Poisson model with log link and \code{"nb"} assumes a negative-binomial model
#' with log link and overdispersion parameter \eqn{\phi}.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`.
#'
#'
#' @references
#'
#' Cairns, A. J. G., Blake, D., Dowd, K., Coughlan, G. D., Epstein, D.,
#' Ong, A., & Balevich, I. (2009). A quantitative comparison of stochastic
#' mortality models using data from England and Wales and the United States.
#' North American Actuarial Journal, 13(1), 1-35.
#'
#'
#' @examples
#'
#' \dontrun{
#' #10-year forecasts for French data for ages 50-90 and years 1970-2017 with a log-Poisson model
#' ages.fit<-50:90
#' years.fit<-1970:2017
#' deathFR<-FRMaleData$Dxt[formatC(ages.fit),formatC(years.fit)]
#' exposureFR<-FRMaleData$Ext[formatC(ages.fit),formatC(years.fit)]
#' fitM6=m6_stan(death = deathFR,exposure=exposureFR, age=ages.fit,forecast = 10, family = "poisson")
#' }
#'
#'
m6_stan <- function(death, exposure,forecast, age, validation=0, family=c("poisson","nb"), ...) {
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
       age=age,
       dval=as.integer(as.vector(death2)),
       eval=as.integer(as.vector(exposure2)),
       Tfor=forecast,Tval=Tval,
       family=family)
  suppressWarnings( {
    out <- rstan::sampling(stanmodels$M6model, data = standata, ...)
      } )
  return(out)
}
