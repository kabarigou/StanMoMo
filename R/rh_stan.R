# Save this file as `R/rh_stan.R`

#' Bayesian Renshaw-Haberman model with Stan
#'
#' Fit and Forecast Bayesian Renshaw-Haberman model (Lee-Carter with cohort effect) introduced in Renshaw and Haberman (2006).
#' The model can be fitted with a Poisson or Negative-Binomial distribution. The function outputs posteriors distributions for each parameter,
#' predicted death rates and log-likelihoods.
#'
#' The created model is either a log-Poisson or a log-Negative-Binomial version of
#' the Renshaw-Haberman model:
#' \deqn{D_{x,t} \sim \mathcal{P}(\mu_{x,t} e_{x,t})}
#' or
#' \deqn{D_{x,t}\sim NB\left(\mu_{x,t} e_{x,t},\phi\right)}
#' with
#' \deqn{\log \mu_{xt} = \alpha_x + \beta_x\kappa_t+\gamma_{t-x}.}
#'
#' To ensure the identifiability of th model, we impose
#' \deqn{\kappa_1=0, \gamma_1=0,\sum gamma_i =0, \gamma_C=0,}
#' where \eqn{C} represents the most recent cohort in the data.
#'
#' For the priors, the model chooses wide priors:
#' \deqn{\alpha_x \sim N(0,100),\beta_{x} \sim Dir(1,\dots,1),\frac{1}{\phi} \sim Half-N(0,1).}
#'
#' For the period term, we consider the standard random walk with drift:
#' \deqn{\kappa_{t}=c+\kappa_{t-1}+\epsilon_{t},\epsilon_{t}\sim N(0,\sigma^2)}
#' with \eqn{c \sim N(0,10),\sigma \sim Exp(0.1)}.
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
#' @param validation Number of years for validation.
#' @param family specifies the random component of the mortality model. \code{"Poisson"} assumes a
#' Poisson model with log link and \code{"nb"} assumes a negative-binomial model
#' with log link and overdispersion parameter \eqn{\phi}.
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`.
#'
#'@references
#' Renshaw, A. E., & Haberman, S. (2006). A cohort-based extension to the
#' Lee-Carter model for mortality reduction factors.
#' Insurance: Mathematics and Economics, 38(3), 556-570.
#'
#' @examples
#'
#'
#' #10-year forecasts for French data for ages 50-90 and years 1970-2017 with a log-Poisson model
#' ages.fit<-70:90
#' years.fit<-1990:2010
#' deathFR<-FRMaleData$Dxt[formatC(ages.fit),formatC(years.fit)]
#' exposureFR<-FRMaleData$Ext[formatC(ages.fit),formatC(years.fit)]
#' iterations<-50 # Toy example, consider at least 2000 iterations
#' fitRH=rh_stan(death = deathFR,exposure=exposureFR, forecast = 5, family = "poisson",
#' iter=iterations,chains=1)
#'
#'
rh_stan <- function(death, exposure,forecast, validation=0, family=c("poisson","nb"), ...) {
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
    out <- rstan::sampling(stanmodels$RHmodel, data = standata, ...)
      } )
  return(out)
}
