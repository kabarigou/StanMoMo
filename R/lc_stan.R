# Save this file as `R/lc_stan.R`

#' Bayesian Lee-Carter with Stan
#'
#' @export
#' @param death Matrix of deaths
#' @param exposure Matrix of exposures
#' @param forecast Number of years to forecast
#' @param ... Arguments passed to `rstan::sampling` (e.g. iter, chains).
#' @return An object of class `stanfit` returned by `rstan::sampling`
#'
lc_stan <- function(death, exposure,forecast, ...) {
  standata <- list(J=nrow(death),T=ncol(death),
                   d=as.integer(as.vector(death)),
                   e=as.integer(as.vector(exposure)),
                   Tfor=forecast)
  out <- rstan::sampling(stanmodels$leecarterfinal, data = standata, ...)
  return(out)
}
