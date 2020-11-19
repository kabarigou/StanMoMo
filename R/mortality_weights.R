# Save this file as `R/mortality_weights.R`

#' Model averaging/weighting via future-out stacking or pseudo-BMA weighting
#'
#' Mortality Model averaging via stacking of predictive distributions and Pseudo-BMA weighting. Based on Yao et al. (2018) but adapted by Barigou et al. (2020)
#' for mortality forecasting.
#'
#' Mortality model averaging via stacking of predictive distributions or pseudo-BMA weighting.
#' Both approaches were proposed in Yao et al. (2018) based leave-one-out cross-validation which is not suited for forecasting.
#' Barigou et al. (2020) adapted both appraoches based on leave-future-out validation which is more appropriate for mortality forecasting.
#'
#'@details
#'
#' The stacking method combines all models by maximizing the leave-future-out
#' predictive density of the combination distribution. That is, it finds the
#' optimal linear combining weights for maximizing the leave-future-out log score.
#'
#' The pseudo-BMA method finds the relative weights proportional to the expected
#' log predictive density of each model.
#'
#' Similar to Yao et al. (2018), we recommend stacking for averaging predictive distributions
#' as pseudo-BMA tends to select only one model.
#'
#'
#' @export
#' @param X A list of stanfit objects.
#' @return A matrix containing one weight for each model and each approach.
#'
#'@references
#' Yao, Y., Vehtari, A., Simpson, D., & Gelman, A. (2018). Using stacking to average Bayesian predictive distributions (with discussion).
#' Bayesian Analysis, 13(3), 917-1007.
#'
#' Barigou K., Goffard P-O., Loisel S., Salhi Y. (2020). Bayesian Model Averaging for mortality forecasting. Working paper.
#'
#'@examples
#'
#' \dontrun{
#' #10-year forecasts for French data for ages 50-90 and years 1970-2017 with a log-Poisson model
#' #where the 10 last years are held out for validation. We search for the model weights between
#' #the Lee-Carter model and the RH model (Lee-Carter with cohort effect).
#' ages.fit<-50:90
#' years.fit<-1970:2017
#' deathFR<-FRMaleData$Dxt[formatC(ages.fit),formatC(years.fit)]
#' exposureFR<-FRMaleData$Ext[formatC(ages.fit),formatC(years.fit)]
#' fitLC=lc_stan(death = deathFR,exposure=exposureFR, forecast = 10, validation=10,family = "poisson")
#' fitRH=rh_stan(death = deathFR,exposure=exposureFR, forecast = 10, validation=10,family = "poisson")
#' model_weights<-mortality_weights(list(fitLC,fitRH))
#' }
mortality_weights <- function(X) {

  log_lik_list <- lapply(X, loo::extract_log_lik,parameter_name = "log_lik2")
  log_sum_exp <- function(u) {
    max_u <- max(u);
    a <- 0;
    for (n in 1:length(u)) {
      a <- a + exp(u[n] - max_u);
    }
    return(max_u + log(a));
  }
  lpd<-lapply(log_lik_list,function(x) {apply(x,2,log_sum_exp)-log(nrow(x))})
  lpd_point<-simplify2array(lpd)

  stacking <- loo::stacking_weights(lpd_point)
  pseudobma <- loo::pseudobma_weights(lpd_point, BB=FALSE)
  return(round(cbind(stacking, pseudobma), 3))
}

#' compute_weights_BMA: compute the model evidence via bridge sampling (via the harmonic mean estimator if bridge sampling fails)
#'
#' @param stan_fits  list of Stan model fits where the marginal likelihood was computed via bridge sampling
#' @param mortality_models vector of mortality models names
#'
#' @return data frame with model evidence for BMA
#' @export
#'
#' @examples
#'
#'
compute_weights_BMA <- function(stan_fits, mortality_models){
  names(stan_fits) <- mortality_models
  log_marg <- sapply(mortality_models, function(mortality_model) stan_fits[[mortality_model]]$logml)
  if(any(is.na(log_marg))){
    log_lik_list <- sapply(mortality_models, function(mortality_model) rowSums(loo::extract_log_lik(stan_fits[[mortality_model]]$stan_output)))
    log_marg <- sapply(mortality_models, function(mortality_model) length(log_lik_list[,mortality_model])-log_sum_exp(-log_lik_list[,mortality_model]))
  }
  res <- data.frame(BMA = exp(log_marg - max(log_marg, na.rm = TRUE))/ sum(exp(log_marg - max(log_marg, na.rm = T)), na.rm = TRUE), fitted_model = mortality_models)
  rownames(res)<-NULL
  return(res)
}
