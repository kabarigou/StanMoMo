# Save this file as `R/mortality_weights.R`

#' Mortality Model averaging via stacking of predictive distributions, pseudo-BMA weighting and pseudo-BMA+ 
#' weighting with the Bayesian bootstrap. Based on Yao et al. (2018) but adapted by Barigou et al. (2020) 
#' for mortality forecasting.
#' 
#' 
#' @export
#' @param X A list of stanfit objects
#' @return A matrix containing one weight for each model and each approach.
#'
#'
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
  pseudobmaplus <- loo::pseudobma_weights(lpd_point) # default is BB=TRUE
  return(round(cbind(stacking, pseudobma, pseudobmaplus), 3))
}
