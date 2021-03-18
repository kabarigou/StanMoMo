#' The 'StanMoMo' package.
#'
#' @description The StanMoMo package performs Bayesian Mortality Modeling with Stan for a variety of popular mortality models.
#' The current package supports the Lee-Carter (LC) model, the Renshaw-Haberman model (LC with cohort effect),
#' the Age-Period-Cohort (APC) model, the Cairns-Blake-Dowd (CBD) model and the M6 model (CBD with cohort effect).
#' By a simple call, the user inputs deaths and exposures and the package outputs the MCMC simulations for each parameter,
#' the log likelihoods and predictions. Moreover, the package includes tools for model selection and Bayesian model averaging
#' by leave-future-out validation.
#'
#' @docType package
#' @name StanMoMo-package
#' @aliases StanMoMo
#' @useDynLib StanMoMo, .registration = TRUE
#' @import methods
#' @import Rcpp
#' @importFrom rstan sampling
#'
#' @references
#' Stan Development Team (2020). RStan: the R interface to Stan. R package version 2.19.3. https://mc-stan.org
#'
NULL
