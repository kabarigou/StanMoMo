#' The 'StanMoMo' package.
#'
#' @description The \pkg{StanMoMo} package performs Bayesian Mortality Modeling with
#' **Stan** for a variety of popular mortality models.
#' The current package supports the Lee-Carter model, the Renshaw-Haberman model,
#' the Age-Period-Cohort model, the Cairns-Blake-Dowd model and the M6 model. By a simple call, the user inputs deaths and exposures and the package outputs the MCMC simulations for each parameter,
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
