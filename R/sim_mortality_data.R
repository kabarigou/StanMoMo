# Save this file as `R/sim_mortality_data.R`

#' sim_death_lc simulate death counts from Lee-Carter mortality model
#'
#' @param a vector of age component
#' @param b vector of age/year component
#' @param k vector of year component
#' @param phi dispersion parameter
#' @param exposure matrix of exposure data
#'
#' @return matrix of death count
#' @export
#'
#' @examples
sim_death_lc <- function(a, b, k, phi, exposure){

  gxt_lc <- exp(sapply(k, function(kt) a + b * kt)) * exposure[, 1:length(k)]
  return(apply(gxt_lc, 1:2, function(gxt) rnbinom(1,size = phi, prob = phi / (phi + gxt)))
         )
}

#' sim_death_apc simulate death counts from Age-Period-Cohort mortality model
#'
#' @param a vector of age component
#' @param k vector of period component
#' @param g vector of cohort component
#' @param phi dispersion parameter
#' @param years vector of calendar years
#' @param ages vectors of ages
#' @param exposure matrix of exposure data
#'
#' @return matrix of death count
#' @export
#'
#' @examples
sim_death_apc <- function(a, k, g, phi, years, ages, exposure){

  cohorts <- sort(unique(as.vector(sapply(years, function(year) year - ages))))

  gxt_apc <- 0 * exposure[,1:length(k)]
  for(i in 1:length(a)){
    for(j in 1:length(k)){
      gxt_apc[i,j] <- exp(a[i] + k[j] + g[match(years[j] - ages[i], cohorts)]) * exposure[i,j]
    }
  }
  return(apply(gxt_apc, 1:2, function(gxt) rnbinom(1,size = phi, prob = phi / (phi + gxt)))
  )
}

sim_mortality_data <- function(a, k, b, g, phi, years, ages, exposure, mortality_model){
  if(mortality_model == "lc"){
    res <- sim_death_lc(a, b, k, phi, exposure)
  }else if(mortality_model == "apc"){
    res <- sim_death_apc(a, k, g, phi, years, ages, exposure)
  }
  return(res)
}

