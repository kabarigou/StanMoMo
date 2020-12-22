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

#' sim_death_rh simulate death counts from Renshaw-Haberman mortality model
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
sim_death_rh <- function(a, b,k, g, phi, years, ages, exposure){

  cohorts <- sort(unique(as.vector(sapply(years, function(year) year - ages))))

  gxt_apc <- 0 * exposure[,1:length(k)]
  for(i in 1:length(a)){
    for(j in 1:length(k)){
      gxt_apc[i,j] <- exp(a[i] + b[i]*k[j] + g[match(years[j] - ages[i], cohorts)]) * exposure[i,j]
    }
  }
  return(apply(gxt_apc, 1:2, function(gxt) rnbinom(1,size = phi, prob = phi / (phi + gxt)))
  )
}

#' sim_death_cbd simulate death counts from the CBD model
#'
#' @param k first vector of period component
#' @param k2 second vector of period component
#' @param phi dispersion parameter
#' @param years vector of calendar years
#' @param ages vectors of ages
#' @param exposure matrix of exposure data
#'
#' @return matrix of death count
#' @export
#'
#' @examples
sim_death_cbd <- function(k, k2, phi, years, ages, exposure){


  gxt_cbd <- 0 * exposure[,1:length(k)]
  for(i in 1:length(ages)){
    for(j in 1:length(k)){
      gxt_cbd[i,j] <- exp(k[j] +(ages[i]-mean(ages))*k2[j]) * exposure[i,j]
    }
  }
  return(apply(gxt_cbd, 1:2, function(gxt) rnbinom(1,size = phi, prob = phi / (phi + gxt)))
  )
}

#' sim_death_m6 simulates death counts from the M6 model
#'
#' @param k first vector of period component
#' @param k2 second vector of period component
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
sim_death_m6 <- function(k, k2, g, phi, years, ages, exposure){
  cohorts <- sort(unique(as.vector(sapply(years, function(year) year - ages))))

  gxt_m6 <- 0 * exposure[,1:length(k)]
  for(i in 1:length(ages)){
    for(j in 1:length(k)){
      gxt_m6[i,j] <- exp(k[j] +(ages[i]-mean(ages))* k2[j] + g[match(years[j] - ages[i], cohorts)]) * exposure[i,j]
    }
  }
  return(apply(gxt_m6, 1:2, function(gxt) rnbinom(1,size = phi, prob = phi / (phi + gxt)))
  )
}


#' sim_mortality_data simulate mortality data from various models
#'
#' @param a vector of age component
#' @param k first vector of time component
#' @param k2 second vector of time component
#' @param b vector of age/time component
#' @param g vector of cohort component
#' @param phi dispersioon parameter
#' @param years vector of calendar year
#' @param ages vector of ages
#' @param exposure matrix of exposure
#' @param mortality_model name of the mortality model that we simulate from
#'
#' @return
#' @export
#'
#' @examples
sim_mortality_data <- function(a, k, k2,b, g, phi, years, ages, exposure, mortality_model){
  if(mortality_model == "lc"){
    res <- sim_death_lc(a, b, k, phi, exposure)
  }else if(mortality_model == "apc"){
    res <- sim_death_apc(a, k, g, phi, years, ages, exposure)
  }else if (mortality_model=="rh"){
    res <- sim_death_rh(a, b,k, g, phi, years, ages, exposure)
  }else if (mortality_model=="cbd"){
    res <- sim_death_cbd(k, k2,phi, years, ages, exposure)
  }else if (mortality_model=="m6"){
    res <- sim_death_m6(k, k2,g,phi, years, ages, exposure)
  }
  return(res)
}

