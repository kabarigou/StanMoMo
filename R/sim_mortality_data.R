# Save this file as `R/sim_mortality_data.R`

#' Simulation of death counts from the Lee-Carter mortality model
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
sim_death_lc <- function(a, b, k, phi, exposure){

  gxt_lc <- exp(sapply(k, function(kt) a + b * kt)) * exposure[, 1:length(k)]
  return(apply(gxt_lc, 1:2, function(gxt) stats::rnbinom(1,size = phi, prob = phi / (phi + gxt)))
  )
}

#' Simulation of death counts from the Age-Period-Cohort mortality model
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
sim_death_apc <- function(a, k, g, phi, years, ages, exposure){

  cohorts <- sort(unique(as.vector(sapply(years, function(year) year - ages))))

  gxt_apc <- 0 * exposure[,1:length(k)]
  for(i in 1:length(a)){
    for(j in 1:length(k)){
      gxt_apc[i,j] <- exp(a[i] + k[j] + g[match(years[j] - ages[i], cohorts)]) * exposure[i,j]
    }
  }
  return(apply(gxt_apc, 1:2, function(gxt) stats::rnbinom(1,size = phi, prob = phi / (phi + gxt)))
  )
}

#' Simulation of death counts from the Renshaw-Haberman mortality model
#'
#' @param a vector of age component
#' @param b vector of age/year component
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
sim_death_rh <- function(a, b,k, g, phi, years, ages, exposure){

  cohorts <- sort(unique(as.vector(sapply(years, function(year) year - ages))))

  gxt_rh <- 0 * exposure[,1:length(k)]
  for(i in 1:length(a)){
    for(j in 1:length(k)){
      gxt_rh[i,j] <- exp(a[i] + b[i]*k[j] + g[match(years[j] - ages[i], cohorts)]) * exposure[i,j]
    }
  }
  return(apply(gxt_rh, 1:2, function(gxt) stats::rnbinom(1,size = phi, prob = phi / (phi + gxt)))
  )
}

#' Simulation of death counts from the CBD model
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
sim_death_cbd <- function(k, k2, phi, years, ages, exposure){


  gxt_cbd <- 0 * exposure[,1:length(k)]
  for(i in 1:length(ages)){
    for(j in 1:length(k)){
      gxt_cbd[i,j] <- exp(k[j] +(ages[i]-mean(ages))*k2[j]) * exposure[i,j]
    }
  }
  return(apply(gxt_cbd, 1:2, function(gxt) stats::rnbinom(1,size = phi, prob = phi / (phi + gxt)))
  )
}

#' Simulation of death counts from the M6 model
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
sim_death_m6 <- function(k, k2, g, phi, years, ages, exposure){
  cohorts <- sort(unique(as.vector(sapply(years, function(year) year - ages))))

  gxt_m6 <- 0 * exposure[,1:length(k)]
  for(i in 1:length(ages)){
    for(j in 1:length(k)){
      gxt_m6[i,j] <- exp(k[j] +(ages[i]-mean(ages))* k2[j] + g[match(years[j] - ages[i], cohorts)]) * exposure[i,j]
    }
  }
  return(apply(gxt_m6, 1:2, function(gxt) stats::rnbinom(1,size = phi, prob = phi / (phi + gxt)))
  )
}

#' Simulation of death counts from a hybrid model that averages the mortality
#' rates from the cbd and rh models
#'
#' @param params_cbd named lsit that contains the parameters of the cbd model
#' @param params_rh named lsit that contains the parameters of the rh model
#' @param years vector of calendar year
#' @param ages vector of ages
#' @param exposure matrix of exposure data
#' @param q mixing parameter (0 <- rh, 1 <- cbd)
#'
#' @return matrix of death count
#' @export
#'
sim_death_mix_cbd_rh <- function(params_cbd, params_rh, years, ages, exposure, q){

  cohorts <- sort(unique(as.vector(sapply(years, function(year) year - ages))))
  gxt_rh <- 0 * exposure[,1:length(params_cbd$k)]
  gxt_cbd <- 0 * exposure[,1:length(params_cbd$k)]
  gxt_mix <- 0 * exposure[,1:length(params_cbd$k)]
  for(i in 1:length(ages)){
    for(j in 1:length(params_cbd$k)){
      gxt_cbd[i,j] <- exp(params_cbd$k[j] +(ages[i]-mean(ages))*params_cbd$k2[j]) * exposure[i,j]
      gxt_rh[i,j] <- exp(params_rh$a[i] + params_rh$b[i] * params_rh$k[j] +
                           params_rh$g[match(years[j] - ages[i], cohorts)]) * exposure[i,j]
      gxt_mix[i,j] <- q * gxt_cbd[i,j] + (1 - q) * gxt_rh[i,j]
    }
  }
  phi <- q * params_cbd$phi + (1 - q) * params_rh$phi
  return(apply(gxt_mix, 1:2, function(gxt) stats::rnbinom(1,size = phi, prob = phi / (phi + gxt)))
  )
}



#' Simulation of mortality data from various models
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
#' @return matrix of death counts
#' @export
#'
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

