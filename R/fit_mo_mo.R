# Save this file as `R/fit_mo_mo.R`

#' fit_mo_mo return the posterior samples as stan object of the fitted mortality model
#'
#' @param mortality_model name of teh mortality model
#' @param death death matrix
#' @param exposure exposure matrix
#' @param ages vector of ages
#' @param validation size of the validation set
#' @param forecast number of calendar year to be forecast
#' @param family underlying count distribution
#' @param chains number of Markov chains
#' @param cores number of cores used
#'
#' @return a stanfit object
#' @export
#'
#' @examples
fit_mo_mo <- function(mortality_model ="lc", death = deathGBR, exposure = exposureGBR, ages = 50:90, validation = 0, forecast = 1, family = "nb",
                      chains=1, cores=4){

  if(mortality_model == "lc"){

    res <- lc_stan(death = death, exposure=exposure, validation=validation, forecast = forecast, family = family ,chains=chains,cores=cores)

  }else if(mortality_model == "apc"){

    res <- apc_stan(death = death,exposure=exposure, validation=validation,forecast = forecast, family = family,chains=chains,cores=cores)

  }else if(mortality_model == "cbd"){

    res <- cbd_stan(death = death,exposure=exposure, age=ages,
                    validation = validation, forecast = forecast, family = family,chains=chains,cores=cores)

  }else if(mortality_model == "rh"){
    res <- rh_stan(death = death, exposure=exposure, validation=validation,forecast = forecast, family = family, chains=chains,cores=cores)
  }else if(mortality_model == "m6"){
    res <- m6_stan(death = death, exposure=exposure,  age=ages, validation=validation,forecast = forecast, family = family, chains=chains,cores=cores)
  }
  return(res)
}

#' extract_map: function to get the mean a posteriori of the parameters based on a stanfit object
#'
#' @param stan_fit a stanfit object
#'
#' @return named list with the point estimates of the parameters
#' @export
#'
#' @examples
extract_map <- function(stan_fit){
  post_mean <-  summarise(dplyr::select(as.data.frame(stan_fit),
                                        starts_with('a['),
                                        starts_with('b['),
                                        starts_with('k['),
                                        starts_with('k2['),
                                        starts_with('g['),
                                        'phi'),
                          across(everything(), mean))
  res <- list(a = as.vector(t(select(post_mean, starts_with('a[')))),
              b = as.vector(t(select(post_mean,starts_with('b[')))),
              k = as.vector(t(select(post_mean, starts_with('k[')))),
              k2 = as.vector(t(select(post_mean, starts_with('k2[')))),
              g = as.vector(t(select(post_mean, starts_with('g[')))),
              phi = as.vector(t(select(post_mean, starts_with('phi'))))
  )
  return(res)
}
