# Save this file as `R/fit_mortality_model.R`

#' Wrapper function to fit and forecast mortality models
#'
#' @param mortality_model name of the mortality model
#' @param death death matrix
#' @param exposure exposure matrix
#' @param ages vector of ages
#' @param validation size of the validation set
#' @param forecast number of calendar years to be forecast
#' @param family underlying count distribution
#' @param chains number of Markov chains
#' @param cores number of cores used
#' @param log_marg Do we compute the marginal likelihood or not?
#' @param iter Length of the Markov chain trajectory
#'
#' @return a stanfit object
#' @export
#'
fit_mo_mo <- function(mortality_model ="lc", death, exposure, ages = 50:90, validation = 0, forecast = 1, family = "nb",
                      chains=1, cores=4, log_marg = F, iter = 2000){
  if(!log_marg){
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

  }else{
    if(mortality_model == "lc"){

      res <- lc_stan(death = death, exposure=exposure, validation=validation, forecast = forecast, family = family ,chains=chains,cores=cores, iter = iter)
      logml <- bridgesampling::bridge_sampler(res, silent = TRUE)$logml
    }else if(mortality_model == "apc"){

      res <- apc_stan(death = death,exposure=exposure, validation=validation,forecast = forecast, family = family,chains=chains,cores=cores, iter = iter)
      logml <- bridgesampling::bridge_sampler(res, silent = TRUE)$logml
    }else if(mortality_model == "cbd"){

      res <- cbd_stan(death = death,exposure=exposure, age=ages,
                      validation = validation, forecast = forecast, family = family,chains=chains,cores=cores, iter = iter)
      logml <- bridgesampling::bridge_sampler(res, silent = TRUE)$logml

    }else if(mortality_model == "rh"){

      res <- rh_stan(death = death, exposure=exposure, validation=validation,forecast = forecast, family = family, chains=chains,cores=cores, iter = iter)
      logml <- bridgesampling::bridge_sampler(res, silent = TRUE)$logml

    }else if(mortality_model == "m6"){

      res <- m6_stan(death = death, exposure=exposure,  age=ages, validation=validation,forecast = forecast, family = family, chains=chains,cores=cores, iter = iter)
      logml <- bridgesampling::bridge_sampler(res, silent = TRUE)$logml


    }
    return(list(stan_output = res, logml = logml))
  }
}



#' Function to get the a posterior means of the parameters based on a stanfit object
#'
#' @param stan_fit a stanfit object
#'
#' @return named list with the point estimates of the parameters
#' @export
#'
extract_map <- function(stan_fit){
  post_mean <-  dplyr::summarise(dplyr::select(as.data.frame(stan_fit),
                                        tidyselect::starts_with('a['),
                                        tidyselect::starts_with('b['),
                                        tidyselect::starts_with('k['),
                                        tidyselect::starts_with('k2['),
                                        tidyselect::starts_with('g['),
                                        'phi'),
                          dplyr::across(tidyselect::everything(), mean))
  res <- list(a = as.vector(t(dplyr::select(post_mean, tidyselect::starts_with('a[')))),
              b = as.vector(t(dplyr::select(post_mean,tidyselect::starts_with('b[')))),
              k = as.vector(t(dplyr::select(post_mean, tidyselect::starts_with('k[')))),
              k2 = as.vector(t(dplyr::select(post_mean, tidyselect::starts_with('k2[')))),
              g = as.vector(t(dplyr::select(post_mean, tidyselect::starts_with('g[')))),
              phi = as.vector(t(dplyr::select(post_mean, tidyselect::starts_with('phi'))))
  )
  return(res)
}
