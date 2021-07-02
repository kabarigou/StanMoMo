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

#' boxplot_post_dist
#'
#' @param stan_fit a stan fit object
#' @param parm_name string to indicate the name of the parameter, to choose from c('a', 'b', 'g', 'k', 'k2', 'phi')
#'
#' @return posterior distribution shown as boxplots
#' @export
#'
#' @examples
#' years <- 1959:2019
#' ages <- 50:90
#' cohorts <- sort(unique(as.vector(sapply(years, function(year) year - ages))))
#' death <- load_HMD_data('BEL', 'Deaths_1x1', years, ages, "Male")$mat
#' exposure <- load_HMD_data('BEL', 'Exposures_1x1', years, ages, "Male")$mat
#' stan_fit <- fit_mo_mo("m6", death , exposure, ages, 0, 5, "nb", 1, 4, log_marg = F)
#' boxplot_post_dist(stan_fit, "k", ages, years)
#' boxplot_post_dist(stan_fit, "a", ages, years)
#'
boxplot_post_dist <- function(stan_fit, parm_name, ages, years){
  map <- extract_map(stan_fit)
  bool <- sapply(map, is.logical)
  parm_names <- names(map)[!bool]
  if(parm_name %in% parm_names){
    if(parm_name %in% c('a', 'b')){
      ages_df <- data.frame("stan" = names(dplyr::select(as.data.frame(stan_fit), starts_with(paste0(parm_name , "[")))),
                            "x" = ages)
      post_df <- dplyr::select(
        dplyr::left_join(
          tidyr::pivot_longer(
            dplyr::select(as.data.frame(stan_fit), starts_with(paste0(parm_name , "["))),
            starts_with(paste0(parm_name , "[")), names_to = "stan", values_to = "parm"),
          ages_df, by = "stan"),
        x, parm)
      ggplot(post_df, aes(x = as.factor(x), y =parm)) +
        geom_boxplot() +
        scale_x_discrete(breaks = post_df$x[seq(1, length(unique(post_df$x)), 5)]) +
        labs(x = "Ages", y = "", title = parm_name) + theme_bw(base_family='sans') +
        theme(
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 20),
          plot.title = element_text(size = 24, hjust = 0.5))
    }
    else if(parm_name %in% c('k', 'k2')){
      years_df <- data.frame("stan" = names(dplyr::select(as.data.frame(stan_fit), starts_with(paste0(parm_name , "[")))),
                             "t" = years[-1])
      post_df <- dplyr::select(
        dplyr::left_join(
          tidyr::pivot_longer(
            dplyr::select(as.data.frame(stan_fit), starts_with(paste0(parm_name , "["))),
            starts_with(paste0(parm_name , "[")), names_to = "stan", values_to = "parm"),
          years_df, by = "stan"),
        t, parm)
      ggplot(post_df, aes(x = as.factor(t), y = parm)) +
        geom_boxplot() +
        scale_x_discrete(breaks = post_df$t[seq(1, length(unique(post_df$t)), 10)]) +
        labs(x = "Calendar years", y = "", title = parm_name) + theme_bw(base_family='sans') +
        theme(
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 20),
          plot.title = element_text(size = 24, hjust = 0.5))
    }
    else if(parm_name ==  "g"){
      cohorts_df <- data.frame(c = sort(unique(as.vector(sapply(years[-1], function(year) year - ages)))),
                               g_stan = names(select(as.data.frame(stan_fit), starts_with("g["))))
      post_df <-select(
        left_join(
          pivot_longer(
            select(as.data.frame(fit_gen_model), starts_with("g[")),
            starts_with("g["),names_to = "g_stan", values_to = "g_t_x"),
          cohorts_df, by = "g_stan"),
        c, g_t_x)
      ggplot(post_df, aes(x = as.factor(c), y = g_t_x)) +
        geom_boxplot() +
        scale_x_discrete(breaks = post_df$c[seq(1, length(unique(post_df$c)), 15)]) +
        labs(x = "Cohorts", y = "", title = parm_name) + theme_bw(base_family='sans') +
        theme(
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 20),
          plot.title = element_text(size = 24, hjust = 0.5))
    }
    else if(parm_name %in% c('phi')){
      post_df <- dplyr::select(as.data.frame(fit_gen_model), starts_with("phi"))
      ggplot(post_df, aes(y = phi)) + geom_boxplot() +
        labs(x = "", y = "", title = parm_name) + theme_bw(base_family='sans') + ggtitle(parm_name)+
        theme(
          axis.title = element_text(size = 20),
          axis.text = element_text(size = 20),
          axis.text.x = element_blank(),
          plot.title = element_text(size = 24, hjust = 0.5))

    }
  }else{
    print(paste('parameter name invalid! Try either ', paste(parm_names, collapse = " ")))
  }
}
