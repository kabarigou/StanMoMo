# Save this file as `R/fit_mo_mo.R`

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
