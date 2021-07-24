#' Fanplot for the mortality predictions
#'
#' The function generates a fan plot representing the forecast death rates based on the \code{ggplot2} and
#' \code{ggfan} R packages.
#'
#' @param stan_fit stanfit object
#' @param ages range of ages
#' @param years range of years
#' @param death matrix of observed deaths
#' @param exposure matrix of observed exposures
#' @param ages.plot ages to be plotted
#'
#'
#' @export
#'
#' @examples
#'
#' \dontrun{
#' years <- 1980:2018
#' ages <- 50:90
#' death <- load_HMD_data('BEL', 'Deaths_1x1', years, ages, "Male")$mat
#' exposure <- load_HMD_data('BEL', 'Exposures_1x1', years, ages, "Male")$mat
#' fitLC<-lc_stan(death = deathFR,exposure=exposureFR, forecast = 10,
#' family = "poisson",chains=1,cores=1)
#' forecasting_plot(fitLC,ages,years,death,exposure,c(65,75,85))
#' }
#'
forecasting_plot <- function(stan_fit, ages, years,death, exposure,ages.plot){
  Age <- Obs <- Year <- NULL
  params<-rstan::extract(stan_fit)
  nforecast<-ncol(params$k_p)
  years.predict<-seq(years[length(years)]+1,years[length(years)]+nforecast,by=1)

  pred<-array(params$mufor,dim=list(nrow(params$mufor),ncol(params$mufor)/nforecast,nforecast),
              dimnames = list(c(1:nrow(params$mufor)),formatC(ages),formatC(years.predict)))

  pred<-as.data.frame.table(pred)
  colnames(pred)<-c("Sim","Age","Year","Obs")
  pred$Age <-as.numeric(as.character(pred$Age))
  pred$Year <-as.numeric(as.character(pred$Year))

  qxt <- death / exposure
  qxt<-as.data.frame.table(qxt)
  colnames(qxt)<-c("Age","Year","Obs")
  qxt$Age <-as.numeric(as.character(qxt$Age))
  qxt$Year <-as.numeric(as.character(qxt$Year))

  years<-c(years,years.predict)

  ggplot2::ggplot()+
    ggplot2::geom_point(data=subset(qxt,Age %in% ages.plot & Year %in% years),ggplot2::aes(x=Year,y=Obs,shape=as.factor(Age)))+
    ggplot2::scale_y_log10()+
    ggfan::geom_fan(data = subset(pred,Age %in% ages.plot &Year %in% years),ggplot2::aes(x=Year,y=Obs,group=Age))+
    ggplot2::theme_bw()+
    ggplot2::scale_fill_distiller(palette="Spectral")+
    ggplot2::labs(y="death rate (log scale)", x="Year",shape="Age")

}
