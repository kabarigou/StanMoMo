# Save this file as `R/load_HMD_data.R`

#' Download mortality data from the Human Mortality Database.
#'
#' Extract different quantities of interest from the HMD website, for instance deaths and exposures of a specific country
#' for a certain age range and year range.
#'
#' @param CNTRY String, HMD country name for instance 'FRATNP'.
#' @param item String, HMD data type for instance 'Deaths_1x1'.
#' @param cal_year_range Array, calendar year of interest, for instance 1970:2017.
#' @param age_range Array, ages of interest, for instance 50:90.
#' @param gender String, gender, for instance "Male".
#'
#' @return A named list containing the data either under a dataframe, a matrix
#' or a vector.
#' @export
#'
#' @examples
#' death_fra <- load_HMD_data('FRATNP', 'Deaths_1x1', 1970:2017, 50:90, "Male")
#' exposure_fra <- load_HMD_data('FRATNP', 'Exposures_1x1', 1970:2017, 50:90, "Male")

load_HMD_data <- function(CNTRY, item, cal_year_range, age_range, gender){
  Age<-NULL
  Year<-NULL
  username <- 'pierre.olivier.goffard@gmail.com'
  password <- 'StanMoMo'
  path <- paste0("https://www.mortality.org/hmd/", CNTRY,
                 "/STATS/", item)
  TEXT <- httr::GET(path, httr::authenticate(username, password),
                    httr::config(ssl_verifypeer = 0L))
  status <- httr::http_status(TEXT)
  DF <- utils::read.table(text = httr::content(TEXT, encoding = "UTF-8"),
                   header = TRUE, skip = 2, na.strings = ".",
                   as.is = TRUE)
  DF <- dplyr::filter(DF, (Age %in% age_range) & (Year %in% cal_year_range))
  DF <-dplyr::select(DF,Year, Age, gender)

  mat <- as.matrix(tibble::column_to_rownames(tidyr::pivot_wider(DF, names_from = Year, values_from = gender), var = "Age"))
  vect <- as.integer(as.vector(mat))
  return(
    list(DF = DF, mat = mat, vect = vect)
  )
}
