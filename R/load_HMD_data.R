# Save this file as `R/load_HMD_data.R`

#' load_HMD_data: a function to download mortality data from the HMD website
#'
#' @param CNTRY String, HMD country name for instance 'FRATNP'
#' @param item String, HMD data type for instance 'Deaths_1x1'
#' @param cal_year_range Array, calendar year of interest, for instance
#' 1970:2017
#' @param age_range Array, ages of interest, for instance 50:90
#' @param gender String, gender, for instance "Male"
#'
#' @return a named list containing the data either under a datafreame, a matrix
#' or a vector
#' @export
#'
#' @examples
#' death_fra <- load_HMD_data('FRATNP', 'Deaths_1x1', 1970:2017, 50:90, "Male")
#' exposure_fra <- load_HMD_data('FRATNP', 'Exposures_1x1', 1970:2017, 50:90, "Male")

load_HMD_data <- function(CNTRY, item, cal_year_range, age_range, gender){
  username <- 'pierre.olivier.goffard@gmail.com'
  password <- 'StanMoMo'
  path <- paste0("https://www.mortality.org/hmd/", CNTRY,
                 "/STATS/", item)
  TEXT <- httr::GET(path, httr::authenticate(username, password),
                    httr::config(ssl_verifypeer = 0L))
  status <- httr::http_status(TEXT)
  DF <- read.table(text = httr::content(TEXT, encoding = "UTF-8"),
                   header = TRUE, skip = 2, na.strings = ".",
                   as.is = TRUE)
  DF <- DF %>% filter( (Age %in% age_range) & (Year %in% cal_year_range) ) %>% select(Year, Age, gender)

  mat <- as.matrix(pivot_wider(DF, names_from = Year, values_from = gender) %>% column_to_rownames(., var = "Age"))
  vect <- as.integer(vec(mat))
  return(
    list(DF = DF, mat = mat, vect = vect)
  )
}
