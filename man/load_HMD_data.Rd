% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load_HMD_data.R
\name{load_HMD_data}
\alias{load_HMD_data}
\title{Download mortality data from the Human Mortality Database.}
\usage{
load_HMD_data(CNTRY, item, cal_year_range, age_range, gender)
}
\arguments{
\item{CNTRY}{String, HMD country name for instance 'FRATNP'.}

\item{item}{String, HMD data type for instance 'Deaths_1x1'.}

\item{cal_year_range}{Array, calendar year of interest, for instance 1970:2017.}

\item{age_range}{Array, ages of interest, for instance 50:90.}

\item{gender}{String, gender, for instance "Male".}
}
\value{
A named list containing the data either under a dataframe, a matrix
or a vector.
}
\description{
Extract different quantities of interest from the HMD website, for instance deaths and exposures of a specific country
for a certain age range and year range.
}
\examples{
death_fra <- load_HMD_data('FRATNP', 'Deaths_1x1', 1970:2017, 50:90, "Male")
exposure_fra <- load_HMD_data('FRATNP', 'Exposures_1x1', 1970:2017, 50:90, "Male")
}
