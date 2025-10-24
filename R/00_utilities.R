#' Utilities Class
#'
#' The `utilities` class is an S7 class designed to group utility functions 
#' that are not directly tied to other specific objects in the package. 
#' These include the functions `gentime`, `paringits`, `harmonicfit`, and `foldlc`.
#'
#' @param times A numeric vector storing the time points.
#' @param series A numeric vector representing the values of the time series.
#' @param series_esd A numeric vector representing the error standard deviations of the time series.
#' @param times_phased A numeric vector of phased times.
#' @param series_phased A numeric vector containing the time series values ordered by phase.
#' @param series_esd_phased A numeric vector containing the error standard deviations of the time series ordered by phase. 
#' @param fitted_values A numeric vector containing the fitted values from the harmonic model.
#' @param residuals A numeric vector containing the residuals from the harmonic model.
#' @param coef A numeric vector containing the estimated coefficients of the harmonic model.
#' @param summary A summary object containing the harmonic model information.
#' @param paired Data Frame with the paired datasets.
#'
#' @section Description:
#' This class acts as a container for standalone methods that perform independent 
#' operations within the package. By grouping them under a single class, the package 
#' achieves better modularity and organization, facilitating maintenance and extensibility.
#'
#' @section Available Methods:
#' - `gentime`: Generates time points based on a specified statistical distribution.
#' - `paringits`: A method for pairing irregular time series.
#' - `harmonicfit`: A method for fitting harmonic models to data.
#' - `foldlc`: A method for folding light curves in time series analysis.
#'
#' @examples
#' # Create a utilities object
#'
#' @export
utilities <- S7::new_class(
  "utilities",
  package = "iAR",
  properties = list(
    times = S7::class_numeric,
    series = S7::class_numeric,
    series_esd = S7::class_numeric,
    times_phased = S7::class_numeric,
    series_phased = S7::class_numeric,
    series_esd_phased = S7::class_numeric,
    fitted_values = S7::class_numeric,
    residuals = S7::class_numeric,
    coef = S7::class_numeric,
    summary = S7::class_list,
    paired = S7::class_numeric
  )
)
