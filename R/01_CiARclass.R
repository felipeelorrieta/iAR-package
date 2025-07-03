#' `CiAR` Class
#'
#' Represents a complex irregular autoregressive (CiAR) time series model. 
#' This class extends the `unidata` class and provides additional properties 
#' for modeling, forecasting, and interpolating irregularly observed time 
#' series data with both negative and positive autocorrelation.
#'
#' @param times A numeric vector representing the time points.
#' @param series A complex vector representing the values of the time series.
#' @param series_esd A numeric vector representing the error standard deviations of the time series.
#' @param series_names An optional character vector of length 1 representing the name of the series. 
#' @param fitted_values A numeric vector containing the fitted values from the model.
#' @param kalmanlik A numeric value representing the Kalman likelihood of the model.
#' @param coef A numeric vector of length 2, containing the coefficients of the model. Each value must lie within [-1, 1]. Defaults to `c(0.9, 0)`.
#' @param tAhead A numeric value specifying the forecast horizon (default: 1).
#' @param forecast A numeric vector containing the forecasted values.
#' @param interpolated_values A numeric vector containing the interpolated values.
#' @param interpolated_times A numeric vector containing the times of the interpolated data points.
#' @param interpolated_series A numeric vector containing the interpolated series.
#' @param zero_mean A logical value indicating if the model assumes a zero-mean process (default: TRUE).
#' @param standardized A logical value indicating if the model assumes a standardized process (default: TRUE).
#'
#' @section Validation:
#' - Inherits all validation rules from the `unidata` class:
#'  - `@times`, `@series`, and `@series_esd` must be numeric vectors.
#'  - `@times` must not contain `NA` values and must be strictly increasing.
#'  - The length of `@series` must match the length of `@times`.
#'  - The length of `@series_esd` must be 0, 1, or equal to the length of `@series`.
#'  - `NA` values in `@series` must correspond exactly (positionally) to `NA` values in `@series_esd`.
#'  - `@series_names`, if provided, must be a character vector of length 1.
#' 
#' - `@coef` must be a numeric vector of length 2 with no dimensions.
#' - Each value in `@coef` must be in the interval [-1, 1].
#' - `@tAhead` must be a strictly positive numeric scalar.
#' 
#' @details
#' The `CiAR` class is designed to handle irregularly observed time series 
#' data with either negative or positive autocorrelation using an autoregressive 
#' approach. It extends the `unidata` class to include functionalities 
#' specific to the `CiAR` model.
#'
#' Key features of the `CiAR` class include:
#' - Support for irregularly observed time series data with negative 
#' or positive autocorrelation.
#' - Forecasting and interpolation functionalities for irregular time points.
#' - Configurable assumptions of zero-mean and standardized processes.
#'
#' @references
#' \insertRef{Elorrieta_2019}{iAR}
#'
#'
#' @examples
#' o=iAR::utilities()
#' o<-gentime(o, n=200, distribution = "expmixture", lambda1 = 130, lambda2 = 6.5,p1 = 0.15, p2 = 0.85)
#' times=o@times
#' my_CiAR <- CiAR(times = times,coef = c(0.9, 0))
#'
#' # Access properties
#' my_CiAR@coef
#'
#' @export
CiAR <- S7::new_class(
  "CiAR",
  parent = unidata,
  package = "iAR",
  properties = list(fitted_values = S7::class_numeric,
                    kalmanlik = S7::class_numeric,
                    coef = S7::new_property(
                      S7::class_numeric,
                      default = c(0.9, 0), # restringir cada parametro a -1,1
                      validator = function(value) {
                        if (!(is.numeric(value) &&
                              length(value) == 2 &&
                              is.null(dim(value)))) {
                          return("'coef' must be a numeric vector of length 2")
                        }
                        if (any(value < -1 | value > 1)) {
                          return("Each value in 'coef' must be between -1 and 1")
                        }
                        NULL
                      }
                    ),
                    tAhead = S7::new_property(
                      S7::class_numeric,
                      default = 1,
                      validator = function(value) {
                        if (!(length(value) == 1 &&
                              is.null(dim(value)) &&
                              value > 0)) {
                          return("'tAhead' must be a strictly positive numeric scalar")
                        }
                        NULL
                      }
                    ),
                    forecast = S7::class_numeric,
                    interpolated_values = S7::class_numeric,
                    interpolated_times = S7::class_numeric,
                    interpolated_series = S7::class_numeric,
                    zero_mean = S7::new_property(class_logical, default = TRUE),
                    standardized = S7::new_property(class_logical, default = TRUE))
)
