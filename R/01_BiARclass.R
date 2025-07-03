#' `BiAR` Class
#'
#' Represents a bivariate irregular autoregressive (BiAR) time series model.
#' This class extends the `multidata` class and provides additional properties
#' for modeling, forecasting, and interpolation of bivariate time series data.
#'
#' @param times A numeric vector representing the time points.
#' @param series A numeric matrix or vector representing the values of the time series.
#' @param series_esd A numeric matrix or vector representing the error standard deviations of the time series.
#' @param series_names An optional character vector representing the name of the series. 
#' @param fitted_values A numeric vector containing the fitted values from the model.
#' @param loglik A numeric value representing the log-likelihood of the model.
#' @param kalmanlik A numeric value representing the Kalman likelihood of the model.
#' @param coef A numeric vector containing the estimated coefficients of the model.
#' @param rho A numeric vector containing the estimated coefficients of the model.
#' @param tAhead A numeric value specifying the forecast horizon (default: 1).
#' @param forecast A numeric vector containing the forecasted values.
#' @param interpolated_values A numeric vector containing the interpolated values.
#' @param interpolated_times A numeric vector containing the times of the interpolated data points.
#' @param interpolated_series A numeric vector containing the interpolated series.
#' @param zero_mean A logical value indicating if the model assumes a zero-mean process (default: TRUE).
#' @param standardized A logical value indicating if the model assumes a standardized process (default: TRUE).
#'
#' @section Validation Rules:
#' - `@times` must be a numeric vector without dimensions and strictly increasing.
#' - `@series` must be a numeric matrix with two columns (bivariate) or be empty.
#' - The number of rows in `@series` must match the length of `@times`.
#' - `@series_esd`, if provided, must be a numeric matrix. Its dimensions must match those of `@series`, or it must have one row and the same number of columns.
#' - If `@series_esd` contains NA values, they must correspond positionally to NA values in `@series`.
#' - `@series_names`, if provided, must be a character vector with length equal to the number of columns in `@series`, and all names must be unique.
#' - `@coef` must be a numeric vector of length 2, with each element strictly between -1 and 1.
#' - `@tAhead` must be a strictly positive numeric scalar.
#'
#' @details
#' The `BiAR` class is designed to handle bivariate irregularly observed time series data
#' using an autoregressive approach. It extends the `multidata` class to include 
#' additional properties for modeling bivariate time series.
#'
#' Key features of the `BiAR` class include:
#' - Support for bivariate time series data.
#' - Forecasting and interpolation functionalities for irregular time points.
#' - Assumptions of zero-mean and standardized processes, configurable by the user.
#' - Estimation of model parameters and likelihoods, including Kalman likelihood.
#'
#' @references
#' \insertRef{Elorrieta_2021}{iAR}
#'
#'
#' @examples
#' o=iAR::utilities()
#' o<-gentime(o, n=200, distribution = "expmixture", lambda1 = 130, lambda2 = 6.5,p1 = 0.15, p2 = 0.85)
#' times=o@times
#' my_BiAR <- BiAR(times = times,coef = c(0.9, 0.3), rho = 0.9)
#'
#' # Access properties
#' my_BiAR@coef
#'
#' @export
BiAR <- S7::new_class(
  "BiAR",
  parent = multidata,
  package = "iAR",
  properties = list(fitted_values = S7::class_numeric,
                    loglik = S7::class_numeric,
                    kalmanlik = S7::class_numeric,
                    coef = S7::new_property(
                      S7::class_numeric,
                      default = c(0.8, 0), # igual que el ciar
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
                    rho = S7::new_property(class_numeric, default = 0),
                    forecast = S7::class_numeric,
                    interpolated_values = S7::class_numeric,
                    interpolated_times = S7::class_numeric,
                    interpolated_series = S7::class_numeric,
                    zero_mean = S7::new_property(class_logical, default = TRUE),
                    standardized = S7::new_property(class_logical, default = TRUE))
)