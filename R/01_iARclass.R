#' `iAR` Class
#'
#' Represents a univariate irregular autoregressive (iAR) time series model. 
#' This class extends the "unidata" class and includes additional properties 
#' for modeling, forecasting, and interpolation.
#'
#' @param times A numeric vector representing the time points.
#' @param series A numeric vector representing the values of the time series.
#' @param series_esd A numeric vector representing the error standard deviations of the time series.
#' @param series_names An optional character vector of length 1 representing the name of the series. 
#' @param family A character string indicating the distribution family of the model. Must be one of `"norm"`, `"t"`, or `"gamma"`. Defaults to `"norm"`.
#' @param fitted_values A numeric vector containing the fitted values from the model.
#' @param loglik A numeric value representing the log-likelihood of the model.
#' @param kalmanlik A numeric value representing the Kalman likelihood of the model.
#' @param coef A numeric vector containing the estimated coefficients of the model (default: 0.9).
#' @param df A numeric value representing the degrees of freedom (`t` distribution) (default: 3).
#' @param sigma A numeric value representing the scale parameter (`t` distribution) (default: 1).
#' @param mean A numeric value representing the estimated mean of the model (`gamma` parameter).
#' @param variance A numeric value representing the estimated variance of the model (`gamma` parameter).
#' @param tAhead A numeric value specifying the forecast horizon (default: 1).
#' @param forecast A numeric vector containing the forecasted values.
#' @param interpolated_values A numeric vector containing the interpolated values.
#' @param interpolated_times A numeric vector containing the times of the interpolated data points.
#' @param interpolated_series A numeric vector containing the interpolated series.
#' @param zero_mean A logical value indicating if the model assumes a zero-mean process (default: TRUE).
#' @param standardized A logical value indicating if the model assumes a standardized process (default: TRUE).
#' @param hessian A logical value indicating whether the Hessian matrix is computed during estimation (default: FALSE).
#' @param summary A list containing the summary of the model fit, including diagnostics and statistical results.
#'
#' @section Validation Rules:
#' - Inherits all validation rules from the `unidata` class:
#'  - `@times`, `@series`, and `@series_esd` must be numeric vectors.
#'  - `@times` must not contain `NA` values and must be strictly increasing.
#'  - The length of `@series` must match the length of `@times`.
#'  - The length of `@series_esd` must be 0, 1, or equal to the length of `@series`.
#'  - `NA` values in `@series` must correspond exactly (positionally) to `NA` values in `@series_esd`.
#'  - `@series_names`, if provided, must be a character vector of length 1.
#' 
#' - `@family` must be one of `"norm"`, `"t"`, or `"gamma"`.
#' - `@coef` must be a numeric scalar strictly between 0 and 1.
#' - `@df` must be a positive integer scalar **only if** `family = "t"`; otherwise, it should not be specified.
#' - `@sigma` must be a strictly positive numeric scalar (used in `"t"` family).
#' - `@mean` and `@variance` must be strictly positive numeric scalars **only if** `family = "gamma"`.
#' - `@tAhead` must be a strictly positive numeric scalar.


#'
#' @details
#' The `iAR` class is designed to handle irregularly observed time series data using an 
#' autoregressive approach. It extends the "unidata" class to include additional 
#' modeling and diagnostic capabilities. Key functionalities include forecasting, 
#' interpolation, and model fitting.
#'
#' The class also supports advanced modeling features, such as:
#' - Different distribution families for the data (e.g., Gaussian, `t`-distribution).
#' - Optional computation of the Hessian matrix for parameter estimation.
#' - Standardized or zero-mean process assumptions.
#'
#' @references
#' \insertRef{Eyheramendy_2018}{iAR}
#'
#' @examples
#' # Create an `iAR` object
#' o=iAR::utilities()
#' o<-gentime(o, n=200, distribution = "expmixture", lambda1 = 130, lambda2 = 6.5,p1 = 0.15, p2 = 0.85)
#' times=o@times
#' my_iAR <- iAR(family = "norm", times = times, coef = 0.9,hessian=TRUE)
#'
#' my_iAR@family
#' my_iAR@coef
#'
#' @export
iAR <- S7::new_class(
  "iAR",
  parent = unidata,
  package = "iAR",
  properties = list(family = S7::new_property(S7::class_character,
                                              validator = function(value) {
                                                if (!value %in% c("norm", "t", "gamma")) {
                                                  return("family must be one of 'norm', 't', or 'gamma'")
                                                }
                                                NULL
                                              },
                                              default = "norm"),
                    fitted_values = S7::class_numeric,
                    loglik = S7::class_numeric,
                    kalmanlik = S7::class_numeric,
                    coef = S7::new_property(
                      S7::class_numeric,
                      default = 0.9,
                      validator = function(value) {
                        if (!(length(value) == 1 && is.null(dim(value)) && value > 0 && value < 1)) {
                          return("'coef' must be a value strictly between 0 and 1")
                        }
                        NULL
                      }
                    ),
                    df = S7::new_property(
                      S7::class_numeric,
                      default = 3,
                      validator = function(value) {
                        if (!(length(value) == 1 &&
                              is.null(dim(value)) &&
                              value > 0 &&
                              value %% 1 == 0)) {
                          return("'df' must be a positive integer value")
                        }
                        NULL
                      }
                    )
                    , # t #colocar df por defecto = 3
                    sigma = S7::new_property(
                      S7::class_numeric,
                      default = 1,
                      validator = function(value) {
                        if (!(length(value) == 1 &&
                              is.null(dim(value)) &&
                              value > 0)) {
                          return("'sigma' must be a positive numeric value")
                        }
                        NULL
                      }
                    ), # t
                    mean = S7::new_property(
                      S7::class_numeric,
                      default = 0), # gamma
                    variance = S7::new_property(
                      S7::class_numeric,
                      default = 1,
                      validator = function(value) {
                        if (!(length(value) == 1 &&
                              is.null(dim(value)) &&
                              value > 0)) {
                          return("'variance' must be a strictly positive numeric scalar")
                        }
                        NULL
                      }
                    ), #gamma
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
                    zero_mean = S7::new_property(class_logical, default = TRUE), # no gamma
                    standardized = S7::new_property(class_logical, default = TRUE), # no gamma # no t
                    hessian = S7::new_property(class_logical, default = FALSE),
                    summary = S7::class_list),
  validator = function(self) {
    fam <- self@family
    if (fam == "t") {
      if (is.na(self@df) || self@df <= 0 || self@df %% 1 != 0) {
        return("For family = 't', df must be a positive integer scalar")
      }
    }
    if (fam == "gamma") {
      if (is.na(self@mean) || is.na(self@variance)) {
        return("For family = 'gamma', both mean and variance must be specified")
      }
    }
    NULL
  }
)