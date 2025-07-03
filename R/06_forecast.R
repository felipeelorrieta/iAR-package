#' Forecast for iAR, CiAR, and BiAR Classes
#'
#' Generates forecasts for the specified time series model. This method is implemented for:
#' 1. Irregular Autoregressive models (`iAR`)
#' 2. Complex Irregular Autoregressive models (`CiAR`)
#' 3. Bivariate Autoregressive models (`BiAR`)
#' 
#' @name forecast
#'
#' @param x An object of class \code{iAR}, \code{CiAR}, or \code{BiAR}, containing the model specification and parameters:
#'   \itemize{
#'     \item For \code{iAR}:
#'       \itemize{
#'         \item \code{family}: The distribution family of the iAR model (one of "norm", "t", or "gamma").
#'         \item \code{series}: A numeric vector representing the time series to forecast.
#'         \item \code{coef}: The coefficient(s) of the iAR model.
#'         \item \code{zero_mean}: Logical, whether the model assumes a zero-mean series.
#'         \item \code{standardized}: Logical, whether the model uses standardized data (only for "norm" family).
#'         \item \code{mean}: The mean parameter (only for "gamma" family).
#'         \item \code{tAhead}: Integer, the number of steps ahead to forecast.
#'       }
#'     \item For \code{CiAR}:
#'       \itemize{
#'         \item \code{coef}: The real and imaginary parts of the CiAR model's coefficients.
#'         \item \code{series}: A numeric vector representing the time series to forecast.
#'         \item \code{times}: A numeric vector specifying the time points of the series.
#'         \item \code{zero_mean}: Logical, whether the model assumes a zero-mean series.
#'         \item \code{standardized}: Logical, whether the model output should be standardized.
#'         \item \code{tAhead}: Integer, the number of steps ahead to forecast.
#'       }
#'     \item For \code{BiAR}:
#'       \itemize{
#'         \item \code{coef}: The coefficients of the BiAR model (real and imaginary parts).
#'         \item \code{series}: A numeric matrix with two columns representing the bivariate time series to forecast.
#'         \item \code{times}: A numeric vector specifying the time points of the series.
#'         \item \code{tAhead}: Integer, the number of steps ahead to forecast.
#'       }
#'   }
#'
#' @param ... Additional arguments (unused).
#'
#' @return An updated object of class \code{iAR}, \code{CiAR}, or \code{BiAR}, where the \code{forecast} property contains the forecasted values.
#'
#' @details
#' This method generates forecasts for the specified time series model. Depending on the class of the input object:
#' \itemize{
#'   \item For \code{iAR}, the function supports three distribution families:
#'     \item "norm" for normal distribution.
#'     \item "t" for t-distribution.
#'     \item "gamma" for gamma distribution.
#'   \item For \code{CiAR}, the function uses complex autoregressive processes.
#'   \item For \code{BiAR}, the function generates forecasts for a bivariate autoregressive process.
#' }
#' All required parameters (e.g., coefficients, time points) must be set before calling this method.
#'
#' @references
#' \insertRef{Eyheramendy_2018}{iAR},\insertRef{Elorrieta_2019}{iAR},\insertRef{Elorrieta_2021}{iAR}
#'
#' @examples
#' # Example 1: Forecasting with a normal iAR model
#' library(iAR)
#' n=100
#' set.seed(6714)
#' o=iAR::utilities()
#' o<-gentime(o, n=n)
#' times=o@times
#' model_norm <- iAR(family = "norm", times = times, coef = 0.9)  
#' model_norm <- sim(model_norm)
#' model_norm <- kalman(model_norm) 
#' model_norm@tAhead=1.3
#' model_norm <- forecast(model_norm)
#' plot(times, model_norm@series, type = "l", main = "Original Series with Forecast")
#' points(max(times)+ model_norm@tAhead, model_norm@forecast, col = "blue", pch = 16)
#' plot_forecast(model_norm)
#' 
#' # Example 2: Forecasting with a CiAR model
#' set.seed(6714)
#' model_CiAR <- CiAR(times = times,coef = c(0.9, 0))
#' model_CiAR <- sim(model_CiAR)
#' y=model_CiAR@series
#' y1=y/sd(y)
#' model_CiAR@series=y1
#' model_CiAR@series_esd=rep(0,n)
#' model_CiAR <- kalman(model_CiAR)
#' print(model_CiAR@coef)
#' model_CiAR@tAhead=1.3
#' model_CiAR <-forecast(model_CiAR)
#' model_CiAR@forecast
#' 
#' # Example 3: Forecasting with a BiAR model
#' n=80
#' set.seed(6714)
#' o=iAR::utilities()
#' o<-gentime(o, n=n)
#' times=o@times
#' model_BiAR <- BiAR(times = times,coef = c(0.9, 0.3), rho = 0.9)
#' model_BiAR <- sim(model_BiAR)
#' y=model_BiAR@series
#' y1=y/apply(y,2,sd)
#' model_BiAR@series=y1
#' model_BiAR@series_esd=matrix(0,n,2)
#' model_BiAR <- kalman(model_BiAR)
#' print(model_BiAR@coef)
#' model_BiAR@tAhead=1.3
#' model_BiAR <-forecast(model_BiAR)
#' model_BiAR@forecast
#'
#' @export
forecast <- S7::new_generic("forecast", "x")
S7::method(forecast, iAR) <- function(x) {
  if(length(x@series)==0) stop("The forecast method needs a time series")
  if(x@family == "norm"){
    if(length(x@coef)==0) stop("The forecast method needs the coefficient of the iAR model")
    res <- iARforecast(coef = x@coef,
                       series = x@series,
                       zero_mean = x@zero_mean,
                       standardized = x@standardized,
                       tAhead = x@tAhead)
    x@forecast <- res
    return(x)
  }
  
  if(x@family == "t"){
    if(length(x@coef)==0) stop("The forecast method needs the coefficient of the iAR-T model")
    res <- iARforecast(coef = x@coef,
                       series = x@series,
                       zero_mean = x@zero_mean,
                       standardized = FALSE, # preguntar a Felipe
                       tAhead = x@tAhead)
    x@forecast <- res
    x@standardized <- FALSE
    return(x)
  }
  
  if(x@family == "gamma"){
    if(length(x@coef)==0) stop("The forecast method needs the coefficients of the iAR-Gamma model")
    res <- iARgforecast(coef = x@coef,
                        series = x@series,
                        mean = x@mean,
                        tAhead = x@tAhead)
    x@forecast <- res
    return(x)
  }
}
S7::method(forecast, CiAR) <- function(x) {
  if(length(x@series)==0) stop("The forecast method needs a time series")
  if(length(x@coef)==0) stop("The forecast method needs the coefficients of the CiAR model")
  res <- CiARforecast(coef = x@coef,
                      series = x@series,
                      times = x@times,
                      zero_mean = x@zero_mean,
                      standardized = x@standardized,
                      tAhead = x@tAhead)$forecast
  x@forecast <- res
  return(x)
}
S7::method(forecast, BiAR) <- function(x) {
  if(length(x@series)==0) stop("The forecast method needs a bivariate time series")
  if(length(x@coef)==0) stop("The forecast method needs the coefficients of the BiAR model")
  res <- t(BiARforecast(coef = x@coef,
                        series1 = x@series[,1],
                        series2 = x@series[,2],
                        times = x@times,
                        # zero_mean = x@zero_mean,
                        # standardized = x@standardized,
                        tAhead = x@tAhead)$forecast)
  x@forecast <- res
  return(x)
}