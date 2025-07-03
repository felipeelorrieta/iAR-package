#' Fitted Values for iAR, CiAR, and BiAR Classes
#'
#' Fitted Values for the provided data. This method is implemented for:
#' 1. Irregular Autoregressive models (`iAR`)
#' 2. Complex Irregular Autoregressive models (`CiAR`)
#' 3. Bivariate Autoregressive models (`BiAR`)
#' 
#' @name fit
#'
#' @param x An object of class \code{iAR}, \code{CiAR}, or \code{BiAR}, containing the model specification and parameters:
#'   \itemize{
#'     \item For \code{iAR}:
#'       \itemize{
#'         \item \code{family}: The distribution family of the iAR model (one of "norm", "t", or "gamma").
#'         \item \code{series}: A numeric vector representing the time series to be fitted.
#'         \item \code{coef}: The coefficient(s) of the iAR model.
#'         \item \code{times}: A numeric vector specifying the time points of the series.
#'         \item \code{zero_mean}: Logical, whether to fit a zero-mean model.
#'         \item \code{standardized}: Logical, whether the model output should be standardized (for "norm" family).
#'         \item \code{mean}: The mean parameter (only for "gamma" family).
#'       }
#'     \item For \code{CiAR}:
#'       \itemize{
#'         \item \code{coef}: The real and imaginary parts of the CiAR model's coefficients.
#'         \item \code{series}: A numeric vector representing the time series to be fitted.
#'         \item \code{times}: A numeric vector specifying the time points of the series.
#'         \item \code{zero_mean}: Logical, whether to fit a zero-mean model.
#'         \item \code{standardized}: Logical, whether the model output should be standardized.
#'         \item \code{c}: A scaling parameter for the CiAR model.
#'       }
#'     \item For \code{BiAR}:
#'       \itemize{
#'         \item \code{coef}: The coefficients of the BiAR model (real and imaginary parts).
#'         \item \code{series}: A numeric matrix with two columns representing the bivariate time series to be fitted.
#'         \item \code{times}: A numeric vector specifying the time points of the series.
#'         \item \code{series_esd}: A numeric matrix for the error structure (optional, used internally).
#'         \item \code{zero_mean}: Logical, whether to fit a zero-mean model.
#'       }
#'   }
#'
#' @param ... Additional arguments (unused).
#'
#' @return An updated object of class \code{iAR}, \code{CiAR}, or \code{BiAR}, where the \code{fitted_values} property contains the fitted time series values.
#'
#' @details
#' This method fits the specified time series model to the data contained in the object. Depending on the class of the input object:
#' \itemize{
#'   \item For \code{iAR}, the function supports three distribution families:
#'     \item "norm" for normal distribution.
#'     \item "t" for t-distribution.
#'     \item "gamma" for gamma distribution.
#'   \item For \code{CiAR}, the function uses complex autoregressive processes.
#'   \item For \code{BiAR}, the function fits a bivariate autoregressive process.
#' }
#' All required parameters (e.g., coefficients, time points) must be set before calling this method.
#'
#' @references
#' \insertRef{Eyheramendy_2018}{iAR},\insertRef{Elorrieta_2019}{iAR},\insertRef{Elorrieta_2021}{iAR}
#'
#' @examples
#' # Example 1: Fitting a normal iAR model
#' library(iAR)
#' n=100
#' set.seed(6714)
#' o=iAR::utilities()
#' o<-gentime(o, n=n)
#' times=o@times
#' model_norm <- iAR(family = "norm", times = times, coef = 0.9)  
#' model_norm <- sim(model_norm)
#' model_norm <- kalman(model_norm) 
#' model_norm <- fit(model_norm)
#' plot(model_norm@times, model_norm@series, type = "l", main = "Original Series")
#' lines(model_norm@times, model_norm@fitted_values, col = "red", lwd = 2)
#' plot_fit(model_norm)
#'
#' # Example 2: Fitting a CiAR model
#' set.seed(6714)
#' model_CiAR <- CiAR(times = times,coef = c(0.9, 0))
#' model_CiAR <- sim(model_CiAR)
#' y=model_CiAR@series
#' y1=y/sd(y)
#' model_CiAR@series=y1
#' model_CiAR@series_esd=rep(0,n)
#' model_CiAR <- kalman(model_CiAR)
#' print(model_CiAR@coef)
#' model_CiAR <- fit(model_CiAR)
#' yhat=model_CiAR@fitted_values
#'
#' # Example 3: Fitting a BiAR model
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
#' model_BiAR <- fit(model_BiAR)
#' print(model_BiAR@rho)
#' yhat=model_BiAR@fitted_values
#'
#' @export
fit <- S7::new_generic("fit", "x")
S7::method(generic = fit, signature = iAR) <- function(x) {
  if(length(x@series)==0) stop("The fit method needs a time series")
  if(x@family == "norm"){
    if(length(x@coef)==0) stop("The fit method needs the coefficient of the iAR model")
    res <- iARfit(coef = x@coef,
                  series = x@series,
                  times = x@times,
                  zero_mean = x@zero_mean,
                  standardized = x@standardized)
    x@fitted_values <- res
    return(x)
  }
  if(x@family == "t"){
    if(length(x@coef)==0) stop("The fit method needs the coefficient of the iAR-T model")
    res <- iARfit(coef = x@coef,
                  series = x@series,
                  times = x@times,
                  zero_mean = x@zero_mean,
                  standardized = FALSE) # preguntar a felipe
    x@fitted_values <- res
    x@standardized <- FALSE
    return(x)
  }
  if(x@family == "gamma"){
    if(length(x@coef)==0) stop("The fit method needs the coefficients of the iAR-Gamma model")
    res <- iARgfit(coef = x@coef,
                   series = x@series,
                   times = x@times,
                   mean = x@mean)
    x@fitted_values <- res
    return(x)
  }
}
S7::method(generic = fit, signature = CiAR)  <- function(x, c = 1) {
  if(length(x@series)==0) stop("The fit method needs a time series")
  if(length(x@coef)==0) stop("The fit method needs the coefficients of the CiAR model")
  res <- CiARfit(coef = x@coef,
                 series = x@series,
                 times = x@times,
                 zero_mean = x@zero_mean,
                 standardized = x@standardized,
                 c = c)$fitted
  x@fitted_values <- res
  return(x)
}
S7::method(generic = fit, signature = BiAR) <- function(x) {
  if(length(x@series)==0) stop("The fit method needs a bivariate time series")
  if(length(x@coef)==0) stop("The fit method needs the coefficients of the BiAR model")
  no_series_esd <- is.integer(x@series_esd)
  if(no_series_esd) x@series_esd <- matrix(0, ncol = 2)
  res <- BiARfit(coef = x@coef,
                   series1 = x@series[,1],
                   series2 = x@series[,2],
                   times = x@times,
                   series_esd1 = x@series_esd[,1],
                   series_esd2 = x@series_esd[,2],
                   zero_mean = x@zero_mean
                   # standardized = x@standardized,
                   )
  x@fitted_values <- t(res$fitted)
  x@rho <- res$rho
  if(no_series_esd) x@series_esd <- integer(0)
  return(x)
}
