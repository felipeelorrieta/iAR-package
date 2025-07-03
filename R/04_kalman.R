#' Maximum Likelihood Estimation of Parameters for iAR, CiAR, and BiAR Models using the Kalman Filter
#'
#' Performs Maximum Likelihood Estimation (MLE) of the parameters of the iAR, CiAR, and BiAR models by maximizing the likelihood function using the Kalman filter.
#' This method applies the Kalman filter to compute the likelihood and estimate the model parameters that maximize the likelihood for each model type.
#' 
#' @name kalman
#' 
#' @param x An object of class \code{iAR}, \code{CiAR}, or \code{BiAR} containing the model specification and parameters:
#'   \itemize{
#'     \item For \code{iAR} models, the object should contain:
#'       \itemize{
#'         \item \code{series}: A numeric vector of the time series data.
#'         \item \code{times}: A numeric vector specifying the time points of the series.
#'         \item \code{series_esd}: A numeric vector specifying the error structure.
#'         \item \code{zero_mean}: A logical value indicating if the series should be zero-centered.
#'         \item \code{standardized}: A logical value indicating if the series should be standardized.
#'       }
#'     \item For \code{CiAR} models, the object should contain:
#'       \itemize{
#'         \item \code{series}: A numeric vector of the time series data.
#'         \item \code{times}: A numeric vector specifying the time points of the series.
#'         \item \code{series_esd}: A numeric vector specifying the error structure.
#'         \item \code{zero_mean}: A logical value indicating if the series should be zero-centered.
#'         \item \code{standardized}: A logical value indicating if the series should be standardized.
#'         \item \code{c}: A numeric value specifying the scale parameter for the CiAR model.
#'         \item \code{niter}: An integer specifying the number of iterations for the Kalman filter.
#'         \item \code{seed}: An integer seed for random number generation (optional).
#'       }
#'     \item For \code{BiAR} models, the object should contain:
#'       \itemize{
#'         \item \code{series}: A numeric matrix containing two columns for bivariate time series data.
#'         \item \code{times}: A numeric vector specifying the time points of the series.
#'         \item \code{series_esd}: A numeric matrix specifying the error structure for both series.
#'         \item \code{zero_mean}: A logical value indicating if the series should be zero-centered.
#'         \item \code{niter}: An integer specifying the number of iterations for the Kalman filter.
#'         \item \code{seed}: An integer seed for random number generation (optional).
#'       }
#'   }
#'
#' @param ... Additional arguments (unused).
#'
#' @return An updated object of class \code{iAR}, \code{CiAR}, or \code{BiAR}, where the \code{coef} property is updated with the estimated model parameters (using MLE) and the \code{kalmanlik} property contains the log-likelihood value of the model.
#'
#' @details
#' This function applies the Kalman filter to perform Maximum Likelihood Estimation (MLE) of the parameters of the autoregressive models (iAR, CiAR, BiAR). The Kalman filter is used to maximize the likelihood function based on the given time series data, and the parameters that maximize the likelihood are estimated.
#' 
#' - For \code{iAR}, the Kalman filter is applied to estimate the model parameters by maximizing the likelihood.
#' - For \code{CiAR}, the Kalman filter is applied to estimate the parameters of the complex autoregressive model by maximizing the likelihood.
#' - For \code{BiAR}, the Kalman filter is applied to estimate the parameters of the bivariate autoregressive model by maximizing the likelihood.
#'
#' The method returns the updated model object, including the estimated parameters and the log-likelihood value.
#'
#' @references
#' \insertRef{Eyheramendy_2018}{iAR},\insertRef{Elorrieta_2019}{iAR},\insertRef{Elorrieta_2021}{iAR}
#'
#' @examples
#' # Example 1: Applying Kalman filter for MLE of iAR model parameters
#' library(iAR)
#' n=100
#' set.seed(6714)
#' o=iAR::utilities()
#' o<-gentime(o, n=n)
#' times=o@times
#' model_norm <- iAR(family = "norm", times = times, coef = 0.9,hessian=TRUE)
#' model_norm <- sim(model_norm)
#' model_norm <- kalman(model_norm)  
#' print(model_norm@coef)  # Access the estimated coefficients
#' print(model_norm@kalmanlik)  # Access the Kalman likelihood value
#' 
#' # Example 2: Applying Kalman filter for MLE of CiAR model parameters
#' set.seed(6714)
#' model_CiAR <- CiAR(times = times,coef = c(0.9, 0))
#' model_CiAR <- sim(model_CiAR)
#' y=model_CiAR@series
#' y1=y/sd(y)
#' model_CiAR@series=y1
#' model_CiAR@series_esd=rep(0,n)
#' model_CiAR <- kalman(model_CiAR)
#' print(model_CiAR@coef)
#' print(model_CiAR@kalmanlik)  
#' 
#' # Example 3: Applying Kalman filter for MLE of BiAR model parameters
#' set.seed(6714)
#' model_BiAR <- BiAR(times = times,coef = c(0.9, 0.3), rho = 0.9)
#' model_BiAR <- sim(model_BiAR)
#' y=model_BiAR@series
#' y1=y/apply(y,2,sd)
#' model_BiAR@series=y1
#' model_BiAR@series_esd=matrix(0,n,2)
#' model_BiAR <- kalman(model_BiAR)
#' print(model_BiAR@coef) 
#' print(model_BiAR@kalmanlik)  
#' 
#' @export
kalman <- S7::new_generic("kalman", "x")
S7::method(kalman, iAR) <- function(x) {
  if(length(x@series) == 0) stop("The Kalman method needs a time series")
  
  if(x@family == "norm") {
    res <- iARkalman(series = x@series,
                     times = x@times,
                     series_esd = x@series_esd,
                     zero_mean = x@zero_mean,
                     standardized = x@standardized)
    x@coef <- res$coef
    x@kalmanlik <- res$kalman
    return(x)
  }
}

S7::method(kalman, CiAR) <- function(x, c = 1, niter = 10, seed = NULL) {
  if(length(x@series) == 0) stop("The Kalman method needs a time series")
  
  res <- CiARkalman(series = x@series,
                    times = x@times,
                    series_esd = x@series_esd,
                    zero_mean = x@zero_mean,
                    standardized = x@standardized,
                    c = c,
                    niter = niter,
                    seed = seed)
  x@coef <- c(res$phiR, res$phiI)
  x@kalmanlik <- res$ll
  return(x)
}

S7::method(kalman, BiAR) <- function(x, niter = 10, seed = NULL) {
  if(length(x@series) == 0) stop("The Kalman method needs a bivariate time series")
  
  no_series_esd <- is.integer(x@series_esd)
  
  if(no_series_esd) x@series_esd <- matrix(0, ncol = 2)
  
  res <- BiARkalman(series1 = x@series[, 1],
                    series2 = x@series[, 2],
                    times = x@times,
                    series_esd1 = x@series_esd[, 1],
                    series_esd2 = x@series_esd[, 2],
                    zero_mean = x@zero_mean,
                    niter = niter,
                    seed = seed)
  x@coef <- c(res$phiR, res$phiI)
  x@kalmanlik <- res$ll
  
  if(no_series_esd) x@series_esd <- integer(0)
  
  return(x)
}
 