#' Simulate Time Series for Multiple iAR Models
#'
#' Simulates a time series for irregular autoregressive (iAR) models, including:
#' 1. Normal iAR model (`iAR`)
#' 2. T-distribution iAR model (`iAR-T`)
#' 3. Gamma-distribution iAR model (`iAR-Gamma`)
#'
#' @name sim
#'
#' @param x An object of class \code{iAR}, \code{CiAR}, or \code{BiAR}, containing the model specification and parameters:
#'   \itemize{
#'     \item For \code{iAR} (irregular AR models), the model family could be "norm", "t", or "gamma", where:
#'       \itemize{
#'         \item \code{family}: The distribution family of the iAR model (one of "norm", "t", or "gamma").
#'         \item \code{coef}: The coefficient(s) of the iAR model.
#'         \item \code{times}: A numeric vector specifying the time points of the series.
#'         \item \code{df}: Degrees of freedom for the t-distribution (only for \code{family = "t"}).
#'         \item \code{sigma}: The scale parameter for the t-distribution (only for \code{family = "t"}).
#'         \item \code{mean}: The mean parameter for the gamma distribution (only for \code{family = "gamma"}).
#'         \item \code{variance}: The variance parameter for the gamma distribution (only for \code{family = "gamma"}).
#'       }
#'     \item For \code{CiAR} (complex irregular autoregressive models):
#'       \itemize{
#'         \item \code{coef}: The real and imaginary parts of the CiAR model's coefficients.
#'         \item \code{times}: A numeric vector specifying the time points of the series.
#'         \item \code{rho}: The correlation parameter for the CiAR model.
#'         \item \code{c}: The scale parameter for the CiAR model.
#'       }
#'     \item For \code{BiAR} (BiAR models):
#'       \itemize{
#'         \item \code{coef}: The coefficients of the BiAR model (real and imaginary).
#'         \item \code{times}: A numeric vector specifying the time points of the series.
#'         \item \code{rho}: The correlation parameter for the BiAR model.
#'         \item \code{series_esd}: The series for the error structure (optional, used internally).
#'       }
#'   }
#' @param ... Additional arguments for generating irregular times. This is used only if no time points have been provided in the \code{times} argument:
#'   \describe{
#'   \item{n}{A positive integer. Length of observation times. Default is 100. This is used only if no time points are provided in \code{times}.}
#'   \item{distribution}{A character string specifying the distribution of the observation times. Default is `"expmixture"`. Available options are:
#'   \describe{
#'     \item{expmixture}{A mixture of two exponential distributions.}
#'     \item{uniform}{A uniform distribution.}
#'     \item{exponential}{A single exponential distribution.}
#'     \item{gamma}{A gamma distribution.}
#'     }
#'     }
#'     \item{lambda1}{Mean (1/rate) of the exponential distribution or the first exponential distribution in a mixture of exponential distributions. Default is `130`.}
#'     \item{lambda2}{Mean (1/rate) of the second exponential distribution in a mixture of exponential distributions. Default is `6.5`.}
#'     \item{p1}{Weight of the first exponential distribution in a mixture of exponential distributions. Default is `0.15`.}
#'     \item{p2}{Weight of the second exponential distribution in a mixture of exponential distributions. Default is `0.85`.}
#'     \item{a}{Shape parameter of a gamma distribution or lower limit of the uniform distribution. Default is `0`.}
#'     \item{b}{Scale parameter of a gamma distribution or upper limit of the uniform distribution. Default is `1`.}
#'   }
#'
#' 
#' @return An updated object of class \code{iAR}, \code{CiAR}, or \code{BiAR}, where the \code{series} property contains the simulated time series.
#'
#' @details
#' This function simulates time series based on the specified model and its parameters. Depending on the class of the input object:
#' \itemize{
#'   \item For \code{iAR} models, it supports three distribution families:
#'     \item "norm" for normal distribution.
#'     \item "t" for t-distribution.
#'     \item "gamma" for gamma distribution.
#'   \item For \code{CiAR} models, it uses complex autoregressive processes to generate the time series.
#'   \item For \code{BiAR} models, it simulates a BiAR process using specified coefficients and correlation.
#' }
#' The coefficients and any family-specific parameters must be set before calling this function.
#'
#' @references
#' \insertRef{Eyheramendy_2018}{iAR},\insertRef{Elorrieta_2019}{iAR},\insertRef{Elorrieta_2021}{iAR}
#'
#' @examples
#' # Example 1: Simulating a normal iAR model
#' library(iAR)
#' n=100
#' set.seed(6714)
#' o=iAR::utilities()
#' o<-gentime(o, n=n)
#' times=o@times
#' model_norm <- iAR(family = "norm", times = times, coef = 0.9,hessian=TRUE)
#' model_norm <- sim(model_norm)
#' plot(model_norm, type = "l", main = "Simulated iAR-Norm Series")
#'
#' # Example 2: Simulating a CiAR model
#' set.seed(6714)
#' model_CiAR <- CiAR(times = times,coef = c(0.9, 0))
#' model_CiAR <- sim(model_CiAR)
#' plot(model_CiAR , type = "l", main = "Simulated CiAR Series")
#'
#' # Example 3: Simulating a BiAR model
#' set.seed(6714)
#' model_BiAR <- BiAR(times = times,coef = c(0.9, 0.3), rho = 0.9)
#' model_BiAR <- sim(model_BiAR)
#' plot(times, model_BiAR@series[,1], type = "l", main = "Simulated BiAR Series")
#'
#' @export
sim <- S7::new_generic("sim", "x")
S7::method(generic = sim, signature = iAR) <- function(x, n = 100, distribution = "expmixture", lambda1 = 130, 
                                                       lambda2 = 6.5, p1 = 0.15, p2 = 0.85, a = 0, b = 1) {
    if (is.null(x@times) || length(x@times) == 0) {
      message("Generating irregular time points using the gentime method of the class utilities...")
      utils <- utilities()
      utils <- gentime(utils, n = n, distribution = distribution,lambda1 = lambda1, 
                       lambda2 = lambda2, p1 = p1, p2 = p2, a = a, b = b)
      x@times <- utils@times
    }
  if(x@family == "norm") {
    if(length(x@coef) == 0) stop("The sim method needs the coefficient of the iAR model")
    res <- iARsample(coef = x@coef, times = x@times)$series
    x@series <- res
    return(x)
  }
  if(x@family == "t") {
    if(length(x@coef) == 0) stop("The sim method needs the coefficient of the iAR-T model")
    res <- iARtsample(coef = x@coef, times = x@times, df = x@df, sigma = x@sigma)$series
    x@series <- res
    return(x)
  }
  if(x@family == "gamma") {
    if(length(x@coef) == 0) stop("The sim method needs the coefficient of the iAR-Gamma model")
    res <- iARgsample(coef = x@coef, times = x@times, mean = x@mean, sigma = x@variance)$series
    x@series <- res
    return(x)
  }
}
S7::method(generic = sim, signature = CiAR) <- function(x, rho = 0, c = 1, n = 100, distribution = "expmixture", lambda1 = 130, 
                                                        lambda2 = 6.5, p1 = 0.15, p2 = 0.85, a = 0, b = 1) {
  if (is.null(x@times) || length(x@times) == 0) {
    message("Generating irregular time points using the gentime method of the class utilities...")
    utils <- utilities()
    utils <- gentime(utils, n = n, distribution = distribution,lambda1 = lambda1, 
                     lambda2 = lambda2, p1 = p1, p2 = p2, a = a, b = b)
    x@times <- utils@times
  }
  if(length(x@coef) == 0) stop("The sim method needs the coefficients of the CiAR model")
  res <- CiARsample(phiR = x@coef[1], phiI = x@coef[2], times = x@times, rho = rho, c = c)$series
  x@series <- res
  return(x)
}
S7::method(generic = sim, signature = BiAR) <- function(x, n = 100, distribution = "expmixture", lambda1 = 130, 
                                                        lambda2 = 6.5, p1 = 0.15, p2 = 0.85, a = 0, b = 1) {
  if (is.null(x@times) || length(x@times) == 0) {
    message("Generating irregular time points using the gentime method of the class utilities...")
    utils <- utilities()
    utils <- gentime(utils, n = n, distribution = distribution,lambda1 = lambda1, 
                     lambda2 = lambda2, p1 = p1, p2 = p2, a = a, b = b)
    x@times <- utils@times
  }
  if(length(x@coef) == 0) stop("The sim method needs the coefficients of the BiAR model")
  no_series_esd <- is.integer(x@series_esd)
  
  if(no_series_esd) x@series_esd <- matrix(0, ncol = 2)
  res <- BiARsample(times = x@times, phiR = x@coef[1], phiI = x@coef[2], series_esd1 = x@series_esd[,1], 
                    series_esd2 = x@series_esd[,2], rho = x@rho)$series
  x@series <- cbind(res[1,],res[2,])
  
  if(no_series_esd) x@series_esd <- integer(0)
  
  return(x)
}