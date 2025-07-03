#' Maximum Likelihood Estimation for iAR Models
#'
#' Maximum Likelihood Estimation for irregular autoregressive (iAR) models, supporting different distribution families:
#' normal (`iAR`), t (`iAR-T`), and gamma (`iAR-Gamma`).
#' 
#' @name loglik
#'
#' @param x An object of class \code{iAR}, containing the model specification, parameters, and the time series to be evaluated:
#'   \itemize{
#'     \item \code{series}: The observed time series.
#'     \item \code{times}: A numeric vector specifying the time points of the series.
#'     \item \code{series_esd}: (Optional) A standardized version of the series.
#'     \item \code{zero_mean}: Logical. Indicates whether the series should be mean-centered.
#'     \item \code{standardized}: Logical. Indicates whether the series is standardized.
#'     \item \code{hessian}: Logical. If \code{TRUE}, the function computes the Hessian matrix for parameter estimation.
#'     \item \code{family}: The distribution family of the iAR model (one of "norm", "t", or "gamma").
#'     \item \code{df}: Degrees of freedom for the t-distribution (only for \code{family = "t"}).
#'     \item \code{sigma}: The scale parameter for the t-distribution (only for \code{family = "t"}).
#'     \item \code{mean}: The mean parameter for the gamma distribution (only for \code{family = "gamma"}).
#'     \item \code{variance}: The variance parameter for the gamma distribution (only for \code{family = "gamma"}).
#'   }
#'
#' @return An updated \code{iAR} object with the following additional attributes:
#'   \itemize{
#'     \item \code{coef}: Estimated model coefficients.
#'     \item \code{loglik}: Log-likelihood value of the model.
#'     \item \code{summary}: A summary table containing parameter estimates, standard errors, and p-values.
#'     \item \code{sigma}: For t and gamma families, the estimated scale parameter.
#'     \item \code{mean}: For the gamma family, the estimated mean parameter.
#'     \item \code{variance}: For the gamma family, the estimated variance parameter.
#'   }
#'
#' @param ... Additional arguments (unused).
#'
#' @details
#' This method estimates the parameters of an iAR model using the Maximum Likelihood Estimation (MLE) approach. 
#' Depending on the chosen distribution family, the corresponding likelihood function is maximized:
#' \itemize{
#'   \item "norm" maximizes the likelihood for a normally-distributed series.
#'   \item "t" maximizes the likelihood for a t-distributed series.
#'   \item "gamma" maximizes the likelihood for a gamma-distributed series.
#' }
#' The function updates the \code{iAR} object with the estimated parameters, the log-likelihood value, and a summary 
#' table that includes standard errors and p-values.
#'
#' @references
#' \insertRef{Eyheramendy_2018}{iAR}
#'
#' @examples
#' # Example: Estimating parameters for a normal iAR model
#' library(iAR)
#' times <- 1:100
#' model <- iAR(family = "norm", times = times, coef = 0.9, hessian = TRUE)
#' model <- sim(model)  # Simulate the series
#' model <- loglik(model)  # Estimate parameters using MLE
#' print(model@coef)  # Access the estimated coefficients
#' print(model@loglik)  # Access the computed log-likelihood
#'
#' @export
loglik <- S7::new_generic("loglik","x")
S7::method(generic=loglik,signature=iAR) <- function(x) {
  if(length(x@series) == 0) stop("The loglik method needs a time series")
  if(x@family == "norm"){
    res <- iARloglik(series = x@series,
                     times = x@times,
                     series_esd = x@series_esd,
                     zero_mean = x@zero_mean,
                     standardized = x@standardized,
                     hessian = x@hessian)
    x@coef <- res$coef
    x@loglik <- res$loglik
    x@summary <- res$summary
    return(x)
  }
  
  if(x@family == "t"){
    res <- iARt(series = x@series,
                times = x@times,
                df = x@df,
                hessian = x@hessian)
    x@coef <- res$coef
    x@sigma <- res$sigma
    x@loglik <- res$ll
    x@summary <- res$summary
    return(x)
  }
  
  if(x@family == "gamma"){
    res <- iARgamma(series = x@series,
                    times = x@times,
                    hessian = x@hessian)
    x@coef <- res$coef
    x@mean <- res$mean
    x@variance <- res$sigma
    x@loglik <- res$ll
    x@summary <- res$summary
    return(x)
  }
}
