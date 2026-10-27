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
#'     \item \code{init}: Initial values for parameter estimation. For \code{family = "norm"}, a single numeric value (default 0.1). For \code{family = "t"}, a numeric vector of length two. For \code{family = "gamma"}, a numeric vector of length three.
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
S7::method(generic=loglik,signature=iAR) <- function(x, init=NULL) {
  if(length(x@series) == 0) stop("The loglik method needs a time series")
  if(any(is.na(x@series))) {
    warning("NA values were removed for parameter estimation.",call. = FALSE)
  }
  keep <- is.finite(x@series)
  series_aux <- x@series[keep]
  times_aux  <- x@times[keep]
  if(length(series_aux) < 3) stop("The Kalman method requires at least 3 observations")
  if(sd(series_aux, na.rm = TRUE) == 0)
    stop("The series is constant (zero variance). Likelihood estimation requires variability in the data.")
  if(length(x@series_esd) != 0) {
    series_esd_aux <- x@series_esd[keep]
  }
  else {
    series_esd_aux <- x@series_esd
  }
  if(x@family == "norm"){
    if(is.null(init))
      init=0.1
    res <- iARloglik(series = series_aux,
                     times = times_aux,
                     series_esd = series_esd_aux,
                     zero_mean = x@zero_mean,
                     standardized = x@standardized,
                     hessian = x@hessian,init=init)
    x@coef <- .clamp_iar_coef(res$coef)
    x@loglik <- res$loglik
    x@summary <- res$summary
    if(x@hessian == TRUE)
    {
      if (!isTRUE(tryCatch({ chol(res$summary$hess); TRUE }, error = function(e) FALSE))) {
        warning("Hessian matrix is not positive definite or not invertible. ",call. = FALSE)
      }
      if (any(!is.finite(res$summary$stderr))) {
        warning("Some standard errors are non-finite.",call. = FALSE)
      }
      if (res$summary$convergence != 0) {
        warning(paste0("Optimization may not have converged (code = ", res$summary$convergence, ")."), call. = FALSE)
      }
    }
    return(x)
  }
  if(x@family == "t"){
    if(is.null(init))
      init=c(0,1)
    res <- iARt(series = series_aux,
                times = times_aux,
                df = x@df,
                hessian = x@hessian, init=init)
    x@coef <- .clamp_iar_coef(res$coef)
    x@sigma <- res$sigma
    x@loglik <- res$ll
    x@summary <- res$summary
    if(x@hessian == TRUE)
    {
      if (!isTRUE(tryCatch({ chol(res$summary$hess); TRUE }, error = function(e) FALSE))) {
        warning("Hessian matrix is not positive definite or not invertible. ",call. = FALSE)
      }
      if (any(!is.finite(res$summary$stderr))) {
        warning("Some standard errors are non-finite.",call. = FALSE)
      }
      if (res$summary$convergence != 0) {
        warning(paste0("Optimization may not have converged (code = ", res$summary$convergence, ")."), call. = FALSE)
      }
    }
    return(x)
  }

  if(x@family == "gamma"){
    if(is.null(init))
      init=c(0,1,1)
    res <- iARgamma(series = series_aux,
                    times = times_aux,
                    hessian = x@hessian,init=init)
    x@coef <- .clamp_iar_coef(res$coef)
    x@mean <- res$mean
    x@variance <- res$sigma
    x@loglik <- res$ll
    x@summary <- res$summary
    if(x@hessian == TRUE)
    {
      if (!isTRUE(tryCatch({ chol(res$summary$hess); TRUE }, error = function(e) FALSE))) {
        warning("Hessian matrix is not positive definite or not invertible. ",call. = FALSE)
      }
      if (any(!is.finite(res$summary$stderr))) {
        warning("Some standard errors are non-finite.",call. = FALSE)
      }
      if (res$summary$convergence != 0) {
        warning(paste0("Optimization may not have converged (code = ", res$summary$convergence, ")."), call. = FALSE)
      }
    }
    return(x)
  }
}
