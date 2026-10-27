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
#'         \item \code{init}: A numeric value specifying the initial value for parameter estimation. Default is 0.1.
#'       }
#'     \item For \code{CiAR} models, the object should contain:
#'       \itemize{
#'         \item \code{series}: A numeric vector of the time series data.
#'         \item \code{times}: A numeric vector specifying the time points of the series.
#'         \item \code{series_esd}: A numeric vector specifying the error structure.
#'         \item \code{zero_mean}: A logical value indicating if the series should be zero-centered.
#'         \item \code{standardized}: A logical value indicating if the series should be standardized.
#'         \item \code{c}: A numeric value specifying the scale parameter for the CiAR model.
#'         \item \code{init}: A numeric vector of length two specifying the initial values for parameter estimation. When provided, these values are included as an additional starting point.
#'       }
#'     \item For \code{BiAR} models, the object should contain:
#'       \itemize{
#'         \item \code{series}: A numeric matrix containing two columns for bivariate time series data.
#'         \item \code{times}: A numeric vector specifying the time points of the series.
#'         \item \code{series_esd}: A numeric matrix specifying the error structure for both series.
#'         \item \code{zero_mean}: A logical value indicating if the series should be zero-centered.
#'         \item \code{init}: A numeric vector of length two specifying the initial values for parameter estimation. When provided, these values are included as an additional starting point.
#'       }
#'   }
#'
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
#' times=gentime(n=n)@times
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
S7::method(kalman, iAR) <- function(x, init = 0.1) {
  if(length(x@series) == 0) stop("The Kalman method needs a time series")
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
  if(x@family == "norm") {
    res <- iARkalman(series = series_aux,
                     times = times_aux,
                     series_esd = series_esd_aux,
                     zero_mean = x@zero_mean,
                     standardized = x@standardized,
                     hessian = x@hessian,init=init)
    x@coef <- .clamp_iar_coef(res$coef)
    x@kalmanlik <- res$kalman
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

S7::method(kalman, CiAR) <- function(x, c = 1, init = rep(0,2)) {
  if(length(x@series) == 0) stop("The Kalman method needs a time series")
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
  res <- CiARkalman(series = series_aux,
                    times = times_aux,
                    series_esd = series_esd_aux,
                    zero_mean = x@zero_mean,
                    standardized = x@standardized,
                    hessian = x@hessian,
                    c = c,init=init)  
  x@coef <- c(res$phiR, res$phiI)
  if (sqrt(sum(x@coef^2)) > 0.995)
    warning(paste0(
      "Estimated modulus sqrt(phiR^2 + phiI^2) = ", round(sqrt(sum(x@coef^2)), 4),
      " is near 1. The model may be near non-stationary; results should be interpreted with caution."
    ), call. = FALSE)
  x@kalmanlik <- res$ll
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

S7::method(kalman, BiAR) <- function(x, init = rep(0,2)) {
  if(length(x@series) == 0) stop("The Kalman method needs a bivariate time series")
  if(any(!is.finite(x@series))) {
    warning("Rows with NA values were removed for parameter estimation.", call. = FALSE)
  }
  keep <- apply(is.finite(x@series), 1, all)
  series_aux <- x@series[keep, , drop = FALSE]
  times_aux  <- x@times[keep]
  if(nrow(series_aux) < 3) stop("The Kalman method requires at least 3 observations")
  if(any(apply(series_aux, 2, sd, na.rm = TRUE) == 0))
    stop("One or both series are constant (zero variance). Likelihood estimation requires variability in the data.")

  no_series_esd <- length(x@series_esd) == 0
  if(no_series_esd) {
    esd1 <- 0
    esd2 <- 0
  } else {
    series_esd_aux <- x@series_esd[keep, , drop = FALSE]
    esd1 <- series_esd_aux[, 1]
    esd2 <- series_esd_aux[, 2]
  }
  res <- BiARkalman(series1 = series_aux[, 1],
                    series2 = series_aux[, 2],
                    times = times_aux,
                    series_esd1 = esd1,
                    series_esd2 = esd2,
                    zero_mean = x@zero_mean,
                    hessian = x@hessian,init=init)
  x@coef <- c(res$phiR, res$phiI)
  if (sqrt(sum(x@coef^2)) > 0.995)
    warning(paste0(
      "Estimated modulus sqrt(phiR^2 + phiI^2) = ", round(sqrt(sum(x@coef^2)), 4),
      " is near 1. The model may be near non-stationary; results should be interpreted with caution."
    ), call. = FALSE)
  x@kalmanlik <- res$ll
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
