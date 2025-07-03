#' Fitted Values of iAR model
#'
#' Fit an iAR model to an irregularly observed time series.
#'
#' @param coef Estimated phi parameter by the iAR model.
#' @param series Array with the time series observations.
#' @param times Array with the irregular observational times.
#' @param standardized logical; if TRUE, the array series is standardized; if FALSE, series contains the raw time series
#' @param zero_mean logical; if TRUE, the array series has zero mean; if FALSE, series has a mean different from zero.
#'
#' @return Fitted values of the iAR model
#' @references
#' \insertRef{Eyheramendy_2018}{iAR}
#'
#' @seealso
#' \code{\link{gentime}}, \code{\link{iARsample}}, \code{\link{iARloglik}}, \code{\link{iARkalman}}
#' @keywords internal
#' @examples
#' \dontshow{
#' set.seed(6714)
#' st<-gentime(n=100)
#' y<-iARsample(coef=0.99,times=st)
#' y<-y$series
#' phi=iARloglik(series=y,times=st)$coef
#' fit=iARfit(coef=phi,series=y,times=st)
#' }
#' @noRd
iARfit<-function(coef,series,times,standardized=TRUE,zero_mean=TRUE)
{
  sigma = 1
  mu = 0
  if (standardized == FALSE)
    sigma = var(series)
  if (zero_mean == TRUE)
    mu = mean(series)
  if (zero_mean == TRUE)
    series = series-mu
  if (standardized == FALSE)
    series=series/sqrt(sigma)
  delta=diff(times)
  y1=series[-c(length(series))]
  fit=c(0,(coef**delta)*y1)
  fit=(fit*sqrt(sigma)+mu)
  return(fitted=fit)
}

