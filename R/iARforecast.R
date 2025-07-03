#' Forecast from iAR model
#'
#' Forecast from models fitted by \code{\link{iARloglik}}
#'
#' @param coef Estimated phi parameter by the iAR model.
#' @param series Array with the time series observations.
#' @param standardized logical; if TRUE, the array series is standardized; if FALSE, series contains the raw time series
#' @param zero_mean logical; if TRUE, the array series has zero mean; if FALSE, series has a mean different from zero.
#' @param tAhead The time ahead for forecast is required.
#'
#' @return Forecasted value from the iAR model
#' @references
#' \insertRef{Eyheramendy_2018}{iAR}
#'
#' @seealso
#' \code{\link{gentime}}, \code{\link{iARsample}}, \code{\link{iARloglik}}, \code{\link{iARkalman}}, \code{\link{iARfit}}
#' @keywords internal
#' @examples
#' \dontshow{
#' set.seed(6714)
#' st<-gentime(n=100)
#' y<-iARsample(coef=0.99,times=st)
#' y<-y$series
#' n=length(y)
#' p=trunc(n*0.99)
#' ytr=y[1:p]
#' yte=y[(p+1):n]
#' str=st[1:p]
#' ste=st[(p+1):n]
#' tahead=ste-str[p]
#' phi=iARloglik(times=ytr,series=str)$coef
#' forIAR=iARforecast(coef=phi,series=ytr,tAhead=tahead)
#' }
#' @noRd
iARforecast<-function(coef,series,standardized=TRUE,zero_mean=TRUE,tAhead)
{
  sigma = 1
  mu = 0
  if (zero_mean == FALSE)
  {
    mu = mean(series)
    series = series-mu
  }
  if (standardized == FALSE)
  {
    sigma = var(series)
    series=series/sqrt(sigma)
  }
  y1=series[length(series)]
  fit=(coef**(tAhead)*y1)
  fit=(fit*sqrt(sigma)+mu)
  return(fitted=fit)
}

