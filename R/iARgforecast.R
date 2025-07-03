#' Forecast from iAR-Gamma model
#'
#' Forecast from models fitted by \code{\link{iARgamma}}
#'
#' @param coef Estimated phi parameter by the iAR-Gamma model.
#' @param mean Estimated mean parameter by the iAR-Gamma model.
#' @param series Array with the time series observations.
#' @param tAhead The time ahead for forecast is required.
#'
#' @return Forecasted value from the iAR-Gamma model
#' @references
#' \insertRef{Eyheramendy_2018}{iAR}
#'
#' @seealso
#' \code{\link{gentime}}, \code{\link{iARgsample}}, \code{\link{iARgamma}}, \code{\link{iARgfit}}
#' @keywords internal
#' @examples
#' \dontshow{
#' n=100
#' set.seed(6714)
#' st<-gentime(n)
#' y<-iARgsample(coef=0.9,times=st,sigma=1,mean=1)
#' y<-y$series
#' n=length(y)
#' p=trunc(n*0.99)
#' ytr=y[1:p]
#' yte=y[(p+1):n]
#' str=st[1:p]
#' ste=st[(p+1):n]
#' tahead=ste-str[p]
#' model<-iARgamma(series=ytr, times=str)
#' phi=model$coef
#' muest=model$mean
#' sigmaest=model$sigma
#' fit=iARgforecast(coef=phi,mean=muest,series=ytr,tAhead=tahead)
#' }
#' @noRd
iARgforecast<-function(coef,mean,series,tAhead)
{
  y1=series[c(length(series))]
  xd=coef**(tAhead)
  yhat = mean + xd * y1
  fit=c(yhat)
  return(fitted=fit)
}
