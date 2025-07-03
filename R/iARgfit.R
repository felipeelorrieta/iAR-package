#' Fitted Values of iAR-Gamma model
#'
#' Fit an iAR-Gamma model to an irregularly observed time series.
#'
#' @param coef Estimated phi parameter by the iAR-Gamma model.
#' @param mean Estimated mean parameter by the iAR-Gamma model.
#' @param series Array with the time series observations.
#' @param times Array with the irregular observational times.
#'
#' @return Fitted values of the iAR-Gamma model
#' @references
#' \insertRef{Eyheramendy_2018}{iAR}
#'
#' @seealso
#' \code{\link{gentime}}, \code{\link{iARgsample}}, \code{\link{iARgamma}}
#' @keywords internal
#' @examples
#' \dontshow{
#' n=300
#' set.seed(6714)
#' st<-gentime(n)
#' y<-iARgsample(coef=0.9,times=st,sigma=1,mean=1)
#' model<-iARgamma(series=y$series, times=st)
#' phi=model$coef
#' muest=model$mean
#' sigmaest=model$sigma
#' fit=iARgfit(coef=phi,mean=muest,series=y$series,times=st)
#' }
#' @noRd
iARgfit<-function(coef,mean,series,times)
{
  delta=diff(times)
  y1=series[-c(length(series))]
  xd=coef**delta
  yhat = mean + xd * y1
  fit=c(mean,yhat)
  return(fitted=fit)
}
