#' Interpolation from iAR-Gamma model
#'
#' Interpolation of missing values from models fitted by \code{\link{iARgamma}}
#'
#' @param coef A given array with the parameters of the iAR-Gamma model. The first element of the array corresponding to the phi parameter, the second to the mean parameter, and the last one to the scale parameter sigma.
#' @param series Array with the time series observations.
#' @param times Array with the irregular observational times.
#' @param yini a single value, initial value for the estimation of the missing value of the time series.
#'
#'
#' @return A list with the following components:
#' \itemize{
#' \item{fitted}{ Estimation of a missing value of the iAR-Gamma process.}
#' \item{ll}{ Value of the negative log likelihood evaluated in the fitted missing values.}
#' }
#' @references
#' \insertRef{Eyheramendy_2018}{iAR}
#'
#' @seealso
#'
#' \code{\link{gentime}}, \code{\link{iARgsample}}, \code{\link{iARgamma}}
#'
#'
#' @keywords internal
#' @examples
#' \dontshow{
#' set.seed(6714)
#' n<-100
#' st<-gentime(n)
#' y<-iARgsample(coef=0.9,times=st,sigma=1,mean=1)
#' model<-iARgamma(series=y$series, times=st)
#' y<-y$series
#' napos=10
#' y0=y
#' y[napos]=NA
#' xest=c(model$coef,model$mean,model$sigma)
#' yest=iARginterpolation(coef=xest,series=y,times=st)
#' yest$fitted
#' mse=(y0[napos]-yest$fitted)^2
#' print(mse)
#' plot(st,y,type='l',xlim=c(st[napos-5],st[napos+5]))
#' points(st,y,pch=20)
#' points(st[napos],yest$fitted,col="red",pch=20)
#' }
#' @noRd
iARginterpolation<-function(coef,series,times,yini=1)
{
  aux<-1e10
  value<-1e10
  br<-0
  countNA = length(which(is.na(series)))
  if(countNA==0){
    stop("This function requires a NA value at the point to be interpolated.")}
  if(countNA > 1){
    stop("This function requires only one NA value.")}
  if (yini==1)
    yini = rgamma(1,1,1)
  out = optim(yini, iARphigamma, lower = 0.0001, upper = Inf,
              coef=coef,series = series, times = times, method="L-BFGS-B")
  par = out$par
  aux = out$value
  return(list(fitted = par, ll = aux))
}
