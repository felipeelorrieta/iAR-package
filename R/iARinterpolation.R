#' Interpolation from iAR model
#'
#' Interpolation of missing values from models fitted by \code{\link{iARkalman}}
#'
#' @param coef A given phi coefficient of the iAR model.
#' @param series Array with the time series observations.
#' @param times Array with the irregular observational times.
#' @param series_esd Array with the measurements error standard deviations.
#' @param yini a single value, initial value for the estimation of the missing value of the time series.
#' @param zero_mean logical; if TRUE, the array y has zero mean; if FALSE, y has a mean different from zero.
#' @param standardized logical; if TRUE, the array y is standardized; if FALSE, y contains the raw time series.
#'
#'
#' @return A list with the following components:
#' \itemize{
#' \item{fitted}{ Estimation of a missing value of the iAR process.}
#' \item{ll}{ Value of the negative log likelihood evaluated in the fitted missing values.}
#' }
#' @references
#' \insertRef{Eyheramendy_2018}{iAR}
#'
#' @seealso
#'
#' \code{\link{gentime}}, \code{\link{iARsample}}, \code{\link{iARkalman}}
#'
#'
#' @keywords internal
#' @examples
#' \dontshow{
#' set.seed(6714)
#' st<-gentime(n=100)
#' y<-iARsample(coef=0.99,times=st)
#' y<-y$series
#' phi=iARkalman(series=y,times=st)$coef
#' print(phi)
#' napos=10
#' y0=y
#' y[napos]=NA
#' xest=phi
#' yest=iARinterpolation(coef=xest,series=y,times=st)
#' yest$fitted
#' mse=(y0[napos]-yest$fitted)^2
#' print(mse)
#' plot(st,y,type='l',xlim=c(st[napos-5],st[napos+5]))
#' points(st,y,pch=20)
#' points(st[napos],yest$fitted,col="red",pch=20)
#' }
#' @noRd
iARinterpolation<-function(coef,series,times,series_esd=0,yini=0,zero_mean=TRUE,standardized=TRUE)
{
  aux<-1e10
  value<-1e10
  br<-0
  countNA = length(which(is.na(series)))
  if(countNA==0){
    stop("This function requires a NA value at the point to be interpolated.")}
  if(countNA > 1){
    stop("This function requires only one NA value.")}
  if(sum(series_esd)==0){
    series_esd=rep(0,length(series))}
  if (yini==0)
    yini = rnorm(1)
  out = optim(c(yini), iARphikalman, lower = -Inf, upper = Inf,
              coef=coef,series = series, times = times, series_esd=series_esd,zero_mean = zero_mean,standardized=standardized, method="BFGS")
  par = out$par
  aux = out$value
  return(list(fitted = par, ll = aux))
}
