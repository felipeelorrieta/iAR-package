#' Interpolation from CiAR model
#'
#' Interpolation of missing values from models fitted by \code{\link{CiARkalman}}
#'
#' @param coef An array with the parameters of the CiAR model. The elements of the array are, in order, the real (phiR) and the imaginary (phiI) part of the coefficient of CiAR model.
#' @param series Array with the time series observations.
#' @param times Array with the irregular observational times.
#' @param series_esd Array with the measurements error standard deviations.
#' @param yini a single value, initial value for the estimation of the missing value of the time series.
#' @param zero_mean logical; if TRUE, the array series has zero mean; if FALSE, series has a mean different from zero.
#' @param standardized logical; if TRUE, the array series is standardized; if FALSE, series contains the raw time series.
#' @param c Nuisance parameter corresponding to the variance of the imaginary part.
#' @param seed a single value, interpreted as the seed of the random process.
#'
#' @return A list with the following components:
#' \itemize{
#' \item{fitted}{ Estimation of a missing value of the CiAR process.}
#' \item{ll}{ Value of the negative log likelihood evaluated in the fitted missing values.}
#' }
#' @references
#' \insertRef{Elorrieta_2019}{iAR}
#'
#' @seealso
#'
#' \code{\link{gentime}}, \code{\link{CiARsample}}, \code{\link{CiARkalman}}
#'
#'
#' @keywords internal
#' @examples
#' \dontshow{
#' n=100
#' set.seed(6714)
#' st<-gentime(n)
#' x=CiARsample(phiR=0.9,phiI=0,times=st,c=1)
#' y=x$series
#' y1=y/sd(y)
#' ciar=CiARkalman(series=y1,times=st)
#' ciar
#' napos=10
#' y0=y1
#' y1[napos]=NA
#' xest=c(ciar$phiR,ciar$phiI)
#' yest=CiARinterpolation(coef=xest,series=y1,times=st)
#' yest$fitted
#' mse=(y0[napos]-yest$fitted)^2
#' print(mse)
#' plot(st,y,type='l',xlim=c(st[napos-5],st[napos+5]))
#' points(st,y,pch=20)
#' points(st[napos],yest$fitted*sd(y),col="red",pch=20)
#' }
#' @noRd
CiARinterpolation<-function(coef,series,times,series_esd=0,yini=0,zero_mean=TRUE,standardized=TRUE,c=1,seed=NULL)
{
  set.seed(seed)
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
  out = optim(c(yini), CiARphikalman, lower = -Inf, upper = Inf,
              coef=coef, series=series, times=times,series_esd=series_esd, zero_mean=zero_mean, standardized=standardized,
              c=c, method="BFGS")
  par = out$par
  aux = out$value
  return(list(fitted = par, ll = aux))
}
