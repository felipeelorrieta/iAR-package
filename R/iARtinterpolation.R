#' Interpolation from iAR-T model
#'
#' Interpolation of missing values from models fitted by \code{\link{iARt}}
#'
#' @param coef A given array with the parameters of the iAR-T model. The first element of the array corresponding to the phi parameter and the second element to the scale parameter sigma
#' @param series Array with the time series observations.
#' @param times Array with the irregular observational times.
#' @param df degrees of freedom
#' @param yini a single value, initial value for the estimation of the missing value of the time series.
#'
#'
#' @return A list with the following components:
#' \itemize{
#' \item{fitted}{ Estimation of a missing value of the iAR-T process.}
#' \item{ll}{ Value of the negative log likelihood evaluated in the fitted missing values.}
#' }
#' @references
#' \insertRef{Eyheramendy_2018}{iAR}
#'
#' @seealso
#'
#' \code{\link{gentime}}, \code{\link{iARtsample}}, \code{\link{iARt}}
#'
#'
#' @keywords internal
#' @examples
#' \dontshow{
#' set.seed(6714)
#' n<-100
#' st<-gentime(n)
#' y<-iARtsample(coef=0.9,times=st,sigma=1,df=3)
#' model<-iARt(series=y$series, times=st)
#' napos=10
#' y0=y$series
#' y=y$series
#' y[napos]=NA
#' xest=c(model$coef,model$sigma)
#' yest=iARtinterpolation(coef=xest,series=y,times=st)
#' yest$fitted
#' mse=(y0[napos]-yest$fitted)^2
#' print(mse)
#' plot(st,y,type='l',xlim=c(st[napos-5],st[napos+5]))
#' points(st,y,pch=20)
#' points(st[napos],yest$fitted,col="red",pch=20)
#' }
#' @noRd
iARtinterpolation<-function(coef,series,times,df=3,yini=0)
{
  aux<-1e10
  value<-1e10
  br<-0
  countNA = length(which(is.na(series)))
  if(countNA==0){
    stop("This function requires a NA value at the point to be interpolated.")}
  if(countNA > 1){
    stop("This function requires only one NA value.")}
  if (yini==0)
    yini = rnorm(1)
  out = optim(yini, iARphit, lower = -Inf, upper = Inf,
              coef=coef,series=series,times=times,df=df, method="L-BFGS-B")
  par = out$par
  aux = out$value
  return(list(fitted = par, ll = aux))
}
