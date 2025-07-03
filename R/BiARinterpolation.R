#' Interpolation from BiAR model
#'
#' Interpolation of missing values from models fitted by \code{\link{BiARkalman}}
#'
#' @param coef An array with the parameters of the BiAR model. The elements of the array are, in order, the real (phiR) and the imaginary (phiI) part of the coefficient of BiAR model.
#' @param series1 Array with the observations of the first time series of the BiAR process.
#' @param series2 Array with the observations of the second time series of the BiAR process.
#' @param times Array with the irregular observational times.
#' @param series_esd1 Array with the measurements error standard deviations of the first time series of the BiAR process.
#' @param series_esd2 Array with the measurements error standard deviations of the second time series of the BiAR process.
#' @param yini1 a single value, initial value of the estimation of the missing value of the first time series of the BiAR process.
#' @param yini2 a single value, initial value of the estimation of the missing value of the second time series of the BiAR process.
#' @param zero_mean logical; if TRUE, the array y has zero mean; if FALSE, y has a mean different from zero.
#' @param niter Number of iterations in which the function nlminb will be repeated.
#' @param seed a single value, interpreted as the seed of the random process.
#' @param nmiss a single value; If 1, only one time series of the BiAR process has a missing value. If 2, both time series of the BiAR process have a missing value.
#'
#' @return A list with the following components:
#' \itemize{
#' \item{fitted}{ Estimation of the missing values of the BiAR process.}
#' \item{ll}{ Value of the negative log likelihood evaluated in the fitted missing values.}
#' }
#' @references
#' \insertRef{Elorrieta_2021}{iAR}
#'
#' @seealso
#'
#' \code{\link{gentime}}, \code{\link{BiARsample}}, \code{\link{BiARphikalman}}
#' @keywords internal
#' @examples
#' \dontshow{
#' set.seed(6713)
#' n=100
#' st<-gentime(n)
#' x=BiARsample(phiR=0.9,phiI=0.3,times=st,rho=0.9)
#' y=x$series
#' y1=y/apply(y,1,sd)
#' yerr1=rep(0,n)
#' yerr2=rep(0,n)
#' biar=BiARkalman(series1=y1[1,],series2=y1[2,],times=st,series_esd1 = yerr1,series_esd2=yerr2)
#' biar
#' napos=10
#' y0=y1
#' y1[1,napos]=NA
#' xest=c(biar$phiR,biar$phiI)
#' yest=BiARinterpolation(coef=xest,series1=y1[1,],series2=y1[2,],times=st,series_esd1=yerr1,
#' series_esd2=yerr2,nmiss=1)
#' yest$fitted
#' mse=(y0[1,napos]-yest$fitted)^2
#' print(mse)
#' par(mfrow=c(2,1))
#' plot(st,x$series[1,],type='l',xlim=c(st[napos-5],st[napos+5]))
#' points(st,x$series[1,],pch=20)
#' points(st[napos],yest$fitted*apply(y,1,sd)[1],col="red",pch=20)
#' plot(st,x$series[2,],type='l',xlim=c(st[napos-5],st[napos+5]))
#' points(st,x$series[2,],pch=20)
#' }
#' @noRd
BiARinterpolation <- function (coef,series1, series2, times, series_esd1 = 0, series_esd2 = 0, yini1=0,yini2=0,zero_mean = TRUE, niter = 10, seed = 1234,nmiss=1) {
  set.seed(seed)
  aux <- 1e+10
  value <- 1e+10
  br <- 0

  if (sum(series_esd1) == 0) {
    series_esd1 = rep(0, length(series1))
  }

  if (sum(series_esd2) == 0) {
    series_esd2 = rep(0, length(series2))
  }

  if (yini1==0)
    yini1 = rnorm(1)

  if (yini2==0)
    yini2 = rnorm(1)

  if (nmiss==2){
    for (i in 1:niter) {
      optim <- nlminb(start = c(yini1, yini2), objective = BiARphikalman,
                      coef=coef,series1 = series1,series2 = series2, times = times, series_esd1 = series_esd1,
                      series_esd2=series_esd2, zero_mean = zero_mean,
                      lower = c(-Inf, -Inf), upper = c(Inf, Inf))

      value <- optim$objective

      if (aux > value) {
        par <- optim$par
        aux <- value
        br <- br + 1
      }

      if (aux <= value & br > 1 & i > trunc(niter/2))
        break
    }
  }

  if (aux == 1e+10)
    par <- c(0, 0)

  if (nmiss==1){
    out = optim(c(yini1), BiARphikalman, lower = -Inf, upper = Inf,
                coef=coef, series1=series1, series2=series2, times=times, series_esd1 = series_esd1,
                series_esd2=series_esd2, zero_mean = zero_mean, method="BFGS")
    par = out$par
    aux = out$value
  }

  return(list(fitted = par, ll = aux))
}
