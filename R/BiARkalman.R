#' Maximum Likelihood Estimation of the BiAR Model via Kalman Recursions
#'
#' Maximum Likelihood Estimation of the BiAR model parameters phiR and phiI. The estimation procedure uses the Kalman Filter to find the maximum of the likelihood.
#'
#'
#' @param series1 Array with the observations of the first time series of the BiAR process.
#' @param series2 Array with the observations of the second time series of the BiAR process.
#' @param times Array with the irregular observational times.
#' @param series_esd1 Array with the measurements error standard deviations of the first time series of the BiAR process.
#' @param series_esd2 Array with the measurements error standard deviations of the second time series of the BiAR process.
#' @param zero_mean logical; if true, the array y has zero mean; if false, y has a mean different from zero.
#' @param hessian logical. Indicates whether to compute the Hessian matrix.; if TRUE, the Hessian matrix will be calculated; if FALSE, the Hessian will not be computed.
#'
#' @return A list with the following components:
#' \itemize{
#' \item{phiR}{ MLE of the autocorrelation coefficient of BiAR model (phiR).}
#' \item{phiI}{ MLE of the cross-correlation coefficient of the BiAR model (phiI).}
#' \item{ll}{ Value of the negative log likelihood evaluated in phiR and phiI.}
#' \item{summary}{ If hessian=T, returns the standard error and p-value of the estimated parameter phi.}
#' }
#' @references
#' \insertRef{Elorrieta_2021}{iAR}
#'
#' @seealso
#'
#' \code{\link{gentime}}, \code{\link{BiARsample}}, \code{\link{BiARphikalman}}
#'
#' @keywords internal
#' @examples
#' \dontshow{
#' n=80
#' set.seed(6714)
#' st<-gentime(n)
#' x=BiARsample(phiR=0.9,phiI=0,times=st,rho=0)
#' y=x$series
#' y1=y/apply(y,1,sd)
#' biar=BiARkalman(series1=y1[1,],series2=y1[2,],times=st,series_esd1 = rep(0,length(y[1,])),
#' series_esd2=rep(0,length(y[1,])))
#' biar
#' }
#' @noRd
BiARkalman<-function (series1, series2, times, series_esd1 = 0, series_esd2 = 0, zero_mean = TRUE, hessian=FALSE,init=rep(0,2))
{
  aux <- 1e+10
  value <- 1e+10
  if (sum(series_esd1) == 0) {
    series_esd1 = rep(0, length(series1))
  }
  if (sum(series_esd2) == 0) {
    series_esd2 = rep(0, length(series2))
  }
  niter <- if(all(init == c(0,0))) 5 else 6
  summary=list()
  for (i in 1:niter) {
    if(i == 1)
      phi=c(0.1,-0.00001)
    if(i == 2)
      phi=c(0.5,0.00001)
    if(i == 3)
      phi=c(-0.00001,-0.5)
    if(i == 4)
      phi=c(0.5,-0.5)
    if(i == 5)
      phi=c(-0.5,0.5)
    if(!all(init == c(0,0)) && i == 6)
      phi = init
    if (Mod(complex(1, real = phi[1], imaginary = phi[2])) < 1) {
      if(hessian==F)
      {
        optim <- nlminb(start = c(phi[1], phi[2]), objective = BiARphikalman,
                      series1 = series1,series2 = series2, times = times, series_esd1 = series_esd1, series_esd2=series_esd2,zero_mean = zero_mean, yest=c(0,0),lower	= c(-1, -1), upper = c(1, 1))
        value <- optim$objective
        par0<-optim$par
      }
      if(hessian==T)
      {
        out<-optim(par=c(phi[1], phi[2]),fn=BiARphikalman,series1 = series1,series2 = series2, times = times, series_esd1 = series_esd1, series_esd2=series_esd2,zero_mean = zero_mean, yest=c(0,0),method="L-BFGS-B",lower=c(-1,-1),upper=c(1,1),hessian=hessian)
        value<-out$value
        par0<-out$par
      }
    }
    if (aux > value) {
      par <- par0
      aux <- value
      if(hessian==T)
      {
        conv <- out$convergence
        hess=out$hessian
      }
    }
  }
  if(hessian==T)
  {
    hess_inv <- tryCatch(solve(hess),error = function(e) NULL)
    if(is.null(hess_inv) || any(!is.finite(hess_inv))) {
      stderr <- rep(NA_real_, length(par))
      tvalue <- rep(NA_real_, length(par))
      pvalue <- rep(NA_real_, length(par))
    } else {
      stderr <- as.numeric(sqrt(diag(hess_inv)))
      tvalue <- par/stderr
      pvalue <- 2*pnorm(abs(tvalue), lower.tail=FALSE)
    }
    summary=list(hess=hess,stderr=stderr,tvalue=tvalue,pvalue=pvalue,convergence=conv)
  }
  if (aux == 1e+10)
    par <- c(0, 0)
  return(list(phiR = par[1], phiI = par[2], ll = aux,summary=summary))
}
