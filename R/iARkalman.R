#' Maximum Likelihood Estimation of the iAR Model via Kalman Recursions
#'
#' Maximum Likelihood Estimation of the iAR model parameter phi. The estimation procedure uses the Kalman Filter to find the maximum of the likelihood.
#'
#' @param series Array with the time series observations.
#' @param times Array with the irregular observational times.
#' @param series_esd Array with the measurements error standard deviations.
#' @param zero_mean logical; if TRUE, the array series has zero mean; if FALSE, series has a mean different from zero.
#' @param standardized logical; if TRUE, the array series is standardized; if FALSE, series contains the raw time series.
#' @param hessian logical. Indicates whether to compute the Hessian matrix.; if TRUE, the Hessian matrix will be calculated; if FALSE, the Hessian will not be computed.
#'
#' @return A list with the following components:
#' \itemize{
#' \item{coef}{ MLE of the phi parameter of the iAR model.}
#' \item{ll}{ Value of the negative log likelihood evaluated in phi.}
#' \item{summary}{ If hessian=T, returns the standard error and p-value of the estimated parameter phi.}
#' }
#' @references
#' \insertRef{Eyheramendy_2018}{iAR}
#'
#' @seealso
#'
#' \code{\link{gentime}}, \code{\link{iARsample}}, \code{\link{arima}},\code{\link{iARphikalman}}
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
#' }
#' @noRd
iARkalman<-function (series, times, series_esd=0,zero_mean=TRUE,standardized=TRUE,hessian=FALSE,init=0.1)
{
  if(!is.numeric(init) || length(init)!=1 || !is.finite(init) ||
     init<0 || init>=1)
    stop("initialization value must be a finite number between 0 and 1.",call.=FALSE)
  if(sum(series_esd)==0){
    series_esd=rep(0,length(series))}
  summary=NULL
  if(hessian==F)
  {
    out = optimize(iARphikalman, interval = c(0, 1), series = series, times = times, series_esd=series_esd,zero_mean = zero_mean,standardized=standardized,yest=0)
    phi = out$minimum
    ll = out$objective
    summary=list()
  }
  if(hessian==T)
  {
    out=optim(par=c(x=init),fn=iARphikalman, series=series, times=times, series_esd=series_esd, zero_mean = zero_mean,yest=0, standardized = standardized,method="L-BFGS-B",lower=0,upper=0.9999999,hessian=hessian)
    phi=as.numeric(out$par)
    conv <- out$convergence
    ll=out$value
    hess=as.numeric(out$hessian)
    hess_inv <- tryCatch(solve(hess),error = function(e) NULL)
    if(is.null(hess_inv) || any(!is.finite(hess_inv))) {
      stderr <- NA_real_
      tvalue <- NA_real_
      pvalue <- NA_real_
    } else {
      stderr <- sqrt(diag(hess_inv))
      tvalue <- phi/stderr
      pvalue <- 2*pnorm(abs(tvalue), lower.tail=FALSE)
    }
    summary=list(hess=hess,stderr=stderr,tvalue=tvalue,pvalue=pvalue,convergence=conv)
  }
  return(list(coef=phi,kalman=ll,summary=summary))
}
