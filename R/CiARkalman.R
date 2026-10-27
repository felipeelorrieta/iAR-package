#' Maximum Likelihood Estimation of the CiAR Model via Kalman Recursions
#'
#' Maximum Likelihood Estimation of the CiAR model parameters phiR and phiI. The estimation procedure uses the Kalman Filter to find the maximum of the likelihood.
#'
#' @param series Array with the time series observations.
#' @param times Array with the irregular observational times.
#' @param series_esd Array with the measurements error standard deviations.
#' @param zero_mean logical; if TRUE, the array series has zero mean; if FALSE, series has a mean different from zero.
#' @param standardized logical; if TRUE, the array series is standardized; if FALSE, series contains the raw time series.
#' @param hessian logical. Indicates whether to compute the Hessian matrix.; if TRUE, the Hessian matrix will be calculated; if FALSE, the Hessian will not be computed.
#' @param c Nuisance parameter corresponding to the variance of the imaginary part.
#'
#' @return A list with the following components:
#' \itemize{
#' \item{phiR}{ MLE of the Real part of the coefficient of CiAR model (phiR).}
#' \item{phiI}{ MLE of the Imaginary part of the coefficient of the CiAR model (phiI).}
#' \item{ll}{ Value of the negative log likelihood evaluated in phiR and phiI.}
#' \item{summary}{ If hessian=T, returns the standard error and p-value of the estimated parameter phi.}
#' }
#' @references
#' \insertRef{Elorrieta_2019}{iAR}
#'
#' @seealso
#'
#' \code{\link{gentime}}, \code{\link{CiARsample}}, \code{\link{CiARphikalman}}
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
#' Mod(complex(real=ciar$phiR,imaginary=ciar$phiI))
#' }
#' @noRd
CiARkalman<-function(series,times,series_esd=0,zero_mean=TRUE,standardized=TRUE,hessian=FALSE,c=1,init=rep(0,2))
{
  aux<-1e10
  value<-1e10
  if(sum(series_esd)==0){
    series_esd=rep(0,length(series))}
  summary=list()
  niter <- if(all(init == c(0,0))) 5 else 6
  for(i in 1:niter)
  {
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
    phiR=phi[1]
    phiI=phi[2]
    if(Mod(complex(1,real=phiR,imaginary=phiI))<1)
    {
      if(hessian==F)
      {
        optim <- nlminb(start=c(phiR,phiI), objective=CiARphikalman, series=series, times=times,
                      series_esd=series_esd, zero_mean=zero_mean, standardized=standardized,
                      c=c,yest=0, lower=c(-1,-1), upper=c(1,1))
        value<-optim$objective
        par0<-optim$par
      }
      if(hessian==T)
      {
        out<-optim(par=c(phiR,phiI),fn=CiARphikalman,series=series, times=times,
                   series_esd=series_esd, zero_mean=zero_mean, standardized=standardized,
                   c=c,yest=0,method="L-BFGS-B",lower=c(-1,-1),upper=c(1,1),hessian=hessian)
        value<-out$value
        par0<-out$par
      }
    }
    if(aux>value)
    {
      par<-par0
      aux<-value
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
  if(aux==1e10)
    par<-c(0,0)
  return(list(phiR=par[1],phiI=par[2],ll=aux,summary=summary))
}

