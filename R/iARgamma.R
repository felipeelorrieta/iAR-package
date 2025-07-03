#' Maximum Likelihood Estimation of the iAR-Gamma model
#'
#' Maximum Likelihood Estimation of the iAR-Gamma model.
#'
#' @param series Array with the time series observations
#' @param times Array with the irregular observational times
#' @param hessian logical. Indicates whether to compute the Hessian matrix.; if TRUE, the Hessian matrix will be calculated; if FALSE, the Hessian will not be computed.
#' 
#' @return A list with the following components:
#' \itemize{
#' \item{coef}{ MLE of the phi parameter of the iAR-Gamma model.}
#' \item{mean}{ MLE of the mean parameter of the iAR-Gamma model.}
#' \item{sigma}{ MLE of the sigma parameter of the iAR-Gamma model.}
#' \item{ll}{ Value of the negative log likelihood evaluated in phi, mean and sigma.}
#' \item{summary}{ If hessian=T, returns the standard error and p-value of the estimated parameter phi.}
#' }
#' @references
#' \insertRef{Eyheramendy_2018}{iAR}
#'
#' @seealso
#'
#' \code{\link{gentime}}, \code{\link{iARgsample}}, \code{\link{iARphigamma}}
#'
#'
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
#' }
#' @noRd
iARgamma<-function(series, times, hessian = FALSE)
{
  aux<-1e10
  value<-1e10
  br<-0
  summary=list()
  for(i in 1:20)
  {
    phi=runif(1)
    mu=mean(series)*runif(1)
    sigma=var(series)*runif(1)
    if(hessian==F)
    {
      optim<-nlminb(start=c(phi,mu,sigma),objective=iARphigamma,series=series,times=times,yest=0,lower=c(0,0.0001,0.0001),upper=c(0.9999,mean(series),var(series)))
      value<-optim$objective
      if(aux>value)
      {
        par<-optim$par
        aux<-value
        br<-br+1
      }
    }
    if(hessian==T)
    {
      out<-optim(par=c(phi,mu,sigma),fn=iARphigamma,series=series,times=times,yest=0,method="L-BFGS-B",lower=c(0,0.0001,0.0001),upper=c(0.9999,mean(series),var(series)),hessian=hessian)
      value<-out$value
      if(aux>value)
      {
        par<-out$par
        aux<-value
        hess=out$hessian
        stderr=as.numeric(sqrt(diag(solve(hess))))
        tvalue=par/stderr
        pvalue=2*pnorm(abs(tvalue),lower.tail=F)
        summary=list(hess=hess,stderr=stderr,tvalue=tvalue,pvalue=pvalue)
        br<-br+1
      }
    }
    if(aux<=value & br>5 & i>10)
        break;
    }
  if(aux==1e10)
    par<-c(0,0,0)
  return(list(coef=par[1],mean=par[2],sigma=par[3],ll=aux,summary=summary)) 
}
