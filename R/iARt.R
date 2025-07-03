#' Maximum Likelihood Estimation of the iAR-T model
#'
#' Maximum Likelihood Estimation of the iAR-T model.
#'
#' @param series Array with the time series observations
#' @param times Array with the irregular observational times
#' @param df degrees of freedom 
#' @param hessian logical. Indicates whether to compute the Hessian matrix.; if TRUE, the Hessian matrix will be calculated; if FALSE, the Hessian will not be computed.
#'
#' @return A list with the following components:
#' \itemize{
#' \item{coef}{ MLE of the phi parameter of the iAR-T model.}
#' \item{sigma}{ MLE of the sigma parameter of the iAR-T model.}
#' \item{ll}{ Value of the negative log likelihood evaluated in phi and sigma.}
#' \item{summary}{ If hessian=T, returns the standard error and p-value of the estimated parameter phi.}
#' }
#' @references
#' \insertRef{Eyheramendy_2018}{iAR}
#'
#' @seealso
#'
#' \code{\link{gentime}}, \code{\link{iARtsample}}, \code{\link{iARphit}}
#'
#' @keywords internal
#' @examples
#' \dontshow{
#' n=300
#' set.seed(6714)
#' st<-gentime(n)
#' y<-iARtsample(coef=0.9,times=st,sigma=1,df=3)
#' model<-iARt(series=y$series, times=st)
#' phi=model$coef
#' sigmaest=model$sigma
#' }
#' @noRd
iARt<-function (series, times,df=3, hessian = FALSE)
{
  aux<-1e10
  value<-1e10
  br<-0
  summary=list()
  for(i in 1:20)
  {
    phi=runif(1)
    sigma=var(series)*runif(1)
    if(hessian==F)
    {
      optim<-nlminb(start=c(phi,sigma),objective=iARphit,series=series,times=times,df=df,yest=0,lower=c(0,0.0001),upper=c(0.9999,2*var(series)))
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
      out<-optim(par=c(phi,sigma),fn=iARphit,series=series,times=times,df=df,yest=0,method="L-BFGS-B",lower=c(0,0.0001),upper=c(0.9999,2*var(series)),hessian=hessian)
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
    if(aux<=value & br>10 & i>15)
      break;
    }
    if(aux==1e10)
      par<-c(0,0)
    return(list(coef=par[1],sigma=par[2],ll=aux,summary=summary))
}
