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
iARt<-function (series, times,df=3, hessian = FALSE, init=c(0,1))
{
  if(!is.numeric(init) || length(init)!=2 || any(!is.finite(init)))
    stop("init must be a numeric vector of length 2.",call.=FALSE)
  if(!all(init==c(0,1)) &&
     (init[1]<0 || init[1]>=1 ||
      init[2]<=0 ))
    stop("Invalid initial values for phi and sigma.",call.=FALSE)
  aux<-1e10
  value<-1e10
  summary=list()
  vals=c(0.1,0.5,0.9)
  sigma=var(series)*vals
  niter <- if(all(init == c(0,1))) 3 else 4
  for(i in 1:niter)
  {
    if(i <= 3)
      phi=c(vals[i],sigma[i])
    if(!all(init == c(0,1)) && i == 4)
      phi = init
    if(hessian==F)
    {
      optim<-nlminb(start=phi,objective=iARphit,series=series,times=times,df=df,yest=0,lower=c(0,0.0001),upper=c(0.9999,2*var(series)))
      value<-optim$objective
      par0<-optim$par
    }
    if(hessian==T)
    {
      out<-optim(par=phi,fn=iARphit,series=series,times=times,df=df,yest=0,method="L-BFGS-B",lower=c(0,0.0001),upper=c(0.9999,2*var(series)),hessian=hessian)
      value<-out$value
      par0<-out$par
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
    if(aux==1e10)
      par<-c(0,0)
    return(list(coef=par[1],sigma=par[2],ll=aux,summary=summary))
}
