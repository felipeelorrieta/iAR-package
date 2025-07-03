#' Maximum Likelihood Estimation of the iAR Model
#'
#' Maximum Likelihood Estimation of the iAR Model.
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
#' \code{\link{gentime}}, \code{\link{iARsample}}, \code{\link{arima}}, \code{\link{iARphiloglik}}
#'
#' @keywords internal
#' @examples
#' \dontshow{
#' #Generating iAR sample
#' set.seed(6714)
#' st<-gentime(n=100)
#' y<-iARsample(coef=0.99,times=st)
#' y<-y$series
#' #Compute Phi
#' phi=iARloglik(series=y,times=st)$coef
#' print(phi)
#' #Compute the standard deviation of innovations
#' n=length(y)
#' d=c(0,diff(st))
#' phi1=phi**d
#' yhat=phi1*as.vector(c(0,y[1:(n-1)]))
#' plot(st,y,type='l')
#' lines(st,yhat,col='red')
#' sigma=var(y)
#' nu=c(sigma,sigma*(1-phi1**(2))[-1])
#' tau<-nu/sigma
#' sigmahat<-mean(c((y-yhat)**2/tau))
#' nuhat<-sigmahat*(1-phi1**(2))
#' nuhat2<-sqrt(nuhat)
#' #Equally spaced models
#' require(arfima)
#' fit2<-arfima(y,order=c(1,0,0))
#' fit<-arima(y,order=c(1,0,0),include.mean=FALSE)
#' syarf<-tacvfARFIMA(phi=fit2$modes[[1]]$phi,dfrac=fit2$modes[[1]]$dfrac,
#' sigma2=fit2$modes[[1]]$sigma,maxlag=20)[1]
#' syar<-fit$sigma/(1-fit$coef[1]**2)
#' print(sigmahat)
#' print(syar)
#' print(syarf)
#' carf<-fit2$modes[[1]]$sigma/syarf
#' car<-(1-fit$coef[1]**2)
#' ciar<-(1-phi1**(2))
#' #Compute the standard deviation of innovations (regular case)
#' sigma=var(y)
#' nuhat3=sqrt(sigma*ciar)
#' searf<-sqrt(sigma*carf)
#' sear<-sqrt(sigma*car)
#' #Plot the standard deviation of innovations
#' plot(st[-1], nuhat3[-1], t="n", axes=FALSE,xlab='Time',ylab='Standard Deviation of Innovations')
#' axis(1)
#' axis(2)
#' segments(x0=st[-1], y0=nuhat3[-1], y1=0, col=8)
#' points(st, nuhat3, pch=20, col=1, bg=1)
#' abline(h=sd(y),col='red',lwd=2)
#' abline(h=sear,col='blue',lwd=2)
#' abline(h=searf,col='green',lwd=2)
#' abline(h=mean(nuhat3[-1]),col='black',lwd=2)
#' }
#' @noRd
iARloglik=function(series,times,series_esd=0,zero_mean=TRUE,standardized=TRUE,hessian=FALSE){
  if(sum(series_esd)==0){
    series_esd=rep(0,length(series))}
  summary=list()
  if(hessian==F)
  {
    out=optimize(iARphiloglik, interval=c(0,1), series=series, times=times, series_esd=series_esd, zero_mean = zero_mean, standardized = standardized)
    phi=out$minimum
    ll=out$objective
  }
  if(hessian==T)
  {
    out=optim(par=c(x=0),fn=iARphiloglik, series=series, times=times, series_esd=series_esd, zero_mean = zero_mean, standardized = standardized,method="L-BFGS-B",lower=0,upper=0.9999999,hessian=hessian)
    phi=as.numeric(out$par)
    ll=out$value
    hess=as.numeric(out$hessian)
    stderr=sqrt(diag(solve(hess)))
    tvalue=phi/stderr
    pvalue=2*pnorm(abs(tvalue),lower.tail=F)
    summary=list(hess=hess,stderr=stderr,tvalue=tvalue,pvalue=pvalue)
  }
  return(list(coef=phi,loglik=ll,summary=summary))
}
