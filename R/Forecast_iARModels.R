#' Forecast from iAR package model's
#'
#' Forecast with any of the models available in the iAR package
#'
#' @param coef Autocorrelation coefficient estimated by the method specified.
#' @param series Array with the time series observations.
#' @param times Array with the observational times.
#' @param tAhead The time ahead for which the forecast is required.
#' @param model model to be used for the forecast. The default is to use the iAR model. Other models available are "iAR-T", "iAR-Gamma", "CiAR" and "BiAR".
#' @param mu Level parameter of the iAR-Gamma process. A positive value.
#' @param phiI Imaginary parameter of CiAR model or Cross-correlation parameter of BiAR model.
#' @param df degrees of freedom parameter of iAR-T model.
#' @param level significance level for the confidence interval. The default value is 95.
#'
#' @return A dataframe with the following columns:
#' \itemize{
#' \item{tAhead}{ The time ahead used for the forecast.}
#' \item{forecast}{ Point forecast in the time ahead required.}
#' \item{stderror}{ Standard error of the forecast.}
#' \item{lowerCI}{ Lower limit of the confidence interval.}
#' \item{upperCI}{ Upper limit of the confidence interval.}
#' }
#' @references
#' \insertRef{Eyheramendy_2018}{iAR}
#'
#' @seealso
#'
#' \code{\link{gentime}}, \code{\link{iARforecast}}, \code{\link{iARgforecast}}, \code{\link{iARforecast}}, \code{\link{BiARforecast}}
#'
#' @keywords internal
#' @examples
#' \dontshow{
#' st <- gentime(n=200,lambda1=15,lambda2=2)
#' y  <- iARsample(coef=0.9,times=st)
#' model<-iARloglik(series=y$series,times=st)
#' phi=model$coef
#' foriAR<-iARforecast(coef=phi,series=y$series,tAhead=c(1.3),standardized=FALSE,zero_mean=FALSE)
#' foriAR
#' foriAR<-Forecast_iARModels(coef=phi,series=y$series,times=st,tAhead=c(1.3,2.6))
#' foriAR
#' }
#' @noRd
Forecast_iARModels<-function(coef,series,times,tAhead,model="iAR",mu=NULL,phiI=NULL,df=NULL,level=95)
{
  if(model == "iAR")
  {
    y=as.numeric(series)
    sigma = sqrt(var(series))
    mu = mean(series)
    y1=(series-mu)/sigma
    fore=iARforecast(coef=coef,series=y1,tAhead=tAhead)
    fore=fore*sigma+mu
    error=sigma*(1-coef**(2*tAhead))
  }
  if(model == "iAR-T")
  {
    fore=iARforecast(coef=coef,series=series,tAhead=tAhead,standardized=FALSE,zero_mean=FALSE)
    error=((df-2)/df)*sigma*(1-coef**(2*tAhead))
  }
  if(model == "iAR-Gamma")
  {
    y1=series
    if(min(series)<0)
      y1=series-min(series)
    fore=iARgforecast(coef=coef,mean=mu,series=series,tAhead=tAhead)
    if(min(series)<0)
      fore=fore+min(series)
    error=sigma*(1-coef**(2*tAhead))
  }
  if(model == "CiAR")
  {
    if(is.null(phiI))
      phiI=0
    sigma = sqrt(var(as.numeric(series)))
    mu = mean(as.numeric(series))
    y1=series
    fore=CiARforecast(coef=c(coef,phiI),series=y1,times=times,tAhead=tAhead)$forecast
    Mod=sqrt(coef^2+phiI^2)
    error=sigma*(1-Mod**(2*tAhead))
  }
  if(model == "BiAR")
  {
    mu = apply(y,1,mean)
    y1=series
    y2=y1[2,]
    y1=y1[1,]
    fore=BiARforecast(coef=c(coef,phiI),series1=y1,series2=y2,times=times,tAhead=tAhead)$forecast
    sigma=matrix(apply(series,1,sd),nrow=2)
    Mod=sqrt(coef^2+phiI^2)
    error=sigma*(1-Mod**(2*tAhead))
  }
  upper=fore+qnorm(level/100)*error
  lower=fore-qnorm(level/100)*error
  fin=data.frame(tAhead=tAhead,forecast=fore,stderror=error,lowerCI=lower,upperCI=upper)
  return(fin)
}
