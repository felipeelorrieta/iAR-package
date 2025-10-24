#' Harmonic Fit to Time Series
#'
#' This function fit an k-harmonic function to time series data.
#'
#' @name harmonicfit
#' 
#' @usage harmonicfit(x, ...)
#'  
#' @param x An object of class `utilities`.
#' @param ... Additional arguments for pairing time series:
#'   \describe{
#'     \item{data}{A data frame with three columns corresponding to the time, values, and standard errors of the irregularly observed time series.}
#'     \item{f1}{frequency (1 / period) of the time series.}
#'     \item{nham}{Number of harmonic components in the model.}
#'     \item{weighted}{logical; if true, a weighted least squares (WLS) estimation is performed using weights based on the standard deviations of the errors. Default is 'FALSE'.}
#'     \item{remove_trend}{logical; if true, the linear trend of time series will be removed before the the harmonic model is fitted.}
#'    }
#' @return An object of class `utilities` with the slots:
#' \item{fitted_values}{Fitted values from the harmonic model.}
#' \item{residuals}{Residuals from the harmonic model.} 
#' \item{coef}{Estimated coefficients of the harmonic model.}
#' \item{summary}{A summary object containing detailed model information.}
#' 
#' @details The function fits a harmonic regression model with 'nham' components to the input time series, optionally removing a linear trend and allowing for weighted estimation when standard errors are available. 
#' @export
#' @references
#' \insertRef{Elorrieta_2021}{iAR}
#'
#' @examples
#' data(clcep)
#' f1=0.060033386
#' o1=iAR::utilities()
#' o1<-phase(o1,data=clcep,f1=f1)
#' #results$R2
#' #results$MSE
#' #results=harmonicfit(file=clcep[,1:2],f1=f1,nham=3)
#' #results$R2
#' #results$MSE
#' #results=harmonicfit(file=clcep[,1:2],f1=f1,weights=clcep[,3])
#' #results$R2
#' #results$MSE
harmonicfit <- S7::new_generic("harmonicfit", "x")
S7::method(harmonicfit, utilities) <- function (x,data, f1, nham = 4, weighted=FALSE, remove_trend=TRUE)
{
  mycurve = data
  t = mycurve[, 1]
  y = mycurve[, 2]
  merr=rep(0,length(t))
  if(dim(data)[2]>2)
    merr = mycurve[, 3]
  if(remove_trend==TRUE & weighted==FALSE)
    y = residuals(lm(mycurve[, 2] ~ t))
  if(remove_trend==TRUE & weighted==TRUE)
    y = residuals(lm(mycurve[, 2] ~ t,weights=merr))
  nam = NULL
  for (i in 1:nham) {
    y = cbind(y, sin(2 * i * pi * f1 * t), cos(2 * i * pi *
                                                 f1 * t))
    nam=c(nam,paste(c("sin","cos"),"comp",i,sep=""))
  }
  df = data.frame(y)
  names(df)[-1]=nam
  model = lm(y ~ ., data = df)
  if (weighted==TRUE)
    model = lm(y ~ ., data = df, weights = merr)
  #res = residuals(model)
  #R2 <- summary(model)$adj.r.squared
  #MSE <- anova(model)$"Mean Sq"[nham * 2 + 1]
  x@fitted_values <- model$fitted.values
  x@residuals <- residuals(model)
  x@coef <- summary(model)$coefficients
  x@summary <-summary(model)
  return(x)
}
