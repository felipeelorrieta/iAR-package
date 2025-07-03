#' Harmonic Fit to Time Series
#'
#' This function fit an k-harmonic function to time series data.
#'
#' @param file A matrix with two columns. The first column corresponds to the observations times, and the second column corresponds to the measures.
#' @param f1 Frequency (1/Period) of the time series
#' @param nham Number of harmonic components in the model
#' @param weights An array with the weights of each observation
#' @param print logical; if true, the summary of the harmonic fitted model will be printed. The default value is false.
#' @param remove_trend logical; if true, the linear trend of time series will be removed before the the harmonic model is fitted.
#'
#' @return  A list with the following components:
#' \item{res}{ Residuals to the harmonic fit of the time series.}
#' \item{t}{ Observations times.}
#' \item{R2}{ Adjusted R-Squared.}
#' \item{MSE}{ Mean Squared Error.}
#' \item{coef}{ Summary of the coefficients estimated by the harmonic model.}
#' @export
#'
#' @examples
#' data(clcep)
#' f1=0.060033386
#' #results=harmonicfit(file=clcep[,1:2],f1=f1)
#' #results$R2
#' #results$MSE
#' #results=harmonicfit(file=clcep[,1:2],f1=f1,nham=3)
#' #results$R2
#' #results$MSE
#' #results=harmonicfit(file=clcep[,1:2],f1=f1,weights=clcep[,3])
#' #results$R2
#' #results$MSE
harmonicfit<-function (file, f1, nham = 4, weights = NULL, print = FALSE,remove_trend=TRUE)
{
  mycurve = file
  t = mycurve[, 1]
  y = mycurve[, 2]
  if(remove_trend==TRUE)
    y = residuals(lm(mycurve[, 2] ~ t))
  nam = NULL
  for (i in 1:nham) {
    y = cbind(y, sin(2 * i * pi * f1 * t), cos(2 * i * pi *
                                                 f1 * t))
    nam=c(nam,paste(c("sin","cos"),"comp",i,sep=""))
  }
  df = data.frame(y)
  names(df)[-1]=nam
  model = lm(y ~ ., data = df)
  if (length(weights) > 0)
    model = lm(y ~ ., data = df, weights = weights)
  if (print == TRUE)
    print(summary(model))
  res = residuals(model)
  R2 <- summary(model)$adj.r.squared
  MSE <- anova(model)$"Mean Sq"[nham * 2 + 1]
  coef <- summary(model)$coefficients
  return(list(res = res, t = t, R2 = R2, MSE = MSE,coef=coef))
}
