#' Test for the significance of the autocorrelation estimated by the iAR package models
#'
#' This function perform a test for the significance of the autocorrelation estimated by the iAR package models. This test is based in to take N disordered samples of the original data.
#'
#' @param series Array with the time series observations.
#' @param times Array with the irregular observational times.
#' @param series_esd Array with the variance of the measurement errors.
#' @param iter Number of disordered samples of the original data (N).
#' @param coef autocorrelation estimated by one of the iAR package models.
#' @param model model used to estimate the autocorrelation parameter ("iAR", "iAR-Gamma", "iAR-T", "CiAR" or "BiAR").
#' @param plot logical; if true, the function return a density plot of the distribution of the bad fitted examples; if false, this function does not return a plot.
#' @param xlim The x-axis limits (x1, x2) of the plot. Only works if plot='TRUE'. See \code{\link{plot.default}} for more details.
#' @param df degrees of freedom parameter of the iAR-T model.
#'
#'
#' @details The null hypothesis of the test is: The autocorrelation coefficient estimated for the time series belongs to the distribution of the coefficients estimated on the disordered data, which are assumed to be uncorrelated. Therefore, if the hypothesis is accepted, it can be concluded that the observations of the time series are uncorrelated.The statistic of the test is log(phi) which was contrasted with a normal distribution with parameters corresponding to the log of the mean and the variance of the phi computed for the N samples of the disordered data. This test differs for \code{\link{iARTest}} in that to perform this test it is not necessary to know the period of the time series.
#' @return A list with the following components:
#' \item{coef}{ MLE of the autocorrelation parameter of the model.}
#' \item{bad}{ MLEs of the autocorrelation parameters of the models that has been fitted to the disordered samples.}
#' \item{norm}{ Mean and variance of the normal distribution of the disordered data.}
#' \item{z0}{ Statistic of the test (log(abs(phi))).}
#' \item{pvalue}{ P-value computed for the test.}
#' @export
#'
#' @references
#' \insertRef{Eyheramendy_2018}{iAR}
#'
#' @seealso
#'
#' \code{\link{Planets}},\code{\link{iARTest}}
#'
#' @examples
#' data(Planets)
#' t<-Planets[,1]
#' res<-Planets[,2]
#' y=res/sqrt(var(res))
#' #res3=iARloglik(y,t,standardized=TRUE)[1]
#' #res3$coef
#' #set.seed(6713)
#' #require(ggplot2)
#' #test<-iARPermutation(series=y,times=t,coef=res3$coef,model="iAR",plot=TRUE,xlim=c(-9.6,-9.45))
iARPermutation<-function (series, times, series_esd=0, iter = 100, coef,model="iAR", plot = TRUE, xlim = c(-1,0),df=3)
{
  phi2 = rep(0, iter)
  if(sum(series_esd)==0)
    series_esd=rep(0,length(series))
  for (i in 1:iter) {
    ord <- sample(1:length(series))
    y1 <- series[ord]
    merr1 <- series_esd[ord]
    if(model=="iAR")
    {
      phi2[i] = iARloglik(series = y1, times = times,merr1)$coef
    }
    if(model=="iAR-T")
    {
      phi2[i] = iARt(y1, times=times,df=df)$coef
    }
    if(model=="iAR-Gamma")
    {
      if(min(y1)<0)
        y1=y1-min(y1)
      phi2[i] = iARgamma(y1, times=times)$coef
    }
    if(model=="CiAR")
    {
      phi2[i] = CiARkalman(series=y1, times=times,series_esd=merr1)[1]$phiR
    }
    if(model=="BiAR")
    {
      y1=series
      y2=y1[2,]
      y1=y1[1,]
      if(sum(series_esd)==0)
        series_esd=matrix(0,nrow=2,ncol=length(y1))
      merr1=series_esd
      merr2=merr1[2,]
      merr1=merr1[1,]
      ord <- sample(1:length(y1))
      y1 <- y1[ord,]
      y2 <- y2[ord,]
      merr1 <- merr1[ord,]
      merr2 <- merr2[ord,]
      phi2[i] = BiARkalman(series1=y1,series2=y2,times=times,series_esd1 = merr1,series_esd2=merr2)$phiR
    }
  }
  mubf <- mean(log(abs(phi2)))
  sdbf <- sd(log(abs(phi2)))
  z0 <- log(abs(coef))
  pvalue <- pnorm(z0, mubf, sdbf)
  out <- NULL
  BAD =phi2
  if (plot == TRUE) {
    phi2 <- as.data.frame(phi2)
    coef <- as.data.frame(coef)
    out <- ggplot(phi2, aes(log(abs(phi2)))) + geom_density(adjust = 2) +
      xlab("") + ylab("") + theme_bw() + ggtitle("") +
      geom_point(data = coef, aes(log(abs(coef))), y = 0, size = 4,
                 col = "red", shape = 17) + xlim(xlim[1], xlim[2]) +
      theme(plot.title = element_text(face = "bold", size = 20),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank())
    out
  }
  return(list(coef = coef, bad=BAD,norm = c(mubf, sdbf), z0 = z0, pvalue = pvalue,
              out = out))
}
