#' Test for the significance of the autocorrelation estimated by the iAR package models in periodic irregularly observed time series
#'
#' This function perform a test for the significance of the autocorrelation estimated by the iAR package models. This test is based on the residuals of the periodical time series fitted with an harmonic model using an incorrect period.
#'
#' @param series Array with the time series observations.
#' @param times Array with the irregular observational times.
#' @param series_esd Array with the variance of the measurement errors.
#' @param f Frequency (1/Period) of the raw time series.
#' @param coef autocorrelation estimated by one of the iAR package models.
#' @param model model used to estimate the autocorrelation parameter ("iAR", "iAR-Gamma", "iAR-T", "CiAR" or "BiAR").
#' @param plot logical; if true, the function return a density plot of the distribution of the bad fitted examples; if false, this function does not return a plot.
#' @param xlim The x-axis limits (x1, x2) of the plot. Only works if plot='TRUE'. See \code{\link{plot.default}} for more details.
#' @param df degrees of freedom parameter of the iAR-T model.
#'
#' @details The null hypothesis of the test is: The autocorrelation estimated in the time series belongs to the distribution of the coefficients estimated for the residuals of the data fitted using wrong periods. Therefore, if the hypothesis is rejected, it can be concluded that the residuals of the harmonic model do not remain a time dependency structure.The statistic of the test is log(phi) which was contrasted with a normal distribution with parameters corresponding to the log of the mean and the variance of the phi computed for the residuals of the bad fitted light curves.
#' @return A list with the following components:
#' \item{phi}{ MLE of the autocorrelation parameter of the iAR/CiAR model.}
#' \item{bad}{ A matrix with two columns. The first column contains the incorrect frequencies used to fit each harmonic model. The second column has the MLEs of the autocorrelation parameters of the iAR/CiAR model that has been fitted to the residuals of the harmonic model fitted using the frequencies of the first column.}
#' \item{norm}{ Mean and variance of the normal distribution of the bad fitted examples.}
#' \item{z0}{ Statistic of the test (log(abs(phi))).}
#' \item{pvalue}{ P-value computed for the test.}
#' @export
#'
#' @references
#' \insertRef{Eyheramendy_2018}{iAR}
#'
#' @seealso
#'
#' \code{\link{clcep}}, \code{\link{harmonicfit}},\code{\link{iARPermutation}}
#'
#' @examples
#' data(clcep)
#' f1=0.060033386
#' #results=harmonicfit(file=clcep,f1=f1)
#' #y=results$res/sqrt(var(results$res))
#' #st=results$t
#' #res3=iARloglik(y,st,standardized=TRUE)[1]
#' #res3$coef
#' #require(ggplot2)
#' #test<-iARTest(series=clcep[,2],times=clcep[,1],f=f1,coef=res3$coef,
#' #model="iAR",plot=TRUE,xlim=c(-10,0.5))
#' #test
iARTest<-function (series, times, series_esd=0,f, coef,model="iAR",plot = TRUE, xlim = c(-1, 0),df=3)
{
  aux = seq(2.5, 47.5, by = 2.5)
  aux = c(-aux, aux)
  aux = sort(aux)
  f0 = f * (1 + aux/100)
  f0 <- sort(f0)
  l1 <- length(f0)
  bad <- rep(0, l1)
  data <- cbind(times, series)
  if(sum(series_esd)==0)
    series_esd=rep(0,length(series))
  for (j in 1:l1) {
    results = harmonicfit(file = data, f1 = f0[j])
    y = results$res/sqrt(var(results$res))
    st = results$t
    if(model=="iAR")
    {
      res3 = iARloglik(y, st,series_esd)[1]
      bad[j] = res3$coef
    }
    if(model=="iAR-T")
    {
      bad[j] = iARt(y, times=st,df=df)$coef
    }
    if(model=="iAR-Gamma")
    {
      y1=y
      if(min(y1)<0)
        y1=y1-min(y1)
      bad[j] = iARgamma(y1, times=st)$coef
    }
    if(model=="CiAR")
    {
      res3 = CiARkalman(series=y, times=st,series_esd=series_esd)[1]
      bad[j] = res3$phiR
    }
  }
  mubf <- mean(log(abs(bad)))
  sdbf <- sd(log(abs(bad)))
  z0 <- log(abs(coef))
  pvalue <- pnorm(z0, mubf, sdbf)
  out <- NULL
  BAD <-data.frame(f0,bad)
  if (plot == TRUE) {
    phi2 = bad
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
  return(list(phi = coef, bad=BAD, norm = c(mubf, sdbf), z0 = z0, pvalue = pvalue, out = out))
}
