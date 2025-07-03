#' Plotting folded light curves
#'
#' This function plots a time series folded on its period.
#' @param file Matrix with the light curve observations. The first column must have the irregular times, the second column must have the brightness magnitudes and the third column must have the measurement errors.
#' @param f1 Frequency (1/Period) of the light curve.
#' @param plot logical; if TRUE, the function returns the plot of folded time series.
#'
#' @return A matrix whose first column has the folded (phased) observational times.
#' @export
#'
#' @examples
#' data(clcep)
#' f1=0.060033386
#' foldlc(clcep,f1)
foldlc=function (file, f1, plot=TRUE)
{
  mycurve = file
  t = mycurve[, 1]
  m = mycurve[, 2]
  merr=rep(0,length(t))
  if(dim(file)[2]>2)
    merr = mycurve[, 3]
  P <- 1/f1
  fold <- (t - t[1])%%(P)/P
  fold <- c(fold, fold + 1)
  m <- rep(m, 2)
  merr <- rep(merr, 2)
  dat1 <- cbind(fold, m, merr)
  dat1 <- as.data.frame(dat1)
  out=NULL
  if(plot==TRUE)
    out <- ggplot(dat1, aes(x = fold, y = m)) + geom_errorbar(aes(ymin = m -
                                                                    merr, ymax = m + merr), width = 0.01, col = "red") +
    geom_point() + scale_y_reverse() + xlab("") + ylab("") +
    theme_bw() + theme(plot.title = element_text(face = "bold",
                                                 size = 20), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       panel.background = element_blank())
  dat1<-dat1[order(dat1[,1],decreasing=TRUE),]
  return(list(folded=dat1,plot=out))
}
