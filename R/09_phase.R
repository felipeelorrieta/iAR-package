#' Computing phased time series
#'
#' This function computes a time series folded on its period.
#'
#' @name phase
#' 
#' @usage phase(x, ...)
#'  
#' @param x An object of class `utilities`.
#' @param ... Additional arguments for pairing time series:
#'   \describe{
#'     \item{data}{A data frame with three columns corresponding to the time, values, and standard errors of the irregularly observed time series.}
#'     \item{f1}{frequency (1 / period) of the time series.}
#'     \item{twop}{logical; if TRUE, the phased series will be duplicated over two cycles (0–2).}
#'    }
#' @return An object of class `utilities` with the slots:
#' \item{series_phased}{A numeric vector containing the time series values ordered by phase.}
#' \item{series_esd_phased}{A numeric vector containing the error standard deviations of the time series ordered by phase.}
#' \item{times_phased}{A numeric vector of phased times (values between 0 and 1, or 0 and 2 if 'two.cycles = TRUE').}
#' 
#' @details 
#' The phase \eqn{\phi} of an observation is computed as
#' \deqn{\phi = \frac{t - t_0}{p} - \mathrm{E}(t),}
#' where \eqn{t_0} is the reference time (by default the first observation),
#' \eqn{p = 1/f_1} is the period, and \eqn{\mathrm{E}(t)} is the integer part of
#' \eqn{(t - t_0)/p}.
#' 
#' @export
#' @references
#' \insertRef{Elorrieta_2021}{iAR}
#'
#' @examples
#' data(clcep)
#' f1=0.060033386
#' o1=iAR::utilities()
#' o1<-phase(o1,data=clcep,f1=f1,twop=TRUE)
#' plot(o1@times_phased,o1@series_phased,pch=20)
phase <- S7::new_generic("phase", "x")
S7::method(phase, utilities) <- function(x, data, f1, twop=TRUE)
{
  mycurve = data
  t = mycurve[, 1]
  m = mycurve[, 2]
  merr=rep(0,length(t))
  if(dim(data)[2]>2)
    merr = mycurve[, 3]
  P <- 1/f1
  fold <- (t - t[1])%%(P)/P
  if(twop==TRUE)
  {
    fold <- c(fold, fold + 1)
    m <- rep(m, 2)
    merr <- rep(merr, 2)
  }
  dat1 <- cbind(fold, m, merr)
  dat1 <- as.data.frame(dat1)
  dat1<-dat1[order(dat1[,1],decreasing=TRUE),]
  x@series_phased <- dat1[,2]
  x@series_esd_phased <- dat1[,3]
  x@times_phased <- dat1[,1]
  return(x)
}
