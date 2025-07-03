#' Pairing two irregularly observed time series
#'
#' This method pairs the observational times of two irregularly observed time series.
#' @name pairingits
#' 
#' @usage pairingits(x, ...)
#' 
#' @param x An object of class `utilities`.
#' @param ... Additional arguments for pairing time series: 
#'   \describe{
#'     \item{lc1}{A data frame with three columns corresponding to the first irregularly observed time series.}
#'     \item{lc2}{A data frame with three columns corresponding to the second irregularly observed time series.}
#'     \item{tol}{A numeric value indicating the tolerance parameter.}
#'   }
#'
#'
#' @return An object of class `utilities` with two slots:
#' \item{series}{A matrix containing the paired time series, where unmatched measurements are filled with `NA`.}
#' \item{series_esd}{A matrix containing the paired error standard deviations of the time series, where unmatched measurements are filled with `NA`.}
#' \item{times}{A numeric vector with the paired observational times.}
#'
#' @details The method checks the observational times in both input time series and pairs the measurements if they fall within the specified tolerance (`tol`). If a measurement in one series cannot be paired, it is filled with `NA` values for the corresponding columns of the other series.
#'
#' @export
#' @references
#' \insertRef{Elorrieta_2021}{iAR}
#'
#' @examples
#' data(cvnovag)
#' data(cvnovar)
#' datag=cvnovag
#' datar=cvnovar
#' o1=iAR::utilities()
#' o1<-pairingits(o1, datag,datar,tol=0.1)
#' pargr1=na.omit(o1@paired)
#' st=apply(pargr1[,c(1,4)],1,mean)
#' model_BiAR <- BiAR(times = st,series=pargr1[,c(2,5)],series_esd=pargr1[,c(3,6)])
#' model_BiAR <- kalman(model_BiAR)
#' 
pairingits <- S7::new_generic("pairingits", "x")
S7::method(pairingits, utilities) <- function(x, lc1,lc2,tol=0.1)
{
  t1=lc1[,1]
  t2=lc2[,1]
  A=cbind(c(t1,t2),c(rep(1,length(t1)),rep(2,length(t2))),c(1:length(t1),1:length(t2)))
  A=A[order(A[,1]),]
  fin=NULL
  i=2
  while(i<=dim(A)[1])
  {
    if(A[i-1,2]!=A[i,2])
    {
      dt=diff(c(A[i,1],A[i-1,1]))
      if(abs(dt)<tol)
      {
        if(A[i,2]>A[i-1,2])
          par=c(lc1[A[i-1,3],],lc2[A[i,3],])
        if(A[i,2]<A[i-1,2])
          par=c(lc1[A[i,3],],lc2[A[i-1,3],])
        i=i+1
      }
      else
      {
        if(A[i-1,2]==1)
          par=c(lc1[A[i-1,3],],rep(NA,3))
        if(A[i-1,2]==2)
          par=c(rep(NA,3),lc2[A[i-1,3],])
      }
    }
    else
    {
      if(A[i-1,2]==1 & A[i,2]==1)
        par=c(lc1[A[i-1,3],],rep(NA,3))
      if(A[i-1,2]==2)
        par=c(rep(NA,3),lc2[A[i-1,3],])
    }
    fin=rbind(fin,c(unlist(par)))
    i=i+1
  }
  i=i-1
  if(i==dim(A)[1])
  {
    if(A[i-1,2]==A[i,2] | (A[i-1,2]!=A[i,2] & abs(dt)>=tol))
    {
      if(A[i,2]==1)
        par=c(lc1[A[i,3],],rep(NA,3))
      if(A[i,2]==2)
        par=c(rep(NA,3),lc2[A[i,3],])
      fin=rbind(fin,c(unlist(par)))
    }
  }
  prop=dim(na.omit(fin))[1]
  x@series <- fin[,c(2,5)]
  x@series_esd <- fin[,c(3,6)]
  x@times <- apply(fin[,c(1,4)],1,mean)
  x@paired <- fin
  return(x)
}
