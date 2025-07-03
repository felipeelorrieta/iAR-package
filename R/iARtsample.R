#' Simulate from an iAR-T Model
#'
#' Simulates an iAR-T Time Series Model.
#'
#' @param coef A coefficient of iAR-T model. A value between 0 and 1.
#' @param times Array with observational times.
#' @param sigma Scale parameter of the iAR-T process. A positive value.
#' @param df degrees of freedom.
#'
#' @return A list with the following components:
#' \itemize{
#' \item{series}{ Array with simulated iAR-t process.}
#' }
#' @references
#' \insertRef{Eyheramendy_2018}{iAR}
#'
#' @seealso
#'
#' \code{\link{gentime}}
#'
#' @keywords internal
#' @examples
#' \dontshow{
#' n=300
#' set.seed(6714)
#' st<-gentime(n)
#' y<-iARtsample(coef=0.9,times=st,sigma=1,df=3)
#' plot(st,y$series,type='l')
#' hist(y$series,breaks=20)
#' }
#' @noRd
iARtsample<-function(coef,times,sigma=1,df=3)
{
  n<-length(times)
  delta<-diff(times)
  y <- vector("numeric", n)
  y[1] <- rnorm(1)
  for (i in 2:n)
  {
    phid=coef**(delta[i-1]) #Expected Value Conditional Distribution
    yhat = phid * y[i-1]  #Mean of conditional distribution
    gL=sigma*(1-phid**(2)) #Variance Value Conditional Distribution
    y[i] <- rt(1, df=df)*sqrt(gL * (df-2)/df) + yhat #Conditional-t IAR
  }
  return(list(series=y))
}
