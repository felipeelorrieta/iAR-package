#' Generating Irregularly spaced times
#'
#' A method for generating time points based on a statistical distribution.
#' The results are stored in the `times` slot of the `utilities` object.
#'
#' @param x An object of class `utilities`.
#' @param ... Additional arguments for generating time points: 
#'   \describe{
#'      \item{n}{A positive integer. Length of observation times.}
#'      \item{distribution}{A character string specifying the distribution of the observation times. Default is `"expmixture"`.
#'   Available options are:
#'   - `"expmixture"`: A mixture of two exponential distributions.
#'   - `"uniform"`: A uniform distribution.
#'   - `"exponential"`: A single exponential distribution.
#'   - `"gamma"`: A gamma distribution.}
#'      \item{lambda1}{Mean (1/rate) of the exponential distribution or the first exponential distribution in a mixture of exponential distributions. Default is `130`.}
#'      \item{lambda2}{Mean (1/rate) of the second exponential distribution in a mixture of exponential distributions. Default is `6.5`.}
#'      \item{p1}{Weight of the first exponential distribution in a mixture of exponential distributions. Default is `0.15`.}
#'      \item{p2}{Weight of the second exponential distribution in a mixture of exponential distributions. Default is `0.85`.}
#'      \item{a}{Shape parameter of a gamma distribution or lower limit of the uniform distribution. Default is `0`.}
#'      \item{b}{Scale parameter of a gamma distribution or upper limit of the uniform distribution. Default is `1`.}
#'   }
#'
#' @return An updated `utilities` object with the generated observation times stored in the `times` slot.
#'
#' @references
#' \insertRef{Eyheramendy_2018}{iAR}
#'
#' @examples
#' set.seed(12917)
#' o1=iAR::utilities()
#' o1<-gentime(o1, n=200, distribution = "expmixture", lambda1 = 130, 
#' lambda2 = 6.5,p1 = 0.15, p2 = 0.85)
#' st=o1@times
#' mean(diff(st))
#'
#' o1=iAR::utilities()
#' o1<-gentime(o1, n=200, distribution = "expmixture", lambda1 = 15, 
#' lambda2 = 2.5,p1 = 0.15, p2 = 0.85)
#' st=o1@times
#' mean(diff(st))
#'
#' @export
gentime <- S7::new_generic("gentime", "x")
S7::method(gentime, utilities) <- function(x, n, distribution = "expmixture", lambda1 = 130, 
                                           lambda2 = 6.5, p1 = 0.15, p2 = 0.85, a = 0, b = 1) {
  if (!is.numeric(n) || length(n) != 1 || n <= 0) {
    stop("'n' must be a positive integer.")
  }
  if(distribution=="expmixture")
  {
    dT <- rep(0, n)
    a <- sample(c(lambda1, lambda2), size = n, prob = c(p1, p2),
                replace = TRUE)
    dT[which(a == lambda1)] = rexp(length(which(a == lambda1)),
                                   rate = 1/lambda1)
    dT[which(a == lambda2)] = rexp(length(which(a == lambda2)),
                                   rate = 1/lambda2)
    sT = cumsum(dT)
  }
  if(distribution=="uniform")
  {
    sT=cumsum(runif(n,a,b))
  }
  if(distribution=="exponential")
  {
    sT=cumsum(rexp(n,rate=1/lambda1))
  }
  if(distribution=="gamma")
  {
    sT=cumsum(rgamma(n,shape=a,scale=b))
  }
  x@times <- sT
  return(x)
}  