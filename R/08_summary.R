#' Summary Method for unidata Objects
#'
#' 
#' This method provides a summary of the `unidata` object, which represents 
#' an iAR (irregular autoregressive) model. The summary includes information 
#' about the model family (e.g., "normal", "t", "gamma"), the coefficients, 
#' standard errors, t-values, and p-values. The output is formatted to display 
#' the relevant statistics based on the model family.
#' @name summary
#' 
#' @param object An object of class `unidata`. This object contains the fitted 
#'   iAR model, including parameters such as coefficients (`coef`), 
#'   standard errors (`stderr`), t-values (`tvalue`), p-values (`pvalue`), 
#'   family type (`family`), and other model-specific parameters.
#' @param ... Additional arguments (unused). 
#'
#'
#' @export
S7::method(summary, unidata) <- function(object,...) {
  x=object
  cat("iAR model", "\n", sep = "")
  cat("\n", sep = "")
  if(x@family == "norm"){
    cat("Family:", "\n ", sep = "")
    cat("", x@family, "\n", sep = "")
    cat("\n", sep = "")
    cat("Coefficients:", "\n", sep = "")
    cat(paste("    ", "Estimate St. Error t value Pr(>|t|)", "\n", sep = ""))
    cat("phi", "     ", format(round(x@coef,2),nsmall=2), "      ", format(round(x@summary$stderr,2),nsmall=2), "   ", format(round(x@summary$tvalue,2),nsmall=2), "     ", format(round(x@summary$pvalue,2),nsmall=2), "\n", sep = "")
  }
  if(x@family == "t"){
    cat("Family:", "\n ", sep = "")
    cat("", x@family, "\n", sep = "")
    cat("\n", sep = "")
    cat("d.f.:", "\n ", sep = "")
    cat("", x@df, "\n", sep = "")
    cat("\n", sep = "")
    cat("Coefficients:", "\n", sep = "")
    cat(paste("      ", "Estimate St. Error t value Pr(>|t|)", "\n", sep = ""))
    cat("phi", "       ", format(round(x@coef,2),nsmall=2), "      ", format(round(x@summary$stderr[1],2),nsmall=2), "   ", format(round(x@summary$tvalue[1],2),nsmall=2), "     ", format(round(x@summary$pvalue[1],2),nsmall=2), "\n", sep = "")
    cat("sigma", "     ", format(round(x@sigma,2),nsmall=2), "      ", format(round(x@summary$stderr[2],2),nsmall=2), "    ", format(round(x@summary$tvalue[2],2),nsmall=2), "     ", format(round(x@summary$pvalue[2],2),nsmall=2), "\n", sep = "")
  }
  if(x@family == "gamma"){
    cat("Family:", "\n ", sep = "")
    cat("", x@family, "\n", sep = "")
    cat("\n", sep = "")
    cat("Coefficients:", "\n", sep = "")
    cat(paste("         ", "Estimate St. Error t value Pr(>|t|)", "\n", sep = ""))
    cat("phi", "          ", format(round(x@coef,2),nsmall=2), "      ", format(round(x@summary$stderr[1],2),nsmall=2), "   ", format(round(x@summary$tvalue[1],2),nsmall=2), "     ", format(round(x@summary$pvalue[1],2),nsmall=2), "\n", sep = "")
    cat("mean", "         ", format(round(x@mean,2),nsmall=2), "      ", format(round(x@summary$stderr[2],2),nsmall=2), "   ", format(round(x@summary$tvalue[2],2),nsmall=2), "     ", format(round(x@summary$pvalue[2],2),nsmall=2), "\n", sep = "")
    cat("variance", "     ", format(round(x@variance,2),nsmall=2), "      ", format(round(x@summary$stderr[3],2),nsmall=2), "    ", format(round(x@summary$tvalue[3],2),nsmall=2), "     ", format(round(x@summary$pvalue[3],2),nsmall=2), "\n", sep = "")
  }
}
